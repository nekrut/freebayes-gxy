//! Per-read CIGAR walker that emits [`AlleleObservation`]s.
//!
//! Port of upstream `AlleleParser::registerAlignment` (`src/AlleleParser.cpp`
//! lines 1329..1838) — the inner loop that converts one aligned read into a
//! vector of allele observations anchored against the reference.
//!
//! ## Coordinates
//! Positions are 0-based and half-open, matching upstream's internal layout
//! (upstream uses 0-based SAM `POSITION`). VCF-level 1-based emission is M4's
//! job.
//!
//! ## M1 scope
//! SNP / Reference / Insertion / Deletion observations only. MNP and Complex
//! alleles are emitted only after haplotype-aware collapsing (M2). Soft-clip
//! `ALLELE_NULL` observations that upstream records for accounting are also
//! deferred — the downstream Bayesian model (M3) does not consume them.

use crate::allele::{Allele, AlleleObservation, Strand};
use rust_htslib::bam::record::{Cigar, CigarStringView};
use rust_htslib::bam::Record;

// ---------------------------------------------------------------------------
// Read filter
// ---------------------------------------------------------------------------

/// Read-level filter criteria. Defaults mirror the freebayes CLI defaults
/// (`src/Parameters.cpp` `setDefaults`): `--min-mapping-quality 1`, keep all
/// base qualities (`--min-base-quality 0`), exclude duplicates and
/// secondary/supplementary/qcfail alignments.
#[derive(Debug, Clone, Copy)]
pub struct ReadFilter {
    pub min_mapping_quality: u8,
    pub min_base_quality: u8,
    pub exclude_duplicates: bool,
    pub exclude_secondary: bool,
    pub exclude_supplementary: bool,
    pub exclude_qc_fail: bool,
    pub exclude_unmapped: bool,
    pub require_proper_pair: bool,
}

impl Default for ReadFilter {
    fn default() -> Self {
        Self {
            min_mapping_quality: 1,
            min_base_quality: 0,
            exclude_duplicates: true,
            exclude_secondary: true,
            exclude_supplementary: true,
            exclude_qc_fail: true,
            exclude_unmapped: true,
            require_proper_pair: false,
        }
    }
}

// BAM flag bits (SAM spec §1.4.2).
const FLAG_UNMAPPED: u16 = 0x4;
const FLAG_REVERSE: u16 = 0x10;
const FLAG_PROPER_PAIR: u16 = 0x2;
const FLAG_SECONDARY: u16 = 0x100;
const FLAG_QC_FAIL: u16 = 0x200;
const FLAG_DUPLICATE: u16 = 0x400;
const FLAG_SUPPLEMENTARY: u16 = 0x800;

impl ReadFilter {
    /// Returns `true` if the alignment passes every enabled check.
    fn admits(&self, mapq: u8, flags: u16) -> bool {
        if mapq < self.min_mapping_quality {
            return false;
        }
        if self.exclude_unmapped && flags & FLAG_UNMAPPED != 0 {
            return false;
        }
        if self.exclude_duplicates && flags & FLAG_DUPLICATE != 0 {
            return false;
        }
        if self.exclude_secondary && flags & FLAG_SECONDARY != 0 {
            return false;
        }
        if self.exclude_supplementary && flags & FLAG_SUPPLEMENTARY != 0 {
            return false;
        }
        if self.exclude_qc_fail && flags & FLAG_QC_FAIL != 0 {
            return false;
        }
        if self.require_proper_pair && flags & FLAG_PROPER_PAIR == 0 {
            return false;
        }
        true
    }
}

// ---------------------------------------------------------------------------
// Alignment view (test-friendly alias for the bits of a BAM record we need)
// ---------------------------------------------------------------------------

/// A minimal view of an aligned read. Letting the walker operate on this view
/// instead of directly on [`Record`] keeps the unit tests in this file
/// hermetic — we construct `AlignmentView`s from raw bytes without standing
/// up a BAM header or serialising records. `walk_record` is the thin wrapper
/// that extracts these fields from a real [`Record`].
#[derive(Debug, Clone)]
pub struct AlignmentView<'a> {
    pub read_name: &'a str,
    /// Decoded query bases, one byte per base (`b'A'`, `b'C'`, `b'G'`, `b'T'`,
    /// or `b'N'`). Not 4-bit packed.
    pub seq: &'a [u8],
    /// Phred-scaled qualities (raw, not offset by 33).
    pub quals: &'a [u8],
    pub cigar: &'a [Cigar],
    /// 0-based leftmost reference position of the first aligned base.
    pub pos: i64,
    pub mapq: u8,
    pub flags: u16,
}

impl<'a> AlignmentView<'a> {
    #[inline]
    fn is_reverse(&self) -> bool {
        self.flags & FLAG_REVERSE != 0
    }
    #[inline]
    fn is_proper_pair(&self) -> bool {
        self.flags & FLAG_PROPER_PAIR != 0
    }
    #[inline]
    fn strand(&self) -> Strand {
        if self.is_reverse() {
            Strand::Reverse
        } else {
            Strand::Forward
        }
    }
}

// ---------------------------------------------------------------------------
// CIGAR walker
// ---------------------------------------------------------------------------

/// Walk the CIGAR of a real BAM record and emit allele observations.
///
/// `ref_seq` is the reference sequence for the contig (or a window within it);
/// `ref_start` is the 0-based reference coordinate corresponding to
/// `ref_seq[0]`. Returns `None` if the read is filtered out by `filter`.
pub fn walk_record(
    record: &Record,
    ref_seq: &[u8],
    ref_start: i64,
    filter: &ReadFilter,
) -> Option<Vec<AlleleObservation>> {
    let cigar_owned: CigarStringView = record.cigar();
    let cigar: Vec<Cigar> = cigar_owned.iter().copied().collect();
    let read_name = std::str::from_utf8(record.qname()).unwrap_or("");
    let seq = record.seq().as_bytes();
    let quals = record.qual().to_vec();
    let view = AlignmentView {
        read_name,
        seq: &seq,
        quals: &quals,
        cigar: &cigar,
        pos: record.pos(),
        mapq: record.mapq(),
        flags: record.flags(),
    };
    walk_alignment(&view, ref_seq, ref_start, filter)
}

/// Walk an [`AlignmentView`] and emit allele observations. See [`walk_record`]
/// for the BAM-record wrapper.
pub fn walk_alignment(
    read: &AlignmentView<'_>,
    ref_seq: &[u8],
    ref_start: i64,
    filter: &ReadFilter,
) -> Option<Vec<AlleleObservation>> {
    if !filter.admits(read.mapq, read.flags) {
        return None;
    }

    let mut out: Vec<AlleleObservation> = Vec::new();

    // See AlleleParser.cpp:1340-1342 — freebayes' naming:
    //   rp  = read position (query cursor, 0-based)
    //   sp  = sample position (reference cursor, 0-based, absolute)
    //   csp = current sequence position (offset into the reference-seq slice)
    let mut rp: usize = 0;
    let mut sp: i64 = read.pos;
    let mut csp: i64 = read.pos - ref_start;

    let strand = read.strand();
    let is_proper_pair = read.is_proper_pair();
    let mapq = read.mapq;
    let read_name = read.read_name.to_string();

    let n_ops = read.cigar.len();
    for (op_idx, op) in read.cigar.iter().enumerate() {
        match *op {
            // Alignment-match / sequence-match / sequence-mismatch: per-base
            // walk that mirrors upstream AlleleParser.cpp:1430..1625.
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                walk_match_run(
                    &mut out,
                    read,
                    ref_seq,
                    &read_name,
                    mapq,
                    strand,
                    is_proper_pair,
                    len as usize,
                    &mut rp,
                    &mut sp,
                    &mut csp,
                );
            }
            // Insertion: bases present in read but not reference.
            // AlleleParser.cpp:1706..1773. Upstream requires `allATGC(readseq)`
            // and a non-empty quality scale; we mirror both.
            Cigar::Ins(len) => {
                let l = len as usize;
                if rp + l <= read.seq.len() {
                    let alt = read.seq[rp..rp + l].to_vec();
                    if all_atgc(&alt) {
                        let bq_sum = sum_u8(&read.quals[rp..rp + l]);
                        out.push(AlleleObservation {
                            allele: Allele::insertion(sp, alt),
                            read_name: read_name.clone(),
                            mapq,
                            // TODO(M3): upstream applies a harmonic-sum
                            // scaling plus a length-ratio adjustment (see
                            // AlleleParser.cpp:1741..1749). For M1 we keep the
                            // plain sum; the Bayesian model (M3) is where the
                            // scaled value is consumed and we revisit it.
                            base_quality_sum: bq_sum,
                            strand,
                            read_position: rp,
                            is_proper_pair,
                            read_ref_start: read.pos,
                            per_base_quals: Vec::new(),
                        });
                    }
                }
                rp += l;
            }
            // Deletion: bases present in reference but not read.
            // AlleleParser.cpp:1626..1704. Upstream suppresses deletions at
            // the very start or end of the read ("without any sequence in the
            // read to support this, it is hard to believe these deletions are
            // real", l.1677) and requires `allATGC(refseq)`.
            Cigar::Del(len) => {
                let l = len as usize;
                let is_edge = op_idx == 0 || op_idx + 1 == n_ops;
                let usize_csp = csp.max(0) as usize;
                if !is_edge && usize_csp + l <= ref_seq.len() {
                    let r = ref_seq[usize_csp..usize_csp + l].to_vec();
                    if all_atgc(&r) {
                        // Flanking-base BQ as a proxy for deletion quality
                        // (AlleleParser.cpp:1638-1673 uses a harmonic-sum over
                        // an l+2 span; M3 will restore the exact formula).
                        let left = if rp > 0 { read.quals[rp - 1] as u32 } else { 0 };
                        let right = if rp < read.quals.len() {
                            read.quals[rp] as u32
                        } else {
                            0
                        };
                        let bq = if rp > 0 && rp < read.quals.len() {
                            (left + right) / 2
                        } else {
                            left + right
                        };
                        out.push(AlleleObservation {
                            allele: Allele::deletion(sp, r),
                            read_name: read_name.clone(),
                            mapq,
                            base_quality_sum: bq,
                            strand,
                            read_position: rp,
                            is_proper_pair,
                            read_ref_start: read.pos,
                            per_base_quals: Vec::new(),
                        });
                    }
                }
                sp += l as i64;
                csp += l as i64;
            }
            // Soft clip: advance the read cursor only. AlleleParser.cpp:1776.
            // Upstream emits a NULL observation here for accounting; M1 drops
            // it because the Bayesian model (M3) does not read soft-clip
            // observations.
            Cigar::SoftClip(len) => {
                rp += len as usize;
            }
            // Hard clip: nothing to do; clipped sequence is absent from
            // `read.seq`. AlleleParser.cpp:1800.
            Cigar::HardClip(_) => {}
            // Reference skip (splice): advance reference cursor only.
            // AlleleParser.cpp:1805-1825.
            Cigar::RefSkip(len) => {
                sp += len as i64;
                csp += len as i64;
            }
            // Padding: no effect (upstream has the same handling commented
            // out at AlleleParser.cpp:1828).
            Cigar::Pad(_) => {}
        }
    }

    Some(out)
}

/// Walk a contiguous `M`/`=`/`X` run, emitting per-run Reference and per-base
/// SNP observations. Ports the inner `for (int i=0; i<l; i++)` loop plus the
/// `inMismatch` / `firstMatch` state machine at
/// `AlleleParser.cpp:1432..1624`.
///
/// Semantics matched to upstream:
/// - A contiguous run of matching ATGC bases produces **one**
///   `AlleleKind::Reference` observation whose `ref_seq` spans the run
///   (upstream emits `ALLELE_REFERENCE` with `length = run_len`, quality set
///   to `MAPPINGQUALITY` — lines 1472-1487 and 1608-1624).
/// - Every mismatched position produces one SNP observation, *regardless* of
///   its base quality — upstream only uses `BQL2` to gate the `mismatches` /
///   `snpCount` counters, not to suppress the allele (lines 1491-1503).
/// - A reference `N` at a matched position is treated as a forced mismatch
///   (line 1465: `if (b != sb || sb == "N")`); the port emits an SNP with
///   `ref = 'N'` and `alt = read_base` when the read base is ATGC.
/// - Non-ATGC read bases at match positions become `ALLELE_NULL` observations
///   (lines 1532-1548).
#[allow(clippy::too_many_arguments)]
fn walk_match_run(
    out: &mut Vec<AlleleObservation>,
    read: &AlignmentView<'_>,
    ref_seq: &[u8],
    read_name: &str,
    mapq: u8,
    strand: Strand,
    is_proper_pair: bool,
    len: usize,
    rp: &mut usize,
    sp: &mut i64,
    csp: &mut i64,
) {
    // Start position of the current in-progress reference run, expressed as
    // (ref_pos, read_pos). `None` means no run is in progress. Mirrors
    // upstream's `firstMatch` cursor (AlleleParser.cpp:1466-1488).
    let mut run_start: Option<(i64, usize)> = None;

    let flush_ref_run =
        |out: &mut Vec<AlleleObservation>, run_start: &mut Option<(i64, usize)>, end_sp: i64| {
            if let Some((start_sp, start_rp)) = run_start.take() {
                let length = (end_sp - start_sp) as usize;
                if length == 0 {
                    return;
                }
                let run_bases = read.seq[start_rp..start_rp + length].to_vec();
                let per_base_quals = read.quals[start_rp..start_rp + length].to_vec();
                debug_assert_eq!(run_bases.len(), length);
                debug_assert_eq!(per_base_quals.len(), length);
                // Sum of per-base BQ across the run — VCF `QR` field.
                // Upstream stores `MAPPINGQUALITY` as the *scalar*
                // quality on the Allele (AlleleParser.cpp:1484, 1621)
                // but accumulates per-base BQ into QR. We capture both:
                // the per-base vec drives downstream per-position BQ
                // after the fb-cli pileup decomposes the run, and the
                // sum here matches what upstream's VCF reports.
                let bq_sum: u32 = per_base_quals.iter().map(|&b| b as u32).sum();
                out.push(AlleleObservation {
                    allele: Allele::reference(start_sp, run_bases),
                    read_name: read_name.to_string(),
                    mapq,
                    base_quality_sum: bq_sum,
                    strand,
                    read_position: start_rp,
                    is_proper_pair,
                    read_ref_start: read.pos,
                    per_base_quals,
                });
            }
        };

    for _ in 0..len {
        if *rp >= read.seq.len() || *csp < 0 || *csp as usize >= ref_seq.len() {
            // Truncated read / ref; flush in-progress run and bail.
            flush_ref_run(out, &mut run_start, *sp);
            return;
        }
        let read_base = read.seq[*rp];
        let ref_base = ref_seq[*csp as usize];
        let bq = read.quals[*rp];

        let ref_is_n = !is_atgc(ref_base);
        let read_is_atgc = is_atgc(read_base);
        let is_match = read_is_atgc && !ref_is_n && ascii_eq_ignore_case(read_base, ref_base);

        if is_match {
            // Extend (or begin) the current reference run.
            run_start.get_or_insert((*sp, *rp));
        } else {
            // Flush any pending reference run first — upstream does the
            // same at AlleleParser.cpp:1466-1488.
            flush_ref_run(out, &mut run_start, *sp);

            // Classify the mismatch. Upstream (AlleleParser.cpp:1515-1548)
            // emits either an ALLELE_SNP (for ATGC read bases) or an
            // ALLELE_NULL (for non-ATGC reads inside a mismatch region);
            // `ref == 'N'` routes through the SNP branch because the read
            // base, when ATGC, is still an informative observation.
            if read_is_atgc {
                out.push(AlleleObservation {
                    allele: Allele::snp(*sp, ref_base, read_base),
                    read_name: read_name.to_string(),
                    mapq,
                    base_quality_sum: bq as u32,
                    strand,
                    read_position: *rp,
                    is_proper_pair,
                    read_ref_start: read.pos,
                    per_base_quals: Vec::new(),
                });
            } else {
                // Non-ATGC read base (typically 'N'). Upstream emits an
                // ALLELE_NULL observation; the Bayesian model (M3) ignores
                // these. We preserve the record for future use.
                out.push(AlleleObservation {
                    allele: Allele::null(*sp, 1),
                    read_name: read_name.to_string(),
                    mapq,
                    base_quality_sum: bq as u32,
                    strand,
                    read_position: *rp,
                    is_proper_pair,
                    read_ref_start: read.pos,
                    per_base_quals: Vec::new(),
                });
            }
        }

        *rp += 1;
        *sp += 1;
        *csp += 1;
    }

    // Flush any trailing in-progress reference run at end of the CIGAR op.
    // Upstream does this at AlleleParser.cpp:1602-1624.
    flush_ref_run(out, &mut run_start, *sp);
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

#[inline]
fn is_atgc(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')
}

#[inline]
fn all_atgc(s: &[u8]) -> bool {
    s.iter().all(|&b| is_atgc(b))
}

#[inline]
fn ascii_eq_ignore_case(a: u8, b: u8) -> bool {
    a.eq_ignore_ascii_case(&b)
}

#[inline]
fn sum_u8(xs: &[u8]) -> u32 {
    xs.iter().map(|&x| x as u32).sum()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::AlleleKind;

    fn view<'a>(
        seq: &'a [u8],
        quals: &'a [u8],
        cigar: &'a [Cigar],
        pos: i64,
        mapq: u8,
        flags: u16,
    ) -> AlignmentView<'a> {
        AlignmentView {
            read_name: "r1",
            seq,
            quals,
            cigar,
            pos,
            mapq,
            flags,
        }
    }

    #[test]
    fn pure_match_all_reference() {
        let seq = b"ACGTACGTAC";
        let quals = vec![30u8; 10];
        let cigar = [Cigar::Match(10)];
        let v = view(seq, &quals, &cigar, 100, 60, 0);
        let obs = walk_alignment(&v, b"ACGTACGTAC", 100, &ReadFilter::default()).unwrap();
        // Per-run emission: one Reference observation covering all 10 bases.
        assert_eq!(obs.len(), 1);
        assert_eq!(obs[0].allele.kind, AlleleKind::Reference);
        assert_eq!(obs[0].allele.position, 100);
        assert_eq!(obs[0].allele.length, 10);
        assert_eq!(obs[0].allele.ref_seq, b"ACGTACGTAC");
        assert_eq!(obs[0].read_position, 0);
        // Per-base BQ fix (M4-E): the walker now sums per-base Phred
        // scores across the run rather than storing MAPQ as a scalar.
        // 10 bases × Q30 = 300. `per_base_quals` carries the breakdown.
        assert_eq!(obs[0].base_quality_sum, 300);
        assert_eq!(obs[0].per_base_quals, vec![30u8; 10]);
    }

    #[test]
    fn pure_match_with_three_snps() {
        // ref  A C G T A C G T A C
        // read A A G T T C G C A C  (positions 1, 4, 7 differ)
        let seq = b"AAGTTCGCAC";
        let quals = vec![30u8; 10];
        let cigar = [Cigar::Match(10)];
        let v = view(seq, &quals, &cigar, 100, 60, 0);
        let obs = walk_alignment(&v, b"ACGTACGTAC", 100, &ReadFilter::default()).unwrap();

        let snps: Vec<_> = obs
            .iter()
            .filter(|o| o.allele.kind == AlleleKind::Snp)
            .collect();
        let refs: Vec<_> = obs
            .iter()
            .filter(|o| o.allele.kind == AlleleKind::Reference)
            .collect();
        // 3 SNPs at 101, 104, 107.
        assert_eq!(snps.len(), 3);
        assert_eq!(snps[0].allele.position, 101);
        assert_eq!(snps[0].allele.ref_seq, b"C");
        assert_eq!(snps[0].allele.alt_seq, b"A");
        assert_eq!(snps[1].allele.position, 104);
        assert_eq!(snps[2].allele.position, 107);
        // 4 Reference runs: [100..101]=1, [102..104]=2, [105..107]=2, [108..110]=2
        assert_eq!(refs.len(), 4);
        let lens: Vec<usize> = refs.iter().map(|o| o.allele.length).collect();
        assert_eq!(lens, vec![1, 2, 2, 2]);
        let positions: Vec<i64> = refs.iter().map(|o| o.allele.position).collect();
        assert_eq!(positions, vec![100, 102, 105, 108]);
    }

    #[test]
    fn single_insertion() {
        // ref: AAAAAAAAAA             (len 10, positions 100..110)
        // read: AAAAA GG AAAAA         (CIGAR 5M 2I 5M)
        let seq = b"AAAAAGGAAAAA";
        let quals = vec![40u8; 12];
        let cigar = [Cigar::Match(5), Cigar::Ins(2), Cigar::Match(5)];
        let v = view(seq, &quals, &cigar, 100, 60, 0);
        let obs = walk_alignment(&v, b"AAAAAAAAAA", 100, &ReadFilter::default()).unwrap();

        // Expected: Ref[100..105] len=5, Ins @105 alt="GG", Ref[105..110] len=5.
        let ins: Vec<_> = obs
            .iter()
            .filter(|o| o.allele.kind == AlleleKind::Insertion)
            .collect();
        assert_eq!(ins.len(), 1);
        assert_eq!(ins[0].allele.position, 105);
        assert!(ins[0].allele.ref_seq.is_empty());
        assert_eq!(ins[0].allele.alt_seq, b"GG");
        assert_eq!(ins[0].base_quality_sum, 80); // 40 + 40
        assert_eq!(ins[0].read_position, 5);

        let refs: Vec<_> = obs
            .iter()
            .filter(|o| o.allele.kind == AlleleKind::Reference)
            .collect();
        assert_eq!(refs.len(), 2);
        assert_eq!(refs[0].allele.position, 100);
        assert_eq!(refs[0].allele.length, 5);
        assert_eq!(refs[1].allele.position, 105);
        assert_eq!(refs[1].allele.length, 5);
        assert_eq!(refs[1].read_position, 7); // after 5M + 2I
    }

    #[test]
    fn single_deletion() {
        // ref: AAAAA CCC AAAAA           (positions 100..113)
        // read: AAAAAAAAAA                (CIGAR 5M 3D 5M)
        let seq = b"AAAAAAAAAA";
        let quals = vec![40u8; 10];
        let cigar = [Cigar::Match(5), Cigar::Del(3), Cigar::Match(5)];
        let v = view(seq, &quals, &cigar, 100, 60, 0);
        let obs = walk_alignment(&v, b"AAAAACCCAAAAA", 100, &ReadFilter::default()).unwrap();

        let dels: Vec<_> = obs
            .iter()
            .filter(|o| o.allele.kind == AlleleKind::Deletion)
            .collect();
        assert_eq!(dels.len(), 1);
        assert_eq!(dels[0].allele.position, 105);
        assert_eq!(dels[0].allele.ref_seq, b"CCC");
        assert!(dels[0].allele.alt_seq.is_empty());
        assert_eq!(dels[0].base_quality_sum, 40); // (40 + 40) / 2

        // Two reference runs (5M + 5M), one Deletion between them.
        let refs: Vec<_> = obs
            .iter()
            .filter(|o| o.allele.kind == AlleleKind::Reference)
            .collect();
        assert_eq!(refs.len(), 2);
        assert_eq!(refs[0].allele.position, 100);
        assert_eq!(refs[0].allele.length, 5);
        assert_eq!(refs[1].allele.position, 108);
        assert_eq!(refs[1].allele.length, 5);
    }

    #[test]
    fn soft_clip_advances_read_only() {
        // 3S 10M 3S, seq length 16; aligned portion matches ref.
        let seq = b"NNNACGTACGTACNNN";
        let mut quals = vec![30u8; 16];
        quals[0] = 5;
        quals[15] = 5;
        let cigar = [Cigar::SoftClip(3), Cigar::Match(10), Cigar::SoftClip(3)];
        let v = view(seq, &quals, &cigar, 100, 60, 0);
        let obs = walk_alignment(&v, b"ACGTACGTAC", 100, &ReadFilter::default()).unwrap();
        // One reference run of length 10, read offset starts after the leading 3S.
        assert_eq!(obs.len(), 1);
        assert_eq!(obs[0].allele.kind, AlleleKind::Reference);
        assert_eq!(obs[0].allele.length, 10);
        assert_eq!(obs[0].allele.position, 100);
        assert_eq!(obs[0].read_position, 3);
    }

    #[test]
    fn ref_skip_advances_reference_only() {
        // 5M 10N 5M — RNA-seq-like split read.
        let seq = b"ACGTAACGTA";
        let quals = vec![30u8; 10];
        let cigar = [Cigar::Match(5), Cigar::RefSkip(10), Cigar::Match(5)];
        let v = view(seq, &quals, &cigar, 0, 60, 0);
        let ref_seq = b"ACGTANNNNNNNNNNACGTA";
        let obs = walk_alignment(&v, ref_seq, 0, &ReadFilter::default()).unwrap();
        // Two reference runs, one either side of the splice junction.
        assert_eq!(obs.len(), 2);
        assert_eq!(obs[0].allele.position, 0);
        assert_eq!(obs[0].allele.length, 5);
        assert_eq!(obs[1].allele.position, 15);
        assert_eq!(obs[1].allele.length, 5);
    }

    #[test]
    fn multi_event_read_preserves_order() {
        // 5M 2I 5M 3D 5M — seq length 17, ref span 18.
        let seq = b"AAAAAGGAAAAAAAAAA";
        let quals = vec![40u8; 17];
        let cigar = [
            Cigar::Match(5),
            Cigar::Ins(2),
            Cigar::Match(5),
            Cigar::Del(3),
            Cigar::Match(5),
        ];
        let v = view(seq, &quals, &cigar, 100, 60, 0);
        let ref_seq = b"AAAAAAAAAACCCAAAAA";
        let obs = walk_alignment(&v, ref_seq, 100, &ReadFilter::default()).unwrap();
        let kinds: Vec<AlleleKind> = obs.iter().map(|o| o.allele.kind).collect();
        // Ref(5), Ins, Ref(5), Del, Ref(5) — 5 observations total.
        let expected = [
            AlleleKind::Reference,
            AlleleKind::Insertion,
            AlleleKind::Reference,
            AlleleKind::Deletion,
            AlleleKind::Reference,
        ];
        assert_eq!(kinds, expected);
        // Ref run positions/lengths: 100..105, 105..110, 113..118.
        assert_eq!(obs[0].allele.position, 100);
        assert_eq!(obs[0].allele.length, 5);
        assert_eq!(obs[1].allele.position, 105); // insertion
        assert_eq!(obs[1].allele.alt_seq, b"GG");
        assert_eq!(obs[2].allele.position, 105);
        assert_eq!(obs[2].allele.length, 5);
        assert_eq!(obs[3].allele.position, 110); // deletion
        assert_eq!(obs[3].allele.ref_seq, b"CCC");
        assert_eq!(obs[4].allele.position, 113);
        assert_eq!(obs[4].allele.length, 5);
    }

    #[test]
    fn edge_deletion_at_read_start_is_dropped() {
        // Upstream (AlleleParser.cpp:1681-1683) refuses to call deletions
        // that sit against the first or last CIGAR op.
        let seq = b"AAAAA";
        let quals = vec![40u8; 5];
        let cigar = [Cigar::Del(3), Cigar::Match(5)];
        let v = view(seq, &quals, &cigar, 100, 60, 0);
        let obs = walk_alignment(&v, b"CCCAAAAA", 100, &ReadFilter::default()).unwrap();
        assert!(obs.iter().all(|o| o.allele.kind != AlleleKind::Deletion));
    }

    #[test]
    fn low_mapq_read_is_filtered() {
        let seq = b"ACGT";
        let quals = vec![30u8; 4];
        let cigar = [Cigar::Match(4)];
        let v = view(seq, &quals, &cigar, 0, 0, 0); // mapq 0
        let filter = ReadFilter::default(); // min_mapping_quality = 1
        assert!(walk_alignment(&v, b"ACGT", 0, &filter).is_none());
    }

    #[test]
    fn duplicate_flag_is_filtered() {
        let seq = b"ACGT";
        let quals = vec![30u8; 4];
        let cigar = [Cigar::Match(4)];
        let v = view(seq, &quals, &cigar, 0, 60, FLAG_DUPLICATE);
        assert!(walk_alignment(&v, b"ACGT", 0, &ReadFilter::default()).is_none());
    }

    #[test]
    fn secondary_alignment_is_filtered() {
        let seq = b"ACGT";
        let quals = vec![30u8; 4];
        let cigar = [Cigar::Match(4)];
        let v = view(seq, &quals, &cigar, 0, 60, FLAG_SECONDARY);
        assert!(walk_alignment(&v, b"ACGT", 0, &ReadFilter::default()).is_none());
    }

    #[test]
    fn low_bq_snp_is_still_emitted() {
        // Upstream AlleleParser.cpp:1491-1503: BQL2 only gates counters
        // (mismatches/snpCount), it never suppresses the SNP observation.
        // The port must mirror that — the per-SNP observation goes out
        // regardless of base quality. `min_base_quality` on `ReadFilter`
        // is retained for future counter-gating (M3).
        let ref_seq = b"ACGTA";
        let seq = b"ACATA"; // G->A at position 2
        let quals = [30u8, 30, 5, 30, 30]; // Q5 mismatch
        let cigar = [Cigar::Match(5)];
        let filter = ReadFilter {
            min_base_quality: 20,
            ..ReadFilter::default()
        };
        let v = view(seq, &quals, &cigar, 0, 60, 0);
        let obs = walk_alignment(&v, ref_seq, 0, &filter).unwrap();
        let snps: Vec<_> = obs
            .iter()
            .filter(|o| o.allele.kind == AlleleKind::Snp)
            .collect();
        assert_eq!(snps.len(), 1, "low-BQ SNP must still be emitted");
        assert_eq!(snps[0].allele.position, 2);
        assert_eq!(snps[0].allele.ref_seq, b"G");
        assert_eq!(snps[0].allele.alt_seq, b"A");
        assert_eq!(snps[0].base_quality_sum, 5);
    }

    #[test]
    fn ref_n_is_forced_mismatch() {
        // Upstream AlleleParser.cpp:1465: `if (b != sb || sb == "N")` — an
        // 'N' in the reference is ALWAYS classified as a mismatch, even if
        // the read base is ATGC. The SNP carries ref='N'/alt=read_base.
        let ref_seq = b"ACNTA";
        let seq = b"ACGTA";
        let quals = vec![30u8; 5];
        let cigar = [Cigar::Match(5)];
        let v = view(seq, &quals, &cigar, 0, 60, 0);
        let obs = walk_alignment(&v, ref_seq, 0, &ReadFilter::default()).unwrap();
        let snps: Vec<_> = obs
            .iter()
            .filter(|o| o.allele.kind == AlleleKind::Snp)
            .collect();
        assert_eq!(snps.len(), 1);
        assert_eq!(snps[0].allele.position, 2);
        assert_eq!(snps[0].allele.ref_seq, b"N");
        assert_eq!(snps[0].allele.alt_seq, b"G");
    }

    #[test]
    fn read_n_at_match_emits_null() {
        // Upstream AlleleParser.cpp:1532-1548: non-ATGC read base inside a
        // match region becomes an ALLELE_NULL observation, not an SNP.
        let ref_seq = b"ACGTA";
        let seq = b"ACNTA"; // read has N at position 2
        let quals = vec![30u8; 5];
        let cigar = [Cigar::Match(5)];
        let v = view(seq, &quals, &cigar, 0, 60, 0);
        let obs = walk_alignment(&v, ref_seq, 0, &ReadFilter::default()).unwrap();
        let nulls: Vec<_> = obs
            .iter()
            .filter(|o| o.allele.kind == AlleleKind::Null)
            .collect();
        assert_eq!(nulls.len(), 1);
        assert_eq!(nulls[0].allele.position, 2);
        // Two reference runs flanking the N.
        let refs: Vec<_> = obs
            .iter()
            .filter(|o| o.allele.kind == AlleleKind::Reference)
            .collect();
        assert_eq!(refs.len(), 2);
        assert_eq!(refs[0].allele.length, 2);
        assert_eq!(refs[1].allele.length, 2);
    }

    #[test]
    fn edge_deletion_at_read_end_is_dropped() {
        // Symmetrical to `edge_deletion_at_read_start_is_dropped`: upstream
        // also refuses deletions that sit against the final CIGAR op.
        let seq = b"AAAAA";
        let quals = vec![40u8; 5];
        let cigar = [Cigar::Match(5), Cigar::Del(3)];
        let v = view(seq, &quals, &cigar, 100, 60, 0);
        let obs = walk_alignment(&v, b"AAAAACCC", 100, &ReadFilter::default()).unwrap();
        assert!(obs.iter().all(|o| o.allele.kind != AlleleKind::Deletion));
    }

    #[test]
    fn hard_clip_is_silent() {
        // 2H 5M 2H: hard clips aren't in the read sequence — the walker
        // should emit only the 5M reference run.
        let seq = b"ACGTA";
        let quals = vec![30u8; 5];
        let cigar = [Cigar::HardClip(2), Cigar::Match(5), Cigar::HardClip(2)];
        let v = view(seq, &quals, &cigar, 10, 60, 0);
        let obs = walk_alignment(&v, b"ACGTA", 10, &ReadFilter::default()).unwrap();
        assert_eq!(obs.len(), 1);
        assert_eq!(obs[0].allele.kind, AlleleKind::Reference);
        assert_eq!(obs[0].allele.position, 10);
        assert_eq!(obs[0].allele.length, 5);
    }

    #[test]
    fn reverse_strand_is_recorded() {
        let seq = b"ACGT";
        let quals = vec![30u8; 4];
        let cigar = [Cigar::Match(4)];
        let v = view(seq, &quals, &cigar, 0, 60, FLAG_REVERSE);
        let obs = walk_alignment(&v, b"ACGT", 0, &ReadFilter::default()).unwrap();
        assert!(obs.iter().all(|o| o.strand == Strand::Reverse));
    }
}
