//! VCF 4.2 writer for freebayes-gxy.
//!
//! M4 expands this beyond the M0 header helper into a real record
//! writer that mirrors a subset of upstream freebayes' VCF schema — the
//! subset `bcftools norm` / `bcftools isec` actually care about:
//!
//! - `##fileformat`, `##source`, `##reference`, `##contig`
//! - `##INFO` declarations for `DP`, `AC`, `AN`, `AF`, `AO`, `RO`,
//!   `TYPE`, `NS`
//! - `##FORMAT` declarations for `GT`, `DP`, `AD`, `RO`, `AO`, `GQ`
//! - One `##FILTER=<ID=PASS,...>` line
//! - Records with proper anchor-base synthesis for INS/DEL and VCF-
//!   style `0/1` / `1/1` GT strings (not the internal `REF/SNP` tag).
//!
//! The upstream schema has ~30 INFO fields (many are diagnostic strand
//! / placement statistics derived from the per-allele observations).
//! M4 Phase A ports only the `bcftools isec`-visible fields; the
//! diagnostic fields land in Phase B as needed. Records are emitted as
//! plain text — `bcftools view -` reads them losslessly.

use std::fmt::Write as _;

/// A minimal contig entry for VCF `##contig` header lines.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Contig {
    pub name: String,
    pub length: u64,
}

/// Build a freebayes-compatible VCF 4.2 header. Includes every INFO /
/// FORMAT / FILTER declaration used by [`write_record`], plus the
/// `##reference` and `##contig` lines.
pub fn build_header(reference: &str, contigs: &[Contig], sample: &str) -> String {
    let mut out = String::new();
    out.push_str("##fileformat=VCFv4.2\n");
    let _ = writeln!(out, "##source=freebayes-gxy-{}", env!("CARGO_PKG_VERSION"));
    let _ = writeln!(out, "##reference={reference}");
    for c in contigs {
        let _ = writeln!(out, "##contig=<ID={},length={}>", c.name, c.length);
    }
    out.push_str("##FILTER=<ID=PASS,Description=\"All filters passed\">\n");
    // INFO declarations — M4 Phase A subset.
    out.push_str(
        "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with data\">\n",
    );
    out.push_str(
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth at the locus\">\n",
    );
    out.push_str("##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Total number of alternate alleles in called genotypes\">\n");
    out.push_str("##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">\n");
    out.push_str("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Estimated allele frequency in the range (0,1]\">\n");
    out.push_str("##INFO=<ID=RO,Number=1,Type=Integer,Description=\"Count of full observations of the reference haplotype.\">\n");
    out.push_str("##INFO=<ID=AO,Number=A,Type=Integer,Description=\"Count of full observations of this alternate haplotype.\">\n");
    out.push_str("##INFO=<ID=TYPE,Number=A,Type=String,Description=\"The type of allele, either snp, mnp, ins, del, or complex.\">\n");
    out.push_str("##INFO=<ID=CIGAR,Number=A,Type=String,Description=\"The extended CIGAR representation of each alternate allele relative to REF.\">\n");
    out.push_str("##INFO=<ID=QR,Number=1,Type=Integer,Description=\"Reference allele quality sum in phred\">\n");
    out.push_str("##INFO=<ID=QA,Number=A,Type=Integer,Description=\"Alternate allele quality sum in phred\">\n");
    // FORMAT declarations.
    out.push_str("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    out.push_str("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">\n");
    out.push_str("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele depth — one value per allele (REF first)\">\n");
    out.push_str("##FORMAT=<ID=RO,Number=1,Type=Integer,Description=\"Reference-allele observation count\">\n");
    out.push_str("##FORMAT=<ID=QR,Number=1,Type=Integer,Description=\"Reference-allele quality sum in Phred\">\n");
    out.push_str("##FORMAT=<ID=AO,Number=A,Type=Integer,Description=\"Alternate-allele observation count\">\n");
    out.push_str("##FORMAT=<ID=QA,Number=A,Type=Integer,Description=\"Alternate-allele quality sum in Phred\">\n");
    out.push_str(
        "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality (Phred)\">\n",
    );
    out.push_str("##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Genotype log-likelihoods, log10-scaled, in VCF-spec F(j/k) order\">\n");
    out.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t");
    out.push_str(sample);
    out.push('\n');
    out
}

/// Classification of a VCF record's primary allele event. Mirrors
/// upstream's `TYPE` INFO string values.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RecordKind {
    Snp,
    Mnp,
    Ins,
    Del,
    Complex,
}

impl RecordKind {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Snp => "snp",
            Self::Mnp => "mnp",
            Self::Ins => "ins",
            Self::Del => "del",
            Self::Complex => "complex",
        }
    }
}

/// One VCF record ready for emission. Sequences are raw ATGC bytes; no
/// anchor-base synthesis happens here — call [`synthesize_anchored`]
/// on your internal (kind, position, ref_seq, alt_seq) tuple first.
///
/// [`Record::default`] gives a skeleton callers can pattern-update:
/// empty chrom, pos 0, no alts, none of the optional fields set.
#[derive(Debug, Clone, Default)]
pub struct Record {
    pub chrom: String,
    /// 1-based position of the leftmost base in `ref_seq`.
    pub pos: i64,
    /// Reference allele bases (already anchored where applicable).
    pub ref_seq: Vec<u8>,
    /// Alternate alleles, in the order they appear in GT indices.
    pub alts: Vec<Vec<u8>>,
    /// Site-level quality (Phred). `None` emits `.`.
    pub qual: Option<f64>,
    /// Total depth at the site (INFO.DP and FORMAT.DP).
    pub depth: u32,
    /// Reference observation count (INFO.RO / FORMAT.RO).
    pub ref_obs: u32,
    /// Alternate observation counts, one per alt (INFO.AO / FORMAT.AO).
    pub alt_obs: Vec<u32>,
    /// GT index vector: `[0, 1]` for het REF/ALT1, `[1, 1]` for hom-alt1, etc.
    pub gt_indices: Vec<u8>,
    /// Phred-scaled genotype quality.
    pub gq: f64,
    /// Per-allele type tag (one per alt, in the same order as `alts`).
    pub alt_kinds: Vec<RecordKind>,
    /// Per-allele CIGAR string relative to REF (one per alt, in the
    /// same order as `alts`). Empty → `CIGAR` INFO omitted.
    pub cigars: Vec<String>,
    /// Reference-allele base-quality sum (Phred). Emitted as `QR` in
    /// INFO and FORMAT. `None` → field omitted.
    pub qual_ref: Option<u32>,
    /// Alternate-allele base-quality sums, one per alt (Phred).
    /// Emitted as `QA`. Empty → `QA` INFO omitted.
    pub qual_alt: Vec<u32>,
    /// Genotype log-likelihoods in VCF spec order (log10). For biallelic
    /// diploid, this is `[P(0/0), P(0/1), P(1/1)]`; for triallelic
    /// diploid it is `[P(0/0), P(0/1), P(1/1), P(0/2), P(1/2), P(2/2)]`
    /// per VCF 4.2 `F(j/k) = k*(k+1)/2 + j`. Emitted as `GL` in FORMAT.
    /// `None` → field omitted.
    pub genotype_log10_likelihoods: Option<Vec<f64>>,
}

/// VCF 4.2 genotype index for a sorted diploid (j, k) with j ≤ k.
///
/// `F(j/k) = k*(k+1)/2 + j`, with positions matching the VCF spec
/// example: `0/0→0, 0/1→1, 1/1→2, 0/2→3, 1/2→4, 2/2→5`.
///
/// Ploidy > 2 requires the general VCF combinatorial formula which is
/// not yet ported — we return `None` there so callers can decide
/// whether to emit `GL` at all for higher-ploidy genotypes.
pub fn vcf_gl_index(gt_indices: &[u8]) -> Option<usize> {
    if gt_indices.len() != 2 {
        return None;
    }
    let mut xs = [gt_indices[0] as usize, gt_indices[1] as usize];
    xs.sort_unstable();
    let (j, k) = (xs[0], xs[1]);
    Some(k * (k + 1) / 2 + j)
}

/// Build a per-alt CIGAR string against the record's REF sequence.
///
/// Mirrors upstream's `CIGAR` INFO emission which packs each alt's
/// edit path relative to REF as a minimal-character string:
/// - `1X` for a single-base SNP
/// - `<N>X` for an MNP of length N
/// - `1M<N>I` for a plain insertion of N bases anchored on the ref base
/// - `1M<N>D` for a plain deletion of N bases anchored on the ref base
/// - `<N>M<I>I<K>X` etc. for complex events — we approximate with a
///   single `<MAX(len_ref,len_alt)>M` for complex kind for now; the
///   full alignment port lands with [`RecordKind::Complex`] refinement.
pub fn alt_cigar(kind: RecordKind, ref_seq: &[u8], alt_seq: &[u8]) -> String {
    match kind {
        RecordKind::Snp => "1X".to_string(),
        RecordKind::Mnp => format!("{}X", ref_seq.len()),
        RecordKind::Ins => {
            // Anchored: REF = anchor, ALT = anchor + N inserted.
            let n_ins = alt_seq.len().saturating_sub(ref_seq.len());
            format!("1M{n_ins}I")
        }
        RecordKind::Del => {
            // Anchored: REF = anchor + N deleted, ALT = anchor.
            let n_del = ref_seq.len().saturating_sub(alt_seq.len());
            format!("1M{n_del}D")
        }
        RecordKind::Complex => {
            // Placeholder: use the longer of the two as the M count
            // with any length delta declared as an indel tail. This
            // is not a minimal-edit CIGAR and will differ from
            // upstream on complex events — TODO for M4 Phase C.
            let m = ref_seq.len().min(alt_seq.len());
            let diff = ref_seq.len() as i64 - alt_seq.len() as i64;
            match diff.cmp(&0) {
                std::cmp::Ordering::Equal => format!("{m}X"),
                std::cmp::Ordering::Greater => format!("{m}M{}D", diff.unsigned_abs()),
                std::cmp::Ordering::Less => format!("{m}M{}I", diff.unsigned_abs()),
            }
        }
    }
}

/// Convert the freebayes-gxy internal `(kind, position, ref_seq, alt_seq)`
/// tuple into a VCF-anchored `(pos, ref_seq, alt_seq, kind)` tuple
/// suitable for populating a [`Record`]. Upstream freebayes:
///
/// - **SNP**: emits as-is — `pos`, single REF base, single ALT base.
/// - **MNP**: emits as-is — N-base REF, N-base ALT of equal length.
/// - **INS**: prepends the anchor reference base at `pos - 1` to both
///   REF and ALT, so `REF=anchor, ALT=anchor+inserted`, and reports the
///   anchor position as `pos - 1` (1-based VCF = 0-based internal).
/// - **DEL**: prepends the anchor base at `pos - 1` to ALT, keeps REF
///   as `anchor + deleted_bases`, reports anchor position.
/// - **Complex**: ships the internal REF/ALT as-is but prepends the
///   anchor when either side would be empty.
///
/// `anchor_base` is the ref byte at `position_zero_based - 1`. The
/// caller is responsible for fetching it from the reference fasta
/// (INS / DEL require it; SNP / MNP ignore it).
///
/// Returns the VCF-anchored 1-based position, REF bytes, ALT bytes.
pub fn synthesize_anchored(
    kind: RecordKind,
    position_zero_based: i64,
    ref_seq: &[u8],
    alt_seq: &[u8],
    anchor_base: Option<u8>,
) -> (i64, Vec<u8>, Vec<u8>) {
    match kind {
        RecordKind::Snp | RecordKind::Mnp => {
            (position_zero_based + 1, ref_seq.to_vec(), alt_seq.to_vec())
        }
        RecordKind::Ins => {
            // Internal: ref_seq empty, alt_seq = inserted bases, position
            // = ref pos of the base immediately after the insertion.
            // VCF: anchor at position - 1 = ref pos of base PRECEDING the
            // insertion; REF = anchor; ALT = anchor + inserted.
            let anchor = anchor_base.expect("INS requires an anchor base");
            let anchor_pos = position_zero_based; // internal pos is anchor+1 in VCF
            let mut alt = Vec::with_capacity(1 + alt_seq.len());
            alt.push(anchor);
            alt.extend_from_slice(alt_seq);
            (anchor_pos, vec![anchor], alt)
        }
        RecordKind::Del => {
            // Internal: ref_seq = deleted bases, alt_seq empty,
            // position = ref pos of first deleted base.
            // VCF: anchor at position - 1; REF = anchor + deleted;
            // ALT = anchor.
            let anchor = anchor_base.expect("DEL requires an anchor base");
            let mut refs = Vec::with_capacity(1 + ref_seq.len());
            refs.push(anchor);
            refs.extend_from_slice(ref_seq);
            (position_zero_based, refs, vec![anchor])
        }
        RecordKind::Complex => {
            // Complex alleles are already (ref_seq, alt_seq) of some
            // length from the M2 clumping pass. If either side is
            // empty (pure INS or DEL composite), prepend anchor; else
            // emit as-is at 1-based pos.
            if ref_seq.is_empty() || alt_seq.is_empty() {
                let anchor = anchor_base.expect("empty-side complex requires an anchor base");
                let mut refs = Vec::with_capacity(1 + ref_seq.len());
                refs.push(anchor);
                refs.extend_from_slice(ref_seq);
                let mut alt = Vec::with_capacity(1 + alt_seq.len());
                alt.push(anchor);
                alt.extend_from_slice(alt_seq);
                (position_zero_based, refs, alt)
            } else {
                (position_zero_based + 1, ref_seq.to_vec(), alt_seq.to_vec())
            }
        }
    }
}

/// Serialise a record to a VCF line (including trailing newline). Uses
/// the header's INFO / FORMAT declarations.
pub fn write_record(rec: &Record) -> String {
    let mut out = String::new();
    // CHROM, POS, ID, REF, ALT.
    out.push_str(&rec.chrom);
    let _ = write!(out, "\t{}\t.\t", rec.pos);
    out.push_str(std::str::from_utf8(&rec.ref_seq).unwrap_or("N"));
    out.push('\t');
    let alt_str = rec
        .alts
        .iter()
        .map(|a| std::str::from_utf8(a).unwrap_or("N").to_string())
        .collect::<Vec<_>>()
        .join(",");
    out.push_str(&alt_str);
    // QUAL, FILTER.
    out.push('\t');
    match rec.qual {
        Some(q) => {
            let _ = write!(out, "{q:.3}");
        }
        None => out.push('.'),
    }
    out.push_str("\tPASS\t");
    // INFO.
    let an = rec.gt_indices.len() as u32;
    let ac: Vec<u32> = (0..rec.alts.len())
        .map(|alt_idx| {
            rec.gt_indices
                .iter()
                .filter(|&&g| g as usize == alt_idx + 1)
                .count() as u32
        })
        .collect();
    let af: Vec<f64> = ac
        .iter()
        .map(|&c| if an == 0 { 0.0 } else { c as f64 / an as f64 })
        .collect();
    let _ = write!(out, "NS=1;DP={};", rec.depth);
    let _ = write!(
        out,
        "AC={};AN={};AF={};",
        comma_join_u32(&ac),
        an,
        comma_join_f64(&af, 6),
    );
    let _ = write!(
        out,
        "RO={};AO={};",
        rec.ref_obs,
        comma_join_u32(&rec.alt_obs)
    );
    let type_str = rec
        .alt_kinds
        .iter()
        .map(|k| k.as_str())
        .collect::<Vec<_>>()
        .join(",");
    let _ = write!(out, "TYPE={type_str}");
    if !rec.cigars.is_empty() {
        let _ = write!(out, ";CIGAR={}", rec.cigars.join(","));
    }
    if let Some(qr) = rec.qual_ref {
        let _ = write!(out, ";QR={qr}");
    }
    if !rec.qual_alt.is_empty() {
        let _ = write!(out, ";QA={}", comma_join_u32(&rec.qual_alt));
    }
    // FORMAT: decide columns based on which optional fields are present.
    let mut fmt_fields = vec!["GT", "DP", "AD", "RO"];
    if rec.qual_ref.is_some() {
        fmt_fields.push("QR");
    }
    fmt_fields.push("AO");
    if !rec.qual_alt.is_empty() {
        fmt_fields.push("QA");
    }
    fmt_fields.push("GQ");
    if rec.genotype_log10_likelihoods.is_some() {
        fmt_fields.push("GL");
    }
    let _ = write!(out, "\t{}\t", fmt_fields.join(":"));
    // Sample fields in the same order.
    let gt_str = rec
        .gt_indices
        .iter()
        .map(|g| g.to_string())
        .collect::<Vec<_>>()
        .join("/");
    let mut ad = Vec::with_capacity(1 + rec.alt_obs.len());
    ad.push(rec.ref_obs);
    ad.extend_from_slice(&rec.alt_obs);
    let gq_int = rec.gq.round() as i64;

    let mut sample_values: Vec<String> = Vec::with_capacity(fmt_fields.len());
    sample_values.push(gt_str);
    sample_values.push(rec.depth.to_string());
    sample_values.push(comma_join_u32(&ad));
    sample_values.push(rec.ref_obs.to_string());
    if let Some(qr) = rec.qual_ref {
        sample_values.push(qr.to_string());
    }
    sample_values.push(comma_join_u32(&rec.alt_obs));
    if !rec.qual_alt.is_empty() {
        sample_values.push(comma_join_u32(&rec.qual_alt));
    }
    sample_values.push(gq_int.to_string());
    if let Some(gl) = &rec.genotype_log10_likelihoods {
        sample_values.push(comma_join_f64(gl, 4));
    }
    out.push_str(&sample_values.join(":"));
    out.push('\n');
    out
}

fn comma_join_u32(xs: &[u32]) -> String {
    xs.iter()
        .map(|x| x.to_string())
        .collect::<Vec<_>>()
        .join(",")
}

fn comma_join_f64(xs: &[f64], precision: usize) -> String {
    xs.iter()
        .map(|x| format!("{x:.*}", precision))
        .collect::<Vec<_>>()
        .join(",")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn header_contains_required_lines() {
        let contigs = vec![Contig {
            name: "chr1".to_string(),
            length: 248_956_422,
        }];
        let h = build_header("/tmp/ref.fa", &contigs, "SAMPLE1");
        assert!(h.starts_with("##fileformat=VCFv4.2\n"));
        assert!(h.contains("##reference=/tmp/ref.fa\n"));
        assert!(h.contains("##contig=<ID=chr1,length=248956422>\n"));
        assert!(h.contains("##INFO=<ID=DP,"));
        assert!(h.contains("##FORMAT=<ID=GT,"));
        assert!(h.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n"));
    }

    #[test]
    fn anchor_synthesis_snp_is_identity() {
        let (pos, r, a) = synthesize_anchored(RecordKind::Snp, 100, b"A", b"G", Some(b'C'));
        // 0-based 100 -> 1-based 101, no anchor prepended.
        assert_eq!(pos, 101);
        assert_eq!(r, b"A");
        assert_eq!(a, b"G");
    }

    #[test]
    fn anchor_synthesis_ins_prepends_anchor() {
        // Internal: INS at 0-based pos 105, alt "GG". Anchor base at
        // 0-based 104 is 'A'. VCF: pos = 105 (1-based), REF=A, ALT=AGG.
        let (pos, r, a) = synthesize_anchored(RecordKind::Ins, 105, b"", b"GG", Some(b'A'));
        assert_eq!(pos, 105);
        assert_eq!(r, b"A");
        assert_eq!(a, b"AGG");
    }

    #[test]
    fn anchor_synthesis_del_prepends_anchor() {
        // Internal: DEL at 0-based pos 200, ref "CCC". Anchor at 199 = 'T'.
        // VCF: pos = 200 (1-based), REF=TCCC, ALT=T.
        let (pos, r, a) = synthesize_anchored(RecordKind::Del, 200, b"CCC", b"", Some(b'T'));
        assert_eq!(pos, 200);
        assert_eq!(r, b"TCCC");
        assert_eq!(a, b"T");
    }

    #[test]
    fn write_record_hom_ref_snp() {
        let r = Record {
            chrom: "chr1".into(),
            pos: 101,
            ref_seq: b"A".to_vec(),
            alts: vec![b"G".to_vec()],
            qual: Some(321.056),
            depth: 11,
            ref_obs: 0,
            alt_obs: vec![11],
            gt_indices: vec![1, 1],
            gq: 39.0,
            alt_kinds: vec![RecordKind::Snp],
            ..Record::default()
        };
        let line = write_record(&r);
        assert!(
            line.starts_with("chr1\t101\t.\tA\tG\t321.056\tPASS\t"),
            "got: {line}"
        );
        assert!(line.contains("NS=1;DP=11;"));
        assert!(line.contains("AC=2;AN=2;AF=1.000000;"));
        assert!(line.contains("RO=0;AO=11;"));
        assert!(line.contains("TYPE=snp"));
        assert!(line.contains("\tGT:DP:AD:RO:AO:GQ\t1/1:11:0,11:0:11:39\n"));
    }

    #[test]
    fn write_record_het_snp_with_gl() {
        let r = Record {
            chrom: "chr1".into(),
            pos: 201,
            ref_seq: b"A".to_vec(),
            alts: vec![b"G".to_vec()],
            qual: Some(233.58),
            depth: 20,
            ref_obs: 10,
            alt_obs: vec![10],
            gt_indices: vec![0, 1],
            gq: 60.0,
            alt_kinds: vec![RecordKind::Snp],
            cigars: vec!["1X".into()],
            qual_ref: Some(300),
            qual_alt: vec![300],
            genotype_log10_likelihoods: Some(vec![-13.1715, 0.0, -13.1715]),
        };
        let line = write_record(&r);
        assert!(line.contains("AC=1;AN=2;AF=0.500000;"));
        assert!(line.contains("RO=10;AO=10;"));
        assert!(line.contains("CIGAR=1X"));
        assert!(line.contains("QR=300"));
        assert!(line.contains("QA=300"));
        assert!(line.contains("\tGT:DP:AD:RO:QR:AO:QA:GQ:GL\t"));
        assert!(line.contains("0/1:20:10,10:10:300:10:300:60:-13.1715,0.0000,-13.1715"));
    }

    #[test]
    fn write_record_hom_ins_anchor_synthesised() {
        let (pos, refs, alts) = synthesize_anchored(RecordKind::Ins, 105, b"", b"TT", Some(b'A'));
        let r = Record {
            chrom: "chr1".into(),
            pos,
            ref_seq: refs,
            alts: vec![alts.clone()],
            qual: Some(100.0),
            depth: 20,
            ref_obs: 0,
            alt_obs: vec![20],
            gt_indices: vec![1, 1],
            gq: 60.0,
            alt_kinds: vec![RecordKind::Ins],
            cigars: vec![alt_cigar(RecordKind::Ins, b"A", &alts)],
            ..Record::default()
        };
        let line = write_record(&r);
        assert!(line.starts_with("chr1\t105\t.\tA\tATT\t100.000\tPASS\t"));
        assert!(line.contains("TYPE=ins"));
        assert!(line.contains("CIGAR=1M2I"));
        assert!(line.contains("\t1/1:20:0,20:0:20:60\n"));
    }

    #[test]
    fn write_record_het_del() {
        let (pos, refs, alts) = synthesize_anchored(RecordKind::Del, 200, b"CC", b"", Some(b'T'));
        let r = Record {
            chrom: "chr1".into(),
            pos,
            ref_seq: refs.clone(),
            alts: vec![alts],
            qual: Some(50.0),
            depth: 20,
            ref_obs: 10,
            alt_obs: vec![10],
            gt_indices: vec![0, 1],
            gq: 30.0,
            alt_kinds: vec![RecordKind::Del],
            cigars: vec![alt_cigar(RecordKind::Del, &refs, b"T")],
            ..Record::default()
        };
        let line = write_record(&r);
        assert!(line.starts_with("chr1\t200\t.\tTCC\tT\t50.000\tPASS\t"));
        assert!(line.contains("TYPE=del"));
        assert!(line.contains("CIGAR=1M2D"));
        assert!(line.contains("\t0/1:20:10,10:10:10:30\n"));
    }

    #[test]
    fn vcf_gl_index_biallelic_diploid() {
        // 0/0 → 0, 0/1 → 1, 1/1 → 2.
        assert_eq!(vcf_gl_index(&[0, 0]), Some(0));
        assert_eq!(vcf_gl_index(&[0, 1]), Some(1));
        assert_eq!(vcf_gl_index(&[1, 0]), Some(1)); // unsorted ok
        assert_eq!(vcf_gl_index(&[1, 1]), Some(2));
    }

    #[test]
    fn vcf_gl_index_triallelic_diploid_matches_spec_order() {
        // Spec: 0/0=0, 0/1=1, 1/1=2, 0/2=3, 1/2=4, 2/2=5.
        assert_eq!(vcf_gl_index(&[0, 0]), Some(0));
        assert_eq!(vcf_gl_index(&[0, 1]), Some(1));
        assert_eq!(vcf_gl_index(&[1, 1]), Some(2));
        assert_eq!(vcf_gl_index(&[0, 2]), Some(3));
        assert_eq!(vcf_gl_index(&[1, 2]), Some(4));
        assert_eq!(vcf_gl_index(&[2, 2]), Some(5));
    }

    #[test]
    fn vcf_gl_index_non_diploid_returns_none() {
        assert_eq!(vcf_gl_index(&[0]), None);
        assert_eq!(vcf_gl_index(&[0, 0, 0]), None);
    }

    #[test]
    fn alt_cigar_kinds() {
        assert_eq!(alt_cigar(RecordKind::Snp, b"A", b"G"), "1X");
        assert_eq!(alt_cigar(RecordKind::Mnp, b"AC", b"GT"), "2X");
        assert_eq!(alt_cigar(RecordKind::Ins, b"A", b"ATT"), "1M2I");
        assert_eq!(alt_cigar(RecordKind::Del, b"ACC", b"A"), "1M2D");
    }

    #[test]
    fn header_preserves_contig_order() {
        let contigs = vec![
            Contig {
                name: "chrZ".to_string(),
                length: 100,
            },
            Contig {
                name: "chrA".to_string(),
                length: 200,
            },
        ];
        let h = build_header("ref.fa", &contigs, "S");
        let z_idx = h.find("chrZ").unwrap();
        let a_idx = h.find("chrA").unwrap();
        assert!(z_idx < a_idx, "contig order should be preserved");
    }
}
