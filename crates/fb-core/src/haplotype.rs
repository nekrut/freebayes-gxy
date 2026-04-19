//! Per-read haplotype clumping.
//!
//! Port of upstream freebayes v1.3.10's
//! `RegisteredAlignment::clumpAlleles` (`src/AlleleParser.cpp:1013-1062`)
//! and `Allele::mergeAllele` (`src/Allele.cpp:1454-1470`). The pass takes
//! a read's ordered sequence of [`AlleleObservation`]s as produced by the
//! M1 CIGAR walker and collapses adjacent non-reference events (with
//! short intervening reference runs, up to `max_complex_gap` bases) into
//! single [`AlleleKind::Complex`] observations.
//!
//! ## Why this lives at the per-read level
//! Upstream does the clumping on each read's `RegisteredAlignment` before
//! the Bayesian model ever sees it, and the genotype likelihood code
//! consumes composite observations per read. Matching that is necessary
//! for byte-level VCF parity — a pair of adjacent SNPs reported as one
//! MNP-style complex allele supports a different haplotype hypothesis
//! than two independent SNPs.
//!
//! ## M2 scope (this module)
//! The three-way-window algorithm + sweep-and-merge from upstream are
//! ported verbatim. The exact upstream CIGAR synthesis on merged alleles
//! (see `mergeCigar` in `Allele.cpp:1468`) and left-alignment of indel
//! composites (`LeftAlign.cpp`) are deferred — those are needed only for
//! M4's VCF emission and can be layered on without reshaping the data
//! types below. TODOs point at the upstream line ranges.
//!
//! ## Non-goals
//! Cross-read haplotype assembly (`AlleleParser::buildHaplotypeAlleles`,
//! `src/AlleleParser.cpp:3214..`) stays as M3+ work — it is the allele
//! *set* at a site, not the per-read observation shape.

use crate::allele::{AlleleKind, AlleleObservation};

/// Collapse adjacent non-reference observations in `observations` into
/// composite `Complex` observations.
///
/// `observations` must be the ordered output of the M1 CIGAR walker for a
/// **single read**, sorted by reference position. `max_complex_gap` is
/// the maximum reference run length that is allowed to bridge two
/// flanking non-reference events into a single clump (freebayes CLI
/// default: `3`, set via `--haplotype-length`). A negative value
/// disables clumping entirely, mirroring upstream's
/// `if (maxComplexGap >= 0)` guard at `AlleleParser.cpp:1018`.
///
/// Returns a newly-allocated vector; the input is not mutated.
pub fn clump_observations(
    observations: &[AlleleObservation],
    max_complex_gap: i64,
) -> Vec<AlleleObservation> {
    if max_complex_gap < 0 || observations.len() < 3 {
        // Upstream's algorithm only fires when at least three alleles are
        // in the list (the sliding 3-way window needs a left and a right
        // flank) — fewer than that means nothing can be merged.
        return observations.to_vec();
    }

    let n = observations.len();
    let mut to_merge = vec![false; n];

    // Three-way sliding window, mirroring AlleleParser.cpp:1019-1040. The
    // inner index runs over the **middle** of each triple, so the valid
    // range is 1..n-1. Upstream uses a signed loop bound `i < n-1`; we
    // skip the edges the same way.
    for i in 1..n - 1 {
        let last = &observations[i - 1].allele;
        let curr = &observations[i].allele;
        let next = &observations[i + 1].allele;

        // Null observations (soft clips, N bases) never participate in
        // clumping (AlleleParser.cpp:1023).
        if last.kind == AlleleKind::Null
            || curr.kind == AlleleKind::Null
            || next.kind == AlleleKind::Null
        {
            continue;
        }

        let last_is_ref = last.kind == AlleleKind::Reference;
        let curr_is_ref = curr.kind == AlleleKind::Reference;
        let next_is_ref = next.kind == AlleleKind::Reference;
        // Upstream's `referenceLength` is the reference-side length of
        // the allele. For our REF observations that equals `length`.
        let curr_ref_len = if curr_is_ref { curr.length as i64 } else { 0 };

        if !last_is_ref && !next_is_ref {
            // Two non-reference flanks with the middle either
            //  - a non-reference event (case !curr_is_ref), or
            //  - a short reference gap (curr_is_ref && ref_len <= gap).
            // Matches AlleleParser.cpp:1024-1028.
            // Equivalent to upstream's
            //   (curr_is_ref && curr_ref_len <= max_complex_gap) || !curr_is_ref
            // collapsed via De Morgan: when `curr_is_ref` is false the
            // length check is immaterial; when it's true, require the gap
            // to fit.
            if !curr_is_ref || curr_ref_len <= max_complex_gap {
                to_merge[i - 1] = true;
                to_merge[i] = true;
                to_merge[i + 1] = true;
            }
        } else if !last_is_ref && !curr_is_ref {
            // Left flank + middle are both non-reference (next can be
            // anything — reference or otherwise). Upstream line 1033-1036.
            to_merge[i - 1] = true;
            to_merge[i] = true;
        } else if !next_is_ref && !curr_is_ref {
            // Mirror of the above: middle + right flank non-reference.
            // Upstream line 1036-1039.
            to_merge[i] = true;
            to_merge[i + 1] = true;
        }
    }

    // Sweep: contiguous runs of `true` in `to_merge` collapse into one
    // composite; everything else is emitted verbatim. Mirrors
    // AlleleParser.cpp:1044-1058, simplified because our while-loop
    // preserves the breaker automatically on the next outer iteration.
    let mut out = Vec::with_capacity(n);
    let mut i = 0;
    while i < n {
        if to_merge[i] {
            let mut merged = observations[i].clone();
            let mut j = i + 1;
            while j < n && to_merge[j] {
                merged.allele.merge_with(&observations[j].allele);
                // Sum base-quality contributions across the clump; the
                // Bayesian model (M3) consumes a per-observation BQ, and
                // upstream uses `averageQuality(baseQualities)` after
                // concatenation — we preserve the sum here and leave the
                // final scaling to M3.
                merged.base_quality_sum = merged
                    .base_quality_sum
                    .saturating_add(observations[j].base_quality_sum);
                j += 1;
            }
            out.push(merged);
            i = j;
        } else {
            out.push(observations[i].clone());
            i += 1;
        }
    }

    out
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::allele::{Allele, Strand};

    fn obs(allele: Allele, bq: u32) -> AlleleObservation {
        AlleleObservation {
            allele,
            read_name: "r1".into(),
            mapq: 60,
            base_quality_sum: bq,
            strand: Strand::Forward,
            read_position: 0,
            is_proper_pair: true,
            read_ref_start: 0,
            per_base_quals: Vec::new(),
        }
    }

    #[test]
    fn too_short_to_clump() {
        // With fewer than three observations, clump_observations must be
        // a no-op — upstream's 3-way window never fires.
        let input = vec![
            obs(Allele::snp(100, b'A', b'G'), 30),
            obs(Allele::snp(101, b'C', b'T'), 30),
        ];
        let out = clump_observations(&input, 3);
        assert_eq!(out.len(), 2);
        assert!(out.iter().all(|o| o.allele.kind == AlleleKind::Snp));
    }

    #[test]
    fn negative_gap_disables_clumping() {
        let input = vec![
            obs(Allele::snp(100, b'A', b'G'), 30),
            obs(Allele::snp(101, b'C', b'T'), 30),
            obs(Allele::snp(102, b'G', b'A'), 30),
        ];
        let out = clump_observations(&input, -1);
        assert_eq!(out.len(), 3);
        assert!(out.iter().all(|o| o.allele.kind == AlleleKind::Snp));
    }

    #[test]
    fn isolated_snp_stays_isolated() {
        // REF-run flanks on both sides → nothing to merge.
        let input = vec![
            obs(Allele::reference(95, b"ACGTA".to_vec()), 60),
            obs(Allele::snp(100, b'A', b'G'), 30),
            obs(Allele::reference(101, b"ACGTA".to_vec()), 60),
        ];
        let out = clump_observations(&input, 3);
        assert_eq!(out.len(), 3);
        assert_eq!(out[0].allele.kind, AlleleKind::Reference);
        assert_eq!(out[1].allele.kind, AlleleKind::Snp);
        assert_eq!(out[2].allele.kind, AlleleKind::Reference);
    }

    #[test]
    fn snp_adjacent_to_insertion_becomes_complex() {
        // REF, SNP, INS, REF — per upstream, SNP+INS clumps into COMPLEX.
        let input = vec![
            obs(Allele::reference(95, b"ACGTA".to_vec()), 60),
            obs(Allele::snp(100, b'A', b'G'), 30),
            obs(Allele::insertion(101, b"TT".to_vec()), 50),
            obs(Allele::reference(101, b"ACGTA".to_vec()), 60),
        ];
        let out = clump_observations(&input, 3);
        // REF, COMPLEX(SNP+INS), REF — 3 observations total.
        assert_eq!(out.len(), 3);
        assert_eq!(out[0].allele.kind, AlleleKind::Reference);
        assert_eq!(out[1].allele.kind, AlleleKind::Complex);
        assert_eq!(out[2].allele.kind, AlleleKind::Reference);
        // Composite carries full ref + alt spans (upstream Allele.cpp:1457-1458).
        assert_eq!(out[1].allele.ref_seq, b"A"); // SNP ref A + INS empty ref
        assert_eq!(out[1].allele.alt_seq, b"GTT"); // SNP alt G + INS alt TT
        assert_eq!(out[1].allele.length, 3); // 1 (SNP) + 2 (INS)
                                             // BQ sums are additive across the clump.
        assert_eq!(out[1].base_quality_sum, 30 + 50);
    }

    #[test]
    fn two_snps_across_short_ref_gap_merge() {
        // SNP, REF(len=2), SNP at gap=3 → clumps to one COMPLEX.
        let input = vec![
            obs(Allele::snp(100, b'A', b'G'), 30),
            obs(Allele::reference(101, b"CC".to_vec()), 60),
            obs(Allele::snp(103, b'T', b'A'), 30),
        ];
        let out = clump_observations(&input, 3);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].allele.kind, AlleleKind::Complex);
        assert_eq!(out[0].allele.position, 100);
        assert_eq!(out[0].allele.ref_seq, b"ACCT"); // A (SNP) + CC (REF) + T (SNP)
        assert_eq!(out[0].allele.alt_seq, b"GCCA"); // G + CC + A
        assert_eq!(out[0].base_quality_sum, 30 + 60 + 30);
    }

    #[test]
    fn two_snps_across_long_ref_gap_do_not_merge() {
        // Same shape but the REF gap is longer than max_complex_gap.
        let input = vec![
            obs(Allele::snp(100, b'A', b'G'), 30),
            obs(Allele::reference(101, b"CCCCC".to_vec()), 60),
            obs(Allele::snp(106, b'T', b'A'), 30),
        ];
        let out = clump_observations(&input, 3);
        assert_eq!(out.len(), 3);
        assert!(out.iter().all(|o| o.allele.kind != AlleleKind::Complex));
    }

    #[test]
    fn null_observation_breaks_clumping() {
        // NULL in the middle of an otherwise clumpable triple blocks
        // merging per upstream AlleleParser.cpp:1023.
        let input = vec![
            obs(Allele::snp(100, b'A', b'G'), 30),
            obs(Allele::null(101, 1), 30),
            obs(Allele::snp(102, b'T', b'A'), 30),
        ];
        let out = clump_observations(&input, 3);
        assert_eq!(out.len(), 3);
        assert!(out.iter().all(|o| o.allele.kind != AlleleKind::Complex));
    }

    #[test]
    fn three_adjacent_non_reference_all_merge() {
        // SNP, INS, SNP with REF flanks — upstream's `!last_ref && !next_ref
        // && !curr_is_ref` branch catches this.
        let input = vec![
            obs(Allele::reference(95, b"AAAAA".to_vec()), 60),
            obs(Allele::snp(100, b'A', b'G'), 30),
            obs(Allele::insertion(101, b"T".to_vec()), 30),
            obs(Allele::snp(101, b'C', b'A'), 30),
            obs(Allele::reference(102, b"AAAAA".to_vec()), 60),
        ];
        let out = clump_observations(&input, 3);
        // REF, COMPLEX(SNP+INS+SNP), REF.
        assert_eq!(out.len(), 3);
        assert_eq!(out[1].allele.kind, AlleleKind::Complex);
        assert_eq!(out[1].allele.ref_seq, b"AC"); // SNP A + INS empty + SNP C
        assert_eq!(out[1].allele.alt_seq, b"GTA"); // SNP G + INS T + SNP A
        assert_eq!(out[1].allele.length, 3); // 1 + 1 + 1
    }

    #[test]
    fn deletion_plus_snp_becomes_complex() {
        let input = vec![
            obs(Allele::reference(95, b"AAAAA".to_vec()), 60),
            obs(Allele::deletion(100, b"TT".to_vec()), 40),
            obs(Allele::snp(102, b'A', b'G'), 30),
            obs(Allele::reference(103, b"AAAAA".to_vec()), 60),
        ];
        let out = clump_observations(&input, 3);
        assert_eq!(out.len(), 3);
        assert_eq!(out[1].allele.kind, AlleleKind::Complex);
        // Composite ref: deletion's "TT" + SNP's "A" = "TTA".
        assert_eq!(out[1].allele.ref_seq, b"TTA");
        // Composite alt: deletion's empty + SNP's "G" = "G".
        assert_eq!(out[1].allele.alt_seq, b"G");
    }

    #[test]
    fn clumping_preserves_non_adjacent_runs() {
        // Two independent clumps separated by a long REF should both
        // survive as Complex but stay separate.
        let input = vec![
            obs(Allele::reference(90, b"A".to_vec()), 60),
            obs(Allele::snp(100, b'A', b'G'), 30),
            obs(Allele::insertion(101, b"T".to_vec()), 30),
            obs(Allele::reference(101, b"A".repeat(20)), 60), // 20 bp REF
            obs(Allele::snp(121, b'C', b'T'), 30),
            obs(Allele::insertion(122, b"G".to_vec()), 30),
            obs(Allele::reference(122, b"A".to_vec()), 60),
        ];
        let out = clump_observations(&input, 3);
        let complexes: Vec<_> = out
            .iter()
            .filter(|o| o.allele.kind == AlleleKind::Complex)
            .collect();
        assert_eq!(complexes.len(), 2);
        assert_eq!(complexes[0].allele.position, 100);
        assert_eq!(complexes[1].allele.position, 121);
    }

    #[test]
    fn gap_equal_to_threshold_merges() {
        // Upstream uses `<=` so an exact-threshold REF run should clump.
        let input = vec![
            obs(Allele::snp(100, b'A', b'G'), 30),
            obs(Allele::reference(101, b"CCC".to_vec()), 60), // len=3
            obs(Allele::snp(104, b'T', b'A'), 30),
        ];
        let out = clump_observations(&input, 3);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].allele.kind, AlleleKind::Complex);
        assert_eq!(out[0].allele.ref_seq, b"ACCCT");
        assert_eq!(out[0].allele.alt_seq, b"GCCCA");
    }

    #[test]
    fn gap_one_past_threshold_does_not_merge() {
        // Boundary check mirror of the above — len > gap blocks clumping.
        let input = vec![
            obs(Allele::snp(100, b'A', b'G'), 30),
            obs(Allele::reference(101, b"CCCC".to_vec()), 60), // len=4
            obs(Allele::snp(105, b'T', b'A'), 30),
        ];
        let out = clump_observations(&input, 3);
        assert_eq!(out.len(), 3);
        assert!(out.iter().all(|o| o.allele.kind != AlleleKind::Complex));
    }

    #[test]
    fn run_of_four_non_reference_events_all_merge() {
        // SNP SNP INS DEL in a row (with REF flanks) — the sweep should
        // produce a single COMPLEX composite spanning all four.
        let input = vec![
            obs(Allele::reference(90, b"A".to_vec()), 60),
            obs(Allele::snp(100, b'A', b'G'), 30),
            obs(Allele::snp(101, b'C', b'T'), 30),
            obs(Allele::insertion(102, b"A".to_vec()), 30),
            obs(Allele::deletion(102, b"G".to_vec()), 30),
            obs(Allele::reference(103, b"A".to_vec()), 60),
        ];
        let out = clump_observations(&input, 3);
        let complexes: Vec<_> = out
            .iter()
            .filter(|o| o.allele.kind == AlleleKind::Complex)
            .collect();
        assert_eq!(complexes.len(), 1);
        // Composite ref: SNP A + SNP C + INS "" + DEL G = "ACG"
        // Composite alt: SNP G + SNP T + INS A + DEL "" = "GTA"
        assert_eq!(complexes[0].allele.ref_seq, b"ACG");
        assert_eq!(complexes[0].allele.alt_seq, b"GTA");
        assert_eq!(complexes[0].allele.length, 4); // 1+1+1+1
    }

    #[test]
    fn insertion_then_deletion_preserves_concatenation_order() {
        // INS + DEL: ref comes from DEL (left-concat order matters), alt
        // from INS. Upstream concatenates in insertion order; verify we
        // do the same.
        let input = vec![
            obs(Allele::reference(90, b"A".to_vec()), 60),
            obs(Allele::insertion(100, b"TT".to_vec()), 40),
            obs(Allele::deletion(100, b"CC".to_vec()), 40),
            obs(Allele::reference(102, b"A".to_vec()), 60),
        ];
        let out = clump_observations(&input, 3);
        assert_eq!(out.len(), 3);
        assert_eq!(out[1].allele.kind, AlleleKind::Complex);
        // INS ref="" + DEL ref="CC" → "CC"
        assert_eq!(out[1].allele.ref_seq, b"CC");
        // INS alt="TT" + DEL alt="" → "TT"
        assert_eq!(out[1].allele.alt_seq, b"TT");
    }

    #[test]
    fn empty_input_is_passthrough() {
        let out = clump_observations(&[], 3);
        assert!(out.is_empty());
    }
}
