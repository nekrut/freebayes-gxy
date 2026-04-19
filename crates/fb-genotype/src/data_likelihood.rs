//! Per-sample data likelihood: `P(observations | genotype)`.
//!
//! Port of the `standardGLs` path of upstream freebayes v1.3.10's
//! `probObservedAllelesGivenGenotype` (`src/DataLikelihood.cpp:5-166`).
//! Given a candidate [`Genotype`] and the set of
//! [`AlleleObservation`]s at a site for one sample, returns the log of
//! the probability that those observations were produced under the
//! genotype.
//!
//! ## Model
//! The standard-GL likelihood for a sample at one site decomposes into:
//!
//! 1. **Error term** (`prod_q_out`). Every observation whose allele is
//!    not one of the genotype's elements contributes `ln(err_prob)`,
//!    where `err_prob = 10^(-Q/10)` for base quality `Q`. Upstream's
//!    [`phred2ln`](../../../../../tmp/freebayes-upstream/src/Utility.cpp)
//!    (line 35-37) does the same conversion.
//!
//! 2. **Read-dependence factor** (`RDF`, default 0.9). If more than one
//!    observation is "outside" the genotype, the error sum is scaled by
//!    `(1 + (n-1)*RDF) / n` to downweight successive reads (upstream
//!    `src/DataLikelihood.cpp:147-149`). A single outside observation
//!    is unaffected; `RDF = 0.0` disables the downweight.
//!
//! 3. **Multinomial sampling term**. Observations whose allele is inside
//!    the genotype are scored via
//!    [`crate::multinomial::multinomial_sampling_prob_ln`] against the
//!    genotype's allele probabilities (upstream line 155).
//!
//! ## Deliberate scope limits
//! The experimental non-standard-GL path (`DataLikelihood.cpp:44-140`)
//! handles partial observations, contamination estimation, and
//! reversible-partials bookkeeping. Those pieces land in a later M3
//! sub-phase and depend on `Contamination.cpp` + `Bias.cpp` ports that
//! this milestone does not include. `Bias` correction is also deferred
//! — the prior reviewer flagged that once `Bias` ships, it must not
//! push any genotype sampling probability to exactly 0, or the
//! `multinomial_sampling_prob_ln` `k == 0` short-circuit guarantee
//! weakens.

use fb_core::AlleleObservation;

use crate::genotype::Genotype;
use crate::multinomial::multinomial_sampling_prob_ln;
use crate::params::Parameters;

/// Upstream freebayes' default read-dependence factor
/// (`src/Parameters.cpp:475`). Exposed as a constant so callers can
/// pass it through without rediscovering the magic number.
pub const DEFAULT_READ_DEPENDENCE_FACTOR: f64 = 0.9;

/// Convert a Phred-scaled integer quality to `ln(err_prob)`.
///
/// `phred2ln(q) = ln(10^(-q/10)) = -q * ln(10) / 10`.
/// Mirrors upstream `phred2ln` at `src/Utility.cpp:35-37`.
#[inline]
pub fn phred_to_ln(q: u8) -> f64 {
    // M_LN10 * q * -0.1 — upstream uses long double, we use f64.
    std::f64::consts::LN_10 * (q as f64) * -0.1
}

/// Compute `log P(observations | genotype)` for a single sample, using
/// the standard-GL model. Mirrors the `parameters.standardGLs` branch
/// of upstream `probObservedAllelesGivenGenotype`
/// (`src/DataLikelihood.cpp:26-43, 147-157`).
///
/// `params.read_dependence_factor` drives the RDF scaling;
/// `params.use_mapping_quality` gates the `max(ln_bq, ln_mq)`
/// cap on outside-observation quality (upstream's
/// `useMappingQuality` flag at `DataLikelihood.cpp:31`).
///
/// Observations whose base quality is unknown should carry
/// `base_quality_sum = 0` and will be treated as Phred 0 (error prob
/// 1.0 — i.e. zero information); callers are responsible for filtering
/// such observations if that behaviour is undesired.
///
/// Special cases:
/// - Empty `observations`: returns `0.0` (log of probability 1 — no
///   evidence means no likelihood penalty).
/// - All observations outside the genotype: returns the error term
///   alone, with RDF applied (upstream line 151-152).
/// - All observations inside the genotype: returns the multinomial
///   term alone (upstream line 155 with `prod_q_out = 0`).
pub fn sample_log_likelihood(
    genotype: &Genotype,
    observations: &[AlleleObservation],
    params: &Parameters,
) -> f64 {
    if observations.is_empty() {
        return 0.0;
    }

    let mut in_counts = vec![0i64; genotype.elements.len()];
    let mut prod_q_out = 0.0f64;
    let mut count_out: i64 = 0;

    for obs in observations {
        // Match the observation against the genotype elements by
        // `(position, kind, ref_seq, alt_seq)` equality. Upstream
        // matches on `currentBase` (a stringified base identity); our
        // canonical Allele equality is strictly stronger but produces
        // the same partitioning for M1/M2-shaped observations.
        let mut matched = false;
        for (i, elem) in genotype.elements.iter().enumerate() {
            if elem.allele == obs.allele {
                in_counts[i] += 1;
                matched = true;
                break;
            }
        }
        if !matched {
            // Upstream uses `max(lnquality, lnmapQuality)` when
            // `useMappingQuality` is on (DataLikelihood.cpp:32-35), and
            // `lnquality` alone otherwise (line 38). The gate is
            // exposed via `Parameters::use_mapping_quality`.
            //
            // TODO(M3-indel): `base_quality_sum` is a per-base Q for
            // SNP / reference observations but a sum-across-bases for
            // indels (M1 simplification; upstream harmonic-sum scaling
            // is also deferred). Clamping to `u8::MAX` silently
            // saturates at Q=255, so ≥~8 Q30 inserted bases collapse
            // to the same effective error term. Revisit alongside the
            // indel-BQ scaling port.
            let q = obs.base_quality_sum.min(u8::MAX as u32) as u8;
            let ln_bq = phred_to_ln(q);
            prod_q_out += if params.use_mapping_quality {
                let ln_mq = phred_to_ln(obs.mapq);
                // In log-space, a larger (less negative) value means a
                // higher error probability — taking `max` selects the
                // less-informative quality, matching upstream's intent.
                ln_bq.max(ln_mq)
            } else {
                ln_bq
            };
            count_out += 1;
        }
    }

    // Read-dependence factor (RDF): downweight successive error reads
    // to avoid over-penalising batch-correlated errors. Upstream
    // `DataLikelihood.cpp:147-149`.
    if count_out > 1 {
        let n = count_out as f64;
        prod_q_out *= (1.0 + (n - 1.0) * params.read_dependence_factor) / n;
    }

    let in_total: i64 = in_counts.iter().sum();
    if in_total == 0 {
        // All observations are errors under this genotype — return just
        // the error term (upstream line 151-152).
        return prod_q_out;
    }

    let probs = genotype.allele_probabilities();
    prod_q_out + multinomial_sampling_prob_ln(&probs, &in_counts)
}

/// Compute `log P(observations | genotype)` for each of a candidate
/// genotype set, returning a vector aligned with `genotypes`.
///
/// This is the sample-level building block the M3 Phase C posterior
/// code aggregates. It's a thin wrapper over [`sample_log_likelihood`]
/// that exists so callers don't have to rewrite the same loop at every
/// site.
pub fn sample_log_likelihoods(
    genotypes: &[Genotype],
    observations: &[AlleleObservation],
    params: &Parameters,
) -> Vec<f64> {
    genotypes
        .iter()
        .map(|g| sample_log_likelihood(g, observations, params))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genotype::enumerate_genotypes;
    use fb_core::{Allele, Strand};

    /// Build a `Parameters` with the given RDF and upstream-default
    /// `use_mapping_quality = true`. The tests in this module were
    /// originally written against a bare `rdf: f64` argument; this
    /// helper keeps them terse.
    fn params(rdf: f64) -> Parameters {
        Parameters {
            read_dependence_factor: rdf,
            ..Parameters::default()
        }
    }

    fn ref_allele() -> Allele {
        Allele::reference(100, vec![b'A'])
    }
    fn snp_ag() -> Allele {
        Allele::snp(100, b'A', b'G')
    }

    fn obs(allele: Allele, mapq: u8, bq_sum: u32) -> AlleleObservation {
        AlleleObservation {
            allele,
            read_name: "r".into(),
            mapq,
            base_quality_sum: bq_sum,
            strand: Strand::Forward,
            read_position: 0,
            is_proper_pair: true,
            read_ref_start: 0,
            per_base_quals: Vec::new(),
        }
    }

    #[test]
    fn phred_to_ln_matches_closed_form() {
        // Q=30 → err_prob = 10^-3 → ln = -3 * ln(10).
        let got = phred_to_ln(30);
        let expected = -3.0 * std::f64::consts::LN_10;
        assert!((got - expected).abs() < 1e-12, "got {got}");
        // Q=0 → err_prob = 1 → ln = 0.
        assert!((phred_to_ln(0) - 0.0).abs() < 1e-12);
    }

    #[test]
    fn empty_observations_returns_zero() {
        let gt = Genotype::from_alleles(vec![ref_allele(), ref_allele()]);
        assert_eq!(sample_log_likelihood(&gt, &[], &params(0.9)), 0.0);
    }

    #[test]
    fn all_reference_observations_favour_hom_ref() {
        let obs_set = vec![obs(ref_allele(), 60, 30); 10];
        let alleles = vec![ref_allele(), snp_ag()];
        let gts = enumerate_genotypes(&alleles, 2);
        let lls = sample_log_likelihoods(&gts, &obs_set, &params(0.9));

        // Find the indices of the three genotypes by tag.
        let mut hom_ref_ll = f64::NEG_INFINITY;
        let mut het_ll = f64::NEG_INFINITY;
        let mut hom_alt_ll = f64::NEG_INFINITY;
        for (g, &ll) in gts.iter().zip(lls.iter()) {
            match g.str_tag().as_str() {
                "REF" => hom_ref_ll = ll,
                "SNP" => hom_alt_ll = ll,
                _ => het_ll = ll,
            }
        }
        assert!(
            hom_ref_ll > het_ll && het_ll > hom_alt_ll,
            "hom-ref {hom_ref_ll} > het {het_ll} > hom-alt {hom_alt_ll}"
        );
    }

    #[test]
    fn mixed_observations_favour_het() {
        let mut obs_set: Vec<AlleleObservation> = Vec::new();
        for _ in 0..10 {
            obs_set.push(obs(ref_allele(), 60, 30));
        }
        for _ in 0..10 {
            obs_set.push(obs(snp_ag(), 60, 30));
        }
        let alleles = vec![ref_allele(), snp_ag()];
        let gts = enumerate_genotypes(&alleles, 2);
        let lls = sample_log_likelihoods(&gts, &obs_set, &params(0.9));

        let (mut hom_ref_ll, mut het_ll, mut hom_alt_ll) =
            (f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
        for (g, &ll) in gts.iter().zip(lls.iter()) {
            match g.str_tag().as_str() {
                "REF" => hom_ref_ll = ll,
                "SNP" => hom_alt_ll = ll,
                _ => het_ll = ll,
            }
        }
        assert!(
            het_ll > hom_ref_ll && het_ll > hom_alt_ll,
            "het {het_ll} > hom_ref {hom_ref_ll} and hom_alt {hom_alt_ll}"
        );
    }

    #[test]
    fn all_outside_observations_yields_error_term_only() {
        // Genotype is hom-SNP; all observations are REF → all "outside".
        // With RDF=1 the outside sum is unscaled (`(1 + (n-1)*1)/n = 1`).
        let gt = Genotype::from_alleles(vec![snp_ag(), snp_ag()]);
        let obs_set = vec![obs(ref_allele(), 60, 30); 5];
        let ll = sample_log_likelihood(&gt, &obs_set, &params(1.0));
        // prod_q_out = 5 * max(ln_bq, ln_mq). For Q=30 BQ and Q=60 MQ:
        // ln_bq = -3*ln10 ≈ -6.91; ln_mq = -6*ln10 ≈ -13.82; max = ln_bq.
        let expected = 5.0 * (-3.0 * std::f64::consts::LN_10);
        assert!(
            (ll - expected).abs() < 1e-10,
            "got {ll} expected {expected}"
        );
    }

    #[test]
    fn read_dependence_factor_direction() {
        // Upstream's RDF semantics: prod_q_out *= (1 + (n-1)*RDF) / n.
        //   RDF=0 → factor = 1/n (reads fully correlated; effective n=1;
        //                         most lenient / least negative ll).
        //   RDF=1 → factor = 1   (reads fully independent; no scaling).
        //   RDF=0.9 → factor = (1 + 9*0.9)/10 = 0.91 (upstream default).
        let gt = Genotype::from_alleles(vec![snp_ag(), snp_ag()]);
        let obs_set = vec![obs(ref_allele(), 60, 30); 10];

        let ll_rdf_0 = sample_log_likelihood(&gt, &obs_set, &params(0.0));
        let ll_rdf_default = sample_log_likelihood(&gt, &obs_set, &params(0.9));
        let ll_rdf_1 = sample_log_likelihood(&gt, &obs_set, &params(1.0));

        // RDF=0 is the most lenient (closest to 0), RDF=1 the strictest.
        assert!(
            ll_rdf_0 > ll_rdf_default && ll_rdf_default > ll_rdf_1,
            "expected RDF=0 ({ll_rdf_0}) > RDF=0.9 ({ll_rdf_default}) > RDF=1 ({ll_rdf_1})"
        );

        // Closed-form cross-checks:
        let one_err = -3.0 * std::f64::consts::LN_10;
        //   RDF=1: 10 * err * 1.0 = 10 * err
        let expected_rdf_1 = 10.0 * one_err;
        //   RDF=0.9: 10 * err * 0.91
        let expected_rdf_default = 10.0 * one_err * 0.91;
        //   RDF=0: 10 * err * 0.1 = err
        let expected_rdf_0 = one_err;
        assert!((ll_rdf_1 - expected_rdf_1).abs() < 1e-10);
        assert!((ll_rdf_default - expected_rdf_default).abs() < 1e-10);
        assert!((ll_rdf_0 - expected_rdf_0).abs() < 1e-10);
    }

    #[test]
    fn single_outside_observation_not_rdf_scaled() {
        // Upstream: `if (countOut > 1)` — a single outside obs is not
        // affected by RDF.
        let gt = Genotype::from_alleles(vec![snp_ag(), snp_ag()]);
        let obs_set = vec![obs(ref_allele(), 60, 30)];
        let ll_rdf_0 = sample_log_likelihood(&gt, &obs_set, &params(0.0));
        let ll_rdf_1 = sample_log_likelihood(&gt, &obs_set, &params(1.0));
        assert!((ll_rdf_0 - ll_rdf_1).abs() < 1e-15);
    }

    #[test]
    fn het_allele_probabilities_reduce_to_one_half() {
        // Sanity: under a het genotype, each allele has sampling prob 0.5.
        let gt = Genotype::from_alleles(vec![ref_allele(), snp_ag()]);
        let probs = gt.allele_probabilities();
        assert_eq!(probs.len(), 2);
        assert!((probs[0] - 0.5).abs() < 1e-15);
        assert!((probs[1] - 0.5).abs() < 1e-15);
    }

    #[test]
    fn higher_base_quality_penalises_errors_more() {
        // A Q=40 error is more "certain" than Q=10 → larger penalty
        // (more negative log-likelihood).
        let gt = Genotype::from_alleles(vec![snp_ag(), snp_ag()]);
        let obs_q40 = vec![obs(ref_allele(), 60, 40)];
        let obs_q10 = vec![obs(ref_allele(), 60, 10)];

        let ll_q40 = sample_log_likelihood(&gt, &obs_q40, &params(0.0));
        let ll_q10 = sample_log_likelihood(&gt, &obs_q10, &params(0.0));
        assert!(
            ll_q40 < ll_q10,
            "Q=40 error ({ll_q40}) should be MORE penalising than Q=10 ({ll_q10})"
        );
    }

    #[test]
    fn ml_genotype_selection_matches_truth_hom_ref() {
        let obs_set = vec![obs(ref_allele(), 60, 30); 20];
        let alleles = vec![ref_allele(), snp_ag()];
        let gts = enumerate_genotypes(&alleles, 2);
        let lls = sample_log_likelihoods(&gts, &obs_set, &params(0.9));

        let best_idx = lls
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .unwrap()
            .0;
        assert_eq!(gts[best_idx].str_tag(), "REF");
    }

    #[test]
    fn ml_genotype_selection_matches_truth_het() {
        // 15 REF + 15 SNP at Q30, MQ60 → het should win.
        let mut obs_set = vec![obs(ref_allele(), 60, 30); 15];
        obs_set.extend(vec![obs(snp_ag(), 60, 30); 15]);
        let alleles = vec![ref_allele(), snp_ag()];
        let gts = enumerate_genotypes(&alleles, 2);
        let lls = sample_log_likelihoods(&gts, &obs_set, &params(0.9));

        let best_idx = lls
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .unwrap()
            .0;
        assert!(
            gts[best_idx].str_tag().contains('/'),
            "expected a het genotype, got {}",
            gts[best_idx].str_tag()
        );
    }

    #[test]
    fn tetraploid_genotype_likelihood_sane() {
        // 4 REF observations vs. a tetraploid genotype. The 4/4 hom-ref
        // genotype (allele probs [1.0]) must outscore every other
        // enumerated genotype. Smoke-test the port works for ploidy > 2.
        let obs_set = vec![obs(ref_allele(), 60, 30); 4];
        let alleles = vec![ref_allele(), snp_ag()];
        let gts = enumerate_genotypes(&alleles, 4);
        assert_eq!(gts.len(), 5); // C(2+4-1, 4) = C(5, 4) = 5.
        let lls = sample_log_likelihoods(&gts, &obs_set, &params(0.9));

        // Find hom-ref (all 4 refs) and pick the ML genotype.
        let (best_idx, _) = lls
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .unwrap();
        assert!(
            gts[best_idx].homozygous,
            "expected homozygous genotype, got {}",
            gts[best_idx].str_tag()
        );
        assert_eq!(gts[best_idx].elements[0].allele, ref_allele());
    }

    #[test]
    fn mixed_inside_outside_closed_form() {
        // Het genotype (REF/SNP, allele probs [0.5, 0.5]), observations:
        //   8 REF inside + 8 SNP inside + 2 "unknown" outside
        // (Approximated here by making the "outside" observation an
        // indel not in the genotype so it is matched as outside.)
        let del = Allele::deletion(100, vec![b'A']);
        let mut obs_set: Vec<AlleleObservation> = Vec::new();
        obs_set.extend(vec![obs(ref_allele(), 60, 30); 8]);
        obs_set.extend(vec![obs(snp_ag(), 60, 30); 8]);
        obs_set.extend(vec![obs(del, 60, 30); 2]);

        let gt = Genotype::from_alleles(vec![ref_allele(), snp_ag()]);
        // RDF=1 to skip the scaling and isolate the closed form.
        let ll = sample_log_likelihood(&gt, &obs_set, &params(1.0));

        // Closed form:
        //   prod_q_out = 2 * max(phred2ln(30), phred2ln(60))
        //              = 2 * (-3 * ln10)
        //   multinomial: C(16, 8, 8) * (0.5)^8 * (0.5)^8
        //     log = lgamma(17) - 2*lgamma(9) + 16 * ln(0.5)
        let ln10 = std::f64::consts::LN_10;
        let prod_q_out = 2.0 * (-3.0 * ln10);
        let lgamma17 = libm::lgamma(17.0);
        let lgamma9 = libm::lgamma(9.0);
        let multi = lgamma17 - 2.0 * lgamma9 + 16.0 * (0.5f64).ln();
        let expected = prod_q_out + multi;
        assert!(
            (ll - expected).abs() < 1e-10,
            "got {ll} expected {expected} (diff {})",
            (ll - expected).abs()
        );
    }

    #[test]
    fn ml_genotype_selection_matches_truth_hom_alt() {
        let obs_set = vec![obs(snp_ag(), 60, 30); 20];
        let alleles = vec![ref_allele(), snp_ag()];
        let gts = enumerate_genotypes(&alleles, 2);
        let lls = sample_log_likelihoods(&gts, &obs_set, &params(0.9));

        let best_idx = lls
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .unwrap()
            .0;
        assert_eq!(gts[best_idx].str_tag(), "SNP");
    }
}
