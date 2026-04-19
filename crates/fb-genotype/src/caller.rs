//! Single-sample genotype caller.
//!
//! Combines the M3 Phase A/B primitives (`Genotype`,
//! `sample_log_likelihood`) with the M3 Phase C-1 prior into a
//! Bayesian posterior over candidate genotypes, returns the ML
//! genotype, and produces a Phred-scaled genotype quality (GQ).
//!
//! ## Model
//! For a set of candidate genotypes `G = {g_1, …, g_K}` and a single
//! sample's observations `O`:
//!
//! ```text
//!   log P(g_i | O) = log P(O | g_i) + log P(g_i) − log P(O)
//!   log P(O)       = logsumexp_j [ log P(O | g_j) + log P(g_j) ]
//! ```
//!
//! The returned [`Call`] carries normalised log-posteriors, the
//! maximum-a-posteriori (MAP) index, and the GQ
//! (`−10 * log10(1 − exp(log_posterior[best]))`).
//!
//! ## Scope limits (M3 Phase C-1)
//! - Single sample only — multi-sample joint calls arrive with the
//!   `GenotypeCombo` port in a later sub-phase.
//! - Only the permutations-weighted / flat prior from
//!   [`crate::prior::genotype_log_prior`] is applied — the Ewens and
//!   HWE terms are Phase C-2+.
//! - Variant-quality scoring (upstream's QUAL column) is also
//!   Phase C-2+; this module only emits sample-level GQ.

use std::f64;

use fb_core::AlleleObservation;

use crate::data_likelihood::sample_log_likelihoods;
use crate::genotype::Genotype;
use crate::params::Parameters;
use crate::prior::genotype_log_priors;
use crate::sum::log_sum_exp;

/// Result of a single-sample genotype call.
#[derive(Debug, Clone)]
pub struct Call {
    /// Index into the candidate genotype slice of the
    /// maximum-a-posteriori call.
    pub best_index: usize,
    /// Log-data-likelihood per candidate, aligned with the input slice.
    pub log_likelihoods: Vec<f64>,
    /// Log-prior per candidate, aligned with the input slice.
    pub log_priors: Vec<f64>,
    /// Normalised log-posterior per candidate, aligned with the input
    /// slice. `logsumexp(log_posteriors) == 0` up to rounding.
    pub log_posteriors: Vec<f64>,
    /// Evidence term `log P(O)` — the normalising constant used to
    /// produce `log_posteriors`.
    pub log_evidence: f64,
    /// Phred-scaled genotype quality of the MAP call:
    /// `-10 * log10(1 - exp(log_posteriors[best_index]))`. When the
    /// posterior is ≥ `1 - 1e-300` (effectively certain), the raw
    /// formula overflows; the value is capped at [`MAX_GQ`].
    pub genotype_quality: f64,
}

/// Upstream freebayes caps per-site QUAL at ~`3000` via ttmath; our
/// simpler port caps GQ at `2000` which keeps the VCF column sane and
/// avoids f64-overflow when `log1p(-exp(log_p))` goes to `-inf`.
pub const MAX_GQ: f64 = 2000.0;

/// Score a set of candidate genotypes against a sample's observations
/// and return the MAP call.
///
/// `genotypes` must be non-empty. `observations` may be empty, in
/// which case every genotype scores at `ll == 0` and the prior alone
/// selects the MAP call.
pub fn call_genotype(
    genotypes: &[Genotype],
    observations: &[AlleleObservation],
    params: &Parameters,
) -> Call {
    assert!(
        !genotypes.is_empty(),
        "call_genotype requires at least one candidate"
    );

    let log_likelihoods = sample_log_likelihoods(genotypes, observations, params);
    let log_priors = genotype_log_priors(genotypes, params);

    // Un-normalised log-posterior proportional to prior × likelihood.
    let joint: Vec<f64> = log_likelihoods
        .iter()
        .zip(log_priors.iter())
        .map(|(ll, lp)| ll + lp)
        .collect();

    let log_evidence = log_sum_exp(&joint);
    let log_posteriors: Vec<f64> = joint.iter().map(|&j| j - log_evidence).collect();

    // MAP index — break ties deterministically in index order (upstream
    // sorts its `GenotypeCombo` vector; ties are rare in practice and
    // our canonical enumeration is already deterministic).
    let best_index = log_posteriors
        .iter()
        .enumerate()
        .fold(
            0usize,
            |best, (i, &lp)| {
                if lp > log_posteriors[best] {
                    i
                } else {
                    best
                }
            },
        );

    let genotype_quality = log_posterior_to_gq(log_posteriors[best_index]);

    Call {
        best_index,
        log_likelihoods,
        log_priors,
        log_posteriors,
        log_evidence,
        genotype_quality,
    }
}

/// Convert a log-posterior probability to a Phred-scaled genotype
/// quality — the probability of the MAP call being **wrong**:
///
/// ```text
///   GQ = -10 * log10(1 - exp(log_p))
/// ```
///
/// Uses `ln_1p(-exp(log_p))` for numerical stability when `log_p` is
/// close to `0` (posterior near 1). Capped at [`MAX_GQ`] to avoid
/// overflow when the posterior is effectively 1.
#[inline]
pub fn log_posterior_to_gq(log_p: f64) -> f64 {
    // -exp(log_p).ln_1p() = ln(1 - exp(log_p)) — valid because
    // exp(log_p) ∈ [0, 1].
    let ln_err = (-log_p.exp()).ln_1p();
    if ln_err.is_infinite() {
        return MAX_GQ;
    }
    let gq = ln_err * -10.0 / f64::consts::LN_10;
    gq.clamp(0.0, MAX_GQ)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genotype::enumerate_genotypes;
    use fb_core::{Allele, Strand};

    fn ref_a() -> Allele {
        Allele::reference(100, vec![b'A'])
    }
    fn snp_ag() -> Allele {
        Allele::snp(100, b'A', b'G')
    }

    fn obs(allele: Allele, mapq: u8, bq: u32) -> AlleleObservation {
        AlleleObservation {
            allele,
            read_name: "r".into(),
            mapq,
            base_quality_sum: bq,
            strand: Strand::Forward,
            read_position: 0,
            is_proper_pair: true,
            read_ref_start: 0,
            per_base_quals: Vec::new(),
        }
    }

    #[test]
    fn log_posterior_to_gq_round_numbers() {
        // posterior = 0.9 → 1 - 0.9 = 0.1 → GQ = -10 * log10(0.1) = 10.
        let log_p = 0.9f64.ln();
        let gq = log_posterior_to_gq(log_p);
        assert!((gq - 10.0).abs() < 1e-10);

        // posterior = 0.99 → GQ = 20.
        let log_p = 0.99f64.ln();
        let gq = log_posterior_to_gq(log_p);
        assert!((gq - 20.0).abs() < 1e-10);

        // posterior = 0.999 → GQ = 30.
        let log_p = 0.999f64.ln();
        let gq = log_posterior_to_gq(log_p);
        assert!((gq - 30.0).abs() < 1e-10);
    }

    #[test]
    fn log_posterior_to_gq_clamps_at_extremes() {
        // Posterior = 1 → 1 - 1 = 0 → log = -inf → clamped to MAX_GQ.
        let gq = log_posterior_to_gq(0.0);
        assert_eq!(gq, MAX_GQ);

        // Posterior = 0 → 1 - 0 = 1 → log = 0 → GQ = 0.
        let gq = log_posterior_to_gq(f64::NEG_INFINITY);
        assert_eq!(gq, 0.0);
    }

    #[test]
    fn call_picks_hom_ref_on_clean_ref_data() {
        let obs_set = vec![obs(ref_a(), 60, 30); 20];
        let alleles = vec![ref_a(), snp_ag()];
        let gts = enumerate_genotypes(&alleles, 2);
        let params = Parameters::default();
        let call = call_genotype(&gts, &obs_set, &params);

        assert_eq!(gts[call.best_index].str_tag(), "REF");
        // With 20 matching Q30 REF observations the posterior mass
        // should concentrate on hom-ref → high GQ.
        assert!(
            call.genotype_quality > 50.0,
            "expected GQ > 50, got {}",
            call.genotype_quality
        );
    }

    #[test]
    fn call_picks_het_on_mixed_data() {
        let mut obs_set = vec![obs(ref_a(), 60, 30); 15];
        obs_set.extend(vec![obs(snp_ag(), 60, 30); 15]);
        let alleles = vec![ref_a(), snp_ag()];
        let gts = enumerate_genotypes(&alleles, 2);
        let params = Parameters::default();
        let call = call_genotype(&gts, &obs_set, &params);

        assert!(
            gts[call.best_index].str_tag().contains('/'),
            "expected het, got {}",
            gts[call.best_index].str_tag()
        );
        assert!(call.genotype_quality > 50.0);
    }

    #[test]
    fn call_picks_hom_alt_on_clean_alt_data() {
        let obs_set = vec![obs(snp_ag(), 60, 30); 20];
        let alleles = vec![ref_a(), snp_ag()];
        let gts = enumerate_genotypes(&alleles, 2);
        let params = Parameters::default();
        let call = call_genotype(&gts, &obs_set, &params);

        assert_eq!(gts[call.best_index].str_tag(), "SNP");
        assert!(call.genotype_quality > 50.0);
    }

    #[test]
    fn posteriors_sum_to_one() {
        let obs_set = vec![obs(ref_a(), 60, 30); 5];
        let alleles = vec![ref_a(), snp_ag()];
        let gts = enumerate_genotypes(&alleles, 2);
        let call = call_genotype(&gts, &obs_set, &Parameters::default());

        let total: f64 = call.log_posteriors.iter().map(|&x| x.exp()).sum();
        assert!(
            (total - 1.0).abs() < 1e-10,
            "posteriors should sum to 1, got {total}"
        );
    }

    #[test]
    fn low_evidence_yields_low_gq() {
        // 2 REF + 2 SNP observations at low Q → ambiguous → low GQ.
        let mut obs_set = vec![obs(ref_a(), 20, 10); 2];
        obs_set.extend(vec![obs(snp_ag(), 20, 10); 2]);
        let alleles = vec![ref_a(), snp_ag()];
        let gts = enumerate_genotypes(&alleles, 2);
        let call = call_genotype(&gts, &obs_set, &Parameters::default());
        assert!(
            call.genotype_quality < 30.0,
            "expected low GQ on ambiguous data, got {}",
            call.genotype_quality
        );
    }

    #[test]
    fn empty_observations_default_prior_is_flat() {
        // Under the default (permute=true) single-sample prior is flat,
        // so the MAP on empty observations falls to the deterministic
        // tiebreak — the first genotype in enumeration order (hom-ref
        // for [REF, SNP] at diploid).
        let alleles = vec![ref_a(), snp_ag()];
        let gts = enumerate_genotypes(&alleles, 2);
        let call = call_genotype(&gts, &[], &Parameters::default());
        assert_eq!(gts[call.best_index].str_tag(), "REF");
        // All three posteriors are equal (1/3 each) with no evidence
        // and flat prior → GQ ≈ -10 * log10(2/3) ≈ 1.76.
        assert!(call.genotype_quality < 5.0);
    }

    #[test]
    fn empty_observations_permute_off_disfavours_het() {
        // With permute=false the prior becomes -permutations_ln, so
        // heterozygous genotypes get -ln(2) and homs get 0. MAP on
        // empty obs is a hom (tiebreak to index 0 = hom-ref).
        let alleles = vec![ref_a(), snp_ag()];
        let gts = enumerate_genotypes(&alleles, 2);
        let params = Parameters {
            use_permutations_prior: false,
            ..Parameters::default()
        };
        let call = call_genotype(&gts, &[], &params);
        assert!(
            !gts[call.best_index].str_tag().contains('/'),
            "expected homozygous MAP under permute=false prior, got {}",
            gts[call.best_index].str_tag()
        );
    }

    #[test]
    fn flat_prior_still_selects_ml_genotype() {
        // With permutations disabled the prior is -permutations_ln → not
        // flat — but 10 Q30 REF observations swamp the prior: MAP = hom-ref.
        let obs_set = vec![obs(ref_a(), 60, 30); 10];
        let alleles = vec![ref_a(), snp_ag()];
        let gts = enumerate_genotypes(&alleles, 2);
        let params = Parameters {
            use_permutations_prior: false,
            ..Parameters::default()
        };
        let call = call_genotype(&gts, &obs_set, &params);
        assert_eq!(gts[call.best_index].str_tag(), "REF");
    }

    #[test]
    fn log_evidence_equals_logsumexp_of_joints() {
        let obs_set = vec![obs(ref_a(), 60, 30); 5];
        let alleles = vec![ref_a(), snp_ag()];
        let gts = enumerate_genotypes(&alleles, 2);
        let call = call_genotype(&gts, &obs_set, &Parameters::default());

        let manual_evidence = {
            let joint: Vec<f64> = call
                .log_likelihoods
                .iter()
                .zip(call.log_priors.iter())
                .map(|(l, p)| l + p)
                .collect();
            log_sum_exp(&joint)
        };
        assert!(
            (call.log_evidence - manual_evidence).abs() < 1e-12,
            "got {} expected {}",
            call.log_evidence,
            manual_evidence
        );
    }

    #[test]
    #[should_panic(expected = "at least one candidate")]
    fn empty_genotypes_panics() {
        call_genotype(&[], &[obs(ref_a(), 60, 30)], &Parameters::default());
    }
}
