//! Log-priors on a single-sample genotype.
//!
//! Port of the subset of upstream freebayes v1.3.10's
//! `GenotypeCombo::calculatePosteriors` prior terms
//! (`src/Genotype.cpp:1495-1596`) that apply to a single sample. The
//! full upstream machinery sums four prior components
//! (`priorProbG_Af`, `priorProbAf`, `priorProbObservations`,
//! `priorProbGenotypesGivenHWE`); this M3 Phase C-1 port covers just
//! `priorProbG_Af` — the probability of the genotype combination given
//! a uniform allele frequency — because that is the only term the
//! single-sample caller needs and it reduces to a closed-form scalar.
//!
//! ## Math
//! Upstream's `GenotypeCombo::probabilityGivenAlleleFrequencyln(permute)`
//! (`src/Genotype.cpp:1391-1405`) returns
//!
//! ```text
//!   (permute ? permutationsln : 0) - multinomialCoefficientLn(n, counts())
//! ```
//!
//! where `n = numberOfAlleles()` is the total ploidy across the combo,
//! `counts()` aggregates allele counts across the combo, and
//! `permutationsln` is the sum of `permutationsln` across the combo's
//! constituent genotypes. For a **single-sample** combo those two terms
//! collapse: `permutationsln` equals `multinomialCoefficientLn(ploidy,
//! counts())` by definition (`src/Genotype.h:53, 68`). That gives
//!
//! - `permute = true`  → `0`
//! - `permute = false` → `-genotype.permutations_ln`
//!
//! So under the default (permute=true) **all single-sample genotypes
//! receive the same prior** — which is the correct Bayesian stance
//! without cross-sample allele-frequency information: the posterior is
//! just the normalised data likelihood. Multi-sample combos (Phase C-2)
//! will break the symmetry via the cross-sample `counts()`
//! aggregation.
//!
//! ## Deferred to later M3 sub-phases
//! - Ewens sampling formula prior on allele frequency
//!   (`alleleFrequencyProbabilityln`, upstream `Ewens.cpp`).
//! - HWE genotype-frequency prior (`hweProbGenotypeFrequencyln`,
//!   upstream `src/Genotype.cpp:1524-1529`).
//! - Allele-balance / strand-balance / placement binomial priors
//!   (upstream `src/Genotype.cpp:1531-1570`).
//! - Contamination-aware priors (`Contamination.cpp`).

use crate::genotype::Genotype;
use crate::params::Parameters;

/// Log-prior for a single-sample genotype, mirroring upstream's
/// `probabilityGivenAlleleFrequencyln` (`src/Genotype.cpp:1391-1405`)
/// specialised to the single-sample case.
///
/// Returns:
/// - `0.0` when `params.use_permutations_prior` is `true` (upstream
///   `permute == true`) — the per-genotype permutations exactly cancel
///   the per-combo multinomial coefficient.
/// - `-genotype.permutations_ln` otherwise — the coefficient term
///   survives, disfavouring heterozygotes.
#[inline]
pub fn genotype_log_prior(genotype: &Genotype, params: &Parameters) -> f64 {
    if params.use_permutations_prior {
        0.0
    } else {
        -genotype.permutations_ln
    }
}

/// Compute the log-prior for each of a candidate genotype set,
/// returning a vector aligned with `genotypes`.
pub fn genotype_log_priors(genotypes: &[Genotype], params: &Parameters) -> Vec<f64> {
    genotypes
        .iter()
        .map(|g| genotype_log_prior(g, params))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genotype::enumerate_genotypes;
    use fb_core::Allele;

    fn ref_a() -> Allele {
        Allele::reference(100, vec![b'A'])
    }
    fn snp_ag() -> Allele {
        Allele::snp(100, b'A', b'G')
    }

    #[test]
    fn default_prior_is_flat_single_sample() {
        // Upstream's probabilityGivenAlleleFrequencyln(permute=true) for
        // a single-sample combo collapses to 0 for every genotype —
        // permutations_ln exactly cancels multinomialCoefficientLn.
        let alleles = vec![ref_a(), snp_ag()];
        let gts = enumerate_genotypes(&alleles, 2);
        let p = Parameters::default();
        let priors = genotype_log_priors(&gts, &p);
        assert_eq!(priors.len(), 3);
        assert!(
            priors.iter().all(|&x| x == 0.0),
            "default single-sample prior must be flat, got {priors:?}"
        );
    }

    #[test]
    fn homozygous_prior_is_zero_in_both_modes() {
        // For homs permutations_ln = 0, so neither arm changes the value.
        let gt = Genotype::from_alleles(vec![ref_a(), ref_a()]);
        let p_perm = Parameters::default();
        let p_flat = Parameters {
            use_permutations_prior: false,
            ..Parameters::default()
        };
        assert_eq!(genotype_log_prior(&gt, &p_perm), 0.0);
        assert_eq!(genotype_log_prior(&gt, &p_flat), 0.0);
    }

    #[test]
    fn heterozygous_prior_without_permute_disfavours_het() {
        // permute=false: prior = -permutations_ln = -ln(2) ≈ -0.693.
        let gt = Genotype::from_alleles(vec![ref_a(), snp_ag()]);
        let p = Parameters {
            use_permutations_prior: false,
            ..Parameters::default()
        };
        let got = genotype_log_prior(&gt, &p);
        assert!((got - (-2.0f64.ln())).abs() < 1e-12, "got {got}");
    }

    #[test]
    fn default_prior_matches_negated_permutations_ln_with_permute_off() {
        // With permute off, prior is the negative of the cached value.
        let alleles = vec![ref_a(), snp_ag()];
        let gts = enumerate_genotypes(&alleles, 2);
        let p = Parameters {
            use_permutations_prior: false,
            ..Parameters::default()
        };
        let priors = genotype_log_priors(&gts, &p);
        for (g, prior) in gts.iter().zip(priors.iter()) {
            assert_eq!(*prior, -g.permutations_ln);
        }
    }

    #[test]
    fn triploid_three_distinct_prior_without_permute_is_neg_ln_6() {
        // 3 distinct alleles in a triploid, permute=false → -ln(6).
        let gt = Genotype::from_alleles(vec![ref_a(), snp_ag(), Allele::snp(100, b'A', b'T')]);
        let p = Parameters {
            use_permutations_prior: false,
            ..Parameters::default()
        };
        let got = genotype_log_prior(&gt, &p);
        assert!((got - (-6.0f64.ln())).abs() < 1e-12);
    }
}
