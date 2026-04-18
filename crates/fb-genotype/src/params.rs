//! Global configuration for the Bayesian genotype model.
//!
//! Collects the tunables that the data-likelihood, prior, and posterior
//! paths read — mirroring the subset of upstream freebayes'
//! `Parameters` struct (`src/Parameters.h`) that our port currently
//! consumes. Adding a field here and threading it into a specific
//! module is the expected way to extend the model in Phase C and
//! beyond; keep the default values locked to upstream's for parity.

/// Configuration for the Bayesian model.
///
/// Defaults mirror upstream's CLI defaults where they exist, and lean
/// conservative otherwise. Set fields explicitly to diverge —
/// round-tripping through [`Parameters::default`] must always produce
/// upstream-compatible behaviour.
#[derive(Debug, Clone, Copy)]
pub struct Parameters {
    /// Read-dependence factor (upstream `--read-dependence-factor`,
    /// `Parameters.cpp:475`). `0.0` = fully correlated reads (effective
    /// sample size 1); `1.0` = fully independent; upstream default
    /// `0.9`.
    pub read_dependence_factor: f64,
    /// Whether to cap per-observation log-error by mapping quality
    /// (`max(lnquality, lnmapQuality)` upstream). Upstream's default is
    /// `false` (`standardGLs = false` & `useMappingQuality = false`);
    /// the M3 Phase B port applied the cap unconditionally. Phase C
    /// exposes the gate.
    pub use_mapping_quality: bool,
    /// Whether the prior should weight genotypes by the number of
    /// ordered realisations (upstream `permute == true`). Mirrors
    /// `multinomialCoefficientLn(ploidy, counts)` — equivalent to
    /// `Genotype::permutations_ln`. Upstream's default is `true`.
    pub use_permutations_prior: bool,
    /// Target ploidy. The caller uses this to drive
    /// [`enumerate_genotypes`](crate::genotype::enumerate_genotypes).
    /// Upstream default: 2 (diploid).
    pub ploidy: u32,
}

impl Default for Parameters {
    fn default() -> Self {
        Self {
            read_dependence_factor: crate::data_likelihood::DEFAULT_READ_DEPENDENCE_FACTOR,
            // Upstream ships with mapping-quality gating off by default
            // (`standardGLs = false`), but the common-use-case
            // `--legacy-gls` path turns it on. Pick the legacy path as
            // our default because it matches what most freebayes users
            // actually run; flip to `false` for strict
            // bit-compatibility with `--standard-gls`.
            use_mapping_quality: true,
            use_permutations_prior: true,
            ploidy: 2,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn defaults_are_upstream_compatible() {
        let p = Parameters::default();
        assert_eq!(p.read_dependence_factor, 0.9);
        assert_eq!(p.ploidy, 2);
        assert!(p.use_permutations_prior);
        assert!(p.use_mapping_quality);
    }
}
