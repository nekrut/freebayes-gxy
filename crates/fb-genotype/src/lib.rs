//! Bayesian genotype model for freebayes-gxy.
//!
//! Ports the numerical core of freebayes' `Genotype.cpp` /
//! `Multinomial.cpp` / `Sum.h` (and the log-space helpers from
//! `Utility.cpp` that those depend on). M3 Phase A — this milestone —
//! lands the reusable primitives:
//!
//! - [`mod@sum`] — linear and log-domain reductions.
//! - [`mod@multinomial`] — factorial-ln, multinomial coefficient, sampling
//!   probability.
//! - [`mod@genotype`] — the [`Genotype`] data type and the
//!   [`enumerate_genotypes`] multichoose port.
//!
//! Phase B (next) layers the per-read data likelihood
//! (`DataLikelihood.cpp`) and the single-sample genotype posterior on
//! top. Phase C wires it into the per-site pipeline so `fb-cli` can
//! emit real VCF genotype calls.

pub mod genotype;
pub mod multinomial;
pub mod sum;

pub use genotype::{enumerate_genotypes, Genotype, GenotypeElement};
pub use multinomial::{
    factorial_ln, multinomial_coefficient_ln, multinomial_sampling_prob_ln, pow_ln,
    sampling_prob_ln,
};
pub use sum::{log_add, log_sum_exp, sum, sum_i64};

/// A numerically stable log-sum-exp over a pair of log-space values.
///
/// Retained as a deprecated alias for the [`log_add`] helper in
/// [`mod@sum`] — removes in M3 Phase B.
#[deprecated(since = "0.0.1", note = "use fb_genotype::log_add instead")]
pub fn logsumexp2(a: f64, b: f64) -> f64 {
    log_add(a, b)
}
