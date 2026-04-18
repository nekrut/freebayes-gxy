//! Bayesian genotype model for freebayes-gxy.
//!
//! Ports the numerical core of freebayes' `Genotype.cpp` /
//! `Multinomial.cpp` / `DataLikelihood.cpp` / `Sum.h` (and the
//! log-space helpers from `Utility.cpp` that those depend on).
//!
//! - [`mod@sum`] — linear and log-domain reductions.
//! - [`mod@multinomial`] — factorial-ln, multinomial coefficient,
//!   sampling probability.
//! - [`mod@genotype`] — the [`Genotype`] data type and the
//!   [`enumerate_genotypes`] multichoose port.
//! - [`mod@data_likelihood`] — `P(observations | genotype)` for a
//!   single sample. The standard-GL path of upstream
//!   `probObservedAllelesGivenGenotype`.
//!
//! ## Roadmap
//! - **M3 Phase A** (shipped): primitives + `Genotype` + enumeration.
//! - **M3 Phase B** (this commit): single-sample data likelihood.
//! - **M3 Phase C** (next): per-site posterior + ML genotype call +
//!   pipeline integration in `fb-cli`. Priors (Ewens sampling,
//!   `Bias.cpp` correction, `Contamination.cpp`) layer on top.

pub mod data_likelihood;
pub mod genotype;
pub mod multinomial;
pub mod sum;

pub use data_likelihood::{
    phred_to_ln, sample_log_likelihood, sample_log_likelihoods, DEFAULT_READ_DEPENDENCE_FACTOR,
};
pub use genotype::{enumerate_genotypes, Genotype, GenotypeElement};
pub use multinomial::{
    factorial_ln, multinomial_coefficient_ln, multinomial_sampling_prob_ln, pow_ln,
    sampling_prob_ln,
};
pub use sum::{log_add, log_sum_exp, sum, sum_i64};
