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
//!   single sample.
//! - [`mod@params`] — the [`Parameters`] configuration struct.
//! - [`mod@prior`] — single-sample genotype log-prior.
//! - [`mod@caller`] — single-sample posterior + MAP call + GQ.
//!
//! ## Roadmap
//! - **M3 Phase A** (shipped): primitives + `Genotype` + enumeration.
//! - **M3 Phase B** (shipped): single-sample data likelihood.
//! - **M3 Phase C-1** (this commit): `Parameters` struct, flat /
//!   permutations-weighted prior, single-sample caller with
//!   log-posterior + GQ.
//! - **M3 Phase C-2** (next): fb-cli pipeline integration — pileup →
//!   candidate-allele set → call loop → TSV output. Plus the
//!   deferred prior terms (Ewens, HWE, allele balance), contamination,
//!   and `Bias.cpp` correction.

pub mod caller;
pub mod data_likelihood;
pub mod genotype;
pub mod multinomial;
pub mod params;
pub mod prior;
pub mod sum;

pub use caller::{call_genotype, log_posterior_to_gq, Call, MAX_GQ};
pub use data_likelihood::{
    phred_to_ln, sample_log_likelihood, sample_log_likelihoods, DEFAULT_READ_DEPENDENCE_FACTOR,
};
pub use genotype::{enumerate_genotypes, Genotype, GenotypeElement};
pub use multinomial::{
    factorial_ln, multinomial_coefficient_ln, multinomial_sampling_prob_ln, pow_ln,
    sampling_prob_ln,
};
pub use params::Parameters;
pub use prior::{genotype_log_prior, genotype_log_priors};
pub use sum::{log_add, log_sum_exp, sum, sum_i64};
