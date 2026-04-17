//! Multinomial probability kernels.
//!
//! Port of upstream freebayes v1.3.10's `Multinomial.cpp` plus the
//! supporting `factorialln` / `powln` helpers from `Utility.cpp`.
//!
//! The Bayesian genotype likelihood calls `multinomial_sampling_prob_ln`
//! to score a vector of allele observation counts against a vector of
//! expected sampling probabilities derived from the genotype. The
//! coefficient helper is exposed because upstream caches it on each
//! `Genotype` instance (`Genotype::permutationsln`,
//! `src/Genotype.h:53`).

use crate::sum::{sum, sum_i64};

/// `log(n!)` via `lgamma(n+1)`. Mirrors upstream `factorialln` at
/// `src/Utility.cpp:302-319`, which delegates to a cache + `lgamma`. We
/// skip the caching layer — `libm::lgamma` is already cheap enough that
/// the cache isn't a measurable win at realistic input sizes, and
/// removing it keeps the primitive stateless (important for threading,
/// which arrives in M5).
///
/// Negative `n` is a caller error; we return `NaN` to surface it loudly.
/// `n == 0` returns `0.0` (upstream: `impl_factorialln` line 313-315).
#[inline]
pub fn factorial_ln(n: i64) -> f64 {
    if n < 0 {
        return f64::NAN;
    }
    if n == 0 {
        return 0.0;
    }
    libm::lgamma((n + 1) as f64)
}

/// `n * ln_x` in log-space. Mirrors upstream `powln(ln_m, n)` at
/// `src/Utility.cpp:69-71`. Splitting out the helper lets readers track
/// calls 1:1 against the upstream multinomial code.
#[inline]
pub fn pow_ln(ln_x: f64, n: i64) -> f64 {
    ln_x * (n as f64)
}

/// Log multinomial coefficient — the `ln(n! / prod(k_i!))` term.
/// Mirrors `multinomialCoefficientLn(int n, const vector<int>& counts)`
/// at `src/Multinomial.cpp:34-39`.
///
/// `counts` must be non-negative and (by convention) sum to `n`;
/// upstream does not enforce the constraint and we follow suit —
/// `n != sum(counts)` simply yields the corresponding generalised
/// coefficient, which the genotype model never calls for.
pub fn multinomial_coefficient_ln(n: i64, counts: &[i64]) -> f64 {
    let ln_count_factorials: Vec<f64> = counts.iter().map(|&c| factorial_ln(c)).collect();
    factorial_ln(n) - sum(&ln_count_factorials)
}

/// Log multinomial sampling probability, i.e.
/// `log(n! / prod(k_i!) * prod(p_i^k_i))` for `k = obs`, `n = sum(obs)`.
/// Mirrors `multinomialSamplingProbLn(const vector<long double>& probs,
/// const vector<int>& obs)` at `src/Multinomial.cpp:21-32`.
///
/// `probs` are linear-space probabilities (we take `ln` here, matching
/// upstream's `powln(log(*p), *o)` call). Zero-probability events with
/// a nonzero observation count produce `-∞`, and events with zero
/// observations contribute `0` to the sum regardless of `prob[i]`.
///
/// **Intentional divergence from upstream:** when `obs[i] == 0` the port
/// short-circuits to `0.0` without evaluating `ln(prob[i])`. Upstream's
/// `powln(log(0.0), 0)` evaluates to `NaN` (`-inf * 0`), but the sole
/// caller (`DataLikelihood.cpp` via `alleleProbabilities()`) never hits
/// that input because its `probs` vector always has `prob[i] > 0` for
/// every allele in the current genotype. The short-circuit is
/// parity-equivalent on all physically reachable inputs and avoids the
/// NaN booby-trap on synthetic test data.
///
/// Panics in debug mode if `probs.len() != obs.len()`.
pub fn multinomial_sampling_prob_ln(probs: &[f64], obs: &[i64]) -> f64 {
    debug_assert_eq!(
        probs.len(),
        obs.len(),
        "probs and obs must have matching lengths"
    );
    let n = sum_i64(obs);

    let ln_count_factorials: Vec<f64> = obs.iter().map(|&c| factorial_ln(c)).collect();
    let pow_terms: Vec<f64> = probs
        .iter()
        .zip(obs.iter())
        .map(|(&p, &k)| {
            if k == 0 {
                // 0 * ln(0) == 0 by convention — avoid the -inf.
                0.0
            } else {
                pow_ln(p.ln(), k)
            }
        })
        .collect();

    factorial_ln(n) - sum(&ln_count_factorials) + sum(&pow_terms)
}

/// `sum_i obs[i] * ln(probs[i])` — just the sampling-probability core
/// without the multinomial coefficient. Mirrors `samplingProbLn` at
/// `src/Multinomial.cpp:41-49`.
pub fn sampling_prob_ln(probs: &[f64], obs: &[i64]) -> f64 {
    debug_assert_eq!(probs.len(), obs.len());
    let mut acc = 0.0f64;
    for (&p, &k) in probs.iter().zip(obs.iter()) {
        if k == 0 {
            continue;
        }
        acc += pow_ln(p.ln(), k);
    }
    acc
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn factorial_ln_small_values() {
        assert_eq!(factorial_ln(0), 0.0);
        assert!((factorial_ln(1) - 0.0).abs() < 1e-15);
        assert!((factorial_ln(5) - (120.0f64).ln()).abs() < 1e-12);
        assert!((factorial_ln(10) - (3_628_800.0f64).ln()).abs() < 1e-10);
    }

    #[test]
    fn factorial_ln_rejects_negative() {
        assert!(factorial_ln(-1).is_nan());
    }

    #[test]
    fn pow_ln_is_multiplication() {
        assert_eq!(pow_ln(2.0, 3), 6.0);
        assert_eq!(pow_ln(-1.5, 4), -6.0);
    }

    #[test]
    fn multinomial_coefficient_trivial() {
        // n=3 distributed as (3,0,0) — coefficient = 3!/3!0!0! = 1.
        let got = multinomial_coefficient_ln(3, &[3, 0, 0]);
        assert!(got.abs() < 1e-12, "got {got}");
    }

    #[test]
    fn multinomial_coefficient_diploid_het() {
        // Diploid heterozygous call: n=2, counts=[1,1] → 2!/1!1! = 2.
        let got = multinomial_coefficient_ln(2, &[1, 1]);
        assert!((got - 2.0f64.ln()).abs() < 1e-12, "got {got}");
    }

    #[test]
    fn multinomial_coefficient_matches_closed_form() {
        // n=10, counts=[5,3,2] → 10!/(5!3!2!) = 2520.
        let got = multinomial_coefficient_ln(10, &[5, 3, 2]);
        assert!((got - 2520.0f64.ln()).abs() < 1e-10, "got {got}");
    }

    #[test]
    fn multinomial_sampling_prob_fair_coin() {
        // Tossing a fair coin 10 times, expecting 5 heads and 5 tails.
        //   C(10,5) = 252
        //   p^5 * (1-p)^5 with p = 0.5 → (0.5)^10 = 1/1024
        // Sampling prob = 252/1024 = 0.24609375
        let got = multinomial_sampling_prob_ln(&[0.5, 0.5], &[5, 5]);
        assert!(
            (got.exp() - 0.24609375).abs() < 1e-12,
            "got exp={} expected 0.24609375",
            got.exp()
        );
    }

    #[test]
    fn multinomial_sampling_prob_all_one_category() {
        // All observations land in category 0 with prob 1.0.
        //   n=3, obs=[3,0,0], probs=[1.0, 0.0, 0.0] → prob = 1.0.
        let got = multinomial_sampling_prob_ln(&[1.0, 0.0, 0.0], &[3, 0, 0]);
        assert!(got.abs() < 1e-12, "got {got}");
    }

    #[test]
    fn multinomial_sampling_prob_zero_prob_with_observations_is_neg_inf() {
        // Observing an event with probability 0 → -inf log-likelihood.
        let got = multinomial_sampling_prob_ln(&[0.5, 0.0], &[1, 1]);
        assert_eq!(got, f64::NEG_INFINITY);
    }

    #[test]
    fn sampling_prob_ln_skips_zero_counts() {
        // obs=[0,0,0] → 0, regardless of probs.
        let got = sampling_prob_ln(&[0.1, 0.2, 0.7], &[0, 0, 0]);
        assert_eq!(got, 0.0);
    }

    #[test]
    fn sampling_prob_ln_mixed_zero_and_nonzero_counts() {
        // obs=[3, 0, 2], probs=[0.5, 0.25, 0.25]
        // Expected: 3*ln(0.5) + 0 + 2*ln(0.25) = -3*ln(2) + 2*(-2*ln(2))
        //         = -3*ln(2) - 4*ln(2) = -7*ln(2) ≈ -4.852030...
        let got = sampling_prob_ln(&[0.5, 0.25, 0.25], &[3, 0, 2]);
        let expected = -7.0 * 2.0f64.ln();
        assert!(
            (got - expected).abs() < 1e-12,
            "got {got}, expected {expected}"
        );
    }
}
