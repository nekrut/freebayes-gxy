//! Log-domain reduction primitives.
//!
//! Ports the subset of upstream freebayes v1.3.10's `Utility.cpp`
//! log-space arithmetic that the Bayesian genotype model depends on.
//! Upstream splits these across `Sum.h` (linear sum), `Utility.cpp`
//! (`logsumexp`, `logsumexp_probs`) and `Multinomial.cpp` (which calls
//! `sum(obs)` as well). The Rust port keeps them together in one module.
//!
//! ## Precision
//! Upstream uses `long double` — typically 80-bit extended precision on
//! x86 — and a `BigFloat` (ttmath) variant for the `logsumexp_probs`
//! hot path. We standardise on `f64` so the port is platform-portable.
//! For the parity target of 1e-9 on per-site likelihoods this is
//! enough in practice; if a specific site trips the bound a
//! `BigFloat`-equivalent can be added behind a feature flag.

/// Standard linear-space sum over a slice. Mirrors the templated
/// `sum<T>(vector<T>&)` at `src/Sum.h:6-13`.
#[inline]
pub fn sum(xs: &[f64]) -> f64 {
    let mut acc = 0.0f64;
    for &x in xs {
        acc += x;
    }
    acc
}

/// Integer sum variant — used by the multinomial code to total observation
/// counts. Mirrors upstream's `sum(vector<int>&)` template instantiation.
#[inline]
pub fn sum_i64(xs: &[i64]) -> i64 {
    xs.iter().copied().sum()
}

/// Numerically-stable log-sum-exp over a slice of log-space values.
///
/// Returns `log(sum_i exp(xs[i]))` without overflow when the `xs` are
/// large-magnitude. The standard max-shift stabilisation is used:
///
/// ```text
///   logsumexp(x) = M + log(sum_i exp(xs[i] - M))     where M = max(xs)
/// ```
///
/// Semantics match upstream `logsumexp` (`src/Utility.cpp:387-412`); we
/// do not mirror the `BigFloat` extended-precision variant
/// (`logsumexp_probs`, `Utility.cpp:368-384`) — that is a TODO pending
/// a demonstrated precision need.
///
/// **Note on upstream's shift-selection branch.** Upstream's
/// `logsumexp` shifts by `minN` rather than `maxN` when
/// `max |x_i| > max x_i` — a guard for inputs where a very-negative
/// value dominates magnitude. For log-probability inputs (all values
/// `<= 0`) that branch never fires, so the port unconditionally
/// max-shifts. If M3 Phase B or later ever feeds this function mixed-
/// sign values with large-negative magnitude, revisit the shift pick.
///
/// Edge cases:
/// - Empty input returns `f64::NEG_INFINITY` (log of an empty sum == 0).
/// - All inputs `-∞` returns `-∞` (no probability mass).
/// - A single finite input returns itself.
pub fn log_sum_exp(xs: &[f64]) -> f64 {
    if xs.is_empty() {
        return f64::NEG_INFINITY;
    }
    let mut max = f64::NEG_INFINITY;
    for &x in xs {
        if x > max {
            max = x;
        }
    }
    if max == f64::NEG_INFINITY {
        return f64::NEG_INFINITY;
    }
    let mut acc = 0.0f64;
    for &x in xs {
        acc += (x - max).exp();
    }
    max + acc.ln()
}

/// Two-argument log-sum-exp, equivalent to `log_sum_exp(&[a, b])`. Kept
/// separate as a hot-path helper — the Bayesian model calls this inside
/// tight loops where allocating a slice would be wasteful.
#[inline]
pub fn log_add(a: f64, b: f64) -> f64 {
    let (hi, lo) = if a >= b { (a, b) } else { (b, a) };
    if hi == f64::NEG_INFINITY {
        return hi;
    }
    hi + (-((hi - lo).abs())).exp().ln_1p()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sum_empty_is_zero() {
        assert_eq!(sum(&[]), 0.0);
    }

    #[test]
    fn sum_matches_naive() {
        assert!((sum(&[1.0, 2.5, -0.5]) - 3.0).abs() < 1e-15);
    }

    #[test]
    fn log_sum_exp_empty_is_neg_infinity() {
        assert_eq!(log_sum_exp(&[]), f64::NEG_INFINITY);
    }

    #[test]
    fn log_sum_exp_single_returns_itself() {
        assert!((log_sum_exp(&[-3.5]) - (-3.5)).abs() < 1e-15);
    }

    #[test]
    fn log_sum_exp_two_matches_naive() {
        let (a, b) = (-1.5f64, -0.25f64);
        let naive = (a.exp() + b.exp()).ln();
        let got = log_sum_exp(&[a, b]);
        assert!((got - naive).abs() < 1e-12, "got {got}, naive {naive}");
    }

    #[test]
    fn log_sum_exp_large_magnitude_is_stable() {
        // Naive computation would overflow; max-shift keeps it finite.
        let xs = [1000.0f64, 1001.0, 1002.0];
        let got = log_sum_exp(&xs);
        // Analytically: log(e^1000 + e^1001 + e^1002) = 1002 + log(1 + e^-1 + e^-2).
        let expected = 1002.0 + (1.0f64 + (-1.0f64).exp() + (-2.0f64).exp()).ln();
        assert!(
            (got - expected).abs() < 1e-12,
            "got {got}, expected {expected}"
        );
    }

    #[test]
    fn log_sum_exp_all_neg_infinity() {
        assert_eq!(
            log_sum_exp(&[f64::NEG_INFINITY, f64::NEG_INFINITY]),
            f64::NEG_INFINITY
        );
    }

    #[test]
    fn log_add_matches_slice_variant() {
        for (a, b) in [(-1.5, -0.25), (0.0, -10.0), (-10.0, 0.0), (100.0, 99.9)] {
            let got = log_add(a, b);
            let via_slice = log_sum_exp(&[a, b]);
            assert!(
                (got - via_slice).abs() < 1e-12,
                "a={a} b={b} got={got} slice={via_slice}"
            );
        }
    }

    #[test]
    fn log_add_neg_infinity() {
        assert_eq!(
            log_add(f64::NEG_INFINITY, f64::NEG_INFINITY),
            f64::NEG_INFINITY
        );
        assert!((log_add(f64::NEG_INFINITY, 0.0) - 0.0).abs() < 1e-15);
    }
}
