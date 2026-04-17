//! Bayesian genotype model for freebayes-gxy.
//!
//! In later milestones this crate hosts the port of freebayes' `Genotype.cpp`,
//! `Multinomial.cpp`, and `Sum.cpp`, including the SIMD log-sum-exp kernel. M0
//! ships only a placeholder so the workspace compiles.

/// A numerically stable log-sum-exp over a pair of log-space values.
///
/// This is the minimum viable primitive we need in the genotype model; it is
/// exercised by the unit test below so the crate contributes to `cargo test`
/// from M0 onward. Later milestones replace callers of this with a vectorised
/// implementation.
pub fn logsumexp2(a: f64, b: f64) -> f64 {
    let (hi, lo) = if a >= b { (a, b) } else { (b, a) };
    if hi == f64::NEG_INFINITY {
        return hi;
    }
    hi + (1.0 + (lo - hi).exp()).ln()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn logsumexp2_matches_naive_for_normal_values() {
        let a: f64 = -1.5;
        let b: f64 = -0.25;
        let naive = (a.exp() + b.exp()).ln();
        let got = logsumexp2(a, b);
        assert!((got - naive).abs() < 1e-12, "got {got}, naive {naive}");
    }

    #[test]
    fn logsumexp2_handles_neg_infinity() {
        assert_eq!(
            logsumexp2(f64::NEG_INFINITY, f64::NEG_INFINITY),
            f64::NEG_INFINITY
        );
        assert!((logsumexp2(f64::NEG_INFINITY, 0.0) - 0.0).abs() < 1e-12);
    }
}
