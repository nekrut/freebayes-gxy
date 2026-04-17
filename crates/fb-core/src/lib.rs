//! Core primitives for freebayes-gxy.
//!
//! This crate will host the pileup engine, allele observation model, haplotype
//! construction, and the building blocks of the Bayesian likelihood pipeline
//! (see `PLAN.md` §4). M0 ships only a placeholder so the workspace compiles
//! and `cargo test` has something to run.

/// Milestone marker for the current build.
///
/// M1 replaces this with real allele/pileup types. The struct is deliberately
/// tiny so downstream crates can depend on `fb-core` from day one.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Milestone {
    pub name: &'static str,
}

impl Milestone {
    /// Current milestone for the workspace scaffolding.
    pub const CURRENT: Self = Self { name: "M0" };
}

/// Returns the current milestone tag.
pub fn current_milestone() -> Milestone {
    Milestone::CURRENT
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn current_milestone_is_m0() {
        assert_eq!(current_milestone().name, "M0");
    }
}
