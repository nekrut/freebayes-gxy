//! Core primitives for freebayes-gxy.
//!
//! M1 hosts the allele observation model and the CIGAR walker that turns one
//! aligned read into per-position [`AlleleObservation`]s. Haplotype
//! construction (M2), the Bayesian genotype model (M3) and downstream crates
//! depend on the types here.

pub mod allele;
pub mod pileup;

pub use allele::{Allele, AlleleKind, AlleleObservation, Strand};
pub use pileup::{walk_alignment, walk_record, AlignmentView, ReadFilter};

/// Milestone marker for the current build.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Milestone {
    pub name: &'static str,
}

impl Milestone {
    /// Current milestone — bumped as work lands.
    pub const CURRENT: Self = Self { name: "M1" };
}

/// Returns the current milestone tag.
pub fn current_milestone() -> Milestone {
    Milestone::CURRENT
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn current_milestone_is_m1() {
        assert_eq!(current_milestone().name, "M1");
    }
}
