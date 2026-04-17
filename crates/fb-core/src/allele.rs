//! Per-position allele types and read-level observations.
//!
//! Port of upstream freebayes v1.3.10's internal allele model. The canonical
//! references are:
//! - `src/Allele.h:60-69`          — the `AlleleType` enum.
//! - `src/Allele.h:104-142`        — the `Allele` class field inventory.
//! - `src/AlleleParser.cpp:1329..` — `registerAlignment`, which walks the
//!   CIGAR and emits the observation records we mirror here.
//!
//! **Coordinate conventions.** All positions in this module are **0-based,
//! half-open** against the reference, matching upstream's internal layout.
//! VCF emission (M4) converts to 1-based at the boundary.

use std::hash::{Hash, Hasher};

/// Allele kind. Mirrors upstream's `enum AlleleType` (`src/Allele.h:60`).
///
/// The `Genotype` variant from upstream is omitted deliberately — it is an
/// abstract tag used by the Bayesian genotype model (M3), not produced by the
/// pileup walker. The numeric values match upstream's bitmask so that future
/// ported code performing bitset filtering stays semantically equivalent.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum AlleleKind {
    /// Matches the reference base(s) — `ALLELE_REFERENCE = 2`.
    Reference = 2,
    /// Multi-nucleotide polymorphism — `ALLELE_MNP = 4`. Emitted only after
    /// M2's collapsing pass; M1's pileup walker never emits this directly.
    Mnp = 4,
    /// Single-nucleotide polymorphism — `ALLELE_SNP = 8`.
    Snp = 8,
    /// Inserted bases in the read relative to the reference — `ALLELE_INSERTION = 16`.
    Insertion = 16,
    /// Deleted bases in the read relative to the reference — `ALLELE_DELETION = 32`.
    Deletion = 32,
    /// Complex event (MNP + indel in close proximity) — `ALLELE_COMPLEX = 64`.
    /// M2's job to synthesise; M1 never emits this directly.
    Complex = 64,
    /// Non-informative observation (e.g. N base, soft-clip) — `ALLELE_NULL = 128`.
    Null = 128,
}

impl AlleleKind {
    /// One-letter tag used for TSV / debug output.
    pub fn as_tag(self) -> &'static str {
        match self {
            AlleleKind::Reference => "REF",
            AlleleKind::Mnp => "MNP",
            AlleleKind::Snp => "SNP",
            AlleleKind::Insertion => "INS",
            AlleleKind::Deletion => "DEL",
            AlleleKind::Complex => "CPX",
            AlleleKind::Null => "NUL",
        }
    }
}

/// Which strand the observation came from.
///
/// Mirrors upstream `enum AlleleStrand { STRAND_FORWARD, STRAND_REVERSE }`
/// (`src/Allele.h:74-77`). In BAM terms, this is `!record.is_reverse()`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Strand {
    Forward,
    Reverse,
}

/// A distinct allele at a reference position.
///
/// We store the reference and alternate sequences verbatim (matching upstream
/// `referenceSequence` / `alternateSequence`, `src/Allele.h:108-109`). `length`
/// is the upstream `Allele::length` field (`src/Allele.h:114`): deletion
/// implies 0 inserted-length? — upstream actually stores the *event* length
/// there (deletion: deleted-base count; insertion: inserted-base count; SNP:
/// 1; reference: run length). We follow that convention.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Allele {
    pub kind: AlleleKind,
    /// Reference bases this allele spans on the reference (0-based).
    pub ref_seq: Vec<u8>,
    /// Alternate bases as they appear in the read. Empty for deletions.
    pub alt_seq: Vec<u8>,
    /// 0-based start of the allele on the reference (half-open interval
    /// `[position, position + ref_seq.len())`). Matches upstream
    /// `Allele::position` (`src/Allele.h:111`).
    pub position: i64,
    /// Upstream's `length` field — event length in the read, not reference
    /// length. For a deletion it is the deleted-base count (same as
    /// `ref_seq.len()`); for a reference run it is the run length. Stored
    /// explicitly so consumers don't need to branch on `kind`.
    pub length: usize,
}

impl Allele {
    /// A reference-matching run `[position, position + length)`.
    pub fn reference(position: i64, ref_seq: Vec<u8>) -> Self {
        let length = ref_seq.len();
        Self {
            kind: AlleleKind::Reference,
            alt_seq: ref_seq.clone(),
            ref_seq,
            position,
            length,
        }
    }

    /// A single-base substitution at `position`. Matches upstream SNP
    /// emission in `AlleleParser.cpp:1516..1530`.
    pub fn snp(position: i64, ref_base: u8, alt_base: u8) -> Self {
        Self {
            kind: AlleleKind::Snp,
            ref_seq: vec![ref_base],
            alt_seq: vec![alt_base],
            position,
            length: 1,
        }
    }

    /// An insertion of `alt_seq` anchored at `position`. Upstream keeps
    /// insertions at the reference coordinate immediately preceding the
    /// inserted bases (`AlleleParser.cpp:1755..1769`); `ref_seq` is left
    /// empty at the internal layer and synthesised (with an anchor base)
    /// only at VCF emission time. We keep the same convention for byte
    /// parity with upstream.
    pub fn insertion(position: i64, alt_seq: Vec<u8>) -> Self {
        let length = alt_seq.len();
        Self {
            kind: AlleleKind::Insertion,
            ref_seq: Vec::new(),
            alt_seq,
            position,
            length,
        }
    }

    /// A deletion of `ref_seq` starting at `position`. Upstream leaves
    /// `alternateSequence` empty (`AlleleParser.cpp:1684-1698`, passes
    /// `nullstr`); the VCF anchor base is synthesised at emission time.
    pub fn deletion(position: i64, ref_seq: Vec<u8>) -> Self {
        let length = ref_seq.len();
        Self {
            kind: AlleleKind::Deletion,
            ref_seq,
            alt_seq: Vec::new(),
            position,
            length,
        }
    }

    /// A null / uninformative observation (soft clip or N base). See
    /// `AlleleParser.cpp:1532-1548` (N base in a match) and `1783-1797`
    /// (soft clip).
    pub fn null(position: i64, length: usize) -> Self {
        Self {
            kind: AlleleKind::Null,
            ref_seq: Vec::new(),
            alt_seq: Vec::new(),
            position,
            length,
        }
    }

    /// Fold `other` into `self`, producing an [`AlleleKind::Complex`]
    /// composite allele. Mirrors upstream `Allele::mergeAllele`
    /// (`src/Allele.cpp:1454-1470`).
    ///
    /// Upstream concatenates `alternateSequence` and sums `length`; we do
    /// the same and additionally concatenate `ref_seq` so the resulting
    /// allele carries the full reference span. That matters because our
    /// internal model keeps `ref_seq` empty for insertions and `alt_seq`
    /// empty for deletions; the composite must express both sides of the
    /// event so downstream consumers (M4's VCF anchor synthesis, M3's
    /// likelihood model) see a complete REF/ALT pair.
    ///
    /// Position invariant: `self.position` is preserved — upstream
    /// anchors the composite at the leftmost event, and
    /// [`haplotype::clump_observations`](crate::haplotype::clump_observations)
    /// builds composites by cloning the leftmost observation and folding
    /// later ones in.
    ///
    /// Visibility is deliberately `pub(crate)`: this is an internal
    /// helper for the clumping pipeline only; constructing Complex
    /// alleles outside that context would skip upstream's bookkeeping
    /// (CIGAR synthesis, flanking-base maintenance, left-alignment)
    /// that M4 will need.
    ///
    /// **TODO(M4):** this does not synthesise a merged CIGAR
    /// (upstream `mergeCigar`, `Allele.cpp:1468`) nor left-align the
    /// resulting indel composite (`LeftAlign.cpp`). Both are required
    /// before VCF emission.
    pub(crate) fn merge_with(&mut self, other: &Allele) {
        self.kind = AlleleKind::Complex;
        self.ref_seq.extend_from_slice(&other.ref_seq);
        self.alt_seq.extend_from_slice(&other.alt_seq);
        self.length += other.length;
    }
}

// Manual `Hash` so that (kind, position, ref, alt) uniquely key an allele —
// `length` is derivable from those and would otherwise perturb the hash.
impl Hash for Allele {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.kind.hash(state);
        self.position.hash(state);
        self.ref_seq.hash(state);
        self.alt_seq.hash(state);
    }
}

/// One read's observation of one allele.
///
/// Mirrors the per-allele fields upstream stores on `Allele` itself
/// (`src/Allele.h:120-138`). We split them out of `Allele` because the
/// variant caller aggregates many observations of the same `Allele`.
#[derive(Debug, Clone)]
pub struct AlleleObservation {
    pub allele: Allele,
    pub read_name: String,
    pub mapq: u8,
    /// Sum of base qualities supporting this observation.
    ///
    /// For SNPs and single-base reference observations this is just the
    /// underlying base Q. For insertions it is the sum of inserted-base Qs.
    /// For deletions upstream averages flanking base Qs
    /// (`AlleleParser.cpp:1658-1674`); we store the (rounded) average here
    /// so that aggregation across observations stays additive.
    pub base_quality_sum: u32,
    pub strand: Strand,
    /// 0-based offset into the read (query sequence) of the first base of
    /// this allele event. For deletions this is the read position *between*
    /// the flanking bases (same convention upstream uses at
    /// `AlleleParser.cpp:1691` — `rp` before the deletion consumes it).
    pub read_position: usize,
    pub is_proper_pair: bool,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn allele_kind_discriminants_match_upstream() {
        // Values taken verbatim from src/Allele.h:60-69.
        assert_eq!(AlleleKind::Reference as u8, 2);
        assert_eq!(AlleleKind::Mnp as u8, 4);
        assert_eq!(AlleleKind::Snp as u8, 8);
        assert_eq!(AlleleKind::Insertion as u8, 16);
        assert_eq!(AlleleKind::Deletion as u8, 32);
        assert_eq!(AlleleKind::Complex as u8, 64);
        assert_eq!(AlleleKind::Null as u8, 128);
    }

    #[test]
    fn snp_constructor_sets_fields() {
        let a = Allele::snp(42, b'A', b'G');
        assert_eq!(a.kind, AlleleKind::Snp);
        assert_eq!(a.position, 42);
        assert_eq!(a.ref_seq, b"A");
        assert_eq!(a.alt_seq, b"G");
        assert_eq!(a.length, 1);
    }

    #[test]
    fn reference_run_preserves_length() {
        let a = Allele::reference(100, b"ACGT".to_vec());
        assert_eq!(a.kind, AlleleKind::Reference);
        assert_eq!(a.length, 4);
    }

    #[test]
    fn insertion_has_empty_ref_seq() {
        let a = Allele::insertion(5, b"CC".to_vec());
        assert!(a.ref_seq.is_empty());
        assert_eq!(a.alt_seq, b"CC");
        assert_eq!(a.length, 2);
    }

    #[test]
    fn deletion_has_empty_alt_seq() {
        let a = Allele::deletion(5, b"AAA".to_vec());
        assert_eq!(a.ref_seq, b"AAA");
        assert!(a.alt_seq.is_empty());
        assert_eq!(a.length, 3);
    }

    #[test]
    fn equality_ignores_derivable_length() {
        let a = Allele::snp(7, b'C', b'T');
        let b = Allele::snp(7, b'C', b'T');
        assert_eq!(a, b);
        let c = Allele::snp(7, b'C', b'A');
        assert_ne!(a, c);
    }
}
