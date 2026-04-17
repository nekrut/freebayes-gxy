//! Genotype type and enumeration.
//!
//! Port of the subset of upstream freebayes v1.3.10's `Genotype.h` /
//! `Genotype.cpp` that the Bayesian model needs at this milestone:
//! a `Genotype` data type (a multiset of alleles with associated counts)
//! and the `all_possible_genotypes` enumeration used to build candidate
//! sets at each variant site.
//!
//! Upstream references:
//! - `src/Genotype.h:29-38` — `GenotypeElement` (allele + count).
//! - `src/Genotype.h:41-110` — `Genotype` class.
//! - `src/Genotype.cpp:354-361` — `allPossibleGenotypes(ploidy, alleles)`,
//!   which delegates to `multichoose` from `vcflib/multichoose.h`.
//!
//! ## M3 Phase A scope (this module)
//! - [`Genotype`] struct with deterministic ordering, homozygosity probe,
//!   and cached `permutations_ln` (upstream `permutationsln`).
//! - [`enumerate_genotypes`] — the `multichoose(ploidy, alleles)` port,
//!   producing every distinct ploidy-sized multiset of the input alleles.
//!
//! Deferred to M3 Phases B/C (documented at call sites):
//! - Data-likelihood per genotype (`DataLikelihood.cpp`).
//! - Prior (`Genotype::priorProbability`).
//! - `GenotypeCombo` (cross-sample genotype joint).
//! - `alleleSamplingProb` with observation bias (`Bias.cpp`).

use std::cmp::Ordering;

use fb_core::Allele;

use crate::multinomial::multinomial_coefficient_ln;

/// One "element" of a genotype — an allele and the number of times it
/// appears in the multiset. Mirrors upstream `GenotypeElement`
/// (`src/Genotype.h:29-38`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GenotypeElement {
    pub allele: Allele,
    pub count: u32,
}

// `f64` has no `Eq`, so `Genotype` derives only `PartialEq`. That is fine
// for the Bayesian model; a canonical `Eq` would require an exact-bit
// equality on `permutations_ln` which we don't want.

/// A genotype — a size-`ploidy` multiset of alleles. Stored as a
/// canonical sorted list of [`GenotypeElement`]s so equality and hashing
/// are order-independent.
///
/// Mirrors upstream `Genotype` (`src/Genotype.h:41`), which inherits
/// from `vector<GenotypeElement>` and carries auxiliary fields like
/// `ploidy`, `homozygous` and `permutationsln` (the
/// `multinomialCoefficientLn(ploidy, counts)` cache at `src/Genotype.h:53`).
#[derive(Debug, Clone, PartialEq)]
pub struct Genotype {
    /// Sorted, deduplicated element list.
    pub elements: Vec<GenotypeElement>,
    /// Ploidy — sum of counts across elements.
    pub ploidy: u32,
    /// `true` when a single element accounts for the whole genotype.
    pub homozygous: bool,
    /// Cached `multinomialCoefficientLn(ploidy, counts)`. Zero when the
    /// genotype is homozygous — upstream leaves it at 0 in that case
    /// (`src/Genotype.h:65-69`).
    pub permutations_ln: f64,
}

impl Genotype {
    /// Build a canonical [`Genotype`] from a raw vector of alleles.
    /// Upstream does the same at `src/Genotype.h:55-71`: sort, group by
    /// base, record counts, stash `permutationsln`.
    ///
    /// `alleles` is consumed to emphasise that the input is semantically
    /// a multiset — the caller loses any notion of insertion order.
    pub fn from_alleles(mut alleles: Vec<Allele>) -> Self {
        // Sort by the allele's natural ordering so adjacent runs can be
        // collapsed into GenotypeElements. `Allele` doesn't derive `Ord`
        // today, so we sort by (position, kind, ref, alt) explicitly.
        alleles.sort_by(allele_cmp);
        let ploidy = alleles.len() as u32;
        let mut elements: Vec<GenotypeElement> = Vec::new();
        for a in alleles.into_iter() {
            if let Some(last) = elements.last_mut() {
                if allele_cmp(&last.allele, &a) == Ordering::Equal {
                    last.count += 1;
                    continue;
                }
            }
            elements.push(GenotypeElement {
                allele: a,
                count: 1,
            });
        }
        let homozygous = elements.len() == 1;
        let permutations_ln = if homozygous {
            0.0
        } else {
            let counts: Vec<i64> = elements.iter().map(|e| e.count as i64).collect();
            multinomial_coefficient_ln(ploidy as i64, &counts)
        };
        Self {
            elements,
            ploidy,
            homozygous,
            permutations_ln,
        }
    }

    /// Sampling probability of each allele under the genotype, i.e.
    /// `count_i / ploidy`. Mirrors upstream `alleleProbabilities`
    /// (`src/Genotype.h:85`) without the `Bias` correction — that lands
    /// with M3 Phase B when the data-likelihood model comes online.
    pub fn allele_probabilities(&self) -> Vec<f64> {
        let ploidy = self.ploidy as f64;
        self.elements
            .iter()
            .map(|e| e.count as f64 / ploidy)
            .collect()
    }

    /// Total number of distinct alleles in the genotype (e.g. 1 for
    /// homozygous, 2 for heterozygous).
    #[inline]
    pub fn unique_allele_count(&self) -> usize {
        self.elements.len()
    }

    /// Count of the specific allele within the genotype, 0 if absent.
    pub fn allele_count(&self, allele: &Allele) -> u32 {
        for e in &self.elements {
            if allele_cmp(&e.allele, allele) == Ordering::Equal {
                return e.count;
            }
        }
        0
    }

    /// Human-readable summary — the allele `kind` tags separated by `/`,
    /// with homozygous genotypes collapsed (e.g. `REF/REF` → `REF`).
    pub fn str_tag(&self) -> String {
        if self.homozygous {
            return self.elements[0].allele.kind.as_tag().to_string();
        }
        let mut parts: Vec<&str> = Vec::new();
        for e in &self.elements {
            for _ in 0..e.count {
                parts.push(e.allele.kind.as_tag());
            }
        }
        parts.join("/")
    }
}

/// Enumerate every ploidy-sized multiset of `alleles`, returning one
/// canonical [`Genotype`] per distinct multiset. Mirrors upstream
/// `allPossibleGenotypes` (`src/Genotype.cpp:354-361`), which calls
/// `multichoose(ploidy, alleles)` from `vcflib/multichoose.h`.
///
/// The combinatorial count is `C(n + k - 1, k)` where `n = alleles.len()`
/// and `k = ploidy`. For diploid (`ploidy=2`) the caller gets
/// `n(n+1)/2` genotypes.
///
/// Duplicate input alleles are deduplicated before enumeration — the
/// caller is typically handing in a site's candidate-allele set which
/// may have been produced by aggregating per-read observations.
pub fn enumerate_genotypes(alleles: &[Allele], ploidy: u32) -> Vec<Genotype> {
    // Deduplicate the input while preserving position/kind ordering.
    let mut unique: Vec<Allele> = alleles.to_vec();
    unique.sort_by(allele_cmp);
    unique.dedup_by(|a, b| allele_cmp(a, b) == Ordering::Equal);

    if ploidy == 0 || unique.is_empty() {
        return Vec::new();
    }

    // Iterative multichoose: generate index multisets of size `ploidy`
    // over the range 0..unique.len(), in lexicographic order. Each
    // multiset becomes one Genotype.
    let n = unique.len();
    let k = ploidy as usize;
    let mut out = Vec::new();
    let mut idx = vec![0usize; k];
    loop {
        let picks: Vec<Allele> = idx.iter().map(|&i| unique[i].clone()).collect();
        out.push(Genotype::from_alleles(picks));

        // Increment the rightmost index, carrying when it rolls over.
        // Multichoose keeps the indices monotonically non-decreasing so
        // that `(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)` is the
        // sequence for n=3, k=2.
        let mut j = k;
        while j > 0 {
            j -= 1;
            idx[j] += 1;
            if idx[j] < n {
                // Reset everything to the right of `j` to match, keeping
                // the non-decreasing invariant.
                for m in j + 1..k {
                    idx[m] = idx[j];
                }
                break;
            }
            if j == 0 {
                return out;
            }
        }
    }
}

/// Total-order comparison of two alleles used for [`Genotype`]
/// canonicalisation. Sorts by (position, kind, ref_seq, alt_seq) so
/// equivalent alleles collide regardless of insertion order.
fn allele_cmp(a: &Allele, b: &Allele) -> Ordering {
    a.position
        .cmp(&b.position)
        .then_with(|| (a.kind as u8).cmp(&(b.kind as u8)))
        .then_with(|| a.ref_seq.cmp(&b.ref_seq))
        .then_with(|| a.alt_seq.cmp(&b.alt_seq))
}

#[cfg(test)]
mod tests {
    use super::*;
    use fb_core::AlleleKind;

    fn ref_a(pos: i64) -> Allele {
        Allele::reference(pos, vec![b'A'])
    }
    fn snp_ag(pos: i64) -> Allele {
        Allele::snp(pos, b'A', b'G')
    }
    fn snp_at(pos: i64) -> Allele {
        Allele::snp(pos, b'A', b'T')
    }

    #[test]
    fn homozygous_diploid_reference() {
        let gt = Genotype::from_alleles(vec![ref_a(100), ref_a(100)]);
        assert_eq!(gt.ploidy, 2);
        assert_eq!(gt.unique_allele_count(), 1);
        assert!(gt.homozygous);
        assert_eq!(gt.permutations_ln, 0.0);
        assert_eq!(gt.elements[0].count, 2);
    }

    #[test]
    fn heterozygous_diploid_snp() {
        let gt = Genotype::from_alleles(vec![ref_a(100), snp_ag(100)]);
        assert_eq!(gt.ploidy, 2);
        assert_eq!(gt.unique_allele_count(), 2);
        assert!(!gt.homozygous);
        // 2!/1!1! = 2.
        assert!((gt.permutations_ln - 2.0f64.ln()).abs() < 1e-12);
        // Allele probabilities = [0.5, 0.5].
        let probs = gt.allele_probabilities();
        assert_eq!(probs.len(), 2);
        assert!((probs[0] - 0.5).abs() < 1e-12);
        assert!((probs[1] - 0.5).abs() < 1e-12);
    }

    #[test]
    fn triploid_with_three_distinct_alleles() {
        let gt = Genotype::from_alleles(vec![ref_a(100), snp_ag(100), snp_at(100)]);
        assert_eq!(gt.ploidy, 3);
        assert_eq!(gt.unique_allele_count(), 3);
        // 3!/1!1!1! = 6.
        assert!((gt.permutations_ln - 6.0f64.ln()).abs() < 1e-12);
    }

    #[test]
    fn from_alleles_collapses_duplicate_input() {
        // Two distinct SNPs AND a dup SNP → should group into {SNP:2, REF:1}
        let gt = Genotype::from_alleles(vec![snp_ag(100), ref_a(100), snp_ag(100)]);
        assert_eq!(gt.ploidy, 3);
        assert_eq!(gt.unique_allele_count(), 2);
        assert_eq!(gt.allele_count(&snp_ag(100)), 2);
        assert_eq!(gt.allele_count(&ref_a(100)), 1);
        // 3!/(2!1!) = 3.
        assert!((gt.permutations_ln - 3.0f64.ln()).abs() < 1e-12);
    }

    #[test]
    fn str_tag_is_slash_separated() {
        let gt = Genotype::from_alleles(vec![ref_a(100), snp_ag(100)]);
        let tag = gt.str_tag();
        // REF and SNP are present; the tag reflects canonical order.
        assert!(tag.contains('/'));
        assert!(tag.contains("REF"));
        assert!(tag.contains("SNP"));
    }

    #[test]
    fn str_tag_homozygous_collapses() {
        let gt = Genotype::from_alleles(vec![ref_a(100), ref_a(100)]);
        assert_eq!(gt.str_tag(), "REF");
    }

    #[test]
    fn enumerate_genotypes_diploid_two_alleles() {
        let alleles = vec![ref_a(100), snp_ag(100)];
        let gts = enumerate_genotypes(&alleles, 2);
        // C(2+2-1, 2) = 3: REF/REF, REF/SNP, SNP/SNP.
        assert_eq!(gts.len(), 3);
        let tags: Vec<_> = gts.iter().map(|g| g.str_tag()).collect();
        assert!(tags.iter().any(|t| t == "REF"));
        assert!(tags.iter().any(|t| t == "SNP"));
        assert!(tags.iter().any(|t| t == "REF/SNP" || t == "SNP/REF"));
    }

    #[test]
    fn enumerate_genotypes_diploid_three_alleles() {
        let alleles = vec![ref_a(100), snp_ag(100), snp_at(100)];
        let gts = enumerate_genotypes(&alleles, 2);
        // C(3+2-1, 2) = 6 distinct diploid genotypes.
        assert_eq!(gts.len(), 6);
    }

    #[test]
    fn enumerate_genotypes_triploid_three_alleles() {
        let alleles = vec![ref_a(100), snp_ag(100), snp_at(100)];
        let gts = enumerate_genotypes(&alleles, 3);
        // C(3+3-1, 3) = 10.
        assert_eq!(gts.len(), 10);
        // Each genotype must have ploidy == 3.
        for g in &gts {
            assert_eq!(g.ploidy, 3);
        }
    }

    #[test]
    fn enumerate_genotypes_deduplicates_input() {
        let alleles = vec![ref_a(100), ref_a(100), snp_ag(100)];
        let gts = enumerate_genotypes(&alleles, 2);
        // Effective unique alleles = 2, so 3 genotypes.
        assert_eq!(gts.len(), 3);
    }

    #[test]
    fn enumerate_genotypes_single_allele_homozygous_only() {
        let alleles = vec![ref_a(100)];
        let gts = enumerate_genotypes(&alleles, 2);
        assert_eq!(gts.len(), 1);
        assert!(gts[0].homozygous);
    }

    #[test]
    fn enumerate_genotypes_ploidy_zero_or_empty_alleles_is_empty() {
        assert!(enumerate_genotypes(&[ref_a(100)], 0).is_empty());
        assert!(enumerate_genotypes(&[], 2).is_empty());
    }

    #[test]
    fn enumerate_genotypes_canonical_ordering_is_stable() {
        // Two calls with the same input must yield identical vectors.
        let alleles = vec![snp_ag(100), ref_a(100), snp_at(100)];
        let a = enumerate_genotypes(&alleles, 2);
        let b = enumerate_genotypes(&alleles, 2);
        assert_eq!(a.len(), b.len());
        for (ga, gb) in a.iter().zip(b.iter()) {
            assert_eq!(ga, gb);
        }
    }

    #[test]
    fn allele_count_returns_zero_for_absent() {
        let gt = Genotype::from_alleles(vec![ref_a(100), ref_a(100)]);
        assert_eq!(gt.allele_count(&snp_ag(100)), 0);
    }

    #[test]
    fn tetraploid_heterozygous_permutations() {
        // Tetraploid with 2 ref + 2 snp → 4!/(2!2!) = 6.
        let gt = Genotype::from_alleles(vec![ref_a(100), ref_a(100), snp_ag(100), snp_ag(100)]);
        assert_eq!(gt.ploidy, 4);
        assert!(!gt.homozygous);
        assert!((gt.permutations_ln - 6.0f64.ln()).abs() < 1e-12);
    }

    #[test]
    fn two_snps_same_position_different_alt_are_distinguished() {
        // Confirm allele_cmp treats (pos, kind, ref, alt) as the identity —
        // same position + same kind + different alt must yield two
        // distinct alleles and a 3-genotype diploid enumeration.
        let alleles = vec![ref_a(100), snp_ag(100), snp_at(100)];
        let gts = enumerate_genotypes(&alleles, 2);
        // 3 input alleles, diploid → 6 genotypes.
        assert_eq!(gts.len(), 6);
        // There should be a het genotype with both A→G and A→T.
        let mixed_het = gts.iter().find(|g| {
            g.elements.len() == 2
                && g.allele_count(&snp_ag(100)) == 1
                && g.allele_count(&snp_at(100)) == 1
        });
        assert!(
            mixed_het.is_some(),
            "expected a het genotype containing A->G and A->T at pos 100"
        );
    }

    #[test]
    fn enumerate_genotypes_n4_diploid_is_ten() {
        // Explicit (n=4, k=2) count check — the reviewer flagged it.
        let alleles = vec![
            ref_a(100),
            snp_ag(100),
            snp_at(100),
            Allele::snp(100, b'A', b'C'),
        ];
        let gts = enumerate_genotypes(&alleles, 2);
        assert_eq!(gts.len(), 10); // C(n+k-1, k) = C(5, 2) = 10
    }

    #[test]
    fn genotype_types_cover_kind_variety() {
        // Spot check: enumeration works with a variety of AlleleKind values.
        let alleles = vec![
            ref_a(100),
            Allele::insertion(100, vec![b'T']),
            Allele::deletion(100, vec![b'C']),
        ];
        let gts = enumerate_genotypes(&alleles, 2);
        assert_eq!(gts.len(), 6); // C(3+2-1, 2)
                                  // Confirm we see INS and DEL represented.
        let kinds: std::collections::HashSet<AlleleKind> = gts
            .iter()
            .flat_map(|g| g.elements.iter().map(|e| e.allele.kind))
            .collect();
        assert!(kinds.contains(&AlleleKind::Insertion));
        assert!(kinds.contains(&AlleleKind::Deletion));
    }
}
