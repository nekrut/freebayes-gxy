//! Work-stealing window scheduler for freebayes-gxy.
//!
//! In later milestones this crate hosts the Rayon-backed tile scheduler that
//! drives parallel variant calling (see `PLAN.md` §4 and M5). M0 ships a
//! placeholder `Window` type and a tiny helper that splits a contig length
//! into fixed-size tiles so downstream crates can depend on a stable surface
//! from day one.

/// A half-open genomic window `[start, end)` on a single contig.
///
/// Positions are zero-based to match `rust_htslib` conventions; the VCF
/// writer adapts to one-based coordinates at emit time.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Window {
    pub contig: String,
    pub start: u64,
    pub end: u64,
}

/// Split a contig of `length` bases into tiles of at most `tile_size` bases.
///
/// Used by the M5 scheduler; exposed now so other crates can unit-test
/// against a stable window layout.
pub fn tile_contig(contig: &str, length: u64, tile_size: u64) -> Vec<Window> {
    assert!(tile_size > 0, "tile_size must be positive");
    let mut out = Vec::new();
    let mut start = 0u64;
    while start < length {
        let end = (start + tile_size).min(length);
        out.push(Window {
            contig: contig.to_string(),
            start,
            end,
        });
        start = end;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tile_contig_splits_evenly() {
        let tiles = tile_contig("chr1", 250, 100);
        assert_eq!(tiles.len(), 3);
        assert_eq!(
            tiles[0],
            Window {
                contig: "chr1".into(),
                start: 0,
                end: 100
            }
        );
        assert_eq!(
            tiles[1],
            Window {
                contig: "chr1".into(),
                start: 100,
                end: 200
            }
        );
        assert_eq!(
            tiles[2],
            Window {
                contig: "chr1".into(),
                start: 200,
                end: 250
            }
        );
    }

    #[test]
    fn tile_contig_handles_exact_multiple() {
        let tiles = tile_contig("chr1", 200, 100);
        assert_eq!(tiles.len(), 2);
        assert_eq!(tiles.last().unwrap().end, 200);
    }

    #[test]
    fn tile_contig_empty_contig() {
        assert!(tile_contig("chr1", 0, 100).is_empty());
    }
}
