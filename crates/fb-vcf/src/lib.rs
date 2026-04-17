//! VCF writer for freebayes-gxy.
//!
//! Later milestones will wrap `rust_htslib::bcf` with freebayes-specific
//! INFO/FORMAT field ordering and the reorder buffer needed for deterministic
//! output across threads. M0 ships the header-building helper used by
//! `fb-cli` so we emit a valid VCF 4.2 header with `##reference` and
//! `##contig` lines.
//!
//! This module intentionally builds the header as plain text rather than via
//! `rust_htslib::bcf::Header` — for M0 we only need to stream the header out,
//! and the text path is easier to snapshot in unit tests without requiring a
//! full BCF writer setup.

/// A minimal contig entry for VCF `##contig` header lines.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Contig {
    pub name: String,
    pub length: u64,
}

/// Build a minimal VCF 4.2 header.
///
/// Produces:
/// - `##fileformat=VCFv4.2`
/// - `##fileDate=YYYYMMDD` (omitted in M0 — date is a side effect we don't
///   want in unit tests; we'll add it in `fb-cli` at emission time once M1
///   starts producing records)
/// - `##source=freebayes-gxy-<version>`
/// - `##reference=<fasta path>`
/// - one `##contig=<ID=name,length=...>` line per contig, in input order
/// - the column header line, with a single sample column
///
/// The trailing newline is included.
pub fn build_header(reference: &str, contigs: &[Contig], sample: &str) -> String {
    let mut out = String::new();
    out.push_str("##fileformat=VCFv4.2\n");
    out.push_str(&format!(
        "##source=freebayes-gxy-{}\n",
        env!("CARGO_PKG_VERSION")
    ));
    out.push_str(&format!("##reference={reference}\n"));
    for c in contigs {
        out.push_str(&format!("##contig=<ID={},length={}>\n", c.name, c.length));
    }
    out.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t");
    out.push_str(sample);
    out.push('\n');
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn header_contains_required_lines() {
        let contigs = vec![
            Contig {
                name: "chr1".to_string(),
                length: 248_956_422,
            },
            Contig {
                name: "chr2".to_string(),
                length: 242_193_529,
            },
        ];
        let h = build_header("/tmp/ref.fa", &contigs, "SAMPLE1");
        assert!(h.starts_with("##fileformat=VCFv4.2\n"));
        assert!(h.contains("##reference=/tmp/ref.fa\n"));
        assert!(h.contains("##contig=<ID=chr1,length=248956422>\n"));
        assert!(h.contains("##contig=<ID=chr2,length=242193529>\n"));
        assert!(h.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n"));
    }

    #[test]
    fn header_preserves_contig_order() {
        let contigs = vec![
            Contig {
                name: "chrZ".to_string(),
                length: 100,
            },
            Contig {
                name: "chrA".to_string(),
                length: 200,
            },
        ];
        let h = build_header("ref.fa", &contigs, "S");
        let z_idx = h.find("chrZ").unwrap();
        let a_idx = h.find("chrA").unwrap();
        assert!(z_idx < a_idx, "contig order should be preserved");
    }
}
