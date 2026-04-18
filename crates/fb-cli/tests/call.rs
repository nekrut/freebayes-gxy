//! Integration test for `freebayes-gxy --call`.
//!
//! Builds a small multi-read BAM covering a few positions — one
//! heterozygous SNP, one homozygous SNP, and surrounding reference —
//! and checks that the Bayesian caller emits the expected per-site
//! genotype with plausible GQ.

use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::process::Command;

use rust_htslib::bam::{self, header::HeaderRecord, record::CigarString, Header};

const REFERENCE: &[u8] = b"ACGTACGTACGTACGT"; // 16 bp, contig "chr1"

fn write_fasta_with_index(dir: &Path) -> (std::path::PathBuf, std::path::PathBuf) {
    let fa = dir.join("ref.fa");
    let contents = format!(">chr1\n{}\n", std::str::from_utf8(REFERENCE).unwrap());
    {
        let mut f = File::create(&fa).unwrap();
        f.write_all(contents.as_bytes()).unwrap();
    }
    let fai = dir.join("ref.fa.fai");
    {
        let mut f = File::create(&fai).unwrap();
        writeln!(
            f,
            "chr1\t{}\t6\t{}\t{}",
            REFERENCE.len(),
            REFERENCE.len(),
            REFERENCE.len() + 1
        )
        .unwrap();
    }
    (fa, fai)
}

fn build_bam(dir: &Path) -> std::path::PathBuf {
    let bam_path = dir.join("calls.bam");
    let mut header = Header::new();
    let mut sq = HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", "chr1");
    sq.push_tag(b"LN", REFERENCE.len() as u32);
    header.push_record(&sq);

    let mut writer =
        bam::Writer::from_path(&bam_path, &header, bam::Format::Bam).expect("create writer");

    // Site layout:
    //   - Position 5 (0-based): heterozygous SNP C→A. 10 reads total:
    //     5 carry the reference C, 5 carry the alt A.
    //   - Position 10 (0-based): homozygous SNP G→T. All 10 reads carry T.
    // Each read is a 16M covering the whole contig, so every position
    // sees depth=10.
    let push = |writer: &mut bam::Writer, name: &[u8], seq: &[u8]| {
        let mut r = bam::Record::new();
        let cigar = CigarString(vec![bam::record::Cigar::Match(16)]);
        let quals = vec![30u8; 16];
        r.set(name, Some(&cigar), seq, &quals);
        r.set_flags(0);
        r.set_tid(0);
        r.set_pos(0);
        r.set_mapq(60);
        r.set_mpos(-1);
        r.set_mtid(-1);
        r.set_insert_size(0);
        writer.write(&r).unwrap();
    };

    // 5 reads with REF at 5 and SNP(G->T) at 10
    //   ref:     A C G T A C G T A C G T A C G T
    //   seq_ref: A C G T A C G T A C T T A C G T
    //   diff:                            ^---- pos 10: G→T
    for i in 0..5 {
        let name = format!("r_ref_{i}");
        push(&mut writer, name.as_bytes(), b"ACGTACGTACTTACGT");
    }
    // 5 reads with SNP(C->A) at 5 and SNP(G->T) at 10
    //   seq_alt: A C G T A A G T A C T T A C G T
    //   diff:               ^-- pos 5 C→A; pos 10 G→T
    for i in 0..5 {
        let name = format!("r_alt_{i}");
        push(&mut writer, name.as_bytes(), b"ACGTAAGTACTTACGT");
    }

    bam_path
}

fn cli_bin() -> std::path::PathBuf {
    env!("CARGO_BIN_EXE_freebayes-gxy").into()
}

#[test]
fn call_emits_het_and_hom_alt_rows() {
    let tmp = tempfile::tempdir().unwrap();
    let (fa, _fai) = write_fasta_with_index(tmp.path());
    let bam = build_bam(tmp.path());

    let output = Command::new(cli_bin())
        .arg("--call")
        .arg("-f")
        .arg(&fa)
        .arg(&bam)
        .output()
        .expect("spawn cli");
    assert!(
        output.status.success(),
        "CLI failed: stdout={}\nstderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );

    let stdout = String::from_utf8(output.stdout).unwrap();
    // M4 Phase A: output is now a VCF, not a TSV. Strip `##` header
    // lines and the `#CHROM` column header.
    let body: Vec<String> = stdout
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(String::from)
        .collect();

    // Expect exactly two called sites.
    assert_eq!(body.len(), 2, "expected 2 called records, got: {body:#?}");

    // VCF columns: CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE.
    let het_row = body
        .iter()
        .find(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            f[0] == "chr1" && f[1] == "6"
        })
        .expect("het row missing");
    let het_fields: Vec<&str> = het_row.split('\t').collect();
    assert_eq!(het_fields[3], "C", "het ref mismatch");
    assert_eq!(het_fields[4], "A", "het alt mismatch");
    // Sample GT must be heterozygous: 0/1.
    let het_gt = het_fields[9].split(':').next().unwrap();
    assert_eq!(het_gt, "0/1", "expected 0/1 at het site, got {het_gt}");
    // FORMAT DP (second field) must match coverage.
    let het_dp = het_fields[9].split(':').nth(1).unwrap();
    assert_eq!(het_dp, "10", "het DP mismatch");

    let hom_row = body
        .iter()
        .find(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            f[0] == "chr1" && f[1] == "11"
        })
        .expect("hom-alt row missing");
    let hom_fields: Vec<&str> = hom_row.split('\t').collect();
    assert_eq!(hom_fields[3], "G");
    assert_eq!(hom_fields[4], "T");
    let hom_gt = hom_fields[9].split(':').next().unwrap();
    assert_eq!(hom_gt, "1/1", "expected 1/1 at hom-alt site, got {hom_gt}");
    let hom_dp = hom_fields[9].split(':').nth(1).unwrap();
    assert_eq!(hom_dp, "10");
    // INFO field sanity: TYPE=snp.
    assert!(
        hom_fields[7].contains("TYPE=snp"),
        "expected TYPE=snp in INFO, got {}",
        hom_fields[7]
    );
}

#[test]
fn call_filters_sites_below_min_alternate_count() {
    // Build a BAM where position 5 has only 1 alt observation out of 10.
    // With default `--min-alternate-count 2`, the site must be filtered.
    let tmp = tempfile::tempdir().unwrap();
    let (fa, _fai) = write_fasta_with_index(tmp.path());
    let bam_path = tmp.path().join("low_alt.bam");
    {
        let mut header = Header::new();
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr1");
        sq.push_tag(b"LN", REFERENCE.len() as u32);
        header.push_record(&sq);
        let mut writer =
            bam::Writer::from_path(&bam_path, &header, bam::Format::Bam).expect("create writer");
        for i in 0..9 {
            let mut r = bam::Record::new();
            let cigar = CigarString(vec![bam::record::Cigar::Match(16)]);
            let quals = vec![30u8; 16];
            let name = format!("r_ref_{i}");
            r.set(name.as_bytes(), Some(&cigar), REFERENCE, &quals);
            r.set_flags(0);
            r.set_tid(0);
            r.set_pos(0);
            r.set_mapq(60);
            r.set_mpos(-1);
            r.set_mtid(-1);
            writer.write(&r).unwrap();
        }
        // Single alt read: C→A at pos 5.
        let mut r = bam::Record::new();
        let cigar = CigarString(vec![bam::record::Cigar::Match(16)]);
        let quals = vec![30u8; 16];
        r.set(b"r_lone_alt", Some(&cigar), b"ACGTAAGTACGTACGT", &quals);
        r.set_flags(0);
        r.set_tid(0);
        r.set_pos(0);
        r.set_mapq(60);
        writer.write(&r).unwrap();
    }

    let output = Command::new(cli_bin())
        .arg("--call")
        .arg("-f")
        .arg(&fa)
        .arg(&bam_path)
        .output()
        .expect("spawn cli");
    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    // VCF: expect header lines but no variant records.
    let body: Vec<&str> = stdout.lines().filter(|l| !l.starts_with('#')).collect();
    assert!(
        body.is_empty(),
        "expected zero called records with 1/10 alt below threshold, got: {body:#?}"
    );
}
