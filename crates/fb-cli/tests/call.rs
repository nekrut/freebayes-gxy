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
    let lines: Vec<&str> = stdout.lines().collect();
    assert_eq!(
        lines[0], "chrom\tpos\tref\talt\tgenotype\tgq\tdp\tlog_posterior",
        "header mismatch"
    );
    let body: Vec<&str> = lines.iter().skip(1).copied().collect();

    // Expect exactly two called sites: 1-based pos 6 (het) and 11 (hom-alt).
    assert_eq!(body.len(), 2, "expected 2 called sites, got: {body:#?}");

    // Het site at 1-based pos 6: ref C, alt A, genotype REF/SNP, dp=10.
    let het_row = body
        .iter()
        .find(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            f[0] == "chr1" && f[1] == "6"
        })
        .expect("het row missing");
    let het_fields: Vec<&str> = het_row.split('\t').collect();
    assert_eq!(het_fields[2], "C", "het ref mismatch");
    assert_eq!(het_fields[3], "A", "het alt mismatch");
    assert!(
        het_fields[4].contains('/'),
        "expected heterozygous genotype at pos 6, got {}",
        het_fields[4]
    );
    assert_eq!(het_fields[6], "10", "het dp mismatch");
    // GQ should be substantial on 5/5 Q30 at MQ60.
    let gq: f64 = het_fields[5].parse().expect("gq parse");
    assert!(gq > 30.0, "expected GQ > 30 on het site, got {gq}");

    // Hom-alt site at 1-based pos 11: ref G, alt T, genotype SNP, dp=10.
    let hom_row = body
        .iter()
        .find(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            f[0] == "chr1" && f[1] == "11"
        })
        .expect("hom-alt row missing");
    let hom_fields: Vec<&str> = hom_row.split('\t').collect();
    assert_eq!(hom_fields[2], "G");
    assert_eq!(hom_fields[3], "T");
    assert_eq!(
        hom_fields[4], "SNP",
        "expected homozygous SNP tag at pos 11, got {}",
        hom_fields[4]
    );
    assert_eq!(hom_fields[6], "10");
    // At 10x Q30 coverage, the next-best het genotype has P = 0.5^10
    // ≈ 10^-3 under the multinomial, capping GQ around 30. Lower
    // coverage or lower Q would push this down further.
    let gq: f64 = hom_fields[5].parse().expect("gq parse");
    assert!(gq > 25.0, "expected GQ > 25 on clean hom-alt, got {gq}");
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
    let body: Vec<&str> = stdout.lines().skip(1).collect();
    assert!(
        body.is_empty(),
        "expected zero called sites with 1/10 alt below threshold, got: {body:#?}"
    );
}
