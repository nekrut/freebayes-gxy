//! Integration test for `freebayes-gxy --dump-alleles`.
//!
//! Builds a tiny single-contig FASTA + single-read BAM on disk, invokes the
//! binary, and checks that the TSV contains the expected SNP / Insertion /
//! Reference rows.

use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::process::Command;

use rust_htslib::bam::{self, header::HeaderRecord, record::CigarString, Header};

const REFERENCE: &[u8] = b"ACGTACGTACGTACGT"; // 16 bp, all contig "chr1"

fn write_fasta_with_index(dir: &Path) -> (std::path::PathBuf, std::path::PathBuf) {
    let fa = dir.join("ref.fa");
    // Single-line FASTA so the index is trivial to hand-roll.
    let contents = format!(">chr1\n{}\n", std::str::from_utf8(REFERENCE).unwrap());
    {
        let mut f = File::create(&fa).unwrap();
        f.write_all(contents.as_bytes()).unwrap();
    }
    // offset = len(">chr1\n") = 6, linebases = 16, linewidth = 17 (with \n).
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
    let bam_path = dir.join("test.bam");
    let mut header = Header::new();
    let mut sq = HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", "chr1");
    sq.push_tag(b"LN", REFERENCE.len() as u32);
    header.push_record(&sq);

    let mut writer =
        bam::Writer::from_path(&bam_path, &header, bam::Format::Bam).expect("create writer");

    // Read 1: 16 M, perfect match except a SNP at position 5 (ref 'C' -> read 'A').
    // ref  = A C G T A C G T A C G T A C G T
    // read = A C G T A A G T A C G T A C G T  (pos 5 differs)
    {
        let mut r = bam::Record::new();
        let cigar = CigarString(vec![bam::record::Cigar::Match(16)]);
        let seq = b"ACGTAAGTACGTACGT";
        let quals = vec![30u8; 16];
        r.set(b"r_snp", Some(&cigar), seq, &quals);
        r.set_flags(0); // rust-htslib defaults to FUNMAP=0x4 — clear it.
        r.set_tid(0);
        r.set_pos(0);
        r.set_mapq(60);
        r.set_mpos(-1);
        r.set_mtid(-1);
        r.set_insert_size(0);
        writer.write(&r).unwrap();
    }

    // Read 2: 5M 2I 5M. ref covered = first 10 of "ACGTACGTAC"; insertion "GG"
    // lives between ref positions 5 and 6 — anchors at pos 5 (0-based).
    {
        let mut r = bam::Record::new();
        let cigar = CigarString(vec![
            bam::record::Cigar::Match(5),
            bam::record::Cigar::Ins(2),
            bam::record::Cigar::Match(5),
        ]);
        let seq = b"ACGTAGGCGTAC"; // 12 bases: 5M + 2I + 5M
        let quals = vec![40u8; 12];
        r.set(b"r_ins", Some(&cigar), seq, &quals);
        r.set_flags(0); // rust-htslib defaults to FUNMAP=0x4 — clear it.
        r.set_tid(0);
        r.set_pos(0);
        r.set_mapq(60);
        r.set_mpos(-1);
        r.set_mtid(-1);
        r.set_insert_size(0);
        writer.write(&r).unwrap();
    }

    bam_path
}

fn cli_bin() -> std::path::PathBuf {
    // `CARGO_BIN_EXE_<binname>` is set by cargo at integration-test time.
    env!("CARGO_BIN_EXE_freebayes-gxy").into()
}

#[test]
fn dump_alleles_emits_snp_and_insertion_rows() {
    let tmp = tempfile::tempdir().unwrap();
    let (fa, _fai) = write_fasta_with_index(tmp.path());
    let bam = build_bam(tmp.path());

    let output = Command::new(cli_bin())
        .arg("--dump-alleles")
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
        lines[0], "chrom\tpos\tref\talt\tkind\tcount\tforward\treverse\tbq_sum",
        "header mismatch"
    );
    let body: Vec<&str> = lines.iter().skip(1).copied().collect();

    // SNP at 1-based position 6 (0-based 5), ref C -> alt A, 1 forward read.
    assert!(
        body.iter().any(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            f[0] == "chr1"
                && f[1] == "6"
                && f[2] == "C"
                && f[3] == "A"
                && f[4] == "SNP"
                && f[5] == "1"
        }),
        "expected SNP row, got: {:#?}",
        body
    );

    // Insertion anchored at 1-based position 6 (0-based 5), alt "GG", 1 read.
    assert!(
        body.iter().any(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            f[0] == "chr1"
                && f[1] == "6"
                && f[2] == "."
                && f[3] == "GG"
                && f[4] == "INS"
                && f[5] == "1"
        }),
        "expected INS row, got: {:#?}",
        body
    );

    // Both reads share an identical leading reference run (5M = "ACGTA") at
    // pos 0 — per-run emission collapses them into a single row with count=2.
    assert!(
        body.iter().any(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            f[0] == "chr1"
                && f[1] == "1"
                && f[2] == "ACGTA"
                && f[3] == "ACGTA"
                && f[4] == "REF"
                && f[5] == "2"
        }),
        "expected REF row at pos 1 with len 5 count 2, got: {:#?}",
        body
    );
}
