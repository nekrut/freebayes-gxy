//! `freebayes-gxy` CLI entry point.
//!
//! M0 shipped a BAM+FASTA opener that wrote a minimal VCF header. M1 adds
//! `--dump-alleles`, which walks the BAM and emits a per-allele TSV summary
//! using [`fb_core::walk_record`]. Variant calling proper arrives in M3–M4.

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{anyhow, Context, Result};
use clap::Parser;
use fb_core::{clump_observations, walk_record, AlleleKind, ReadFilter, Strand};
use fb_vcf::{build_header, Contig};
use rust_htslib::bam::{self, Read as _};
use rust_htslib::faidx;
use tracing::{debug, info, warn};

/// freebayes-gxy — a Rust reimplementation of freebayes.
#[derive(Debug, Parser)]
#[command(name = "freebayes-gxy", version, about, long_about = None)]
struct Cli {
    /// Reference FASTA (must be indexed, i.e. `.fai` present).
    #[arg(short = 'f', long = "fasta-reference", value_name = "FASTA")]
    fasta: PathBuf,

    /// Output path. When omitted, output is written to stdout.
    #[arg(short = 'v', long = "vcf", value_name = "VCF")]
    output: Option<PathBuf>,

    /// Sample name to use in the `#CHROM` line.
    #[arg(long = "sample", default_value = "SAMPLE")]
    sample: String,

    /// Emit a per-allele TSV instead of a VCF (M1 debug dump).
    #[arg(long = "dump-alleles")]
    dump_alleles: bool,

    /// Maximum reference-run length (bp) allowed to bridge two flanking
    /// non-reference events into a single COMPLEX allele (M2 clumping).
    /// Set to a negative value to disable clumping. `--haplotype-length`
    /// is the upstream freebayes CLI spelling and is accepted as an
    /// alias. Use the `--flag=value` form for negative values, as clap
    /// parses bare `-1` as a short flag.
    #[arg(
        long = "haplotype-length",
        alias = "max-complex-gap",
        default_value_t = 3,
        value_name = "BP"
    )]
    haplotype_length: i64,

    /// Input BAM file (positional, matches upstream freebayes).
    #[arg(value_name = "BAM")]
    bam: PathBuf,
}

fn main() -> Result<()> {
    tracing_subscriber::fmt()
        .with_env_filter(
            tracing_subscriber::EnvFilter::try_from_default_env()
                .unwrap_or_else(|_| tracing_subscriber::EnvFilter::new("info")),
        )
        .with_target(false)
        .with_writer(std::io::stderr)
        .init();

    let cli = Cli::parse();
    info!(
        milestone = fb_core::current_milestone().name,
        "freebayes-gxy starting"
    );

    if cli.dump_alleles {
        run_dump_alleles(&cli)
    } else {
        run_emit_header(&cli)
    }
}

fn run_emit_header(cli: &Cli) -> Result<()> {
    open_fasta(&cli.fasta)
        .with_context(|| format!("failed to open reference FASTA {:?}", cli.fasta))?;
    let contigs = read_bam_contigs(&cli.bam)
        .with_context(|| format!("failed to read BAM header from {:?}", cli.bam))?;
    info!(n_contigs = contigs.len(), "loaded BAM header");

    let reference = cli.fasta.to_string_lossy().into_owned();
    let header = build_header(&reference, &contigs, &cli.sample);
    write_output(cli.output.as_deref(), header.as_bytes()).context("failed to write VCF header")?;
    Ok(())
}

/// Aggregation key for the TSV dump: one row per distinct allele at a
/// reference position on a given contig.
#[derive(PartialEq, Eq, PartialOrd, Ord)]
struct AlleleKey {
    tid: i32,
    pos: i64,
    kind_tag: &'static str,
    ref_seq: Vec<u8>,
    alt_seq: Vec<u8>,
}

#[derive(Default)]
struct AlleleCounts {
    total: u32,
    forward: u32,
    reverse: u32,
    bq_sum: u64,
}

fn run_dump_alleles(cli: &Cli) -> Result<()> {
    let fasta = open_fasta(&cli.fasta)
        .with_context(|| format!("failed to open reference FASTA {:?}", cli.fasta))?;
    let mut bam_reader = bam::Reader::from_path(&cli.bam)
        .with_context(|| format!("failed to open BAM {:?}", cli.bam))?;
    // Build contig name + length table once; we fetch reference sequences
    // lazily per tid as records stream in.
    let contigs = read_bam_contigs(&cli.bam)?;
    let target_names: Vec<String> = contigs.iter().map(|c| c.name.clone()).collect();
    let target_lens: Vec<u64> = contigs.iter().map(|c| c.length).collect();

    info!(
        n_contigs = contigs.len(),
        "dumping allele observations (M1)"
    );

    let filter = ReadFilter::default();
    let mut current_tid: i32 = -1;
    let mut current_ref: Vec<u8> = Vec::new();
    let mut agg: BTreeMap<AlleleKey, AlleleCounts> = BTreeMap::new();
    let mut n_reads: u64 = 0;
    let mut n_filtered: u64 = 0;

    let mut record = bam::Record::new();
    while let Some(result) = bam_reader.read(&mut record) {
        result.context("BAM record read failed")?;
        if record.tid() < 0 {
            continue;
        }
        let tid = record.tid();
        if tid != current_tid {
            current_ref = fetch_contig(&fasta, &target_names, &target_lens, tid as usize)?;
            current_tid = tid;
            debug!(
                contig = target_names[tid as usize].as_str(),
                "fetched reference"
            );
        }
        n_reads += 1;
        debug!(
            qname = %std::str::from_utf8(record.qname()).unwrap_or("?"),
            mapq = record.mapq(),
            flags = record.flags(),
            pos = record.pos(),
            "read record"
        );
        let observations = match walk_record(&record, &current_ref, 0, &filter) {
            Some(v) => v,
            None => {
                n_filtered += 1;
                continue;
            }
        };
        // M2 clumping: collapse adjacent non-reference events (with up to
        // `max_complex_gap` bp of intervening reference) into Complex
        // alleles before aggregating for the TSV.
        let observations = clump_observations(&observations, cli.haplotype_length);
        for obs in observations {
            if obs.allele.kind == AlleleKind::Null {
                continue;
            }
            let key = AlleleKey {
                tid,
                pos: obs.allele.position,
                kind_tag: obs.allele.kind.as_tag(),
                ref_seq: obs.allele.ref_seq.clone(),
                alt_seq: obs.allele.alt_seq.clone(),
            };
            let entry = agg.entry(key).or_default();
            entry.total += 1;
            match obs.strand {
                Strand::Forward => entry.forward += 1,
                Strand::Reverse => entry.reverse += 1,
            }
            entry.bq_sum += u64::from(obs.base_quality_sum);
        }
    }

    info!(
        n_reads = n_reads,
        n_filtered = n_filtered,
        n_alleles = agg.len(),
        "dump complete"
    );

    emit_tsv(cli.output.as_deref(), &target_names, &agg)?;
    Ok(())
}

fn emit_tsv(
    output: Option<&Path>,
    target_names: &[String],
    agg: &BTreeMap<AlleleKey, AlleleCounts>,
) -> Result<()> {
    let mut buf: Vec<u8> = Vec::new();
    writeln!(
        buf,
        "chrom\tpos\tref\talt\tkind\tcount\tforward\treverse\tbq_sum"
    )?;
    for (key, counts) in agg {
        let chrom = target_names
            .get(key.tid as usize)
            .map(|s| s.as_str())
            .unwrap_or("?");
        // Convert 0-based internal to 1-based VCF-style for display.
        let display_pos = key.pos + 1;
        let ref_display = if key.ref_seq.is_empty() {
            "."
        } else {
            std::str::from_utf8(&key.ref_seq).unwrap_or(".")
        };
        let alt_display = if key.alt_seq.is_empty() {
            "."
        } else {
            std::str::from_utf8(&key.alt_seq).unwrap_or(".")
        };
        writeln!(
            buf,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            chrom,
            display_pos,
            ref_display,
            alt_display,
            key.kind_tag,
            counts.total,
            counts.forward,
            counts.reverse,
            counts.bq_sum,
        )?;
    }
    write_output(output, &buf)
}

fn fetch_contig(
    fasta: &faidx::Reader,
    names: &[String],
    lens: &[u64],
    tid: usize,
) -> Result<Vec<u8>> {
    let name = names
        .get(tid)
        .ok_or_else(|| anyhow!("contig tid {tid} not in BAM header"))?;
    let len = *lens
        .get(tid)
        .ok_or_else(|| anyhow!("contig tid {tid} has no recorded length"))?;
    if len == 0 {
        warn!(contig = name.as_str(), "contig length 0 in BAM header");
        return Ok(Vec::new());
    }
    let end = (len - 1) as usize;
    let seq = fasta
        .fetch_seq(name, 0, end)
        .with_context(|| format!("fetch_seq failed for contig {name:?}"))?;
    // Normalise to uppercase ASCII bases so our ATGC checks behave.
    Ok(seq.iter().map(|b| b.to_ascii_uppercase()).collect())
}

fn open_fasta(path: &Path) -> Result<faidx::Reader> {
    faidx::Reader::from_path(path).map_err(anyhow::Error::from)
}

fn read_bam_contigs(path: &Path) -> Result<Vec<Contig>> {
    let reader = bam::Reader::from_path(path)?;
    let header = reader.header();
    let names = header.target_names();
    let mut out = Vec::with_capacity(names.len());
    for (tid, name_bytes) in names.iter().enumerate() {
        let name = std::str::from_utf8(name_bytes)
            .with_context(|| format!("contig name at tid={tid} is not valid UTF-8"))?
            .to_string();
        let length = header.target_len(tid as u32).with_context(|| {
            format!("missing length for contig {name:?} (tid={tid}) in BAM header")
        })?;
        out.push(Contig { name, length });
    }
    Ok(out)
}

fn write_output(path: Option<&Path>, bytes: &[u8]) -> Result<()> {
    match path {
        Some(p) => {
            let file = File::create(p).with_context(|| format!("creating {p:?}"))?;
            let mut w = BufWriter::new(file);
            w.write_all(bytes)?;
            w.flush()?;
        }
        None => {
            let stdout = std::io::stdout();
            let mut w = stdout.lock();
            w.write_all(bytes)?;
            w.flush()?;
        }
    }
    Ok(())
}
