//! `freebayes-gxy` CLI entry point.
//!
//! For M0 this binary only opens the reference FASTA and the input BAM, then
//! writes a minimal VCF 4.2 header derived from the BAM's `@SQ` records to
//! stdout or `-v <path>`. Variant calling arrives in M1+.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use clap::Parser;
use fb_vcf::{build_header, Contig};
use rust_htslib::bam::{self, Read as _};
use rust_htslib::faidx;
use tracing::info;

/// freebayes-gxy — a Rust reimplementation of freebayes (M0 scaffolding).
#[derive(Debug, Parser)]
#[command(name = "freebayes-gxy", version, about, long_about = None)]
struct Cli {
    /// Reference FASTA (must be indexed, i.e. `.fai` present).
    #[arg(short = 'f', long = "fasta-reference", value_name = "FASTA")]
    fasta: PathBuf,

    /// Output VCF path. When omitted, the header is written to stdout.
    #[arg(short = 'v', long = "vcf", value_name = "VCF")]
    vcf: Option<PathBuf>,

    /// Sample name to use in the `#CHROM` line.
    ///
    /// M1+ will read this from the BAM `@RG SM:` tag; for M0 we take it as
    /// an option with a default so the placeholder header is well-formed.
    #[arg(long = "sample", default_value = "SAMPLE")]
    sample: String,

    /// Input BAM file (positional, matches upstream freebayes).
    #[arg(value_name = "BAM")]
    bam: PathBuf,
}

fn main() -> Result<()> {
    // Honour RUST_LOG; default to `info` so users see what the CLI is doing.
    tracing_subscriber::fmt()
        .with_env_filter(
            tracing_subscriber::EnvFilter::try_from_default_env()
                .unwrap_or_else(|_| tracing_subscriber::EnvFilter::new("info")),
        )
        .with_target(false)
        .init();

    let cli = Cli::parse();
    info!(
        milestone = fb_core::current_milestone().name,
        "freebayes-gxy starting"
    );

    // Open FASTA to validate it exists and is indexable. We don't read any
    // sequence in M0; M1 will use the index for reference lookups.
    open_fasta(&cli.fasta)
        .with_context(|| format!("failed to open reference FASTA {:?}", cli.fasta))?;

    // Open BAM and extract contig names + lengths from the header.
    let contigs = read_bam_contigs(&cli.bam)
        .with_context(|| format!("failed to read BAM header from {:?}", cli.bam))?;
    info!(n_contigs = contigs.len(), "loaded BAM header");

    let reference = cli.fasta.to_string_lossy().into_owned();
    let header = build_header(&reference, &contigs, &cli.sample);

    write_output(cli.vcf.as_deref(), header.as_bytes()).context("failed to write VCF header")?;

    Ok(())
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
