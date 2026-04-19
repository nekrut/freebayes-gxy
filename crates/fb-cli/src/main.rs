//! `freebayes-gxy` CLI entry point.
//!
//! Modes:
//! - `--dump-alleles` (M1+M2): per-allele TSV summary.
//! - `--call` (M3 + M4 Phase A): per-position pileup + Bayesian single-
//!   sample genotype call, freebayes-compatible VCF output with proper
//!   anchor-base synthesis for indels.
//! - default: emits the M0 VCF header only.

use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use anyhow::{anyhow, Context, Result};
use clap::Parser;
use fb_core::{
    clump_observations, walk_record, Allele, AlleleKind, AlleleObservation, ReadFilter, Strand,
};
use fb_genotype::{call_genotype, enumerate_genotypes, Parameters};
use fb_vcf::{
    alt_cigar, build_header, synthesize_anchored, vcf_gl_index, write_record, Contig, Record,
    RecordKind,
};
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

    /// Run the single-sample Bayesian caller and emit a per-site TSV
    /// (M3 Phase C-2). Mutually exclusive with `--dump-alleles`.
    #[arg(long = "call", conflicts_with = "dump_alleles")]
    call: bool,

    /// Ploidy for the Bayesian caller (upstream default: 2).
    #[arg(long = "ploidy", default_value_t = 2, value_name = "N")]
    ploidy: u32,

    /// Minimum count of non-reference observations at a site for it to
    /// be considered a candidate. Upstream
    /// `--min-alternate-count` default is 2.
    #[arg(long = "min-alternate-count", default_value_t = 2, value_name = "N")]
    min_alternate_count: u32,

    /// Minimum fraction of non-reference observations at a site for it
    /// to be considered a candidate. Upstream
    /// `--min-alternate-fraction` default is 0.05.
    #[arg(
        long = "min-alternate-fraction",
        default_value_t = 0.05,
        value_name = "F"
    )]
    min_alternate_fraction: f64,

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

    /// Number of worker threads for the per-site caller (M5 Phase A).
    /// `1` (default) uses the original single-threaded path. Higher
    /// values shard the BAM into genomic tiles (size controlled by
    /// `--tile-size`) and run the pileup + call loop in parallel via
    /// Rayon. Output is re-sorted so the VCF is identical across
    /// thread counts.
    #[arg(short = 't', long = "threads", default_value_t = 1, value_name = "N")]
    threads: usize,

    /// Tile size (bp) for the M5 parallel scheduler. The BAM is split
    /// into non-overlapping windows of this size along each contig;
    /// each window is processed by one worker. Smaller values
    /// increase parallelism at the cost of per-tile overhead; larger
    /// values reduce overhead but limit concurrency on short contigs.
    /// Default 100 kb matches upstream `freebayes-parallel`
    /// recommendations for ~30× WGS. Only consulted when
    /// `--threads > 1`.
    #[arg(long = "tile-size", default_value_t = 100_000, value_name = "BP")]
    tile_size: u64,

    /// Number of background threads to hand to htslib for BGZF
    /// decompression of the input BAM (M5 Phase G). `0` (default)
    /// leaves the reader single-threaded, matching M5-A..M5-F
    /// behaviour. Values `> 0` call `bam::Reader::set_threads(n)`;
    /// htslib decodes blocks in a worker pool while the main thread
    /// does CIGAR/pileup work. Independent of `--threads`, which
    /// controls the per-site call stage.
    #[arg(long = "bam-threads", default_value_t = 0, value_name = "N")]
    bam_threads: usize,

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

    if cli.call {
        run_call(&cli)
    } else if cli.dump_alleles {
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

// ---------------------------------------------------------------------------
// M3 Phase C-2: per-position pileup + Bayesian caller
// ---------------------------------------------------------------------------

/// Pileup of allele observations, split into **variant sites** and
/// **reference runs** so we don't pay the O(total aligned bases)
/// allocation cost of decomposing every REF run into per-position
/// observations.
///
/// Before M5 Phase D the pileup materialised a
/// `BTreeMap<(tid, pos), Vec<AlleleObservation>>` at ingest — with a
/// separate `AlleleObservation` at every base of every REF run. On
/// the 100 kb benchmark fixture that was 95 %+ of wall-clock
/// (3M allocations × 3 heap objects each = ~9M allocs) even though
/// ~99 % of those entries were never consumed (pure-REF positions
/// short-circuit in `call_site`). See `docs/profile-report-m5c.md`.
///
/// The Phase D shape:
/// - `variant_sites` holds per-position observation vectors **only**
///   for positions where at least one non-REF event (SNP / INS / DEL
///   / Complex) landed. These are the sites that `call_site` actually
///   visits.
/// - `ref_runs_per_tid` stores REF runs as intervals, sorted by
///   start. At call time, `runs_covering(tid, pos)` binary-searches
///   for the ~depth runs overlapping a candidate site; we
///   materialise their single-base REF observations on demand.
///
/// Semantics are byte-identical to the old design — the materialised
/// per-position REF observation inherits the parent run's per-base
/// Phred (via `ref_seq[offset]` and `per_base_quals[offset]`), the
/// read's `read_ref_start` for the drift-2 event-span filter, and the
/// strand / mapq / name needed for the per-read dedupe. All M3 / M4
/// invariants are preserved.
#[derive(Default)]
struct Pileup {
    variant_sites: BTreeMap<(i32, i64), Vec<AlleleObservation>>,
    ref_runs_per_tid: HashMap<i32, TidRuns>,
}

/// Per-contig collection of REF runs. Kept sorted by `start`, with
/// `max_len` cached so `runs_covering` can bound its search window.
#[derive(Default)]
struct TidRuns {
    runs: Vec<RefRun>,
    max_len: usize,
    // `true` once `runs` has been sorted by start. `add_run` sets
    // `false`; `ensure_sorted` (called before the first query) sorts.
    sorted: bool,
}

/// A single contiguous REF match run produced by the M1 walker.
/// Materialisation of per-position REF observations happens at
/// `call_site` time via `materialize_at`.
#[derive(Clone)]
struct RefRun {
    start: i64,
    length: usize,
    read_name: String,
    mapq: u8,
    strand: Strand,
    read_position: usize,
    is_proper_pair: bool,
    read_ref_start: i64,
    per_base_quals: Vec<u8>,
    ref_seq: Vec<u8>,
    /// Fallback BQ scalar for decomposition when `per_base_quals` is
    /// shorter than `length` (shouldn't happen for walker-emitted
    /// runs, but test observations can construct with empty vec).
    fallback_bq: u32,
}

impl RefRun {
    /// True if this run covers the given 0-based reference position.
    #[inline]
    fn covers(&self, pos: i64) -> bool {
        self.start <= pos && pos < self.start + self.length as i64
    }

    /// Materialise the single-base REF `AlleleObservation` at `pos`
    /// from this run. Caller has already verified `covers(pos)`.
    fn materialize_at(&self, pos: i64) -> AlleleObservation {
        let offset = (pos - self.start) as usize;
        let byte = self.ref_seq[offset];
        let per_pos_bq: u32 = self
            .per_base_quals
            .get(offset)
            .copied()
            .map(u32::from)
            .unwrap_or(self.fallback_bq);
        AlleleObservation {
            allele: Allele::reference(pos, vec![byte]),
            read_name: self.read_name.clone(),
            mapq: self.mapq,
            base_quality_sum: per_pos_bq,
            strand: self.strand,
            read_position: self.read_position + offset,
            is_proper_pair: self.is_proper_pair,
            read_ref_start: self.read_ref_start,
            per_base_quals: Vec::new(),
        }
    }
}

impl TidRuns {
    fn ensure_sorted(&mut self) {
        if !self.sorted {
            self.runs.sort_unstable_by_key(|r| r.start);
            self.sorted = true;
        }
    }

    /// Iterator over runs covering `pos`. Relies on `runs` being
    /// sorted by start. Narrows via two `partition_point` calls
    /// bounded by `max_len`, then filters the candidate window.
    fn covering(&self, pos: i64) -> impl Iterator<Item = &RefRun> {
        debug_assert!(self.sorted, "TidRuns::covering requires ensure_sorted()");
        let min_start = pos.saturating_sub(self.max_len as i64 - 1).max(0);
        let lo = self.runs.partition_point(|r| r.start < min_start);
        let hi = self.runs.partition_point(|r| r.start <= pos);
        self.runs[lo..hi].iter().filter(move |r| r.covers(pos))
    }
}

impl Pileup {
    fn add_read_observations(&mut self, tid: i32, observations: Vec<AlleleObservation>) {
        for obs in observations {
            match obs.allele.kind {
                AlleleKind::Reference => {
                    // Phase D: park REF runs as intervals; no per-
                    // position decomposition.
                    let entry = self.ref_runs_per_tid.entry(tid).or_default();
                    let length = obs.allele.length;
                    if length > entry.max_len {
                        entry.max_len = length;
                    }
                    entry.sorted = false;
                    entry.runs.push(RefRun {
                        start: obs.allele.position,
                        length,
                        read_name: obs.read_name,
                        mapq: obs.mapq,
                        strand: obs.strand,
                        read_position: obs.read_position,
                        is_proper_pair: obs.is_proper_pair,
                        read_ref_start: obs.read_ref_start,
                        per_base_quals: obs.per_base_quals,
                        ref_seq: obs.allele.ref_seq,
                        fallback_bq: obs.base_quality_sum,
                    });
                }
                AlleleKind::Null => {
                    // Soft-clip / N-base — not informative.
                }
                _ => {
                    self.variant_sites
                        .entry((tid, obs.allele.position))
                        .or_default()
                        .push(obs);
                }
            }
        }
    }

    /// Sort every contig's REF run list by start. Must be called
    /// once after all `add_read_observations` and before any
    /// `covering` query.
    fn finalize(&mut self) {
        for runs in self.ref_runs_per_tid.values_mut() {
            runs.ensure_sorted();
        }
    }
}

/// Build the combined observation vector for a candidate site: the
/// site's variant observations (passed in, already held by the
/// pileup) followed by the per-position REF observations materialised
/// from every REF run that covers `pos`.
fn materialize_site_observations(
    pileup: &Pileup,
    tid: i32,
    pos: i64,
    variant_obs: &[AlleleObservation],
) -> Vec<AlleleObservation> {
    let mut combined: Vec<AlleleObservation> = Vec::with_capacity(variant_obs.len() + 32);
    combined.extend_from_slice(variant_obs);
    if let Some(runs) = pileup.ref_runs_per_tid.get(&tid) {
        for run in runs.covering(pos) {
            combined.push(run.materialize_at(pos));
        }
    }
    combined
}

/// Given a position's observation bucket, assemble the candidate
/// allele set (reference + distinct variants meeting the min-count /
/// min-fraction thresholds), run the Bayesian caller, and package the
/// VCF-ready fields.
fn call_site(
    tid: i32,
    pos: i64,
    observations: &[AlleleObservation],
    ref_base: u8,
    cli: &Cli,
    params: &Parameters,
) -> Option<SiteCall> {
    // Collect candidate variant alleles with their observation counts.
    let mut variant_counts: HashMap<Allele, u32> = HashMap::new();
    let mut total: u32 = 0;
    for obs in observations {
        total += 1;
        if obs.allele.kind != AlleleKind::Reference {
            *variant_counts.entry(obs.allele.clone()).or_insert(0) += 1;
        }
    }
    if total == 0 {
        return None;
    }

    // Apply min-alt-count / min-alt-fraction thresholds and sort
    // candidates deterministically so ALT ordering is stable across
    // runs (required for VCF parity).
    let total_f = total as f64;
    let mut surviving: Vec<(Allele, u32)> = variant_counts
        .into_iter()
        .filter(|(_, c)| {
            *c >= cli.min_alternate_count && (*c as f64) / total_f >= cli.min_alternate_fraction
        })
        .collect();
    if surviving.is_empty() {
        return None;
    }
    surviving.sort_by(|(a, _), (b, _)| {
        a.position
            .cmp(&b.position)
            .then_with(|| (a.kind as u8).cmp(&(b.kind as u8)))
            .then_with(|| a.ref_seq.cmp(&b.ref_seq))
            .then_with(|| a.alt_seq.cmp(&b.alt_seq))
    });

    // Candidate allele set = reference + surviving variants (in canonical order).
    let mut alleles: Vec<Allele> = Vec::with_capacity(surviving.len() + 1);
    alleles.push(Allele::reference(pos, vec![ref_base]));
    alleles.extend(surviving.into_iter().map(|(a, _)| a));

    let gts = enumerate_genotypes(&alleles, cli.ploidy);
    if gts.is_empty() {
        return None;
    }

    // Drift 2 fix (event-span admission). At a candidate site with any
    // non-SNP allele, a REF observation is only informative if the
    // source read spans the variant's anchor base (0-based `pos - 1`).
    // Reads starting at `pos` or later have a match op at `pos` but do
    // not see the anchor, so they can't distinguish REF from INS/DEL
    // — upstream's event-span pileup excludes them. Replicate that
    // here by filtering REF observations with `read_ref_start >= pos`
    // whenever the candidate set contains an INS, DEL, or Complex.
    let site_has_indel_candidate = alleles[1..].iter().any(|a| {
        matches!(
            a.kind,
            AlleleKind::Insertion | AlleleKind::Deletion | AlleleKind::Complex
        )
    });
    let filtered_observations: Vec<AlleleObservation>;
    let effective_observations: &[AlleleObservation] = if site_has_indel_candidate {
        // Step 1: drop REF observations from reads whose alignment does
        // not span the anchor (read_ref_start >= pos).
        let span_filtered: Vec<&AlleleObservation> = observations
            .iter()
            .filter(|o| o.allele.kind != AlleleKind::Reference || o.read_ref_start < pos)
            .collect();
        // Step 2: dedupe by read_name, preferring non-REF. Upstream
        // treats a read's event (INS/DEL) as its single observation at
        // the anchor — our walker can emit both the event AND a
        // trailing REF at the same position when the read carries the
        // indel, so without this step the same read double-counts
        // toward both AO and RO.
        use std::collections::HashMap as StdHashMap;
        let mut per_read: StdHashMap<&str, &AlleleObservation> = StdHashMap::new();
        for obs in span_filtered {
            let key = obs.read_name.as_str();
            let promote = match per_read.get(key) {
                Some(existing) => {
                    existing.allele.kind == AlleleKind::Reference
                        && obs.allele.kind != AlleleKind::Reference
                }
                None => true,
            };
            if promote {
                per_read.insert(key, obs);
            }
        }
        filtered_observations = per_read.into_values().cloned().collect();
        &filtered_observations
    } else {
        observations
    };

    // Recount DP / RO / AO from the event-span-filtered set so the
    // VCF sample fields reflect what the caller scored. For SNP-only
    // sites this is a no-op (filtered == observations).
    let effective_total: u32 = effective_observations.len() as u32;
    let effective_ref_obs: u32 = effective_observations
        .iter()
        .filter(|o| o.allele.kind == AlleleKind::Reference)
        .count() as u32;
    let alt_obs: Vec<u32> = alleles
        .iter()
        .skip(1)
        .map(|a| {
            effective_observations
                .iter()
                .filter(|o| &o.allele == a)
                .count() as u32
        })
        .collect();

    let call = call_genotype(&gts, effective_observations, params);

    // VCF-style GT indices: walk the MAP genotype's elements and map
    // each to its index in `alleles` (0 = REF, 1..N = alts in the
    // sorted order above).
    let best = &gts[call.best_index];
    let mut gt_indices: Vec<u8> = Vec::with_capacity(best.ploidy as usize);
    for elem in &best.elements {
        let idx = alleles
            .iter()
            .position(|a| a == &elem.allele)
            .expect("best genotype's element must be in the candidate set") as u8;
        for _ in 0..elem.count {
            gt_indices.push(idx);
        }
    }
    gt_indices.sort(); // canonical GT order (e.g. 0/1 not 1/0).

    // Site-level variant quality: P(site is NOT variant) is the
    // hom-reference genotype's posterior. QUAL is its Phred-scaled
    // complement: -10 * log10(P(hom-ref)), clamped at MAX_GQ.
    let hom_ref_idx = gts
        .iter()
        .position(|g| g.homozygous && g.elements[0].allele == alleles[0]);
    let qual = hom_ref_idx.map(|i| {
        let log_p = call.log_posteriors[i];
        if log_p == f64::NEG_INFINITY {
            fb_genotype::MAX_GQ
        } else {
            (-10.0 * log_p / std::f64::consts::LN_10).clamp(0.0, fb_genotype::MAX_GQ)
        }
    });

    let alt_kinds: Vec<RecordKind> = alleles
        .iter()
        .skip(1)
        .map(|a| match a.kind {
            AlleleKind::Snp => RecordKind::Snp,
            AlleleKind::Mnp => RecordKind::Mnp,
            AlleleKind::Insertion => RecordKind::Ins,
            AlleleKind::Deletion => RecordKind::Del,
            AlleleKind::Complex => RecordKind::Complex,
            // REF/Null shouldn't appear in the alt set.
            _ => RecordKind::Snp,
        })
        .collect();

    // QR / QA: sum of base qualities for observations of each allele.
    // Walk the observation list one more time matching against
    // `alleles` (REF at index 0, alts at 1..). Upstream computes these
    // in the calling context alongside the likelihood; our per-position
    // pileup doesn't thread them through the caller yet, so recompute
    // here. Cheap at realistic coverage.
    let mut qual_ref: u32 = 0;
    let mut qual_alt: Vec<u32> = vec![0; alt_obs.len()];
    for obs in effective_observations {
        if let Some(idx) = alleles.iter().position(|a| a == &obs.allele) {
            if idx == 0 {
                qual_ref = qual_ref.saturating_add(obs.base_quality_sum);
            } else {
                qual_alt[idx - 1] = qual_alt[idx - 1].saturating_add(obs.base_quality_sum);
            }
        }
    }

    // GL: log10-scaled data likelihoods in VCF spec `F(j/k)` order.
    // Our enumerate_genotypes produces non-decreasing index tuples via
    // multichoose; the resulting order agrees with VCF spec for
    // biallelic diploid but diverges for triallelic+. Use `vcf_gl_index`
    // to remap into spec order.
    let genotype_log10_likelihoods = if params.ploidy == 2 {
        let mut slots: Vec<f64> = vec![f64::NAN; gts.len()];
        let mut ok = true;
        for (g_idx, g) in gts.iter().enumerate() {
            let mut gt_indices_for_g: Vec<u8> = Vec::with_capacity(2);
            for elem in &g.elements {
                let allele_idx = alleles.iter().position(|a| a == &elem.allele).unwrap_or(0) as u8;
                for _ in 0..elem.count {
                    gt_indices_for_g.push(allele_idx);
                }
            }
            match vcf_gl_index(&gt_indices_for_g) {
                Some(slot) if slot < slots.len() => {
                    // ln -> log10: divide by ln(10).
                    slots[slot] = call.log_likelihoods[g_idx] / std::f64::consts::LN_10;
                }
                _ => {
                    ok = false;
                    break;
                }
            }
        }
        // Normalise so max = 0 (upstream convention).
        if ok {
            let max = slots.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            if max.is_finite() {
                for x in slots.iter_mut() {
                    if x.is_finite() {
                        *x -= max;
                    }
                }
            }
            Some(slots)
        } else {
            None
        }
    } else {
        None
    };

    // Per-alt CIGAR strings relative to the VCF-anchored REF/ALT pair.
    // We synthesise the anchor here too so the CIGAR matches what
    // `emit_call_vcf` will eventually emit.
    let cigars: Vec<String> = alleles
        .iter()
        .skip(1)
        .zip(alt_kinds.iter())
        .map(|(alt, kind)| {
            let (_pos, vcf_ref, vcf_alt) =
                synthesize_anchored(*kind, pos, &alt.ref_seq, &alt.alt_seq, Some(ref_base));
            alt_cigar(*kind, &vcf_ref, &vcf_alt)
        })
        .collect();

    Some(SiteCall {
        tid,
        pos,
        ref_base,
        depth: effective_total,
        ref_obs: effective_ref_obs,
        alt_obs,
        gt_indices,
        genotype_quality: call.genotype_quality,
        qual,
        qual_ref,
        qual_alt,
        alt_alleles: alleles.into_iter().skip(1).collect(),
        alt_kinds,
        cigars,
        genotype_log10_likelihoods,
    })
}

struct SiteCall {
    tid: i32,
    pos: i64,
    ref_base: u8,
    depth: u32,
    ref_obs: u32,
    alt_obs: Vec<u32>,
    gt_indices: Vec<u8>,
    genotype_quality: f64,
    qual: Option<f64>,
    qual_ref: u32,
    qual_alt: Vec<u32>,
    alt_alleles: Vec<Allele>,
    alt_kinds: Vec<RecordKind>,
    cigars: Vec<String>,
    genotype_log10_likelihoods: Option<Vec<f64>>,
}

/// M5 parallel caller — single-threaded pileup + parallel per-site
/// Bayesian call.
///
/// Phase A tried to parallelise the BAM read loop too, by giving each
/// worker its own [`bam::IndexedReader`] and a non-overlapping genomic
/// window. That segfaults with rust-htslib 0.47 / hts_sys when
/// multiple readers are opened on the same file concurrently (known
/// thread-unsafety in the index-loading path). Phase B swaps the
/// design: do all BAM I/O on the main thread into a shared
/// [`Pileup`], then use Rayon to chunk the ordered `(tid, pos)`
/// positions and score each chunk in parallel via
/// [`call_site`]. The reference is shared via `Arc<Vec<u8>>` per
/// contig — no per-worker FASTA copy.
///
/// Output is sorted by `(tid, pos)` so the VCF is byte-identical to
/// the single-threaded path at any thread count.
fn run_call_parallel(cli: &Cli) -> Result<()> {
    use rayon::prelude::*;

    let fasta = open_fasta(&cli.fasta)
        .with_context(|| format!("failed to open reference FASTA {:?}", cli.fasta))?;
    let contigs = read_bam_contigs(&cli.bam)?;
    let target_names: Vec<String> = contigs.iter().map(|c| c.name.clone()).collect();
    let target_lens: Vec<u64> = contigs.iter().map(|c| c.length).collect();

    info!(
        n_contigs = contigs.len(),
        threads = cli.threads,
        "running Bayesian caller (M5 parallel: shared pileup, per-site compute)"
    );

    // Pre-load every contig's reference into `Arc<Vec<u8>>` once.
    let ref_cache: Vec<Arc<Vec<u8>>> = contigs
        .iter()
        .enumerate()
        .map(|(tid, _)| -> Result<Arc<Vec<u8>>> {
            Ok(Arc::new(fetch_contig(
                &fasta,
                &target_names,
                &target_lens,
                tid,
            )?))
        })
        .collect::<Result<Vec<_>>>()?;

    // --- Single-threaded pileup construction ---
    let mut bam_reader = bam::Reader::from_path(&cli.bam)
        .with_context(|| format!("failed to open BAM {:?}", cli.bam))?;
    if cli.bam_threads > 0 {
        bam_reader
            .set_threads(cli.bam_threads)
            .with_context(|| format!("bam::Reader::set_threads({})", cli.bam_threads))?;
    }
    let filter = ReadFilter::default();
    let mut pileup = Pileup::default();
    let mut n_reads: u64 = 0;
    let mut n_filtered: u64 = 0;
    let mut current_tid: i32 = -1;
    let mut current_ref: &[u8] = &[];
    let mut record = bam::Record::new();
    while let Some(result) = bam_reader.read(&mut record) {
        result.context("BAM record read failed")?;
        if record.tid() < 0 {
            continue;
        }
        let tid = record.tid();
        if tid != current_tid {
            current_ref = &ref_cache[tid as usize];
            current_tid = tid;
        }
        n_reads += 1;
        let observations = match walk_record(&record, current_ref, 0, &filter) {
            Some(v) => v,
            None => {
                n_filtered += 1;
                continue;
            }
        };
        let observations = clump_observations(&observations, cli.haplotype_length);
        pileup.add_read_observations(tid, observations);
    }
    pileup.finalize();

    info!(
        n_reads = n_reads,
        n_filtered = n_filtered,
        n_variant_sites = pileup.variant_sites.len(),
        "pileup complete"
    );

    // --- Parallel per-site Bayesian call ---
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(cli.threads)
        .build()
        .map_err(|e| anyhow!("failed to build rayon pool: {e}"))?;

    let params = Parameters {
        ploidy: cli.ploidy,
        ..Parameters::default()
    };

    // Iterate only the positions that actually carry a variant
    // observation — pure-REF sites never need to be called. For each
    // such site, merge the variant observations with REF obs
    // materialised on-demand from overlapping ref runs.
    let sites: Vec<((i32, i64), &Vec<AlleleObservation>)> =
        pileup.variant_sites.iter().map(|(k, v)| (*k, v)).collect();

    let chunk_size = (sites.len() / (cli.threads.max(1) * 4))
        .max(1)
        .min(cli.tile_size as usize);

    let site_calls_unsorted: Vec<SiteCall> = pool.install(|| {
        sites
            .par_chunks(chunk_size)
            .flat_map(|chunk| -> Vec<SiteCall> {
                let mut local: Vec<SiteCall> = Vec::with_capacity(chunk.len());
                for ((tid, pos), variant_obs) in chunk {
                    let reference = &ref_cache[*tid as usize];
                    let Some(&ref_base) = reference.get(*pos as usize) else {
                        continue;
                    };
                    let combined = materialize_site_observations(&pileup, *tid, *pos, variant_obs);
                    if let Some(call) = call_site(*tid, *pos, &combined, ref_base, cli, &params) {
                        local.push(call);
                    }
                }
                local
            })
            .collect()
    });
    let mut site_calls = site_calls_unsorted;
    site_calls.sort_by_key(|c| (c.tid, c.pos));

    info!(
        n_calls = site_calls.len(),
        chunk_size = chunk_size,
        "parallel calling complete"
    );

    emit_call_vcf(
        cli,
        &fasta,
        &target_names,
        &target_lens,
        &contigs,
        &site_calls,
    )
}

fn run_call(cli: &Cli) -> Result<()> {
    if cli.threads > 1 {
        return run_call_parallel(cli);
    }
    let fasta = open_fasta(&cli.fasta)
        .with_context(|| format!("failed to open reference FASTA {:?}", cli.fasta))?;
    let mut bam_reader = bam::Reader::from_path(&cli.bam)
        .with_context(|| format!("failed to open BAM {:?}", cli.bam))?;
    if cli.bam_threads > 0 {
        bam_reader
            .set_threads(cli.bam_threads)
            .with_context(|| format!("bam::Reader::set_threads({})", cli.bam_threads))?;
    }
    let contigs = read_bam_contigs(&cli.bam)?;
    let target_names: Vec<String> = contigs.iter().map(|c| c.name.clone()).collect();
    let target_lens: Vec<u64> = contigs.iter().map(|c| c.length).collect();

    info!(
        n_contigs = contigs.len(),
        ploidy = cli.ploidy,
        min_alt_count = cli.min_alternate_count,
        min_alt_fraction = cli.min_alternate_fraction,
        "running Bayesian caller (M3 Phase C-2)"
    );

    let filter = ReadFilter::default();
    let params = Parameters {
        ploidy: cli.ploidy,
        ..Parameters::default()
    };

    let mut current_tid: i32 = -1;
    let mut current_ref: Vec<u8> = Vec::new();
    let mut pileup = Pileup::default();
    let mut n_reads: u64 = 0;
    let mut n_filtered: u64 = 0;

    // Optional fine-grained profiling of the pileup stage. Enable
    // with `FBGXY_PROFILE=1` in the environment to get a breakdown of
    // time spent in read-from-BAM vs walk vs clump vs pileup-insert.
    // Costs ~1 ns per call; off by default.
    let profile = std::env::var("FBGXY_PROFILE").is_ok();
    let mut t_read = std::time::Duration::ZERO;
    let mut t_walk = std::time::Duration::ZERO;
    let mut t_clump = std::time::Duration::ZERO;
    let mut t_insert = std::time::Duration::ZERO;

    let mut record = bam::Record::new();
    loop {
        let t0 = std::time::Instant::now();
        let step = bam_reader.read(&mut record);
        if profile {
            t_read += t0.elapsed();
        }
        let Some(result) = step else { break };
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

        let t1 = std::time::Instant::now();
        let observations = match walk_record(&record, &current_ref, 0, &filter) {
            Some(v) => v,
            None => {
                n_filtered += 1;
                if profile {
                    t_walk += t1.elapsed();
                }
                continue;
            }
        };
        if profile {
            t_walk += t1.elapsed();
        }

        let t2 = std::time::Instant::now();
        let observations = clump_observations(&observations, cli.haplotype_length);
        if profile {
            t_clump += t2.elapsed();
        }

        let t3 = std::time::Instant::now();
        pileup.add_read_observations(tid, observations);
        if profile {
            t_insert += t3.elapsed();
        }
    }

    pileup.finalize();
    info!(
        n_reads = n_reads,
        n_filtered = n_filtered,
        n_variant_sites = pileup.variant_sites.len(),
        "pileup complete"
    );
    if profile {
        let total_ms = |d: std::time::Duration| d.as_secs_f64() * 1000.0;
        info!(
            t_read_ms = total_ms(t_read),
            t_walk_ms = total_ms(t_walk),
            t_clump_ms = total_ms(t_clump),
            t_insert_ms = total_ms(t_insert),
            "pileup-stage breakdown (FBGXY_PROFILE=1)"
        );
    }

    // Iterate variant sites only; materialise REF obs on demand from
    // overlapping runs in the Phase D ref_runs_per_tid index.
    let mut site_calls: Vec<SiteCall> = Vec::new();
    let mut current_tid: i32 = -1;
    let mut current_ref: Vec<u8> = Vec::new();
    for ((tid, pos), variant_obs) in pileup.variant_sites.iter() {
        if *tid != current_tid {
            current_ref = fetch_contig(&fasta, &target_names, &target_lens, *tid as usize)?;
            current_tid = *tid;
        }
        let ref_base = match current_ref.get(*pos as usize) {
            Some(&b) => b,
            None => continue,
        };
        let combined = materialize_site_observations(&pileup, *tid, *pos, variant_obs);
        if let Some(call) = call_site(*tid, *pos, &combined, ref_base, cli, &params) {
            site_calls.push(call);
        }
    }

    info!(n_calls = site_calls.len(), "calling complete");
    emit_call_vcf(
        cli,
        &fasta,
        &target_names,
        &target_lens,
        &contigs,
        &site_calls,
    )?;
    Ok(())
}

/// Emit a freebayes-compatible VCF from the per-site calls. For
/// indel records, we look up the anchor base (at `position - 1`) from
/// the reference FASTA to synthesise the upstream VCF shape.
fn emit_call_vcf(
    cli: &Cli,
    fasta: &faidx::Reader,
    target_names: &[String],
    target_lens: &[u64],
    contigs: &[Contig],
    calls: &[SiteCall],
) -> Result<()> {
    let reference = cli.fasta.to_string_lossy().into_owned();
    let mut buf: Vec<u8> = Vec::new();
    buf.extend_from_slice(build_header(&reference, contigs, &cli.sample).as_bytes());

    // Fetch the reference lazily per tid for anchor-base lookups.
    let mut current_tid: i32 = -1;
    let mut current_ref: Vec<u8> = Vec::new();

    for call in calls {
        if call.tid != current_tid {
            current_ref = fetch_contig(fasta, target_names, target_lens, call.tid as usize)?;
            current_tid = call.tid;
        }
        let chrom = target_names
            .get(call.tid as usize)
            .map(|s| s.to_string())
            .unwrap_or_else(|| "?".to_string());

        // For each alt, synthesise the VCF-anchored (pos, ref, alt) tuple.
        // Multi-alt sites pick the first alt to define the anchor layout;
        // multiallelic INS/DEL mixing is an M4-B scenario. For Phase A
        // each call has 1 alt so the "first alt" always agrees with
        // the single alt present.
        let first_kind = call.alt_kinds.first().copied().unwrap_or(RecordKind::Snp);
        let first_alt = call
            .alt_alleles
            .first()
            .expect("call_site guarantees ≥1 alt");

        // Anchor base: the ref byte at (0-based) pos - 1, needed for
        // INS / DEL synthesis. SNPs ignore it.
        let anchor_idx = call.pos - 1;
        let anchor_base = if anchor_idx >= 0 && (anchor_idx as usize) < current_ref.len() {
            Some(current_ref[anchor_idx as usize])
        } else {
            None
        };

        let (vcf_pos, vcf_ref, vcf_alt) = synthesize_anchored(
            first_kind,
            call.pos,
            &first_alt.ref_seq,
            &first_alt.alt_seq,
            anchor_base,
        );

        let record = Record {
            chrom,
            pos: vcf_pos,
            ref_seq: vcf_ref,
            alts: vec![vcf_alt],
            qual: call.qual,
            depth: call.depth,
            ref_obs: call.ref_obs,
            alt_obs: call.alt_obs.clone(),
            gt_indices: call.gt_indices.clone(),
            gq: call.genotype_quality,
            alt_kinds: call.alt_kinds.clone(),
            cigars: call.cigars.clone(),
            qual_ref: Some(call.qual_ref),
            qual_alt: call.qual_alt.clone(),
            genotype_log10_likelihoods: call.genotype_log10_likelihoods.clone(),
        };
        buf.extend_from_slice(write_record(&record).as_bytes());

        // Explicitly drop unused bindings to keep lint clean.
        let _ = &call.ref_base;
    }
    write_output(cli.output.as_deref(), &buf)
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
