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

/// Per-position pileup of allele observations.
///
/// Reference runs from the M1 CIGAR walker cover many positions with a
/// single observation; this aggregator **decomposes** each ref run into
/// per-position single-base Reference observations so the genotype
/// model can count per-site depth correctly. Non-reference observations
/// (SNP, INS, DEL, Complex) are bucketed at their anchor position
/// as-is.
///
/// **Quality scalar caveat:** decomposed per-position Reference
/// observations inherit the parent run's `base_quality_sum`, which the
/// M1 walker sets to the read's `mapq` (upstream's convention — see
/// `fb_core::pileup::walk_match_run`). That is accurate for
/// single-observation math but does not recover per-base BQ for the
/// run's interior. Fixing this is a Phase C-3 / M4 task: either emit
/// per-position Reference observations from the walker directly, or
/// cache per-base BQs on the run and index them here.
///
/// Memory cost: O(total aligned bases). Fine for synthetic test inputs
/// and small regions; a streaming or interval-tree backed aggregator is
/// a Phase C-3 optimisation for full-WGS runs.
#[derive(Default)]
struct Pileup {
    positions: BTreeMap<(i32, i64), Vec<AlleleObservation>>,
}

impl Pileup {
    fn add_read_observations(&mut self, tid: i32, observations: Vec<AlleleObservation>) {
        for obs in observations {
            match obs.allele.kind {
                AlleleKind::Reference => {
                    // Decompose the ref run into per-position
                    // single-base Reference observations. The run's
                    // `ref_seq` is the aligned read bases (which equal
                    // the ref bases by construction — see `Allele::
                    // reference` in fb-core).
                    for offset in 0..obs.allele.length {
                        let pos = obs.allele.position + offset as i64;
                        let byte = obs.allele.ref_seq[offset];
                        // Per-position BQ: prefer the walker's cached
                        // per_base_quals slice (real Phred-scaled BQ),
                        // fall back to the obs-level scalar if somehow
                        // missing (e.g. during test construction).
                        let per_pos_bq: u32 = obs
                            .per_base_quals
                            .get(offset)
                            .copied()
                            .map(u32::from)
                            .unwrap_or(obs.base_quality_sum);
                        let per_pos = AlleleObservation {
                            allele: Allele::reference(pos, vec![byte]),
                            read_name: obs.read_name.clone(),
                            mapq: obs.mapq,
                            base_quality_sum: per_pos_bq,
                            strand: obs.strand,
                            read_position: obs.read_position + offset,
                            is_proper_pair: obs.is_proper_pair,
                            // Preserve the parent read's alignment start so
                            // call_site can filter REF observations at indel
                            // sites (drift 2 fix).
                            read_ref_start: obs.read_ref_start,
                            // Per-position refs have a single BQ which we
                            // store in `base_quality_sum`; leave the vec
                            // empty so downstream code uses the scalar.
                            per_base_quals: Vec::new(),
                        };
                        self.positions.entry((tid, pos)).or_default().push(per_pos);
                    }
                }
                AlleleKind::Null => {
                    // Soft-clip / N-base observations don't inform the
                    // Bayesian model; drop them at the pileup stage.
                }
                _ => {
                    self.positions
                        .entry((tid, obs.allele.position))
                        .or_default()
                        .push(obs);
                }
            }
        }
    }
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

fn run_call(cli: &Cli) -> Result<()> {
    let fasta = open_fasta(&cli.fasta)
        .with_context(|| format!("failed to open reference FASTA {:?}", cli.fasta))?;
    let mut bam_reader = bam::Reader::from_path(&cli.bam)
        .with_context(|| format!("failed to open BAM {:?}", cli.bam))?;
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
        let observations = match walk_record(&record, &current_ref, 0, &filter) {
            Some(v) => v,
            None => {
                n_filtered += 1;
                continue;
            }
        };
        let observations = clump_observations(&observations, cli.haplotype_length);
        pileup.add_read_observations(tid, observations);
    }

    info!(
        n_reads = n_reads,
        n_filtered = n_filtered,
        n_positions = pileup.positions.len(),
        "pileup complete"
    );

    // Re-fetch ref per tid as we iterate — same pattern as the dump
    // path. Positions are in (tid, pos) order thanks to the BTreeMap.
    let mut site_calls: Vec<SiteCall> = Vec::new();
    let mut current_tid: i32 = -1;
    let mut current_ref: Vec<u8> = Vec::new();
    for ((tid, pos), observations) in pileup.positions.iter() {
        if *tid != current_tid {
            current_ref = fetch_contig(&fasta, &target_names, &target_lens, *tid as usize)?;
            current_tid = *tid;
        }
        let ref_base = match current_ref.get(*pos as usize) {
            Some(&b) => b,
            None => continue, // out-of-bounds (shouldn't happen post-walker)
        };
        if let Some(call) = call_site(*tid, *pos, observations, ref_base, cli, &params) {
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
