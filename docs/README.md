# freebayes-gxy — project dossier

A from-scratch Rust rewrite of [freebayes](https://github.com/freebayes/freebayes)
targeting byte-level VCF parity with upstream **v1.3.10** (with
`--legacy-gls`) while delivering native multithreading and a clean
library API for Galaxy / nf-core.

This document synthesises the 14 phase reports under `docs/` into a
single thread of the story so far. For the canonical design, see
[`../PLAN.md`](../PLAN.md).

## Current status (2026-04-19)

| Axis | Result |
|---|---|
| **Records emitted** | Identical to upstream on every parity fixture (100% record parity). |
| **F1 vs synthetic truth** | **1.0000** on 2 kb / 10 kb / 100 kb / 1 Mb fixtures. |
| **GL numeric parity vs upstream `--legacy-gls`** | **Byte-identical** across all SNP sites. |
| **Wall-clock (100 kb, single-threaded)** | **8.6×** faster than upstream serial, **7.0×** faster than `freebayes-parallel -j 8` (post-M5-H). |
| **Wall-clock (1 Mb, single-threaded)** | **12.7×** faster than upstream serial, **3.4×** faster than `freebayes-parallel -j 8` (post-M5-H). |
| **Wall-clock (1 Mb, -t 4 bt=4)** | **17.3×** upstream serial, **4.7×** `freebayes-parallel -j 8` (post-M5-H). |
| **Peak RSS (1 Mb)** | ~230 MB (**16×** upstream, 40× `freebayes-parallel` workers). |
| **Test suite** | 129/129 pass; `clippy`, `fmt`, `cargo doc` clean. |

## Crate architecture

```
crates/
├── fb-core        — allele observations, walker, pileup, haplotype clumping,
│                    lazy REF-run materialisation
├── fb-genotype    — genotype enumeration, data likelihood, prior, single-sample
│                    caller (standardGLs path), GQ
├── fb-vcf        — record writer: anchor-base synthesis, CIGAR, QR, QA, GL
├── fb-scheduler   — rayon-driven window scheduler (used by the --call path)
└── fb-cli         — `freebayes-gxy` binary; wires the others together
```

The call path is roughly:

1. `fb-cli` opens the reference + BAM, builds a `Pileup` by walking reads
   once (single-threaded — rust-htslib's `bam::Reader` is not
   multi-reader safe, see M5-A / M5-C).
2. At ingest: SNP/INS/DEL observations land in
   `Pileup::variant_sites: BTreeMap<(tid, pos), Vec<AlleleObservation>>`.
   REF observations are stored as interval runs
   (`ref_runs_per_tid: HashMap<tid, TidRuns>`) — **not** decomposed into
   per-position observations.
3. `finalize()` sorts each contig's run list by start once.
4. `run_call_parallel` chunks `variant_sites` via `rayon::par_chunks`
   and, per candidate site, materialises covering REF observations
   on demand via a bounded `partition_point` binary search.
5. Per site, `fb-genotype::Caller` computes the Bayesian posterior and
   emits a record through `fb-vcf`.

See [`performance-report-m5d.md`](performance-report-m5d.md) for the
pileup data structure and
[`threading-report-m5b.md`](threading-report-m5b.md) for the threading
model.

## How we got here — parity

The M3–M4 arc closed three distinct "drifts" between gxy and upstream:

| Drift | Symptom | Fix | Report |
|---|---|---|---|
| 1 | Indel representation lacked anchor bases. | `synthesize_anchored` in the VCF writer. | [`parity-report-m4a.md`](parity-report-m4a.md) |
| 2 | Hom→het miscalls on indels from over-eager REF pileup. | Event-span admission + per-read dedupe at indel sites. | [`parity-report-m4d.md`](parity-report-m4d.md) |
| 3 | Harness mis-mapped `GT` tag, inflating FN rate. | `bcftools norm`-based harness + tag-aware comparator. | [`parity-report-m4a.md`](parity-report-m4a.md), [`parity-report-m4c.md`](parity-report-m4c.md) |

After M4-E (per-position REF BQ — QR / QA now track upstream exactly)
and M4-F (upstream run with `--legacy-gls` — matches the `standardGLs`
formula we ported), GL numeric parity became byte-identical. M4-G
deconflicted overlapping truth positions so hap.py scores cleanly on
the 10 kb fixture; F1 = 1.0000 for both callers.

All parity reports: `parity-report-m3.md` through
`parity-report-m4g.md`.

## How we got here — performance

| Phase | Change | Result |
|---|---|---|
| M5-A | Parallelise the I/O + pileup stage across rayon workers. | Segfaults at t ≥ 4 — rust-htslib's multi-reader race. |
| M5-B | Pivot: single-threaded pileup, parallel call stage, shared `Arc<Vec<u8>>` reference. | Correct but flat scaling (the call stage is only ~14% of runtime). |
| M5-C (experiment) | Hypothesis: sequential index load → multi-thread fetch would dodge the race. | Falsified — the race is deeper than `from_path`. |
| M5-C (profile) | Env-gated `FBGXY_PROFILE=1` probe. | Revealed 95.7% of runtime was `Pileup::add_read_observations` decomposing REF runs into ~9M heap observations. |
| **M5-D** | Lazy REF-run materialisation: keep REF obs as intervals, materialise at call sites only. | **19.9× single-threaded speedup.** |
| M5-E | Head-to-head on 100 kb vs upstream + `freebayes-parallel`. | gxy-t1 = 4.7× upstream, 3.9× `freebayes-parallel -j 8`. |
| M5-F | Chromosome-scale (1 Mb) benchmark. | gxy-t1 = 7.9× upstream serial, 2.1× `freebayes-parallel -j 8`. Thread scaling saturates (bam::Reader::read + walk is now the wall). 16× RSS overhead flagged. |
| M5-G | Optional `--bam-threads N` (htslib BGZF decoder pool). | 2.4× on the read stage but only ~13% of wall, so no visible speedup on this fixture. Flag kept as a no-cost knob for network/dense BAMs. |
| M5-H | Producer/consumer pipeline attempted, found glibc-malloc pathological on cross-thread free. Swapped to `mimalloc` — delivered 1.5–1.8× wall reduction with 3 lines of code instead of 100. |

Reports: `threading-report-m5{a,b,c-experiment}.md`,
`profile-report-m5c.md`, `performance-report-m5d.md`,
`bench-report-m5{e,f}.md`, `perf-report-m5{g,h}.md`.

## Moral

> The bottleneck was never where the PLAN (or I) thought it was. A
> single pass of flame-graph-style profiling (the env-var-gated
> `FBGXY_PROFILE=1` probe from Phase C) exposed the real cost and
> steered us away from spending 4–6 hours on a `noodles` migration or
> complex producer-consumer pipelines, toward a 2-hour data-structure
> refactor that delivered 20× speedup. **Profile before optimising.**
> — [`performance-report-m5d.md`](performance-report-m5d.md)

## Build + run

```bash
# Workspace build
cargo build --release --workspace

# Parity / call
./target/release/freebayes-gxy --call \
    --threads 8 \
    -f reference.fa \
    sample.bam > out.vcf

# Single-sample only in this port (upstream multi-sample is out of scope for now).
```

### Parity harness

```bash
cd tests/parity/fixture-large    # or fixture-bench / fixture-chrscale
python3 build_fixture.py         # builds ref.fa, sample.bam, truth.vcf.gz
bash ../harness-large.sh         # runs both callers, bcftools-norm-diffs, reports F1
```

The harness passes `--legacy-gls` to upstream for apples-to-apples GL
comparison (default upstream uses the experimental GL path).

### Benchmarking

```bash
export PATH=/usr/lib/vcflib/bin:/path/to/freebayes-upstream/scripts:/path/to/freebayes-upstream/build:$PATH
python3 tests/parity/fixture-chrscale/build_fixture.py
cd tests/parity/fixture-chrscale
fasta_generate_regions.py ref.fa.fai 100000 > regions.100k.txt

time freebayes --legacy-gls -f ref.fa sample.bam > /dev/null
time freebayes-parallel regions.100k.txt 8 --legacy-gls -f ref.fa sample.bam > /dev/null
time ./target/release/freebayes-gxy --call --threads 8 -f ref.fa sample.bam > /dev/null
```

## Is gxy the right choice?

From [`bench-report-m5f.md`](bench-report-m5f.md):

| Scenario | Winner |
|---|---|
| Short contigs (< 500 kb), low thread count | **gxy** (≫ 2× faster even single-threaded) |
| 1–10 Mb contigs, any thread count | **gxy** (still 2× faster, lower wall) |
| Multi-chromosome (> 100 Mb), ≥ 8 cores | `freebayes-parallel` likely wins on wall |
| Memory-constrained host | **`freebayes-parallel`** (16–40× less RSS) |
| Accuracy / F1 | Tied at 1.0000 on every parity fixture |
| VCF byte parity (with `--legacy-gls`) | Identical across every tool |

The 16× memory cost is the big asterisk. For Galaxy / nf-core /
clinical deployments where RAM is bounded, gxy's current design is
unsuitable at full-WGS scale without streaming pileup.

## Open items / roadmap

In descending order of impact:

1. **Streaming pileup.** The single biggest win available. Process
   reads in coordinate order, flush completed sites once no future
   read can still contribute, free the REF runs and variant
   observations behind them. Projected: chromosome-scale RSS from
   ~50 GB down to ~500 MB. Would also unlock concurrent callers over
   a single BAM cursor.
2. **`noodles` migration for thread-safe BAM I/O.** Pure Rust, no
   hts_sys race. Needed if chromosome-scale single-threaded
   `bam::Reader::read` is the bottleneck (it is, per M5-F). Scoped in
   M5-C-experiment; not yet implemented.
3. **`rust-htslib::set_threads(n)`.** htslib's internal BGZF-decoder
   threading is thread-safe within a single reader. May unlock
   parallel decompression of the read stream without the multi-reader
   race. Low-hanging.
4. **Per-thread chunk tuning.** At 1 Mb the par_chunks size is ~3k
   sites/thread, which may be too coarse after M5-D. Worth probing.
5. **Multi-sample.** Currently single-sample only. Upstream's joint
   calling + `--populations` prior are out of scope for M4/M5 but
   blocker for real pop-gen use.
6. **Structural variants, complex events.** Upstream's haplotype
   reconstruction over long windows is not fully ported (we do M2
   per-read clumping but not cross-read window haplotypes).
7. **Full-WGS parity regression.** hap.py at chr-scale to catch any
   drift hidden behind the synthetic fixtures. M5-F only checked
   record counts at 1 Mb.

## Phase report index

| Phase | Topic | File |
|---|---|---|
| M3 | Initial parity scaffolding (12/12 identical SNP sites). | `parity-report-m3.md` |
| M4-A | VCF emission + anchor-base synthesis (drift 1 closed). | `parity-report-m4a.md` |
| M4-B | GL + CIGAR + QR + QA (VCF schema sufficient for hap.py). | `parity-report-m4b.md` |
| M4-C | hap.py-style F1 scoring on 10 kb fixture. | `parity-report-m4c.md` |
| M4-D | Event-span pileup + per-read dedupe (drift 2 closed). | `parity-report-m4d.md` |
| M4-E | Per-position REF BQ — QR/QA match upstream. | `parity-report-m4e.md` |
| M4-F | GL byte-parity via `--legacy-gls`. | `parity-report-m4f.md` |
| M4-G | Deconflict truth positions; F1 = 1.0000. | `parity-report-m4g.md` |
| M5-A | First threading attempt — segfaults on parallel I/O. | `threading-report-m5a.md` |
| M5-B | Shared ref cache + 100 kb fixture + honest scaling. | `threading-report-m5b.md` |
| M5-C-x | Experiment: sequential index load hypothesis (falsified). | `threading-report-m5c-experiment.md` |
| M5-C-p | Profile — pileup decomposition is 95.7% of wall. | `profile-report-m5c.md` |
| M5-D | Lazy REF-run materialisation — 19.9× speedup. | `performance-report-m5d.md` |
| M5-E | Head-to-head vs upstream + `freebayes-parallel` (100 kb). | `bench-report-m5e.md` |
| M5-F | Chromosome-scale (1 Mb) benchmark + memory caveat. | `bench-report-m5f.md` |
| M5-G | `--bam-threads` (htslib BGZF decoder pool) — works, bounded. | `perf-report-m5g.md` |
| M5-H | mimalloc global allocator — 1.5–1.8× wall, no code-flow change. | `perf-report-m5h.md` |

## Parity fixtures

| Fixture | Size | Variants | Reads | Purpose |
|---|---|---|---|---|
| `tests/parity/fixture-small` | 2 kb | 20 | 200 | Fast sanity + golden VCF byte-diff. |
| `tests/parity/fixture-large` | 10 kb | 113 | 2k | F1 regression against hap.py-style truth. |
| `tests/parity/fixture-bench` | 100 kb | 1163 | ~20k | Head-to-head perf benchmark. |
| `tests/parity/fixture-chrscale` | 1 Mb | 11663 | ~200k | Chromosome-scale scaling + memory. |

All fixtures are deterministically reproducible
(`random.seed(42)`).

## License

MIT — see [`../LICENSE`](../LICENSE).
