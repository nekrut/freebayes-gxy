# Performance report — M5 Phase D (lazy REF-run materialisation)

**Branch:** `claude/rewrite-freebayes-oL5RQ` @ pending commit

## Summary

**19.9× single-threaded speedup.** 3.04 s → 0.15 s on the 100 kb /
1163-variant / 19971-read benchmark fixture. Hits the upper end of
the predicted 10–20× range from the M5-C profile report. Byte-
identical output; all 129 tests pass; F1 = 1.0000 on both parity
fixtures.

## Scoreboard (100 kb benchmark fixture, release mode, 3-run avg)

| build   | threads | wall   | vs. pre-D | vs. pre-anything |
|---      |:---:    | -----: | --------: | ---------------: |
| pre-D   | 1       | 3.041 s | 1.00×    | 1.00×            |
| pre-D   | 8       | 2.995 s | 1.02×    | 1.02×            |
| post-D  | 1       | 0.153 s | **19.9×** | 19.9×            |
| post-D  | 2       | 0.148 s | 20.6×    | 20.6×            |
| post-D  | 4       | 0.143 s | 21.3×    | 21.3×            |
| post-D  | 8       | 0.134 s | **22.7×** | 22.7×            |

Threading now does something (modestly — a further 1.14× from t=1 to
t=8) because the pileup allocation cost that used to dominate is
gone.

## Pileup-stage breakdown (100 kb, single-threaded)

| stage             | pre-D   | post-D  | delta |
|---                |--------:|--------:|------:|
| BAM read          | 29 ms   | 13 ms   | 2.2× faster (better cache locality? jitter) |
| CIGAR walk        | 57 ms   | 26 ms   | 2.2× faster |
| M2 clump          | 22 ms   | 15 ms   | 1.5× faster |
| **pileup insert** | **2476 ms** | **18 ms** | **138× faster** |
| total pileup      | 2584 ms | 72 ms   | 36× faster  |

The other stages also got slightly faster — probably because the
massive heap churn was evicting their working set from L1/L2. Once
we stopped allocating 9M observations, everything else warmed up too.

## What changed

`Pileup::add_read_observations` used to decompose every REF run into
per-position single-base `AlleleObservation`s at ingest time —
allocating 3 heap objects per base (Vec for `ref_seq`, String for
`read_name.clone()`, empty Vec for `per_base_quals`) and inserting
into a `BTreeMap<(i32, i64), Vec<_>>`. On the 100 kb fixture that
was ~3M positions × ~30× coverage × 3 allocs = ~270M allocation
operations overall, and 99 % of the resulting entries were never
read because their positions had no variant observation.

Phase D splits the pileup:

```rust
struct Pileup {
    variant_sites: BTreeMap<(i32, i64), Vec<AlleleObservation>>,
    ref_runs_per_tid: HashMap<i32, TidRuns>,
}

struct TidRuns {
    runs: Vec<RefRun>,
    max_len: usize,
    sorted: bool,
}

struct RefRun { start, length, read_name, mapq, strand, ... }
```

- **At ingest:** REF observations are pushed as intervals (one per
  run); non-REF observations go into `variant_sites[pos]` as before.
- **At `finalize()`:** each contig's run list is sorted by start once
  (O(n log n) one-time cost).
- **At call time:** `run_call` iterates `variant_sites` (≈1 % of
  reference positions), and for each site queries
  `ref_runs_per_tid[tid].covering(pos)` via a bounded binary search
  (`partition_point` on start, windowed by `max_len`). Each
  covering run is materialised into a single-base REF
  `AlleleObservation` on demand via `RefRun::materialize_at(pos)`.

The materialised observation is semantically **identical** to what
the old decomposer produced — same read_name / mapq / strand /
read_position / is_proper_pair / read_ref_start / per-position BQ.
That preserves the drift-2 event-span filter, the per-read dedupe,
QR/QA accounting, and every other M3/M4 invariant.

## Correctness

- `cargo test --workspace` — 129/129 pass.
- `cargo clippy --workspace --all-targets -- -D warnings` — clean.
- `cargo fmt --all --check` — clean.
- Small fixture (2 kb, 20 variants): VCF byte-identical before/after.
  `bcftools norm` diff: 0 lines.
- Large fixture (10 kb, 113 variants): F1 = 1.0000 unchanged.
- Benchmark fixture (100 kb, 1163 variants): F1 = 1.0000 with lazy
  REF materialisation; VCF single-threaded identical to
  multi-threaded output.

## Implementation notes

- `TidRuns::covering` uses `max_len`-bounded windowing so the
  backward-scan cost is O(depth) at each site — ~30 runs at 30×
  coverage, not O(total runs).
- Materialisation cost: one `Vec<u8>` (1 byte) + one `String` clone
  per overlapping run per candidate site. With ~1200 candidate
  sites × ~30 runs = ~36K materialisations, vs. the ~3M
  materialisations in the old decomposer. 83× fewer allocations and
  only at sites we actually call.
- `add_read_observations` takes observations by value so `ref_seq`,
  `read_name`, and `per_base_quals` move into the `RefRun` without
  cloning at ingest.

## Threading, in light of Phase D

At t=1 we're at 150 ms wall, with a breakdown roughly:
- BAM read + walk + clump: ~55 ms
- Pileup insert (now cheap): ~18 ms
- Sort + finalise: tiny
- Per-site call (over 1163 sites): ~50 ms (including REF
  materialisation)
- VCF emit / misc: ~30 ms

The per-site call stage can still be parallelised (already is via
Phase B's design) but the overall wins are capped by the other
stages. On this fixture t=8 → 134 ms = 1.14× further speedup.

Realistic chromosome-scale expectation: single-threaded stays
bottlenecked on `bam::Reader::read` + walk (which we can't
parallelise without the hts_sys workarounds from the M5-B report).
But 1/20th the wall-clock for the compute is already a massive
practical win.

## What's still open

- **`noodles` migration** — no longer needed for speedup. Might
  still be worth it for thread-safety if chromosome-scale
  performance is needed.
- **Chromosome-scale benchmark** — the 100 kb fixture is the
  biggest we have; a real chr-scale test (GIAB HG002 chr22) would
  show how the curve extrapolates.
- **Per-site call stage optimisation** — ~50 ms for 1163 calls = 43
  µs/call. Genotype enumeration + posterior is small; the likely
  hot spot is the `materialize_at` loop. Profiling can confirm.
- **VCF emission** — probably trivial at this scale but worth
  checking.

## Moral of the story

The bottleneck was never where the PLAN (or I) thought it was. A
single pass of flame-graph-style profiling (the env-var-gated
`FBGXY_PROFILE=1` probe from Phase C) exposed the real cost and
steered us away from spending 4–6 hours on a `noodles` migration
or complex producer-consumer pipelines, toward a 2-hour data-
structure refactor that delivered 20× speedup.

**Profile before optimising.**
