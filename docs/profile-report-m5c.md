# Profile report — M5 Phase C (where is the pileup time actually going?)

**Branch:** `claude/rewrite-freebayes-oL5RQ` @ pending commit
**Hot take:** The bottleneck is neither I/O nor the CIGAR walker. It's
our own [`Pileup::add_read_observations`] allocating ~9M heap objects
decomposing Reference runs into per-position observations.

## Methodology

Added four `std::time::Instant` probes to the single-threaded
`run_call` pileup loop, gated behind `FBGXY_PROFILE=1`:

- `t_read`: the `bam::Reader::read(&mut record)` call (BAM I/O +
  bgzf decompress)
- `t_walk`: `walk_record` (CIGAR walker → observations)
- `t_clump`: `clump_observations` (M2 haplotype clumping)
- `t_insert`: `Pileup::add_read_observations` (per-position REF
  decomposition + insertion into `BTreeMap<(tid, pos), Vec<Obs>>`)

Ran on the 100 kb / 19971-read benchmark fixture, release mode,
3 runs.

## Results (per-stage ms, 3-run average)

| stage | time | % of pileup | % of total |
|---|---:|---:|---:|
| t_read | 28.6 ms | 1.1% | 0.9% |
| t_walk | 57.4 ms | 2.2% | 1.9% |
| t_clump | 22.1 ms | 0.9% | 0.7% |
| **t_insert** | **2475.9 ms** | **95.7%** | **82.3%** |
| (pileup total) | 2584.0 ms | 100% | 85.9% |
| call stage (single-threaded) | ~91 ms | — | 3.0% |
| misc | ~335 ms | — | 11.1% |
| **total wall** | **3010 ms** | | **100%** |

## Interpretation

- **BAM I/O is 1% of the job.** `rust-htslib` + `hts_sys` aren't the
  bottleneck. Switching to `noodles` would buy us nothing here (and
  possibly cost us, since its BGZF decoder is typically slower than
  htslib's tuned C path).
- **CIGAR walking + clumping together are 3%.** Even if option #2
  (single-reader producer + parallel-consumer walk) gave perfect
  scaling, the ceiling is a 3% slice of runtime.
- **96% of the pileup stage is our own data-structure hot loop** —
  `Pileup::add_read_observations` iterating over every base of every
  Reference run and allocating a new `AlleleObservation` at each
  position.

## What's actually happening

Looking at the code:

```rust
AlleleKind::Reference => {
    for offset in 0..obs.allele.length {
        let pos = obs.allele.position + offset as i64;
        let byte = obs.allele.ref_seq[offset];
        let per_pos = AlleleObservation {
            allele: Allele::reference(pos, vec![byte]),      // Vec alloc
            read_name: obs.read_name.clone(),                // String alloc
            ...
            per_base_quals: Vec::new(),                      // Vec alloc
        };
        self.positions.entry((tid, pos)).or_default().push(per_pos);
        // ^ BTreeMap entry + Vec push
    }
}
```

For each read's Reference run of ~100–150 bp, we allocate 3 heap
objects per position × 150 positions × 19971 reads ≈ **9 million
allocations** — and each `BTreeMap::entry` does a tree walk.

On top of the allocation cost, the resulting pileup carries
`~3 million` full-fat `AlleleObservation`s (each ~100 bytes) — and
**99% of them are never consumed** because `call_site` short-circuits
on positions with no variant candidate (`candidates.is_empty()`).

## Why this went unnoticed

Development was on the 2 kb / 20-variant fixture where this cost was
~50 ms — not visible. Moving to the 100 kb fixture exposed the
quadratic-ish scaling of the per-position decomposition.

## The fix (M5 Phase D plan)

Replace the eagerly-decomposed `BTreeMap<(tid, pos), Vec<Obs>>` with
a split pileup:

```rust
struct Pileup {
    // Only positions that have ≥1 non-Reference observation. Each
    // entry holds the variant observations at that position.
    variant_sites: BTreeMap<(i32, i64), Vec<AlleleObservation>>,
    // Reference runs stored as intervals, grouped per tid, sorted by
    // start. At call time we binary-search for the runs covering a
    // given position.
    ref_runs_per_tid: HashMap<i32, Vec<RefRun>>,
}
```

At `call_site` we only visit positions in `variant_sites` (≈1% of
the ref length). For each such position, we query
`ref_runs_per_tid[tid]` for overlapping runs via binary search +
linear scan. We reconstruct REF observations on demand from at most
~30 covering runs (at 30× coverage with 150 bp reads staggered 5 bp).

**Expected win:** 10–20× single-threaded speedup on the 100 kb
fixture. The call stage stays the same (already 91 ms); the 2475 ms
collapse since we stop allocating per-position REF observations.

## Threading, in light of this

- Single-threaded wall after Phase D: estimated 200–400 ms.
- Call stage (currently 91 ms, already parallelised): drops to
  ~50 ms with Phase B parallelism.
- Any additional I/O parallelism: 29 ms × N could save ~20 ms.
- Realistic wall-clock: **~150 ms total, 20× from today without
  changing BAM libs.**

Threading was the wrong first optimization. The right first
optimization is to stop decomposing REF runs that nobody reads.

## Action items

1. **Phase D: pileup data-structure redesign.** Keep REF runs as
   intervals; materialise per-position REF observations lazily at
   call-site. Preserves all M3 / M4 semantics (event-span filter,
   per-read dedupe, QR/QA) because we reconstruct on demand.
2. **Revisit threading after Phase D.** With the 2.5 s allocation
   cost gone, the call stage goes from 3% to ~25% of runtime —
   parallelising it has real payoff.
3. **Ship the profiling probe.** `FBGXY_PROFILE=1` is useful — land
   it behind the env var so future profiling is a one-liner.

## Workspace state

- 129 tests pass.
- cargo build / test / clippy / fmt all green.
- Profile instrumentation is env-var gated; default behaviour
  unchanged.
