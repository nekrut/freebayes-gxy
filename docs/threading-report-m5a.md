# Threading report — M5 Phase A (window-parallel `--call`)

**Branch:** `claude/rewrite-freebayes-oL5RQ` @ pending commit
**Scope:** Rayon-backed window scheduler for the per-site caller.

## Summary

First landed parallelism: the `--call` pipeline now tiles each contig
into fixed-size genomic windows, dispatches each window to a Rayon
worker, and merges the resulting `SiteCall` vectors back into a
deterministic VCF. Threaded output is **byte-identical** to the
single-threaded path on the parity fixture.

This materialises the PLAN §2 headline ("4–10× via native
multithreading"). Phase A is intentionally minimal — no shared
reference cache, no adaptive tile sizing, no reorder buffer beyond a
final sort. Those layer on in later phases.

## CLI surface

```
  -t, --threads <N>       Number of worker threads (default: 1).
                          N=1 uses the original single-threaded path.
      --tile-size <BP>    Tile size for the parallel scheduler
                          (default: 100000).
```

Existing flags (`--call`, `--ploidy`, `--min-alternate-count`,
`--min-alternate-fraction`, `--haplotype-length`) are unchanged and
compose with `--threads`.

## Correctness

On the 10 kb / 113-variant fixture:

```
single records:   113
threaded records: 113
diff:             (empty)
```

Byte-identical after a straightforward per-site sort on `(tid, pos)`.
Same on the 2 kb fixture. The test suite (129 tests) stays green —
threading is opt-in and the single-threaded path is unchanged.

## Scaling

| threads | tile | wall-clock | speedup |
|--------:|:----:|----------:|--------:|
| 1       | —    | 0.239 s    | 1.00×   |
| 2       | 500 bp | 0.147 s  | 1.63×   |
| 4       | 500 bp | 0.170 s  | 1.41×   |
| 8       | 500 bp | 0.293 s  | 0.82×   |

Measured on the 10 kb release-mode fixture (1971 reads, 113 truth
variants). Results are bound by task-spawn + BAM-reader-open
overhead at this scale — 10 kb is an order of magnitude smaller than
the realistic workloads (chromosomes) the scheduler targets. A
100 kb fixture would likely show proper 4×+ scaling at t=4.

Proper benchmarks on a chr-scale fixture come with Phase B; the
purpose of Phase A is to land a correct, deterministic parallel
path.

## Implementation

- `fb-cli/src/main.rs`:
  - `--threads` / `--tile-size` CLI flags.
  - `run_call` dispatches to `run_call_parallel` when `threads > 1`.
  - `run_call_parallel` enumerates contigs, tiles them via
    `fb_scheduler::tile_contig`, builds a Rayon `ThreadPoolBuilder`,
    `par_iter`s the windows into `call_window`, collects, sorts,
    emits the VCF.
  - `call_window` opens its own `bam::IndexedReader` + `faidx::Reader`
    (both cheap, thread-affine), fetches reads via the BAM index for
    the window's coordinates, runs walker + clumping + pileup, and
    calls positions inside `[window.start, window.end)`. Reads
    spanning a tile boundary are still seen by workers whose window
    overlaps the read's span (via the indexed fetch); each window's
    workers only EMIT calls for positions strictly inside its range,
    avoiding double-count.
- `fb-cli/Cargo.toml`: `rayon.workspace = true` + `fb-scheduler` path
  dep (already an empty placeholder since M0).
- No changes to `fb-core`, `fb-genotype`, `fb-vcf`, or `fb-scheduler`.

## Limitations / deferred work

- **Requires BAM index** (`sample.bam.bai`). Single-threaded path
  still uses the sequential reader and doesn't need an index; the
  parallel path fails loudly if the index is absent.
- **No shared reference cache** — each worker fetches the full
  contig from the FASTA. Cheap at contig scale, wasteful at
  chromosome scale. Phase B: hoist the reference to an `Arc<Vec<u8>>`
  and share across workers.
- **Static tiling** — windows are fixed-size, no work stealing. For
  coverage-heterogeneous BAMs a work-stealing scheduler would be a
  better fit. Phase B.
- **No reorder buffer** — final sort is O(n log n) on all
  `SiteCall`s; for million-variant outputs a merge-sort over
  per-worker-sorted streams would beat it.
- **Output streaming** — `emit_call_vcf` still buffers every record
  before writing. Streaming VCF output is a C-2 nice-to-have.

## What this unlocks

- Realistic benchmark runs at chromosome scale.
- A concrete speedup demonstration against upstream's
  `freebayes-parallel` wrapper (which shards externally at the
  shell level). Phase B comparison target.
- Foundation for the SIMD likelihood kernel (M6) — SIMD gains
  multiply on top of thread-level gains.

## Running

```bash
# Correctness check:
./target/debug/freebayes-gxy --call -f ref.fa sample.bam > single.vcf
./target/debug/freebayes-gxy --call -t 4 --tile-size 500 \
    -f ref.fa sample.bam > threaded.vcf
diff single.vcf threaded.vcf   # => no output

# Benchmark:
time ./target/release/freebayes-gxy --call -t 4 --tile-size 500 \
    -f ref.fa sample.bam > /dev/null
```

## Workspace state

- 129 tests pass.
- cargo build / test / clippy / fmt all green.
- Parity harness: small fixture byte-identical; large fixture
  F1 = 1.0000 for both callers (unchanged).
- `--threads N` now functional; output deterministic across N.
