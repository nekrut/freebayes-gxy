# Threading report — M5 Phase B (pileup-bound + honest speedup measurement)

**Branch:** `claude/rewrite-freebayes-oL5RQ` @ pending commit

## Summary

Phase B changes the parallel strategy from Phase A's per-window
parallelism (which segfaulted on 100 kb workloads due to `hts_sys`
thread-unsafety) to a **single-threaded pileup + parallel per-site
Bayesian call** design. This is correct and stable, but profiling on
the new 100 kb / 1163-variant fixture reveals a hard truth:

> **The per-site Bayesian call is only ~3% of runtime.** The pileup
> stage (BAM read + CIGAR walk + clump + observation insert) takes
> ~2.6 s out of 3.0 s total. No amount of compute parallelism can
> speed up work that isn't happening.

So Phase B is correct, adds a shared reference cache (nice), and
documents the real performance picture. Phase C has to either make
BAM I/O actually parallel or trade memory for a different pipeline
shape.

## What Phase B delivers

- **Shared reference cache.** Every contig is pre-loaded into
  `Arc<Vec<u8>>` and shared across workers instead of re-fetched per
  tile. Cheap on the current fixture; meaningful on chromosome
  scale.
- **100 kb / 1163-variant benchmark fixture** at
  `tests/parity/fixture-bench/` (10× the `fixture-large/` fixture).
  19 971 single-end 150 bp reads, ~30× coverage, Q40 / MQ60. Same
  deconflicted-residues truth layout as `fixture-large/`.
- **Parallel per-site call** via Rayon `par_chunks` over the
  pileup's ordered position list. 4 chunks per thread for
  rebalancing. Byte-identical output vs. single-threaded at every N.
- **Correctness vs. Phase A.** Phase A segfaulted at t ≥ 4 on the
  100 kb fixture (multiple concurrent `IndexedReader::from_path`
  opens race in hts_sys). Phase B routes all BAM I/O through one
  `bam::Reader` on the main thread and avoids the race.

## Scaling (100 kb fixture, release mode, 3-run avg)

| threads | real | user | user / real | speedup |
|--------:|-----:|-----:|------------:|--------:|
|       1 | 3.04 s | 2.47 s | 0.81× | 1.00× |
|       2 | 3.01 s | — | — | 1.01× |
|       4 | 2.96 s | — | — | 1.03× |
|       8 | 3.00 s | — | — | 1.01× |

**Speedup is flat.** Breakdown on t=1:

| stage | wall |
|---|---|
| BAM read + walk + clump + pileup insert | 2.61 s (86%) |
| Parallel per-site call  | 0.09 s (3%) |
| Misc (config, emit, flush) | 0.34 s (11%) |

On t=4 the call stage drops to ~48 ms (≈2× speedup on the 3% slice)
— real, but invisible at the program level.

## Root cause: hts_sys thread-unsafety

Phase A's per-window `IndexedReader` design is the straightforward
parallel-I/O approach: each worker opens an indexed reader on the
same BAM and fetches its region. **This segfaults with rust-htslib
0.47 / hts_sys when multiple readers open the same file
concurrently.** Reproduced:

```
$ ./target/release/freebayes-gxy --call --threads 4 --tile-size 10000 \
    -f ref.fa sample.bam
Segmentation fault (core dumped)
```

The crash is deterministic; at tile sizes large enough that only one
worker spawns (`--tile-size 100000` = the contig length), t=4
completes successfully. So the issue is multi-reader concurrency,
not anything in our code.

Workarounds that don't involve waiting on a fixed hts_sys:

1. **Sequential index load, parallel record reads** — open all
   readers on the main thread (one per worker), `.fetch` them
   sequentially, then hand each worker its own reader via rayon's
   `install`. Depends on whether the race is at `from_path` (likely)
   or at `fetch` time. Untested.
2. **Single-reader producer + parallel-consumer walk** — one thread
   reads records from the BAM into a channel, N workers consume
   records, run walk + clump, and emit observations back to a
   shared pileup via a lock-free map. Non-trivial.
3. **Re-compress to indexed BAM tiles offline** — pre-shard the BAM
   into per-contig / per-window files, open each with its own
   reader. Expensive for large BAMs.
4. **Use a different htslib crate / version** — `noodles` is pure
   Rust and multi-threadable by construction. Larger dep change.

## CLI

```
  -t, --threads <N>       Worker threads for the per-site call
                          stage (default: 1).
      --tile-size <BP>    Tile-size cap on the per-worker chunk.
                          Default 100000.
```

With Phase B's pipeline, `--tile-size` bounds the number of
positions a single worker handles in one chunk; it no longer shards
the BAM.

## Running

```bash
# Build 100kb fixture:
python3 tests/parity/fixture-bench/build_fixture.py

# Correctness:
./target/release/freebayes-gxy --call -f ref.fa sample.bam > a.vcf
./target/release/freebayes-gxy --call -t 4 -f ref.fa sample.bam > b.vcf
diff a.vcf b.vcf   # => no output

# Benchmark:
for t in 1 2 4 8; do
  { time ./target/release/freebayes-gxy --call -t $t -f ref.fa sample.bam \
      > /dev/null ; } 2>&1 | grep real
done
```

## What this means for the roadmap

The PLAN §2 "4–10× via native multithreading" promise cannot be
delivered by parallel compute alone; it needs parallel I/O. The
cleanest path forward is Phase C = option 2 (single-reader producer
+ parallel walk), which avoids the hts_sys races and unlocks real
parallelism on the work that dominates.

Alternative: accept that our parallel speedup is compute-only and
document that `freebayes-gxy` matches upstream in accuracy (F1 = 1.0)
and numerical likelihoods (GL byte-identical vs `--legacy-gls`) but
matches `freebayes-parallel`'s shell-level sharding in performance
— not the PLAN's in-process vision.

## Workspace state

- 129 tests pass.
- cargo build / test / clippy / fmt all green.
- Parity: small fixture byte-identical; large fixture F1 = 1.0;
  100 kb fixture byte-identical t=1 vs t=N.

## Honest assessment

Phase B ships a correct, stable, documented parallel compute path
with a shared reference cache and a real 100 kb benchmark fixture.
It does **not** deliver wall-clock speedup. Calling that out
explicitly in this report rather than burying it — the next
contributor needs to know exactly where to look.
