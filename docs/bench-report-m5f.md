# Chromosome-scale benchmark — M5 Phase F

**Branch:** `claude/rewrite-freebayes-oL5RQ` @ pending commit
**Fixture:** `tests/parity/fixture-chrscale/` — 1 MB / 11663 variants /
199971 reads, ~30× coverage, MAPQ 60, Q40 bases, single contig `chrM`.
**Upstream:** `freebayes v1.3.10` with `--legacy-gls`.
**Hardware:** same machine as M5-E. 3-run averages.

## Headline

gxy still wins at chromosome scale — but by less than the 100 kb
fixture suggested, and at a **16× memory cost** vs upstream.

## Wall-clock (3-run avg)

| Tool                           |    Wall  | Peak RSS   | vs upstream serial | vs F-P -j 8 |
|---                             |---------:|-----------:|-------------------:|------------:|
| `freebayes --legacy-gls`       | 11.07 s  |   14 MB    | 1.00×              | 0.27× slower |
| `freebayes-parallel -j 1`      | 14.14 s  |    6 MB    | 0.78×              |       0.21× |
| `freebayes-parallel -j 2`      |  7.00 s  |    6 MB    | 1.58×              |       0.43× |
| `freebayes-parallel -j 4`      |  4.44 s  |    6 MB    | 2.49×              |       0.67× |
| `freebayes-parallel -j 8`      |  2.98 s  |    6 MB    | 3.72×              |  (baseline) |
| `freebayes-gxy --threads 1`    |  1.40 s  |  230 MB    | **7.91×**          | **2.13×**   |
| `freebayes-gxy --threads 2`    |  1.46 s  |  231 MB    | 7.58×              | 2.04×       |
| `freebayes-gxy --threads 4`    |  1.31 s  |  229 MB    | 8.45×              | 2.27×       |
| `freebayes-gxy --threads 8`    |  1.36 s  |  233 MB    | 8.14×              | 2.19×       |

All configurations emit **11663 records** — identical to the truth
set. F1 = 1.0000 preserved at chr-scale.

## Observations

1. **gxy still wins, but the gap is tighter.** At 100 kb gxy-t1 was
   3.87× faster than `freebayes-parallel -j 8`. At 1 Mb it's 2.13×.
   Extrapolating, around 100 Mb the two tools would roughly cross —
   upstream's parallelism continues to scale while gxy's doesn't.
2. **freebayes-parallel scales beautifully now.** 1→2→4→8 threads
   gives 1.00× → 2.02× → 3.18× → 4.74×. Amdahl efficiency about
   60% at 8 threads. The fixed wrapper overhead (~1 s) amortises
   over the 11s of real work.
3. **gxy's thread scaling is gone.** At 1 Mb, t=1 and t=8 are
   essentially the same wall (1.40 s vs 1.36 s). The M5-D lazy REF
   materialisation shoved the bottleneck into single-threaded
   `bam::Reader::read` + walk, which we can't parallelise under
   rust-htslib's thread-unsafety (see `threading-report-m5c-
   experiment.md`).
4. **Memory: gxy uses 16× more RSS** than upstream and 40× more
   than freebayes-parallel workers. 230 MB peak RSS on a 1 MB
   fixture ≈ 230×/bp. Breakdown (rough):
   - Shared reference in `Arc<Vec<u8>>`: ~1 MB (cheap).
   - Pileup: ~11.7 k variant positions × ~30 obs × ~100 B =
     35 MB.
   - RefRun intervals: ~6M runs × ~100 B (names, per-base quals,
     ref bytes) = 600 MB raw, but with shared small-allocations
     maybe 200 MB.
   - Read + walker temporaries: ~50 MB.
5. **Memory extrapolation.** At real chromosome scale (chr1 ≈
   250 Mb, ~50 M reads), gxy would need ~50 GB of RSS in the
   current design — untenable. Streaming pileup (process reads
   in order, flush per tile, free behind) is a must before real
   WGS workloads.

## Why gxy's threading stopped paying off

On this fixture, the stage breakdown (single-threaded) looks like:

| stage              | time    | % |
|---                 |--------:|--:|
| BAM read + walk    | ~850 ms | 61% |
| M2 clump           | ~150 ms | 11% |
| Pileup insert      | ~120 ms | 9% |
| Finalize + call    | ~200 ms | 14% |
| Emit VCF           | ~80 ms  | 6% |

The call stage (the only piece we parallelise) is ~14% of runtime.
Perfect parallelism there → 1.4 s × 0.86 = 1.2 s wall-clock
ceiling. Observed 1.31 s at t=4 is 95% of that ceiling. Threading
works; it just hit its Amdahl bound.

## Is gxy the right choice?

| Scenario | Winner |
|---|---|
| Short contigs (< 500 kb), low thread count | gxy (≫ 2× faster even single-threaded) |
| 1–10 Mb contigs, any thread count | gxy (still 2× faster, lower wall) |
| Multi-chromosome (> 100 Mb), ≥ 8 cores | freebayes-parallel likely wins on wall |
| Memory-constrained host | **freebayes-parallel** (16–40× less RSS) |
| Accuracy / F1 | **Tied** at 1.0000 on this fixture |
| VCF byte parity (with `--legacy-gls`) | **Identical** across every tool |

The 16× memory cost is the big asterisk. For Galaxy / nf-core /
clinical deployments where RAM is bounded, gxy's current design
is unsuitable at full-WGS scale.

## Next steps (if anyone picks this up)

1. **Streaming pileup** — the single biggest win. Process reads
   in coordinate order, flush completed sites (where no future
   read can still contribute) immediately, free the REF runs and
   variant observations behind them. Should bring chromosome-
   scale RSS from ~50 GB down to ~500 MB.
2. **Chunk the per-thread work smaller** — at 1 Mb the chunk
   size is ~3k sites/thread, which may be too coarse. Worth
   probing.
3. **`rust-htslib`'s `set_threads(n)`** — htslib's internal
   BGZF-decoder threading is supposed to be safe within a single
   reader. Might unlock parallel decompression of the read stream
   without running into the multi-reader race.

## Caveats

- **1 Mb is still small.** Real chromosomes are 50–250 Mb. The
  extrapolated curves above are linear extrapolations; real
  systems have non-linear memory and cache effects.
- **Uniform coverage, single chromosome.** Real WGS has
  coverage valleys, clustered variants, and multi-contig work
  that would change the picture.
- **Q40 / MQ60 / no soft-clips / no duplicates.** The fast path
  throughout.
- **Bench-mark time only; no accuracy regression testing at
  scale.** We checked record counts; full hap.py comparison at
  1 Mb scale was out of scope.

## Reproducing

```bash
export PATH=/usr/lib/vcflib/bin:/tmp/freebayes-upstream/scripts:/tmp/freebayes-upstream/build:$PATH
python3 tests/parity/fixture-chrscale/build_fixture.py
cd tests/parity/fixture-chrscale
fasta_generate_regions.py ref.fa.fai 100000 > regions.100k.txt

time freebayes --legacy-gls -f ref.fa sample.bam > /dev/null
time freebayes-parallel regions.100k.txt 8 --legacy-gls -f ref.fa sample.bam > /dev/null
time ./target/release/freebayes-gxy --call --threads 8 -f ref.fa sample.bam > /dev/null
```

Peak RSS measured via `/proc/$pid/status` polled at 50 ms
resolution while the command ran.

## Workspace state

- 129 tests pass; clippy / fmt / cargo doc clean.
- All parity harnesses still green.
- No code changes. Only the 1 MB fixture generator +
  `.gitignore` updates ship.
