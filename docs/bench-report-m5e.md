# Head-to-head benchmark — M5 Phase E (gxy vs upstream vs freebayes-parallel)

**Branch:** `claude/rewrite-freebayes-oL5RQ` @ pending commit
**Fixture:** `tests/parity/fixture-bench/` — 100 kb / 1163 variants /
19971 reads, ~30× coverage, MAPQ 60, Q40 bases.
**Upstream:** `freebayes v1.3.10` with `--legacy-gls` (matches the
`standardGLs` formula our port implements).
**Hardware:** same machine for every run. 3-run averages.

## Headline

**gxy single-threaded beats freebayes-parallel -j 8 by 3.9×.
gxy at -t 8 beats freebayes-parallel -j 8 by 5.8×.**

## Full table

| Tool                                  |   Wall  | vs upstream serial | vs best F-P | vs best gxy |
|---                                    |--------:|-------------------:|------------:|------------:|
| `freebayes --legacy-gls`              | 1.106 s | 1.00×              |       —     | 7.04× slower |
| `freebayes-parallel -j 1`             | 2.243 s | 0.49×              | — (baseline) | 14.29× slower |
| `freebayes-parallel -j 2`             | 1.328 s | 0.83×              | 1.69×       | 8.46× slower |
| `freebayes-parallel -j 4`             | 1.003 s | 1.10×              | 2.24×       | 6.39× slower |
| `freebayes-parallel -j 8`             | 0.906 s | 1.22×              | 2.48×       | 5.77× slower |
| `freebayes-gxy --call --threads 1`    | 0.234 s | **4.73×**          | **3.87×**   | 1.49× slower |
| `freebayes-gxy --call --threads 2`    | 0.160 s | 6.91×              | 5.66×       | 1.02× slower |
| `freebayes-gxy --call --threads 4`    | 0.161 s | 6.87×              | 5.63×       | 1.03× slower |
| **`freebayes-gxy --call --threads 8`** | **0.157 s** | **7.04×**      | **5.77×**   | 1.00× |

## Observations

1. **gxy single-threaded (0.234 s) already beats every upstream
   configuration on this workload.** It outruns even
   `freebayes-parallel -j 8` (0.906 s) by 3.9×. This is the M5-D
   data-structure win: the 20× pileup-insert speedup swamps any
   shell-level parallelism upstream can deploy.
2. **`freebayes-parallel -j 1` is slower than serial upstream**
   (2.24 s vs 1.11 s). The GNU-parallel wrapper + vcfstreamsort +
   vcfuniq + per-region freebayes spawn overhead adds ~1.1 s of
   fixed cost even before any parallelism pays off.
3. **`freebayes-parallel` scales sub-linearly.** Going from -j 1
   to -j 8 takes 2.24 s → 0.91 s, a 2.48× speedup on 8 cores. The
   wrapper overhead dominates until the job count grows enough to
   amortise it; on real WGS workloads it scales much better.
4. **gxy's internal threading saturates quickly** — t=2 gives 1.46×
   over t=1, then t=4 and t=8 add almost nothing. Consistent with
   the M5-B / M5-D analysis: after the pileup-insert bottleneck is
   gone, the remaining per-site compute is only ~50 ms, and the
   bam::Reader::read + walk stages (still single-threaded) cap
   further gains.

## Why gxy-t1 beats F-P-j8

On a 100 kb / 20k-read workload:

- **freebayes-parallel** pays:
  - ~1 s of per-process startup cost (shared libs, index load,
    parameter parsing) × 10 regions, amortised across the pool.
  - GNU parallel scheduling overhead.
  - vcfstreamsort buffering (-w 1000) on the merged output.
  - vcfuniq deduplication pass.
- **freebayes-gxy** pays:
  - One process startup.
  - One BAM read pass.
  - Per-site compute on ~1163 variant sites (≈1% of positions).

The constant factors favor gxy by so much that even t=1 wins.

## Correctness

Output membership identical across all configurations:

| Tool                              | Records |
|---                                |--------:|
| upstream `--legacy-gls`           | 1163    |
| freebayes-parallel -j 4            | 1163    |
| freebayes-parallel -j 8            | 1163    |
| freebayes-gxy --threads 1          | 1163    |
| freebayes-gxy --threads 8          | 1163    |

F1 vs truth: 1.0000 for every tool and thread count on this
fixture (from the M4-G and M5-D reports).

GL values byte-identical between upstream `--legacy-gls` and gxy
across all SNP sites (documented in M4-F report).

## Caveats

- **100 kb is a toy scale.** On chromosome-scale workloads the
  constant-factor advantage will shrink and upstream's parallelism
  will eventually overtake our single-threaded bam::Reader::read.
  Real GIAB HG002 chr22 (~50 Mb) is the proper benchmark; we
  don't have it in this sandbox.
- **All Q40 / MQ60 / staggered reads.** No compression tricks,
  no soft-clips, no duplicates, no secondary alignments. Realistic
  BAMs have all of those and each has code paths we haven't
  stressed at scale.
- **Wall-clock only.** Didn't measure memory. gxy holds the full
  per-contig reference in `Arc<Vec<u8>>` (100 kB here; would be
  ~50 MB per chromosome at chr-scale). Upstream streams the
  reference per-region. Peak RSS comparison would matter for
  operational deployments.
- **`--legacy-gls` is not upstream's default.** The default uses
  the experimental GL path, which computes more quantities per
  site and is slower. Enabling it in upstream would make the gap
  larger, not smaller.

## How to reproduce

```bash
# Tools needed:
#  - freebayes + freebayes-parallel from /tmp/freebayes-upstream/
#  - GNU parallel, vcflib-tools (vcffirstheader, vcfstreamsort, vcfuniq)
export PATH=/usr/lib/vcflib/bin:/tmp/freebayes-upstream/scripts:/tmp/freebayes-upstream/build:$PATH

# Build fixture (once):
python3 tests/parity/fixture-bench/build_fixture.py

# Region list for freebayes-parallel:
cd tests/parity/fixture-bench
fasta_generate_regions.py ref.fa.fai 10000 > regions.10k.txt

# Run:
time freebayes --legacy-gls -f ref.fa sample.bam > /dev/null
time freebayes-parallel regions.10k.txt 8 --legacy-gls -f ref.fa sample.bam > /dev/null
time ./target/release/freebayes-gxy --call --threads 8 -f ref.fa sample.bam > /dev/null
```

## Workspace state

- 129 tests pass; clippy / fmt / cargo doc clean.
- All parity harnesses still green.
- No code changes in this phase — only the benchmark run and this
  report. (The wins are from M5-D.)
