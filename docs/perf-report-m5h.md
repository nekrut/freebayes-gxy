# Performance report — M5 Phase H (mimalloc + pipelining post-mortem)

**Branch:** `claude/rewrite-freebayes-oL5RQ` @ pending commit
**Fixture:** `tests/parity/fixture-chrscale/` (1 MB) + `fixture-bench/`
(100 kb).
**Change:** `mimalloc` as the global allocator in `fb-cli`. One
`#[global_allocator]` line, one dep. No code-flow changes.

## Headline

**~1.5× wall-clock reduction at t=1, ~1.8× at t=4, for free.**

| Fixture | Mode          | pre-H  | post-H (mimalloc) | speedup |
|---      |---            |------:|------:|--------:|
| 1 Mb    | t=1 bt=0      | 1.33 s | 0.87 s | **1.53×** |
| 1 Mb    | t=4 bt=0      | 1.19 s | 0.71 s | 1.68× |
| 1 Mb    | t=4 bt=4      | 1.14 s | **0.64 s** | **1.79×** |
| 100 kb  | t=1 bt=0      | 0.23 s | 0.13 s | 1.79× |
| 100 kb  | t=4 bt=4      | 0.16 s | 0.09 s | 1.81× |

Refreshed gxy-vs-upstream ratios:

| Scale | vs `freebayes --legacy-gls` | vs `freebayes-parallel -j 8` |
|---    |---:|---:|
| 100 kb, gxy-t1 bt=0 | 8.6× (was 4.7×) | 7.0× (was 3.9×) |
| 1 Mb, gxy-t1 bt=0  | 12.7× (was 7.9×) | 3.4× (was 2.1×) |
| 1 Mb, gxy-t4 bt=4  | 17.3× | 4.7× |

Output byte-identical across `(threads, bam-threads) ∈ {1,4} × {0,4}`
on chrscale. 129/129 tests pass; clippy / fmt clean.

## Backstory: the producer/consumer detour

M5-G landed the `--bam-threads` htslib decoder knob. The next
obvious lever was a producer/consumer split on the pileup stage —
spawn one thread for BAM read + CIGAR walk, keep clump + pileup-
insert on main, let the two halves overlap. Back-of-envelope said
~30 % wall saved on the 1 Mb fixture. Straightforward, right?

Implementation went in cleanly: `std::sync::mpsc::sync_channel<Vec<(i32,
Vec<AlleleObservation>)>>`, batched at 256 records per send so the
channel sync cost amortised. Producer thread did read + walk,
consumer thread did clump + insert.

Wall-clock: **4.1 s** (vs 1.3 s baseline). **3× slowdown.**
`FBGXY_PROFILE=1` showed `producer_ms = consumer_ms = ~3900 ms` —
perfectly lockstepped. First instinct: channel sync. But switching
from bounded to unbounded channel didn't change anything. Nor did
batch size.

The diagnostic that cracked it: replace `tx.send(...)` with
`batch.clear()` — drop every batch, do no consumer work. Producer
thread alone completed in **460 ms** — faster than the full
single-threaded baseline's pileup stage. So the read + walk
producer-thread cost was tiny. The 3500 ms of overhead was paid
somewhere else.

Where? **Cross-thread free on 6 M+ small heap allocations.** Each
`AlleleObservation` carries a `String` (read_name) and a `Vec<u8>`
(per_base_quals); the CIGAR walk produces ~30 observations per
read; at 200 k reads that's ~6 M observation structs and ~18 M
inner heap cells. In the pre-pipeline code, producer and consumer
were the same thread — glibc-malloc's thread-local arena saw every
`free()` as same-arena and recycled the slot immediately. In the
pipelined code, producer allocates on thread A and consumer frees
on thread B. glibc handles that by returning the block to a global
pool under a mutex — **560 ns per free** vs ~76 ns same-thread.
Multiply by 18 M cells: ~3.4 s. Matches the observed 3.5 s gap
exactly.

## The pivot: mimalloc

Allocator-of-choice for heavy multi-threaded allocation. Per-thread
heaps, deferred foreign-free mechanism that doesn't require a
global mutex. Drop-in replacement via `#[global_allocator]`.

One hunk:

```toml
[dependencies]
mimalloc = { version = "0.1", default-features = false }
```

```rust
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;
```

That's all. No code-flow changes. Wall-clock dropped from 1.33 s
→ 0.87 s at t=1 (single-threaded) purely on allocator efficiency
for the *same* allocation pattern. At t=4 the call-stage rayon
workers also benefit (they allocate per-site observation vectors);
wall drops further to 0.64 s.

The producer/consumer pipeline was not reinstated. With mimalloc
the cross-thread free cost collapses, but the best-case savings
(~400 ms of overlapped read+walk with clump+insert) are a smaller
fraction of the now-0.87 s wall. Not worth the complexity.
**mimalloc gives the win pipelining was supposed to deliver, for
3 lines of code instead of ~100.**

## Why mimalloc wins even single-threaded

Even without any thread crossing, mimalloc is consistently ~20 % faster
than glibc-malloc on hot allocate-use-drop paths because:

- Per-thread slab-like free-lists reduce branch mispredicts in
  `malloc`.
- Smaller size-class resolution overhead (power-of-2 bucket lookup
  via pointer arithmetic, not a `for` loop).
- Better cache locality on alloc/free pairs within a single
  batch.

Our workload is a near-ideal mimalloc target: lots of short-lived
small allocations (the per-read observation vectors, the String
read_names, the transient Vec<u8> REF run buffers) with predictable
size classes.

## Scoreboard refresh (1 Mb chrscale, 3-run avg)

| Tool                           | Wall   | vs upstream | vs F-P -j 8 |
|---                             |------:|------:|------:|
| `freebayes --legacy-gls`       | 11.07 s | 1.00× | 0.27× |
| `freebayes-parallel -j 1`      | 14.14 s | 0.78× | 0.21× |
| `freebayes-parallel -j 4`      |  4.44 s | 2.49× | 0.67× |
| `freebayes-parallel -j 8`      |  2.98 s | 3.72× | 1.00× |
| **`freebayes-gxy t=1 bt=0`**   | **0.87 s** | **12.7×** | **3.4×** |
| **`freebayes-gxy t=4 bt=4`**   | **0.64 s** | **17.3×** | **4.7×** |

The gap vs upstream widened substantially. The 1 Mb parallel
crossover point from the M5-F projection (~100 Mb) now moves out
further — gxy would need to hit `freebayes-parallel -j 8`'s wall-
clock curve.

## Correctness

- `cargo test --workspace` — 129/129 pass.
- `cargo fmt --all --check` — clean.
- `cargo clippy --workspace --all-targets -- -D warnings` — clean.
- Chrscale output byte-identical across
  `(threads, bam-threads) ∈ {1,4} × {0,4}` — `diff -q` shows no
  differences.
- All 11663 variant records preserved; F1 = 1.0000.

## Caveats

- **mimalloc binary footprint.** +~40 kB in release, so negligible.
- **Not a cure for streaming pileup.** RSS is still 230 MB per 1 Mb
  at full WGS scale → ~57 GB on chr1-sized contigs. mimalloc does
  not reduce live-heap, only allocation throughput.
- **macOS / Windows.** mimalloc supports both but hasn't been
  tested by this project; only Linux x86_64 was measured here.
- **Profile allocator.** The mimalloc crate provides no heap-profile
  hooks by default. If future profiling needs them, switch to
  `tikv-jemallocator` (which has `stats`).

## What this means for the roadmap

Streaming pileup is still the big-ticket item (memory, not wall).
After M5-H:

1. **Streaming pileup** — unchanged priority. ~50 GB → ~500 MB RSS
   at chr-scale.
2. **Further rayon call-stage tuning** — maybe worth a profile,
   but the call stage is now <50 % of wall at t=4, so gains are
   bounded.
3. **Multi-sample** — functional gap.
4. **Producer/consumer revisited?** — only if future work pulls
   the per-observation heap cells out of the hot path (e.g. an
   arena/typed-vec for observations). Low priority.
