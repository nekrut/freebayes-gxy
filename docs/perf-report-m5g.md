# Performance report — M5 Phase G (htslib BGZF-decoder threading)

**Branch:** `claude/rewrite-freebayes-oL5RQ` @ pending commit
**Fixture:** `tests/parity/fixture-chrscale/` — 1 MB / 11 663 variants /
199 971 reads (same as M5-F).
**Change:** new `--bam-threads N` CLI flag that calls
`rust_htslib::bam::Reader::set_threads(N)` after opening the reader.
Default `0` → no change (matches M5-F behaviour).

## Summary

**Net result: small but real, mostly hidden by Amdahl.**
`set_threads(4)` gives a **2.4× speedup on the BAM-read stage**
(193 ms → 81 ms at t=1) but wall-clock barely moves
(1.526 s → 1.520 s at t=1, ≈ 0.4 % — within noise) because read is
only ~13 % of total wall after M5-D's pileup-insert fix.

Keeping the flag in, off by default. Real-world BAMs that are more
compressed (GATK recalibrated, deep WGS) may see a larger win. The
knob costs nothing when unused.

## Scoreboard (3-run average wall, chrscale 1 Mb)

| `--bam-threads` | `--threads 1` | `--threads 4` |
|---:|---:|---:|
| 0 (default)     | 1.526 s | 1.399 s |
| 1               | 1.521 s | 1.345 s |
| 2               | 1.519 s | 1.366 s |
| 4               | 1.528 s | 1.347 s |
| 8               | 1.520 s | 1.368 s |

Variance across bam-threads at either call-thread setting is < 4 %
and non-monotonic — indistinguishable from noise.

## Pileup-stage breakdown (t=1, `FBGXY_PROFILE=1`)

| Stage           | bam-threads=0 | bam-threads=4 | delta |
|---              |-------------:|-------------:|------:|
| BAM read        |   193 ms     |    81 ms     | **2.4× faster** |
| CIGAR walk      |   313 ms     |   319 ms     | ~same |
| M2 clump        |   197 ms     |   202 ms     | ~same |
| Pileup insert   |   242 ms     |   245 ms     | ~same |
| **total pileup** | **945 ms**  | **847 ms**   | 1.12× |

The read-stage speedup is real — htslib's BGZF decoder thread pool
is doing useful work. But it saves ~112 ms out of a 1.5 s wall, so
the visible wall-clock benefit is drowned by jitter in the other
stages.

## Why the bounded payoff

After M5-D's lazy REF-run materialisation, the profile changed shape:

| stage            | share of t=1 wall |
|---               |------------------:|
| BAM read (hts_read_record + BGZF decode) | 13 % |
| CIGAR walk       | 21 % |
| M2 clump         | 13 % |
| Pileup insert    | 16 % |
| Call + emit VCF  | ~35 % |
| Misc             | ~2 % |

`set_threads(n)` only parallelises the first line. Even a perfect
decode-to-zero would only save 13 %. The CIGAR walk, clump, insert,
and call stages are all main-thread CPU work that needs a different
lever — streaming pileup, finer chunking, or moving walk into a
producer thread.

## Correctness

- `cargo test --workspace` — 129/129 pass.
- `cargo clippy --workspace --all-targets -- -D warnings` — clean.
- `cargo fmt --all --check` — clean.
- Output byte-identical across `--bam-threads ∈ {0, 2, 4}` at t=1 on
  the chrscale fixture (`diff -q` returned no differences).
- No change to records emitted (11663, matches truth).

## When the flag will matter more

The tiny wall-clock improvement here reflects our synthetic fixture,
not the flag's real ceiling. Cases where `--bam-threads N` should
deliver more:

- **Network-mounted BAMs (NFS, S3/cloud FS)** — read becomes I/O-bound
  rather than decode-bound; extra threads hide latency.
- **GATK-recalibrated BAMs with BQSR tags** — denser blocks, slower
  decode per byte.
- **High-coverage (100–300×) WGS** — read share of total wall grows
  because per-variant call cost doesn't scale with coverage beyond
  the ~30× already sampled.
- **When `--threads` is high and the call stage is already
  saturated** — BAM read becomes the remaining sequential bottleneck.

## Implementation

Single-point change in `crates/fb-cli/src/main.rs`:

```rust
#[arg(long = "bam-threads", default_value_t = 0, value_name = "N")]
bam_threads: usize,
```

Plus two call sites (one in `run_call`, one in the pileup pass of
`run_call_parallel`):

```rust
let mut bam_reader = bam::Reader::from_path(&cli.bam)?;
if cli.bam_threads > 0 {
    bam_reader
        .set_threads(cli.bam_threads)
        .with_context(|| format!("bam::Reader::set_threads({})", cli.bam_threads))?;
}
```

Note: `set_threads` is called on the single `bam::Reader` inside one
worker. It adds *background* decode threads inside htslib's own
pool and is safe within one reader — we're not reopening from
multiple threads (which M5-A / M5-C established is race-prone).

## Next step

`set_threads` has delivered what it's going to deliver on a dense
synthetic fixture. The remaining single-threaded wall (≈1.3 s at 1 Mb,
extrapolated to ~50 s at 50 Mb chr-scale if linear) is in walk +
clump + call, and the dominant memory cost is still the non-streaming
pileup. Streaming pileup (process in coordinate order, flush finished
sites, free behind) is the next real win.
