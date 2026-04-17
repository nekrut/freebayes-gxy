# freebayes-gxy — Rust Rewrite Plan

A from-scratch Rust reimplementation of [freebayes](https://github.com/freebayes/freebayes), targeting byte-level VCF parity with upstream **v1.3.10** while delivering native multithreading, modern I/O, and a clean embeddable library API for Galaxy and nf-core pipelines.

## 1. Goals

- **Parity** with freebayes v1.3.10 on the GIAB HG002 high-confidence intervals: hap.py F1 drift ≤ 0.001 for SNPs and indels
- **5–15× single-machine speedup** on 16+ core hosts vs. the `freebayes-parallel` wrapper
- **2–4× peak RSS reduction** on high-coverage (>100×) inputs
- **MIT-licensed**, no GPL dependencies
- **Library-first**: expose a Rust crate (`freebayes-core`) plus a CLI (`freebayes-gxy`) with stable semver
- **Drop-in CLI compatibility** for the flags most pipelines depend on (see §3)

## 2. Non-goals

- No algorithmic novelty — this is a port, not a new caller. Algorithmic improvements deferred to post-1.0.
- No GPU support (see [design note](docs/no-gpu.md) — kernel too small, warp-divergent, host/device transfer dominates)
- No somatic or joint-genotyping modes beyond what upstream freebayes already supports
- No VCF schema changes — output must be bit-reproducible where feasible, record-equivalent otherwise

## 3. Parity target (v1.3.10)

### VCF output
Record-level identity on GIAB HG002 chr20, measured by `bcftools isec`:
- ≥ 99.9% identical records (POS + REF + ALT + GT)
- GQ/DP/AO/RO numeric drift within floating-point rounding (ULP ≤ 4)
- QUAL drift within 0.01

### CLI surface (v1.0 must support)

```
-f  --fasta-reference
-r  --region               (-L --targets)
-v  --vcf
-@  --variant-input
-t  --threads              [NEW — replaces freebayes-parallel]
    --min-alternate-count
    --min-alternate-fraction
    --min-coverage
    --ploidy
    --pooled-discrete / --pooled-continuous
    --haplotype-length
    --use-best-n-alleles
    --report-monomorphic
    --standard-filters
```

Deferred to v1.1: `--cnv-map`, `--populations`, `--theta` advanced Bayesian priors, `--genotyping-max-iterations`.

## 4. Architecture

```
                 ┌──────────────────────┐
BAM(s) ─────────▶│  window scheduler    │
                 │  (work-stealing)     │
                 └──────────┬───────────┘
                            ▼
         ┌────────────────────────────────┐
         │  per-thread worker pipeline    │
         │                                │
         │  htslib iter → pileup          │
         │    → allele observations       │
         │    → haplotype construction    │
         │    → Bayesian likelihood       │
         │    → genotype call             │
         └────────────────┬───────────────┘
                          ▼
                 ┌──────────────────────┐
                 │  ordered VCF writer  │
                 │  (reorder buffer)    │
                 └──────────────────────┘
```

### Crate layout

```
freebayes-gxy/
├── Cargo.toml
├── crates/
│   ├── fb-core/          # library: pileup, alleles, haplotypes, likelihood
│   ├── fb-genotype/      # Bayesian model, priors, SIMD log-sum-exp
│   ├── fb-vcf/           # VCF writer (on top of rust-htslib)
│   ├── fb-scheduler/     # work-stealing window scheduler
│   └── fb-cli/           # binary freebayes-gxy
├── tests/
│   ├── parity/           # byte-diff vs upstream v1.3.10
│   └── integration/      # GIAB hap.py runs
└── docs/
```

### Key dependencies
- **rust-htslib** — BAM/CRAM/VCF I/O (proven, wraps htslib C)
- **rayon** — work-stealing parallelism for window scheduler
- **bio** — reference sequence / FASTA index
- **wide** or **std::simd** — portable SIMD for likelihood kernel
- **clap** — CLI parsing
- **tracing** — structured logging
- **anyhow** / **thiserror** — error handling

No NIH of htslib. No custom allocators until profiling demands it.

## 5. Milestones

### M0 — Scaffolding (week 1)
- Cargo workspace, CI (GitHub Actions: fmt/clippy/test/bench), hap.py Docker harness
- `fb-cli` reads BAM + reference, prints a placeholder VCF header
- Parity harness skeleton: run upstream v1.3.10 + freebayes-gxy, diff VCFs
- **Exit criteria:** `cargo test` green; CI runs on push; baseline VCF captured for HG002 chr20 10 Mb window

### M1 — Pileup + allele observations (weeks 2–3)
- Port `AlleleParser.cpp` semantics: read filtering, base quality handling, MNP/complex allele detection
- Produce per-position `AlleleObservations` struct identical to upstream's internal model
- **Exit criteria:** allele counts match upstream within tolerance on 10 Mb window; unit tests cover indel anchoring edge cases

### M2 — Haplotype construction (weeks 4–5)
- Port `Haplotype.cpp` + `Allele.cpp` haplotype windowing; sliding window over complex events
- Match upstream's `--haplotype-length` behavior exactly
- **Exit criteria:** haplotype sets match upstream for ≥99% of positions on chr20 test slice

### M3 — Bayesian genotype model (weeks 6–7)
- Port `Genotype.cpp` + `Multinomial.cpp` + `Sum.cpp`
- Scalar implementation first; full-precision log-sum-exp; multinomial PMF
- **Exit criteria:** per-site likelihood values match upstream within 1e-9; genotype calls byte-identical

### M4 — VCF writer parity (week 8)
- Record emission through rust-htslib; INFO/FORMAT field ordering, flag semantics, FILTER column
- Deterministic sort-stable output independent of thread count
- **Exit criteria:** `bcftools isec` on single-threaded chr20 run shows 100% record identity (bar ULP QUAL drift)

### M5 — Native threading (weeks 9–10)
- Work-stealing scheduler over genomic windows (default 100 kb tiles, overlap 1 kb)
- Reorder buffer for VCF writer — output sorted regardless of completion order
- Tile-boundary deduplication for overlapping haplotype windows
- **Exit criteria:** 16-thread output matches 1-thread output record-for-record on HG002 chr20; wall-clock 8–12× faster

### M6 — SIMD likelihood kernel (week 11)
- AVX2 + NEON log-sum-exp and multinomial PMF
- Runtime ISA dispatch; scalar fallback retained
- **Exit criteria:** likelihood microbenchmark 2–3× faster; full-run speedup ≥1.3× vs. M5 scalar

### M7 — Full WGS validation (weeks 12–13)
- Run GIAB HG002 30× WGS autosomes
- hap.py against v4.2.1 truth set
- Memory & wall-clock scaling curves @ 1/4/8/16/32/64 threads
- **Exit criteria:** F1 drift ≤ 0.001 SNPs and indels; end-to-end ≥5× speedup vs. `freebayes-parallel -j 16`

### M8 — Release hardening (weeks 14–16)
- Fuzz harness on BAM parsing + VCF writer (cargo-fuzz)
- CRAM support verification
- Reference-based reproducibility across macOS/Linux/aarch64
- Documentation: crate docs, CLI man page, migration guide from upstream freebayes
- Galaxy wrapper + Conda recipe + nf-core module
- **Exit criteria:** v1.0.0 tagged, published to crates.io + Bioconda

### Total: **~16 weeks wall-clock** with one operator driving Claude Code daily

## 6. Validation strategy

### Per-commit CI
- `cargo test` (unit tests, all crates)
- Parity sweep on HG002 chr20 10 Mb slice (takes ~2 min)
- Clippy, rustfmt, cargo-audit

### Per-milestone integration
- Full HG002 chr20 parity against stock v1.3.10
- Record-level `bcftools isec` report archived as CI artifact

### Pre-release
- GIAB HG002 full WGS hap.py report
- HG003 + HG004 generalization check (not in training flows for model fitting, though freebayes has no ML)
- Amplicon panel stress test (high coverage >1000×)
- Polyploid test (--ploidy 4, Arabidopsis, potato)

### Benchmarks (automated, checked into `benchmarks/`)
- `hyperfine` wall-clock table vs `freebayes-parallel -j N` for N in {1, 4, 8, 16, 32}
- Peak RSS via `/usr/bin/time -v`
- Profile-guided optimization runs quarterly

## 7. Risks

| Risk | Mitigation |
|---|---|
| Floating-point divergence in log-sum-exp → VCF drift | Match upstream's exact scalar order-of-operations in M3; SIMD path validated against scalar reference |
| vcflib INFO/FORMAT field emission quirks | Bit-diff against upstream VCF early (M4); maintain a "quirks" compatibility module |
| Haplotype window boundary effects when tiles are parallel | Overlap windows by `--haplotype-length`, dedupe records at tile joins |
| rust-htslib coverage gaps (rare BAM/CRAM features) | Contribute upstream; keep htslib C fallback behind feature flag if needed |
| Memory regressions on high-coverage data | Streaming pileup; per-thread bounded allocator pools; benchmark at every milestone |
| Upstream freebayes changes during rewrite | Freeze parity target at v1.3.10 tag; track upstream changes in `docs/upstream-sync.md` |

## 8. Deferred to post-1.0

- GPU support (see `docs/no-gpu.md` — rejected for architectural reasons)
- Joint cohort genotyping (gVCF output)
- Streaming BAM ingestion from S3/GCS without full local copy
- Alternative Bayesian priors (Dirichlet process, admixture-aware)
- Long-read mode (currently degrades gracefully; no native tuning)
- WebAssembly build for browser-side micro-callers

## 9. How we'll drive this with Claude Code

- One 2–4 hour session per workday
- Every session begins with loading: this plan, the current milestone design doc, the upstream reference file being ported
- Every session ends with a committed parity-harness run
- Small PRs per subsystem; no omnibus commits
- ThreadSanitizer + Miri runs weekly
- Milestone retros logged in `docs/retros/`

## 10. Non-code deliverables

- `docs/architecture.md` — design rationale
- `docs/parity-report.md` — per-release parity evidence
- `docs/migration.md` — how to swap freebayes → freebayes-gxy in existing pipelines
- `docs/benchmarks.md` — reproducible benchmark methodology + results
- Conference/paper target: Bioinformatics Advances application note, Q3 2026

---

**Stack summary:** Rust 1.80+ · rust-htslib · rayon · std::simd · clap · tracing · MIT license · GitHub Actions CI · Bioconda + crates.io distribution
