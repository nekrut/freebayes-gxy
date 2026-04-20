# freebayes-gxy

[![CI](https://github.com/nekrut/freebayes-gxy/actions/workflows/ci.yml/badge.svg)](https://github.com/nekrut/freebayes-gxy/actions/workflows/ci.yml)

A from-scratch Rust rewrite of [freebayes](https://github.com/freebayes/freebayes),
targeting byte-level VCF parity with upstream **v1.3.10** (with
`--legacy-gls`) while delivering native multithreading, modern I/O, and a
clean embeddable library API for Galaxy and nf-core pipelines.

**Status (2026-04-19):** single-sample calling is end-to-end functional.
F1 = 1.0000 vs synthetic truth on 2 kb / 10 kb / 100 kb / 1 Mb fixtures;
GL values byte-identical to upstream `--legacy-gls`; **8.6–17.3× faster
than upstream serial** (data: M5-E / M5-F / M5-H). 129/129 tests pass;
`clippy` / `fmt` / `cargo doc` clean.

See [`docs/README.md`](docs/README.md) for the full project dossier
(architecture, performance story, parity story, roadmap) and
[`PLAN.md`](PLAN.md) for the original design.

## Crates

| Crate           | Purpose                                              |
|-----------------|------------------------------------------------------|
| `fb-core`       | Pileup, allele observations, haplotypes, likelihood. |
| `fb-genotype`   | Bayesian genotype model and log-sum-exp kernels.     |
| `fb-vcf`        | VCF writer (wraps `rust-htslib`).                    |
| `fb-scheduler`  | Work-stealing window scheduler (Rayon).              |
| `fb-cli`        | Binary `freebayes-gxy`.                              |

## Quickstart

```bash
# Build the workspace
cargo build --workspace

# CLI help
cargo run -- --help

# Call single-sample variants
./target/release/freebayes-gxy --call --threads 8 \
    -f reference.fa sample.bam > out.vcf
```

System requirements: stable Rust (MSRV **1.80**), a C toolchain, and the
usual `rust-htslib` deps (`clang`/`libclang`, `zlib`, `bzip2`, `xz`, `curl`,
`openssl`). On Debian/Ubuntu:

```bash
sudo apt-get install build-essential cmake pkg-config zlib1g-dev \
    libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev libclang-dev clang
```

## Testing

```bash
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
```

The parity harness skeleton lives in [`tests/parity/`](tests/parity/);
point `compare.sh` at an upstream freebayes v1.3.10 VCF and a
freebayes-gxy VCF to get a `bcftools isec` summary.

## License

MIT — see [`LICENSE`](LICENSE).
