# freebayes-gxy

[![CI](https://github.com/nekrut/freebayes-gxy/actions/workflows/ci.yml/badge.svg)](https://github.com/nekrut/freebayes-gxy/actions/workflows/ci.yml)

A from-scratch Rust rewrite of [freebayes](https://github.com/freebayes/freebayes),
targeting byte-level VCF parity with upstream **v1.3.10** while delivering
native multithreading, modern I/O, and a clean embeddable library API for
Galaxy and nf-core pipelines.

See [`PLAN.md`](PLAN.md) for the full design, milestone breakdown, and parity
targets. **Status: M0 scaffolding — no variant calling yet.**

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

# Emit a placeholder VCF header for an indexed BAM
cargo run -- -f reference.fa input.bam
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
