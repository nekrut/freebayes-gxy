# Parity harness

This directory holds the scaffolding for the parity sweep against upstream
freebayes **v1.3.10** (see `PLAN.md` §3). M0 ships only the diff driver; the
baseline pipeline is filled in during M1–M4 as the Rust caller starts
emitting real records.

## Contents

- `compare.sh` — takes two VCF paths, runs `bcftools isec`, prints a
  summary of shared vs unique records. Fails loudly if `bcftools` is not on
  `PATH`.

## Wiring in a real v1.3.10 baseline

The intended workflow (to be automated in CI once M4 lands) is:

1. **Fetch inputs.** Pin a small GIAB HG002 slice — the plan uses a 10 Mb
   window on chr20. Store the BAM, BAI, and reference FASTA (with `.fai`)
   under a cache dir that CI can restore via `actions/cache`.

2. **Run upstream freebayes v1.3.10.** Easiest via the official Docker
   image:

   ```bash
   docker run --rm -v "$PWD:/data" \
       quay.io/biocontainers/freebayes:1.3.10--py310h077b44d_0 \
       freebayes -f /data/ref.fa /data/HG002.chr20.bam \
       > baseline.vcf
   bgzip -f baseline.vcf
   tabix -p vcf baseline.vcf.gz
   ```

3. **Run freebayes-gxy.** Same region, same reference:

   ```bash
   cargo run --release -p fb-cli -- \
       -f ref.fa -v candidate.vcf HG002.chr20.bam
   bgzip -f candidate.vcf
   tabix -p vcf candidate.vcf.gz
   ```

4. **Diff.**

   ```bash
   ./tests/parity/compare.sh baseline.vcf.gz candidate.vcf.gz parity-out
   ```

   The summary should report 100% shared records by M4 (modulo ULP-level
   QUAL drift — see `PLAN.md` §3).

## Exit criteria by milestone

| Milestone | Parity expectation                                 |
|-----------|----------------------------------------------------|
| M0        | `compare.sh` runs, fails cleanly without bcftools  |
| M1        | allele counts match on 10 Mb window                |
| M4        | ≥99.9% shared records on `bcftools isec`           |
| M7        | hap.py F1 drift ≤ 0.001 on full WGS                |
