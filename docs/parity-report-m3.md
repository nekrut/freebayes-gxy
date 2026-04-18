# Parity report — M3 (Phase C-2 HEAD)

**Branch:** `claude/rewrite-freebayes-oL5RQ` @ `a6fc480`
**Upstream:** `freebayes v1.3.10` built from `git clone --branch v1.3.10` at `/tmp/freebayes-upstream/build/freebayes`
**Date:** 2026-04-18

## Summary

**12/12 byte-identical parity** on the M3 synthetic fixture across
`(chrom, pos, ref, alt, gt_tag)`. No site misses, no extra calls, no
genotype flips.

This is the first cross-tool end-to-end sanity check for the Bayesian
call pipeline (`--call` → pileup → posterior → MAP). It covers the
diploid-SNP happy path at 20× mean coverage with Q40 bases and MQ60
reads. It does **not** yet cover indels, haplotype-windowed complex
events, low-coverage / low-BQ edge regimes, multi-sample joint calls,
or the full VCF schema.

## Fixture

- Single 1000-bp contig `chrS`, pseudo-random bases seeded with
  `random.seed(42)`.
- **12 truth variants**: 9 heterozygous SNPs at 1-based positions
  `101, 201, …, 901` (A→G) + 3 homozygous SNPs at `51` (A→T), `551`
  (C→G), `951` (G→A).
- 181 single-end reads, 100 bp, staggered every 5 bp, ~20× coverage
  per site (30× at hom sites with smaller window), `SM:SAMPLE`,
  `RG:synth`, `MAPQ=60`, all bases Q40.
- Coordinate-sorted BAM + `.bai` index; no duplicates, no secondary,
  no supplementary.

Fixture generator: `tests/parity/fixture/build_fixture.py`.
Harness driver: `tests/parity/harness.sh`.

## Commands run

```bash
# Build upstream (once):
cd /tmp/freebayes-upstream && meson setup build --buildtype release && \
    meson compile -C build

# Build fixture:
python3 tests/parity/fixture/build_fixture.py

# Run harness:
tests/parity/harness.sh \
    /tmp/freebayes-upstream/build/freebayes \
    ./target/debug/freebayes-gxy \
    tests/parity/fixture
```

## Result

```
==> running upstream freebayes...
==> running freebayes-gxy --call...

==> upstream calls: 12
==> gxy calls:      12

==> PARITY: byte-identical on (chrom, pos, ref, alt, gt_tag)
```

Full normalised site table (both tools produce it identically):

| chrom | pos | ref | alt | gt_tag  |
|-------|-----|-----|-----|---------|
| chrS  | 51  | A   | T   | SNP     |
| chrS  | 101 | A   | G   | REF/SNP |
| chrS  | 201 | A   | G   | REF/SNP |
| chrS  | 301 | A   | G   | REF/SNP |
| chrS  | 401 | A   | G   | REF/SNP |
| chrS  | 501 | A   | G   | REF/SNP |
| chrS  | 551 | C   | G   | SNP     |
| chrS  | 601 | A   | G   | REF/SNP |
| chrS  | 701 | A   | G   | REF/SNP |
| chrS  | 801 | A   | G   | REF/SNP |
| chrS  | 901 | A   | G   | REF/SNP |
| chrS  | 951 | G   | A   | SNP     |

## What parity does NOT yet assert

The harness collapses both tools into a shared `(chrom, pos, ref, alt,
gt_tag)` schema, so the following do **not** contribute to the
byte-identical verdict and remain open questions:

- **Posterior magnitudes.** Upstream QUAL vs. freebayes-gxy GQ; the
  scales differ (upstream uses a cohort-style variant QUAL, our port
  currently emits per-sample GQ). Numerical agreement within the 1e-9
  tolerance of PLAN §5 M3 has not been measured.
- **Indel calling.** The fixture has none. Phase C-3 should add ≥ 1
  insertion, 1 deletion, and 1 complex event to the fixture.
- **Multi-sample joint calls.** Only one sample in this fixture.
  `GenotypeCombo` posteriors are unported.
- **Haplotype-windowed alleles.** Upstream's `buildHaplotypeAlleles`
  can emit multi-base composite alleles at a site; our per-position
  caller does not. Tested only where upstream and port happen to
  converge (single-SNP sites).
- **Low-coverage / low-BQ regimes.** All bases here are Q40 / MQ60;
  the drift likely to appear at Q10–Q20 (where the base-quality
  `max(ln_bq, ln_mq)` gate matters) is untested.
- **VCF schema.** Our CLI emits TSV; byte-level VCF parity lands with
  M4.

## Parity-harness exit criteria (vs. PLAN §5 M3)

- **"per-site likelihood values match upstream within 1e-9"** — not
  measured; harness only diffs the `(pos, ref, alt, gt_tag)` tuples.
  Blocked on either emitting likelihoods from both tools in a common
  form, or pinning upstream's internal GL values from its `GL` FORMAT
  field and the port's from the `log_posteriors` vector.
- **"genotype calls byte-identical"** — ✅ met on this fixture for
  SNP-only diploid data.

## Next steps

1. **Phase C-3 fixture expansion**: add indels (INS, DEL, COMPLEX),
   a mixed-coverage site, and a low-BQ cluster. Expect drift to
   surface; catalogue and triage.
2. **Likelihood-level parity**: extract `GL` from upstream VCF and
   `log_posteriors` from gxy (needs a `--emit-likelihoods` flag), diff
   per genotype.
3. **M4 (VCF emission)**: promote the TSV to proper VCF with
   deterministic ALT ordering, 0/1 GT syntax, anchor-base synthesis
   for indels, and the full INFO field set. Then re-run parity with
   `bcftools isec` / `hap.py` for the standard F1/precision/recall
   numbers.
4. **Larger fixtures**: a 100 kb and 1 Mb slice of a real reference
   + simulated reads (e.g. `wgsim`) to stress per-position calling
   at scale, then eventually GIAB HG002 chr20.
