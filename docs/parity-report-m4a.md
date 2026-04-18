# Parity report — M4 Phase A (VCF emission)

**Branch:** `claude/rewrite-freebayes-oL5RQ` @ pending commit
**Upstream:** `freebayes v1.3.10` at `/tmp/freebayes-upstream/build/freebayes`
**Fixture:** `tests/parity/fixture/` (2 kb contig, 12 SNPs + 4 INS + 4 DEL)

## Summary

| Metric                                       | M3 C-2 (TSV) | M4 Phase A (VCF + bcftools norm) |
|---                                           |---           |---                                |
| Sites detected                               | 20 / 20      | **20 / 20**                       |
| SNP (CHROM+POS+REF+ALT+GT identity)          | 12 / 12      | **12 / 12** ✅                    |
| INS+DEL (CHROM+POS+REF+ALT identity)         |  0 / 8       | **8 / 8** ✅                      |
| INS+DEL (full CHROM+POS+REF+ALT+GT identity) |  0 / 8       | 4 / 8                             |
| **Overall parity**                           | 12 / 20 = 60% | **16 / 20 = 80%**                |

**Two of the three drifts from the M3 report are now closed.** Only the
hom→het flip on hom indels (drift 2) remains, and it is a caller-layer
issue (event-span vs. per-position pileup) rather than a VCF-emission
one.

## What M4 Phase A delivers

- **Real VCF output from `--call`.** The per-site TSV from M3 Phase C-2
  is replaced with a freebayes-compatible VCF 4.2 record stream:
  - VCF-style GT syntax (`0/1`, `1/1`) instead of the internal
    `REF/SNP` / `SNP` tag.
  - Anchor-base synthesis for INS (prepend ref base at `pos - 1` to REF
    and ALT) and DEL (prepend to REF).
  - QUAL column = `-10 * log10(P(hom-ref | obs))`, capped at
    `MAX_GQ = 2000`.
  - FILTER column: `PASS`.
  - INFO fields: `NS, DP, AC, AN, AF, RO, AO, TYPE` — the
    `bcftools isec`-visible subset. (Upstream has ~30; the rest are
    diagnostic and deferred to Phase B.)
  - FORMAT fields: `GT:DP:AD:RO:AO:GQ` — matches upstream's minus the
    `QR`/`QA` quality sums and the `GL` vector. (Those land with the
    likelihood-parity gate.)
  - Deterministic ALT ordering via `(kind, position, ref_seq, alt_seq)`
    sort before genotype enumeration.

- **`fb-vcf` crate expanded.** New public API:
  `Record`, `RecordKind`, `synthesize_anchored`, `write_record`. 9
  unit tests cover header-builder, anchor synthesis for SNP / INS /
  DEL, and record serialisation for hom-ref / het / hom-alt / het-del.

- **Parity harness uses `bcftools norm`.** Both VCFs are normalised
  (left-aligned, allele-canonicalised) before the diff. Eliminates
  the drift-1 representation differences that dominated the M3
  report.

## What remains (drift 2 only)

Four hom indels are reported as `0/1` by gxy, as `1/1` by upstream:

```
< chrS  1450  C      CG    1/1
< chrS  1550  CG     C     1/1
< chrS  1651  G      GCC   1/1
< chrS  1750  CGTA   C     1/1
---
> chrS  1450  C      CG    0/1
> chrS  1550  CG     C     0/1
> chrS  1651  G      GCC   0/1
> chrS  1750  CGTA   C     0/1
```

**Root cause (unchanged from M3 report):** our caller operates on a
**position-based pileup**. A read whose alignment starts at the
deletion position and carries `100M` contributes a REF observation at
that position, because its match op covers it. Upstream's
**event-span-based pileup** excludes such reads from the site's
observation set, since they do not span the anchor + event window.
With 1 contaminating REF obs + 19 DEL obs, our posterior splits
`REF/DEL` ≈ 0.95 vs. `DEL/DEL` ≈ 0.05, MAP is het (GQ ≈ 13).

**Fix location:** Phase C-3 pileup refactor (respect event-span
semantics when admitting observations at an indel site) OR the full
`AlleleParser::buildHaplotypeAlleles` port (M4 Phase B or later).

## What is still not asserted

- **QUAL magnitude parity.** We emit a site QUAL but no tolerance
  check exists. Numerical 1e-9 gate from PLAN §5 M3 remains open.
- **Complex alleles.** No complex variants in this fixture.
- **Multi-allelic sites.** No triallelic variants in this fixture.
- **Multi-sample joint calling.** Single-sample only.
- **Low-BQ / low-coverage regimes.** All bases Q40, MQ 60, ~20×.
- **Haplotype-windowed composite alleles.** The current per-position
  caller emits each event separately; upstream's `buildHaplotype­Alleles`
  can bundle multiple events into one multi-base composite.

## Commands run

```bash
python3 tests/parity/fixture/build_fixture.py
tests/parity/harness.sh \
    /tmp/freebayes-upstream/build/freebayes \
    ./target/debug/freebayes-gxy \
    tests/parity/fixture

# The harness: --call outputs VCF directly, bcftools norm canonicalises
# both, and a 5-column diff on (CHROM, POS, REF, ALT, GT) surfaces drift.
```

## Harness output (abridged)

```
==> running upstream freebayes...
==> running freebayes-gxy --call...
==> normalising via bcftools norm...

==> upstream calls (post-norm): 20
==> gxy calls      (post-norm): 20

==> PARITY DRIFT (post bcftools norm)
    site-level mismatches (CHROM+POS+REF+ALT): 0
    upstream-only lines (incl. GT diffs):      4
    gxy-only lines      (incl. GT diffs):      4
```

## Next steps

1. **Fix drift 2 (hom→het on hom indels).** Two candidates:
   - Phase C-3 pileup refactor: exclude reads that don't span the full
     indel anchor + event window from that site's observation set.
   - Full `buildHaplotypeAlleles` port: replace per-position pileup
     with per-haplotype-window pileup. Higher-effort, closer to
     upstream's internals.
2. **Expand INFO/FORMAT for real-world pipelines.** Specifically `GL`
   (the genotype log-likelihood triple), `QR`, `QA`, and at minimum
   `CIGAR`. These unlock `hap.py` / `vcfeval` style comparisons.
3. **Close the 1e-9 likelihood gate.** Emit the `log_posteriors`
   vector as `GL` in VCF format and compare to upstream.
4. **Larger fixture.** 10 kb slice of a real reference + simulated
   reads (wgsim) would surface scale / streaming issues.
