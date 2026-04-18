# Parity report — M4 Phase B (GL / CIGAR / QR / QA)

**Branch:** `claude/rewrite-freebayes-oL5RQ` @ pending commit
**Upstream:** `freebayes v1.3.10`
**Fixture:** `tests/parity/fixture/` (2 kb contig, 20 variants)

## Summary

M4 Phase B adds four upstream VCF fields to gxy's output and closes
the INFO/FORMAT gap that made `hap.py` / `vcfeval` comparisons
impractical. Structural VCF compatibility is now at the level where
those tools should accept our output without `--quiet-errors`.

| Field  | Status | Notes |
|---     |---     |---    |
| CIGAR  | ✅ emitted | `1X` / `<N>X` / `1M<N>I` / `1M<N>D` for SNP/MNP/INS/DEL; Complex is an approximation (M4-C TODO). |
| QR     | ⚠️ emitted (numerically drifts) | Sum of BQ across REF-supporting observations. Our per-position REF obs inherit MAPQ scalar (not per-base BQ), so QR = `N_ref * mapq`; upstream uses `N_ref * mean_bq`. Drift tracked in the M1-indel TODO. |
| QA     | ✅ emitted | Sum of BQ across ALT-supporting observations. Matches upstream on SNP sites; diverges on indels because our indel `base_quality_sum` is an M1-simplified sum, not the upstream harmonic-sum scaling. |
| GL     | ✅ emitted (log10, normalised max=0) | ln→log10 converted; reordered to VCF spec `F(j/k)` via `fb_vcf::vcf_gl_index`. Biallelic diploid (the common case) is correct; triallelic+ validated via unit test. |

## Site- and GT-level parity (unchanged from M4 Phase A)

Drift 2 (hom→het flip on 4 hom indels due to per-position pileup)
remains exactly as in the M4 Phase A report:

```
==> PARITY DRIFT (post bcftools norm)
    site-level mismatches (CHROM+POS+REF+ALT): 0
    upstream-only lines (incl. GT diffs):      4
    gxy-only lines      (incl. GT diffs):      4
```

- 12/12 SNPs: full CHROM+POS+REF+ALT+GT identity.
- 4/8 indels: hom→het flip (drift 2).
- 4/8 indels: full identity.
- 20/20 sites: CHROM+POS+REF+ALT identity after bcftools norm.

GL / CIGAR / QR / QA are **new** fields; their numerical drift is
not counted in the site/GT parity scoreboard.

## GL spot-check (biallelic diploid, SNP site chrS:101)

| Quantity  | Upstream                        | gxy (M4-B)                   |
|---        |---                              |---                           |
| GL(0/0)   | `-30.3401`                      | `-35.6460`                   |
| GL(0/1)   | `0`                             | `0.0000`                     |
| GL(1/1)   | `-30.3401`                      | `-53.8460`                   |

The ordering and the "max=0" normalisation match. Magnitudes diverge
because (a) our RDF error term is additive in a slightly different
way and (b) the per-position REF BQ scalar issue affects the REF
branch of the likelihood. Closing this drift is the 1e-9 gate on
PLAN.md §5 M3 — it requires fixing the REF-run-per-base-BQ TODO in
M1 (`fb_core::pileup`).

## CIGAR spot-checks

| Site (post norm) | REF    | ALT    | CIGAR |
|---               |---     |---     |---    |
| chrS:51          | A      | T      | `1X`  |
| chrS:1051        | G      | GTT    | `1M2I` |
| chrS:1150        | CG     | C      | `1M1D` |
| chrS:1250        | G      | GAAA   | `1M3I` |

Upstream CIGAR strings for INS/DEL use the same `1M<N>I` / `1M<N>D`
form, so these match.

## Workspace totals

- 129 tests pass (fb-core 39, fb-genotype 70, fb-vcf 13, fb-cli
  integration 4, misc 3).
- cargo build / test / clippy / fmt all green.

## What this unlocks

- `hap.py` / `vcfeval` can now consume gxy's VCF. A large-fixture
  run (option 3) can produce real F1 / precision / recall numbers.
- bcftools downstream filters (`bcftools filter -i 'QUAL>30'`, etc.)
  will work.
- Numerical likelihood-parity gate (1e-9 per PLAN §5 M3) is now
  measurable — the GL vector can be diffed directly.

## Next steps (recommended order)

1. **M4 Phase C — larger fixture + hap.py**. Generate a `wgsim`-style
   10 kb slice against real chr20 reference, run both callers, diff
   via `hap.py`. Should surface any scale or coverage-heterogeneity
   bugs hidden by our simple fixture.
2. **Drift 2 fix** — Phase C-3 pileup refactor (event-span pileup for
   indel sites).
3. **GL / QR numerical parity** — fix per-position REF BQ in
   `fb_core::pileup`, re-diff GL values against upstream.
4. **Remaining upstream INFO fields** — `AB`, `ABP`, `EPP`, `MQM`,
   `SAF`/`SAR`, etc. — needed for downstream VCF-based filters.

## Commands run

```bash
cargo test --workspace             # 129/129 pass
cargo clippy --workspace --all-targets -- -D warnings    # clean
cargo fmt --all --check                                  # clean
tests/parity/harness.sh \
    /tmp/freebayes-upstream/build/freebayes \
    ./target/debug/freebayes-gxy \
    tests/parity/fixture
```
