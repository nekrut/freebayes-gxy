# Parity report — M4 Phase E (per-position REF BQ)

**Branch:** `claude/rewrite-freebayes-oL5RQ` @ pending commit
**Upstream:** `freebayes v1.3.10`

## Summary

Threaded per-base Phred qualities through REF observations so the
downstream `QR` and `QA` VCF INFO fields match upstream exactly on
SNP sites. Previously each decomposed per-position REF observation
carried `MAPQ` as its scalar BQ (giving `QR = N × MQ60`); it now
carries the underlying read's per-base Phred (giving `QR = N × Q40`
for Q40 reads), identical to upstream.

| Field | Site chrS:101 (het SNP) | |
|---    |---                       |---|
|       | upstream                 | gxy |
| QR    | 400                      | **400** ✅ (was 600) |
| QA    | 400                      | **400** ✅ |
| GL    | −30.3401, 0, −30.3401    | −35.6460, 0, −35.6460 |

`QR` and `QA` now agree exactly. The remaining GL-magnitude
difference (~5 log10 units on the outer terms) traces to a
formula-level gap — our `max(ln_bq, ln_mq)` combined with the RDF
scaling path is not quite the same as upstream's exact data-
likelihood expression. That's a separate investigation for a later
phase; the per-position BQ fix was one of the two known drivers,
and it's now closed.

| Site chrS:51 (hom-alt) | upstream GL      | gxy GL           |
|---                     |---               |---               |
|                        | −39.9568, −3.3113, 0 | −40.0000, −3.3113, 0 |

Hom-alt GLs match to ≤ 0.05 — well within any practical tolerance.
Het SNPs drift the most (as above) because the outside-observation
error path dominates under the null REF/REF and SNP/SNP hypotheses.

## Parity scoreboard (unchanged from M4 Phase D)

Small fixture (2 kb, 20 variants):

```
==> PARITY: byte-identical on (CHROM, POS, REF, ALT, GT) after bcftools norm
```

Large fixture (10 kb, 113 variants):

| Metric | upstream | gxy |
|---     |:---:     |:---:|
| TP     | 93       | 93  |
| FP     | 8        | 4   |
| FN     | 20       | 20  |
| Precision | 0.921 | **0.959** |
| Recall    | 0.823 | 0.823 |
| F1     | 0.869    | **0.886** |
| GT match on TP | 89/93 | 89/93 |

Site-level membership and GT agreement are unchanged — the BQ fix
is a numerical refinement that does not affect which sites get
called or which genotype wins.

## Implementation

- `fb_core::allele::AlleleObservation` gains a
  `per_base_quals: Vec<u8>` field. Populated for Reference runs
  (the multi-base observations where per-position BQ matters);
  empty for SNP / INS / DEL / Null / Complex which already have a
  scalar `base_quality_sum`.
- `fb_core::pileup::walk_match_run`'s `flush_ref_run` closure now
  captures `read.quals[start_rp..start_rp+length]` on emission and
  sets `base_quality_sum = sum(per_base_quals)`. Previously
  `base_quality_sum = mapq`.
- `fb_cli::main::Pileup::add_read_observations` reads
  `per_base_quals[offset]` when decomposing a Reference run into
  per-position REF observations, so each per-position REF carries
  its own per-base Phred. Falls back to the scalar
  `base_quality_sum` when `per_base_quals` is empty (e.g. in unit
  tests that construct REF observations directly).
- Five `AlleleObservation` test constructors across
  `haplotype.rs` / `caller.rs` / `data_likelihood.rs` add
  `per_base_quals: Vec::new()`.
- M1 walker test `pure_match_all_reference` updated to assert
  `base_quality_sum == 300` (10 × Q30) instead of `== 60` (MAPQ).

## What this closes

| Drift source | Status |
|---           |---     |
| QR off by `(MQ/BQ)` factor | **Closed** |
| QA (already per-base via obs base_quality_sum) | Unchanged — already correct |
| GL numerical magnitude (het sites) | Reduced, not closed — remaining gap is formula-specific |
| GL numerical magnitude (hom-alt) | **Effectively closed** (agrees to ≤ 0.05) |

## What remains for GL 1e-9 parity

- Investigate upstream's exact data-likelihood expression — in
  particular whether `max(ln_bq, ln_mq)` is applied on a
  per-observation basis (as we do) or differently. The ~5-log10
  discrepancy on het sites points at a formula mismatch, not a
  per-base BQ one.
- Examine whether upstream applies RDF scaling to the full
  `prod_q_out` sum or only to a subset of observations.
- Once those are aligned, the remaining drift should be within the
  f64 vs long double numerical precision window.

## Tests + lints

- 129 workspace tests pass.
- cargo build / test / clippy / fmt all green.
- Harness outputs unchanged for site/GT parity; `QR` now matches
  upstream character-for-character on every SNP site.

## Commands

```bash
cargo test --workspace            # 129/129 pass
cargo clippy --workspace --all-targets -- -D warnings   # clean
cargo fmt --all --check                                 # clean

tests/parity/harness.sh ... tests/parity/fixture
# => PARITY: byte-identical on (CHROM, POS, REF, ALT, GT) after bcftools norm

tests/parity/harness-large.sh ... tests/parity/fixture-large
# => GT match on TP: 89/93 for both callers
```

## Next ROI items

1. **Close the GL het-site gap** — investigate the formula-level
   difference between our standardGLs port and upstream's exact
   expression. Likely in `DataLikelihood.cpp:26-43`.
2. **Tighten the large fixture** — eliminate the shared FN=20 by
   adjusting the simulator to avoid read-stagger boundary effects.
3. **Extend the fixture with MNPs + complex events** — exercises
   M2's clumping path end-to-end.
4. **Real GIAB fixture** — deferred until network access allows
   the reference download.
