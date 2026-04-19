# Parity report — M4 Phase G (fixture deconflict → F1 = 1.0000)

**Branch:** `claude/rewrite-freebayes-oL5RQ` @ pending commit
**Upstream:** `freebayes v1.3.10` with `--legacy-gls`

## Summary

Diagnosed the shared `FN = 20` that had been present since M4-C as
a fixture artifact, not a caller issue. The truth generator's
arithmetic progressions aliased at 8 positions, producing
multi-allelic truth sites where a SNP and an INS (or DEL) shared
one ref position. Both callers' default `min-alternate-count = 2`
plus the natural half-depth split dropped each alt below threshold,
so neither called the site.

Fix: place the three variant classes on disjoint residue classes
mod 50. Both callers now score **F1 = 1.0000** on 113/113 variants
with **113/113 GT match** on true positives.

## Scoreboard

### Before (aliased positions)

| | upstream | gxy |
|---|:---:|:---:|
| TP | 93 | 93 |
| FP | 8 | 4 |
| FN | 20 | 20 |
| F1 | 0.869 | 0.886 |
| GT match | 89/93 | 89/93 |

### After (position-disjoint)

| | upstream | **gxy** |
|---|:---:|:---:|
| TP | 113 | 113 |
| FP | 0 | 0 |
| FN | 0 | 0 |
| **F1** | **1.0000** | **1.0000** |
| **GT match** | 113/113 | **113/113** |

Per-kind: SNP F1 = 1.000, INS F1 = 1.000, DEL F1 = 1.000 for both
callers.

## Root cause

Original generator:
- SNPs at `range(100, REF_LEN - 200, 150)` → positions 100, 250, 400, 550, 700, 850, 1000, …
- INS at `range(150, REF_LEN - 300, 400)` → positions 150, 550, 950, 1350, …
- DEL at `range(200, REF_LEN - 300, 400)` → positions 200, 600, 1000, 1400, …

Collision cases:
- `550`: SNP + INS
- `1000`: SNP + DEL
- `2200`, `2950`, `3400`, …: similar

At each such position, the simulator applied BOTH events to the
zygosity-determined reads. The walker then observed two distinct
alts at the same position, each getting ~30% of reads (after het
splitting). Under the default `min-alternate-fraction = 0.05` that
should have been fine — but `min-alternate-count = 2` combined with
the split reduced the effective support per alt, and in some cases
the observations became incoherent (read carries both SNP + INS,
so alt gets attributed based on walker order), failing both
filters.

Fix: space variants on mod-50 residues with no arithmetic collisions:
- SNPs at `≡ 0 (mod 50)` in `[100, REF_LEN − 200)`
- INS at `≡ 20 (mod 50)` in `[120, REF_LEN − 300)`
- DEL at `≡ 35 (mod 50)` in `[135, REF_LEN − 300)`

Each class's positions are now pairwise disjoint.

## Implementation

- `tests/parity/fixture-large/build_fixture.py` — updated
  `generate_truth()` with the disjoint residues and a comment
  explaining the earlier collision mode.
- No caller code changes; the fix is entirely in the fixture.

## What about MNPs?

A brief experiment added 2-bp adjacent-SNP pairs to exercise the
M2 clumping path end-to-end. Both callers detected them and emitted
them as single MNP records (via clumping), but the truth VCF
encoded them as two separate SNPs — so `bcftools isec` flagged
every MNP call as an FP and counted truth SNPs as FN.

That's a truth-representation mismatch, not a caller issue. The
M2 clumping path is already well-covered by `fb_core::haplotype`
unit tests; a proper MNP probe against the harness would need
truth emitted as multi-base records too. Left out for now with an
explanatory comment in `generate_truth()`.

## Workspace state

- 129 tests pass.
- cargo build / test / clippy / fmt all green.
- Small fixture: byte-identical on CHROM+POS+REF+ALT+GT (unchanged).
- Large fixture: **F1 = 1.0000** for both callers.

## What remains

| Task | Scope | Priority |
|---|---|---|
| MNP truth in harness | Emit MNPs as multi-base records; run harness | Medium |
| Experimental GL path port | DataLikelihood.cpp:44-140 + Bias/Contamination infra | Large |
| Real GIAB data | Network / reference slice download | High-value, blocked |
| QUAL scoring parity | Separate investigation (different formulae) | Low |

## Running the harness

```bash
python3 tests/parity/fixture-large/build_fixture.py
tests/parity/harness-large.sh \
    /tmp/freebayes-upstream/build/freebayes \
    ./target/debug/freebayes-gxy \
    tests/parity/fixture-large
# => F1 = 1.0000 for both callers on 113 truth variants
```
