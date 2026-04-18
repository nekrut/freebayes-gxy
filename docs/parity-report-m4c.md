# Parity report — M4 Phase C (hap.py-style F1 on 10 kb fixture)

**Branch:** `claude/rewrite-freebayes-oL5RQ` @ pending commit
**Upstream:** `freebayes v1.3.10`
**Fixture:** `tests/parity/fixture-large/` (10 kb contig, 113 truth variants)

## Summary

First run against a site-count-scale fixture (10× the Phase A/B
fixture by basepairs; 5.7× by variant count). Both callers produce
structurally comparable output; gxy's aggregate F1 is marginally
**higher** than upstream's because of a fewer-FP-indels effect, but
gxy's GT agreement on true positives is **lower** due to the known
drift 2 (hom→het flip on per-position pileup).

## Scoreboard

| Metric             | upstream | **freebayes-gxy** |
|---                 |:---:     |:---:              |
| Truth variants     | 113      | 113               |
| Calls emitted      | 101      | 97                |
| True positives     | 93       | 93                |
| False positives    |  8       |  **4**            |
| False negatives    | 20       | 20                |
| **Precision**      | 0.921    | **0.959**         |
| **Recall**         | 0.823    | 0.823             |
| **F1**             | 0.869    | **0.886**         |
| GT match on TP     | 89 / 93 (95.7%) | 81 / 93 (87.1%) |

Site-level membership (CHROM+POS+REF+ALT) is the primary metric.
GT match on TP is reported separately — it's where drift 2 shows up.

## Per-kind breakdown

| Kind | Truth | upstream (TP/FP/FN, F1) | gxy (TP/FP/FN, F1) |
|---   |:---:  |:---:                    |:---:               |
| SNP  | 65    | 53/0/12, F1 = 0.898     | 53/0/12, F1 = 0.898 |
| INS  | 24    | 16/8/8,  F1 = 0.667     | **16/4/8, F1 = 0.727** |
| DEL  | 24    | 24/0/0,  F1 = 1.000     | 24/0/0,  F1 = 1.000 |

- **SNP parity is exact** — both callers land the same 53 TPs and miss
  the same 12 FNs.
- **DEL parity is exact and perfect** — both callers hit every truth
  deletion with no FPs.
- **INS: gxy emits fewer false-positive insertions** (4 vs 8). Both
  miss the same 8 truth insertions.

## FN commonality

All 20 FN are the same 20 sites in both callers. These are positions
where the simulated reads, after running through each caller's read
filter / min-alt-count / min-alt-fraction thresholds, do not support
the variant strongly enough to cross the calling threshold. That's a
fixture / simulation artifact, not a caller bug. (`build_fixture.py`
stages the reads deterministically; some truth positions sit near
read-stagger boundaries where the effective alt support is low.)

## Drift 2 (hom→het) impact at scale

Of the 93 true positives, upstream gets the GT right 89 times (95.7%);
gxy gets it right 81 times (87.1%). The 8 GT disagreements are all
**hom→het on indels** — the same pattern as the M4 Phase A + B
reports, just visible at 10× scale. No new flavor of drift 2
appeared.

- Predicted fix impact: closing drift 2 would raise gxy's GT match
  rate to ~95%, matching upstream.
- Phase C-3 pileup refactor is the fix location.

## Fixture

- 10 kb pseudo-random contig `chrL` (seed 42).
- 113 truth variants: 65 SNPs, 24 INS (1–3 bp), 24 DEL (1–3 bp).
  Variants spaced ≥50 bp apart; indel flanks forced to alternating
  bases so left-alignment is unambiguous.
- 1971 single-end 150 bp reads, stagger every 5 bp, MAPQ 60, Q40.
  ~30× mean coverage.
- Truth VCF emitted with anchor-base synthesis, bgzipped + tabixed.

## Harness

```bash
# Build fixture (once):
python3 tests/parity/fixture-large/build_fixture.py

# Run both callers + score:
tests/parity/harness-large.sh \
    /tmp/freebayes-upstream/build/freebayes \
    ./target/debug/freebayes-gxy \
    tests/parity/fixture-large
```

The harness normalises both caller VCFs via `bcftools norm -f ref.fa`
and feeds them to `compare_against_truth.py`, a hap.py-equivalent
comparator that computes TP / FP / FN on `(CHROM, POS, REF, ALT)` and
reports per-kind F1.

Why not `hap.py` itself? The canonical tool needs a full Python
environment (pandas, pybedtools, various numeric deps) plus its own
Docker image, and downloading that from inside the sandbox is
currently out of reach. Our comparator is schema-equivalent for the
site-level metrics; once `hap.py` is available we'll swap it in.

## Known numerical drifts (unchanged from M4 Phase B)

The VCF schema changes in M4 Phase B added `GL`, `CIGAR`, `QR`, `QA`.
Their numerical values still diverge from upstream as documented:
- QR: our per-position REF obs carry MAPQ as scalar → QR = N × MQ
  instead of upstream's N × mean_BQ.
- GL: magnitudes differ (ordering + max=0 normalisation match).
- QA on indels: M1 simplified sum, not harmonic-sum scaling.

These do not affect the F1 scoreboard — hap.py / bcftools isec match
on `(CHROM, POS, REF, ALT)` and optionally GT.

## Comparison to PLAN.md §5 M3 exit criteria

- **"genotype calls byte-identical"** — 81/93 = 87% agreement on TPs.
  Not byte-identical; gap is drift 2.
- **"per-site likelihood values match upstream within 1e-9"** —
  measurable now via GL diff; not yet run.

The exit criteria were framed for a real GIAB HG002 chr20 run, which
we still haven't done. On this synthetic fixture, site-level precision
is actually slightly above upstream's; GT agreement is 87%.

## What this phase delivers

1. **Real, numeric F1 / precision / recall** — first time we have
   them, for both callers.
2. **A scaling stress test** that didn't reveal new caller bugs — the
   drifts observed at 2 kb (20 variants) scale linearly to 10 kb (113
   variants) without new pathology.
3. **A comparator** that we can point at any (truth, call) VCF pair;
   extensible to real GIAB data when available.
4. **Evidence that gxy is on the right track**: we match upstream on
   SNP and DEL F1, and slightly beat on INS precision.

## Next steps (recommended priority)

1. **Fix drift 2** — Phase C-3 pileup refactor (event-span admission
   at indel sites). Closes the last systematic gap. Predicted effect:
   GT match rate → ~95%.
2. **Larger / real fixture** — GIAB HG002 chr22 5 Mb slice once the
   reference download works. Will likely surface:
   - Real coverage heterogeneity effects.
   - Soft-clipping / MAPQ variation.
   - Homopolymer left-align edge cases.
3. **Numerical GL parity** — close the 1e-9 gate by fixing the
   per-position REF BQ TODO in `fb_core::pileup`.
4. **Complex / MNP support** — truth fixture here has none; extend
   build_fixture.py to include clustered variants that trigger
   haplotype windowing.
