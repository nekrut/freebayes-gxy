# Parity report — M3 (Phase C-2 HEAD) with extended fixture

**Branch:** `claude/rewrite-freebayes-oL5RQ` @ `a6fc480`
**Upstream:** `freebayes v1.3.10` built from `git clone --branch v1.3.10` at `/tmp/freebayes-upstream/build/freebayes`
**Fixture rev:** extended 2026-04-18 to add INS + DEL truth.

## Summary

| Variant class | Truth | Upstream calls | gxy calls | `(chrom, pos, ref, alt, gt_tag)` identical |
|---|---|---|---|---|
| SNP             | 12 (9 het + 3 hom)  | 12  | 12 | **12 / 12 ✅** |
| INS             |  4 (2 het + 2 hom)  |  4  |  4 |   0 /  4 ❌ |
| DEL             |  4 (2 het + 2 hom)  |  4  |  4 |   0 /  4 ❌ |
| **Total**       | **20**              | **20** | **20** | **12 / 20 (60%)** |

Both tools detect **every truth variant** and make **no spurious calls**;
the drift is entirely on indel REF/ALT representation and on genotype
zygosity for hom indels. Three known deferrals (documented in PLAN §5
and the code) explain every line of drift — none are surprises.

## SNP parity

Same as the M3 Phase C-1 report: every SNP site, every REF/ALT, every
diploid GT tag matches upstream on synthetic 20× Q40/MQ60 data. The
harness emits `chrS  101  A  G  REF/SNP` (and friends) identically from
both tools.

## Indel parity — drift catalog

### Drift 1: Anchor-base / position representation

Upstream pads indel REF and ALT with the preceding anchor base and one
trailing base, and reports the 1-based position of the anchor. gxy emits
the raw event at the 1-based position of the first variant base.

| Truth      | Upstream                        | gxy                  |
|---         |---                              |---                   |
| INS TT@1051| `1051  GTA    GTTTA`            | `1052  T     TT`     |
| INS AAA@1251| `1251 GTA    GAAATA`           | `1252  T     AAA`    |
| INS G@1451 | `1450  CGT    CGGT`             | `1452  T     G`      |
| INS CC@1651| `1651  GTA    GCCTA`            | `1652  T     CC`     |
| DEL 1@1150 | `1150  CGT    CT`               | `1151  G     .`      |
| DEL 2@1350 | `1350  CGTA   CA`               | `1351  G     .`      |
| DEL 1@1550 | `1550  CGT    CT`               | `1551  G     .`      |
| DEL 3@1750 | `1750  CGTAC  CC`               | `1751  G     .`      |

**Root cause**: our `--call` TSV emits the internal `Allele` shape
(`ref_seq` empty for INS, `alt_seq` empty for DEL — see
`fb_core::allele`), while upstream normalises to VCF convention. This is
tracked as [VCF anchor-base synthesis](../crates/fb-cli/src/main.rs) in
the M4 deferred list. `bcftools norm -f ref.fa` can round-trip either
form into a canonical representation.

**Fix location**: M4 VCF writer (`fb-vcf`).

### Drift 2: Hom → het flip on hom indels

Upstream calls the four hom indels as `1/1`; gxy calls them as
`REF/INS` / `REF/DEL`. The four het indels agree (both tools emit a
heterozygous call).

**Root cause**: our caller operates a **position-based pileup** — for a
hom DEL at ref pos P, every read that has a match operation at P
contributes a REF observation at P, even if that read started exactly
at P and hence does not actually span the deletion event. Upstream
operates an **event-span-based** pileup — a read must span the full
deletion anchor + event window to contribute.

Concretely in this fixture, the read `r_1550_0` starts exactly at
0-based pos 1550 with CIGAR `100M`. It does not span the deletion
`CGT → CT` at upstream's anchor 1549, so upstream counts 20 AO / 0 RO.
gxy's per-position pileup counts 1 REF obs (from that read's match op
at pos 1550) + 19 DEL obs, and the posterior puts ~95% mass on
`REF/DEL` vs. ~5% on `DEL/DEL` (GQ ≈ 13).

**Fix location**: either Phase C-3 pileup refactor (respect event-span
semantics) or M4 haplotype-window assembly port
(`AlleleParser::buildHaplotypeAlleles`). Both are in the deferred list.

### Drift 3: Harness normalisation

The harness' upstream-VCF normaliser maps any non-zero GT allele to the
tag `SNP`. For indel calls this makes upstream's GT look like `SNP` or
`REF/SNP` even when the underlying variant is INS/DEL. This is a
cosmetic harness quirk — it over-reports drift on the `gt_tag` column
for indels. Worth fixing in the harness before Phase C-3 so the drift
report distinguishes zygosity from variant class.

## Fixture (2026-04-18 extension)

- Contig: 2000 bp pseudo-random `chrS` (seed 42) with enforced
  alternating-base neighbourhoods around every indel anchor (prevents
  left-alignment ambiguity).
- **20 truth variants**: 12 SNPs + 4 INS + 4 DEL. Details in
  `tests/parity/fixture/build_fixture.py`.
- 381 single-end 100 bp reads with 5 bp stagger, MAPQ 60, Q40 bases.
- Known edge case: read `r_1550_0` starts at exactly the DEL position
  and does not carry the DEL. Kept intentionally — surfaces drift 2.

## Commands run

```bash
# Upstream built once at /tmp/freebayes-upstream/build/freebayes
python3 tests/parity/fixture/build_fixture.py
tests/parity/harness.sh \
    /tmp/freebayes-upstream/build/freebayes \
    ./target/debug/freebayes-gxy \
    tests/parity/fixture
```

## Parity-harness exit criteria (vs. PLAN §5 M3)

- **"genotype calls byte-identical"** on (pos, ref, alt, gt_tag) —
  **12 / 20 (60%)** on extended fixture. SNPs are 100%; indels are 0%
  because of drifts 1 + 2 above. Both drifts are understood and have
  assigned M4 fix locations.
- **"per-site likelihood values match upstream within 1e-9"** — still
  not measured. Needs a common log-likelihood export.

## Next steps

1. **Harness improvement** (cheap): infer variant class from REF/ALT
   lengths in the normaliser so drift-3 stops contaminating the gt_tag
   diff. Re-run to get the real indel drift number (spoiler: still
   nonzero, dominated by drifts 1 + 2).
2. **`bcftools norm` on both outputs** to eliminate drift 1 and measure
   how much of the remaining drift is just zygosity flips vs. genuine
   call disagreements.
3. **Phase C-3 or M4** — fix drifts 1 + 2 at the caller layer before
   re-running. The position offset fix is trivial (switch to upstream's
   anchor convention); the hom/het flip needs the event-span pileup
   refactor.
4. **Larger fixtures** (after the caller fixes): GIAB HG002 chr20 10 Mb
   slice against the published truth set.
