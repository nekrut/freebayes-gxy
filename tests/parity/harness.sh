#!/usr/bin/env bash
# Parity harness for freebayes-gxy against upstream freebayes v1.3.10.
#
# Usage:
#   harness.sh <upstream_freebayes> <gxy_binary> <fixture_dir>
#
# Runs both callers against the fixture BAM, then normalises each
# output VCF with `bcftools norm -f ref.fa` (left-alignment + allele
# canonicalisation) and compares (CHROM, POS, REF, ALT, GT) tuples.
#
# Exit codes:
#   0  — byte-identical on (CHROM, POS, REF, ALT, GT) after normalisation
#   3  — parity drift (see $DIR/diff.txt for the full picture)

set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "usage: $0 <upstream_freebayes> <gxy_binary> <fixture_dir>" >&2
    exit 2
fi

UPSTREAM="$1"
GXY="$2"
DIR="$3"

[[ -x "$UPSTREAM" ]] || { echo "error: upstream binary not executable: $UPSTREAM" >&2; exit 1; }
[[ -x "$GXY" ]]      || { echo "error: gxy binary not executable: $GXY" >&2; exit 1; }
[[ -d "$DIR" ]]      || { echo "error: fixture dir not found: $DIR" >&2; exit 1; }

FA="$DIR/ref.fa"
BAM="$DIR/sample.bam"
[[ -f "$FA" && -f "$BAM" ]] || { echo "error: fixture missing ref.fa or sample.bam" >&2; exit 1; }

command -v bcftools >/dev/null 2>&1 || { echo "error: bcftools not on PATH" >&2; exit 1; }

UPSTREAM_VCF="$DIR/upstream.vcf"
GXY_VCF="$DIR/gxy.vcf"

echo "==> running upstream freebayes..."
"$UPSTREAM" -f "$FA" "$BAM" > "$UPSTREAM_VCF" 2> "$DIR/upstream.log"
echo "==> running freebayes-gxy --call..."
"$GXY" --call -f "$FA" "$BAM" > "$GXY_VCF" 2> "$DIR/gxy.log"

# Normalise via bcftools norm (left-align, split multi-allelics,
# canonicalise REF/ALT padding). This removes anchor-base and
# trailing-flanking-base representation differences between callers.
echo "==> normalising via bcftools norm..."
bcftools norm -f "$FA" "$UPSTREAM_VCF" 2> "$DIR/upstream.norm.log" \
    | grep -v '^##' \
    | awk -F'\t' '{print $1"\t"$2"\t"$4"\t"$5"\t"substr($10,1,3)}' \
    > "$DIR/upstream.norm.tsv"
bcftools norm -f "$FA" "$GXY_VCF" 2> "$DIR/gxy.norm.log" \
    | grep -v '^##' \
    | awk -F'\t' '{print $1"\t"$2"\t"$4"\t"$5"\t"substr($10,1,3)}' \
    > "$DIR/gxy.norm.tsv"

n_upstream=$(tail -n +2 "$DIR/upstream.norm.tsv" | wc -l)
n_gxy=$(tail -n +2 "$DIR/gxy.norm.tsv" | wc -l)
echo ""
echo "==> upstream calls (post-norm): $n_upstream"
echo "==> gxy calls      (post-norm): $n_gxy"

# Sort both by (CHROM, POS) for stable diff.
sort -t$'\t' -k1,1 -k2,2n "$DIR/upstream.norm.tsv" > "$DIR/upstream.sorted.tsv"
sort -t$'\t' -k1,1 -k2,2n "$DIR/gxy.norm.tsv"      > "$DIR/gxy.sorted.tsv"

# Site-level parity (CHROM+POS+REF+ALT, ignoring GT).
paste <(cut -f1-4 "$DIR/upstream.sorted.tsv") <(cut -f1-4 "$DIR/gxy.sorted.tsv") \
    | awk -F'\t' 'NR>1 && ($1!=$5 || $2!=$6 || $3!=$7 || $4!=$8) {print}' \
    > "$DIR/site.diff.txt" || true
n_site_diff=$(wc -l < "$DIR/site.diff.txt")

# Full (CHROM+POS+REF+ALT+GT) parity.
if diff "$DIR/upstream.sorted.tsv" "$DIR/gxy.sorted.tsv" > "$DIR/diff.txt"; then
    echo "==> PARITY: byte-identical on (CHROM, POS, REF, ALT, GT) after bcftools norm"
    rm -f "$DIR/diff.txt" "$DIR/site.diff.txt"
    exit 0
fi

n_upstream_only=$(grep -c '^<' "$DIR/diff.txt" || echo 0)
n_gxy_only=$(grep -c '^>' "$DIR/diff.txt" || echo 0)
echo ""
echo "==> PARITY DRIFT (post bcftools norm)"
echo "    site-level mismatches (CHROM+POS+REF+ALT): $n_site_diff"
echo "    upstream-only lines (incl. GT diffs):      $n_upstream_only"
echo "    gxy-only lines      (incl. GT diffs):      $n_gxy_only"
echo "    full diff: $DIR/diff.txt"
echo ""
echo "--- first 40 lines of diff ---"
head -40 "$DIR/diff.txt"
exit 3
