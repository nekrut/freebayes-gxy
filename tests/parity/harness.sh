#!/usr/bin/env bash
# Parity harness for freebayes-gxy against upstream freebayes v1.3.10.
#
# Usage:
#   harness.sh <upstream_freebayes> <gxy_binary> <fixture_dir>
#
# Compares per-site (chrom, pos, ref, alt, GT) tuples between upstream's
# VCF and gxy's TSV. Not byte-level — we normalise both outputs into a
# common (chrom, pos, ref, alt, gt_tag) schema and diff them.

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

UPSTREAM_VCF="$DIR/upstream.vcf"
GXY_TSV="$DIR/gxy.tsv"

echo "==> running upstream freebayes..."
"$UPSTREAM" -f "$FA" "$BAM" > "$UPSTREAM_VCF" 2> "$DIR/upstream.log"
echo "==> running freebayes-gxy --call..."
"$GXY" --call -f "$FA" "$BAM" > "$GXY_TSV" 2> "$DIR/gxy.log"

# Normalise upstream VCF into (chrom, pos, ref, alt, gt_tag) rows.
#   gt_tag: REF (hom-ref), SNP (hom-alt), REF/SNP (het for SNP), etc.
#   GT field 0/0 -> REF, 1/1 -> SNP, 0/1 -> REF/SNP, 1/0 -> REF/SNP.
python3 - "$UPSTREAM_VCF" > "$DIR/upstream.normalized.tsv" <<'PY'
import sys
vcf = sys.argv[1]
out = ["chrom\tpos\tref\talt\tgt_tag"]
with open(vcf) as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        chrom, pos, _id, ref, alt, _qual, _filt, _info, fmt, sample = parts[:10]
        fmt_fields = fmt.split(":")
        sample_fields = sample.split(":")
        gt = sample_fields[fmt_fields.index("GT")]
        # Normalise GT syntax to a gt_tag matching gxy's str_tag output.
        # Upstream emits variant sites only (not hom-ref), so we expect 0/1 or 1/1.
        allele_tag = {"0": "REF", "1": "SNP"}  # multi-allelic sites out of scope here
        a, b = gt.replace("|", "/").split("/")
        tag_a = allele_tag.get(a, "?")
        tag_b = allele_tag.get(b, "?")
        if tag_a == tag_b:
            gt_tag = tag_a
        else:
            # Canonical order matches gxy's sort: REF before SNP.
            gt_tag = "REF/SNP" if "REF" in (tag_a, tag_b) and "SNP" in (tag_a, tag_b) else f"{tag_a}/{tag_b}"
        out.append(f"{chrom}\t{pos}\t{ref}\t{alt}\t{gt_tag}")
print("\n".join(out))
PY

# Normalise gxy TSV into the same schema.
awk -F'\t' 'NR==1 {print "chrom\tpos\tref\talt\tgt_tag"; next} {print $1"\t"$2"\t"$3"\t"$4"\t"$5}' \
    "$GXY_TSV" > "$DIR/gxy.normalized.tsv"

echo ""
echo "==> upstream calls: $(tail -n +2 "$DIR/upstream.normalized.tsv" | wc -l)"
echo "==> gxy calls:      $(tail -n +2 "$DIR/gxy.normalized.tsv" | wc -l)"
echo ""

# Diff the two normalised tables (sorted for stable diff).
sort "$DIR/upstream.normalized.tsv" > "$DIR/upstream.sorted.tsv"
sort "$DIR/gxy.normalized.tsv"      > "$DIR/gxy.sorted.tsv"

if diff "$DIR/upstream.sorted.tsv" "$DIR/gxy.sorted.tsv" > "$DIR/diff.txt"; then
    echo "==> PARITY: byte-identical on (chrom, pos, ref, alt, gt_tag)"
    rm -f "$DIR/diff.txt"
    exit 0
else
    # Count distinct lines that differ: both '<' and '>' are emitted
    # per differing line in plain-diff output.
    n_upstream_only=$(grep -c '^<' "$DIR/diff.txt" || echo 0)
    n_gxy_only=$(grep -c '^>' "$DIR/diff.txt" || echo 0)
    echo "==> PARITY DRIFT"
    echo "    upstream-only lines: $n_upstream_only"
    echo "    gxy-only lines:      $n_gxy_only"
    echo "    full diff at $DIR/diff.txt"
    echo ""
    echo "--- first 40 lines of diff ---"
    head -40 "$DIR/diff.txt"
    exit 3
fi
