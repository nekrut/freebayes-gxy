#!/usr/bin/env bash
# Phase C parity harness — run both callers on the 10 kb fixture and
# compute hap.py-style F1 / precision / recall against the truth VCF.
#
# Usage:
#   harness-large.sh <upstream_freebayes> <gxy_binary> <fixture_dir>

set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "usage: $0 <upstream_freebayes> <gxy_binary> <fixture_dir>" >&2
    exit 2
fi

UPSTREAM="$1"
GXY="$2"
DIR="$3"

[[ -x "$UPSTREAM" ]] || { echo "error: upstream not executable: $UPSTREAM" >&2; exit 1; }
[[ -x "$GXY" ]]      || { echo "error: gxy not executable: $GXY" >&2; exit 1; }
[[ -d "$DIR" ]]      || { echo "error: fixture dir missing: $DIR" >&2; exit 1; }

FA="$DIR/ref.fa"
BAM="$DIR/sample.bam"
TRUTH="$DIR/truth.vcf.gz"
[[ -f "$FA" && -f "$BAM" && -f "$TRUTH" ]] || {
    echo "error: fixture missing ref.fa / sample.bam / truth.vcf.gz" >&2
    exit 1
}

for bin in bcftools bgzip tabix python3; do
    command -v "$bin" >/dev/null 2>&1 || {
        echo "error: $bin not on PATH" >&2
        exit 1
    }
done

UPSTREAM_VCF="$DIR/upstream.vcf"
GXY_VCF="$DIR/gxy.vcf"

echo "==> running upstream (--legacy-gls for apples-to-apples standardGLs comparison)..."
"$UPSTREAM" --legacy-gls -f "$FA" "$BAM" > "$UPSTREAM_VCF" 2> "$DIR/upstream.log"
echo "==> running gxy..."
"$GXY" --call -f "$FA" "$BAM" > "$GXY_VCF" 2> "$DIR/gxy.log"

echo "==> normalising..."
for tag in upstream gxy truth; do
    src="$DIR/${tag}.vcf"
    # Always re-bgzip so stale .gz files don't silently shadow a fresh
    # re-run of the caller. bgzip -k keeps the uncompressed source.
    # truth.vcf.gz is the primary — regenerate only if the .vcf source
    # exists, else trust the existing .gz.
    if [[ -f "$src" ]]; then
        bgzip -kf "$src"
    fi
    tabix -p vcf -f "$src.gz"
    bcftools norm -f "$FA" -Oz "$src.gz" -o "$DIR/${tag}.norm.vcf.gz" 2> "$DIR/${tag}.norm.log"
    tabix -p vcf -f "$DIR/${tag}.norm.vcf.gz"
done

# bcftools norm regenerates truth from .gz; the source stays as is.
# Truth.vcf.gz was already built by build_fixture.py.

echo ""
python3 "$(dirname "$0")/compare_against_truth.py" \
    "$DIR/truth.norm.vcf.gz" "$DIR/upstream.norm.vcf.gz" upstream
python3 "$(dirname "$0")/compare_against_truth.py" \
    "$DIR/truth.norm.vcf.gz" "$DIR/gxy.norm.vcf.gz" freebayes-gxy
