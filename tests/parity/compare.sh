#!/usr/bin/env bash
# compare.sh — diff two VCFs with bcftools isec and print a parity summary.
#
# Usage: compare.sh <baseline.vcf[.gz]> <candidate.vcf[.gz]> [out_dir]
#
# This is the M0 skeleton used by the parity harness. It assumes both inputs
# are already sorted, compressed, and tabix-indexed (bcftools isec insists on
# that). See tests/parity/README.md for how to wire in a real freebayes
# v1.3.10 baseline run.

set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
    echo "usage: $0 <baseline.vcf[.gz]> <candidate.vcf[.gz]> [out_dir]" >&2
    exit 2
fi

BASELINE="$1"
CANDIDATE="$2"
OUTDIR="${3:-parity-out}"

if ! command -v bcftools >/dev/null 2>&1; then
    echo "error: bcftools not on PATH — install bcftools >= 1.17" >&2
    exit 127
fi

for f in "$BASELINE" "$CANDIDATE"; do
    [[ -f "$f" ]] || { echo "error: not a file: $f" >&2; exit 1; }
done

mkdir -p "$OUTDIR"

# isec writes: 0000.vcf (baseline-only), 0001.vcf (candidate-only),
# 0002.vcf (records identical in baseline), 0003.vcf (records identical in candidate).
bcftools isec -p "$OUTDIR" "$BASELINE" "$CANDIDATE"

only_baseline=$(grep -cv '^#' "$OUTDIR/0000.vcf" 2>/dev/null || echo 0)
only_candidate=$(grep -cv '^#' "$OUTDIR/0001.vcf" 2>/dev/null || echo 0)
shared=$(grep -cv '^#' "$OUTDIR/0002.vcf" 2>/dev/null || echo 0)

printf 'parity summary (%s vs %s)\n' "$BASELINE" "$CANDIDATE"
printf '  shared records      : %s\n' "$shared"
printf '  baseline-only       : %s\n' "$only_baseline"
printf '  candidate-only      : %s\n' "$only_candidate"
printf '  isec output dir     : %s\n' "$OUTDIR"
