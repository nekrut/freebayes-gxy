#!/usr/bin/env python3
"""Compute hap.py-style TP / FP / FN / precision / recall / F1 for a
caller's VCF against a truth VCF.

Both VCFs must be bgzipped + indexed, and should be pre-normalised
via `bcftools norm -f ref.fa` for consistent allele representation.

Compares on (CHROM, POS, REF, ALT). GT agreement is reported as a
secondary metric but does not affect the primary precision/recall on
the site-level set.

Usage:
  compare_against_truth.py <truth.vcf.gz> <call.vcf.gz> <label>
"""

from __future__ import annotations

import sys
from pathlib import Path

import pysam


def load_vcf(path: Path):
    """Return {(chrom, pos, ref, alt): gt} keyed by allele identity."""
    out = {}
    with pysam.VariantFile(str(path)) as vcf:
        for rec in vcf:
            for alt in rec.alts or []:
                key = (rec.chrom, rec.pos, rec.ref, alt)
                # Pick the genotype of the first sample, or None if missing.
                gt = None
                if rec.samples:
                    sample = next(iter(rec.samples.values()))
                    if "GT" in sample:
                        gt_tuple = sample["GT"]
                        if all(g is not None for g in gt_tuple):
                            gt = "/".join(str(g) for g in sorted(gt_tuple))
                out[key] = gt
    return out


def kind_of(ref: str, alt: str) -> str:
    if len(ref) == len(alt) == 1:
        return "SNP"
    if len(ref) == len(alt):
        return "MNP"
    if len(ref) < len(alt):
        return "INS"
    return "DEL"


def main():
    if len(sys.argv) != 4:
        print("usage: compare_against_truth.py <truth.vcf.gz> <call.vcf.gz> <label>", file=sys.stderr)
        sys.exit(2)
    truth_path = Path(sys.argv[1])
    call_path = Path(sys.argv[2])
    label = sys.argv[3]

    truth = load_vcf(truth_path)
    calls = load_vcf(call_path)

    tp = set(truth) & set(calls)
    fn = set(truth) - set(calls)
    fp = set(calls) - set(truth)

    # GT agreement on the TP set.
    gt_match = sum(1 for k in tp if truth[k] is not None and truth[k] == calls[k])

    # Per-kind breakdown.
    by_kind = {}
    for kind in ("SNP", "MNP", "INS", "DEL"):
        t_k = {k for k in truth if kind_of(k[2], k[3]) == kind}
        c_k = {k for k in calls if kind_of(k[2], k[3]) == kind}
        tp_k = t_k & c_k
        fp_k = c_k - t_k
        fn_k = t_k - c_k
        by_kind[kind] = (len(t_k), len(c_k), len(tp_k), len(fp_k), len(fn_k))

    precision = len(tp) / (len(tp) + len(fp)) if (len(tp) + len(fp)) else 0.0
    recall = len(tp) / (len(tp) + len(fn)) if (len(tp) + len(fn)) else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0

    print(f"=== {label} ===")
    print(f"  truth  : {len(truth)} variants")
    print(f"  calls  : {len(calls)} variants")
    print(f"  TP     : {len(tp)}")
    print(f"  FP     : {len(fp)}")
    print(f"  FN     : {len(fn)}")
    print(f"  GT match on TP: {gt_match}/{len(tp)}")
    print(f"  precision      = {precision:.4f}")
    print(f"  recall         = {recall:.4f}")
    print(f"  F1             = {f1:.4f}")
    print("  per-kind (truth, called, TP, FP, FN):")
    for kind, (t, c, ntp, nfp, nfn) in by_kind.items():
        f1k = 2 * ntp / (ntp + nfp + ntp + nfn) if (ntp + nfp + ntp + nfn) else 0.0
        print(f"    {kind:4s}: truth={t:3d} called={c:3d} TP={ntp:3d} FP={nfp:3d} FN={nfn:3d} F1={f1k:.3f}")
    print()


if __name__ == "__main__":
    main()
