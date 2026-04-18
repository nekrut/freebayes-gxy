#!/usr/bin/env python3
"""Build a synthetic FASTA + BAM fixture for the parity harness.

Produces:
  ref.fa, ref.fa.fai             — single 1kb contig `chrS`
  truth.tsv                      — ground-truth variants (chrom, pos, ref, alt, gt)
  sample.bam, sample.bam.bai     — 30x short-read simulation covering
                                   the truth variants with Q30 bases
                                   and MAPQ 60

The truth is deterministic and simple:
  - 10 heterozygous SNPs at evenly-spaced positions (100..1000 step 100)
  - 3 homozygous SNPs interleaved (at positions 50, 550, 950)
  - No indels in this first pass — indel parity is a separate, harder probe
    covered by Phase C-3.

Reads are 100bp PE-style singletons (not paired for simplicity), 30x depth.
For each read covering a het site, half the reads carry REF and half ALT.
For hom sites, all reads carry ALT.
"""

from __future__ import annotations

import os
import random
import subprocess
import sys
from pathlib import Path

import pysam

random.seed(42)

CONTIG = "chrS"
REF_LEN = 1000
READ_LEN = 100
DEPTH = 30
OUT_DIR = Path(__file__).resolve().parent

HET_SNPS = [(p, "A", "G") for p in range(100, 1000, 100)]  # 9 sites: 100..900
HOM_SNPS = [(50, "A", "T"), (550, "C", "G"), (950, "G", "A")]  # 3 sites

# All mutated positions (0-based) and their truth
TRUTH = {}
for pos, ref, alt in HET_SNPS:
    TRUTH[pos] = (ref, alt, "het")
for pos, ref, alt in HOM_SNPS:
    TRUTH[pos] = (ref, alt, "hom")


def build_reference() -> str:
    # Deterministic pseudo-random ref so the truth REF bases are set correctly.
    bases = ["A", "C", "G", "T"]
    seq = [random.choice(bases) for _ in range(REF_LEN)]
    # Overwrite positions listed in TRUTH with the expected REF bases so
    # truth REF matches the reference.
    for pos, (ref, _alt, _) in TRUTH.items():
        seq[pos] = ref
    return "".join(seq)


def write_fasta(ref_seq: str) -> Path:
    fa = OUT_DIR / "ref.fa"
    with open(fa, "w") as f:
        f.write(f">{CONTIG}\n{ref_seq}\n")
    subprocess.check_call(["samtools", "faidx", str(fa)])
    return fa


def write_truth():
    with open(OUT_DIR / "truth.tsv", "w") as f:
        f.write("chrom\tpos\tref\talt\tgt\n")
        for pos in sorted(TRUTH):
            ref, alt, gt = TRUTH[pos]
            # 1-based pos for display
            f.write(f"{CONTIG}\t{pos + 1}\t{ref}\t{alt}\t{gt}\n")


def build_reads(ref_seq: str) -> list[tuple[str, int, str]]:
    """Return a list of (qname, 0-based pos, seq) tuples at DEPTH coverage."""
    reads = []
    # For every start position, emit DEPTH reads. Each read is READ_LEN bases.
    for start in range(0, REF_LEN - READ_LEN + 1, 5):  # 5bp stagger
        # Each start gets DEPTH/20 reads = ~1.5 reads per start → 30x overall
        # (every position is covered by READ_LEN/5 = 20 staggered windows,
        # each emitting DEPTH/20 = 1.5 reads; round to 2 for DEPTH=30.)
        n = max(1, DEPTH // (READ_LEN // 5))
        for k in range(n):
            qname = f"r_{start}_{k}"
            # Build the read seq against the ref. Apply variants probabilistically.
            base_seq = list(ref_seq[start : start + READ_LEN])
            for pos, (ref, alt, gt) in TRUTH.items():
                if start <= pos < start + READ_LEN:
                    offset = pos - start
                    if gt == "hom":
                        # Every read carries ALT
                        base_seq[offset] = alt
                    else:  # het
                        # Half the reads carry ALT, deterministic by (start, k, pos).
                        if (start + k + pos) % 2 == 0:
                            base_seq[offset] = alt
                        # else leave as REF
            reads.append((qname, start, "".join(base_seq)))
    return reads


def write_bam(reads: list[tuple[str, int, str]]) -> Path:
    bam_path = OUT_DIR / "sample.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": CONTIG, "LN": REF_LEN}],
        "RG": [{"ID": "synth", "SM": "SAMPLE", "PL": "ILLUMINA", "LB": "lib1"}],
    }
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out:
        # Sort reads by position for coordinate-sorted output.
        for qname, pos, seq in sorted(reads, key=lambda r: r[1]):
            a = pysam.AlignedSegment()
            a.query_name = qname
            a.query_sequence = seq
            a.reference_id = 0
            a.reference_start = pos
            a.mapping_quality = 60
            a.cigar = [(0, len(seq))]  # 0 = BAM_CMATCH
            a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))  # Q40
            a.flag = 0
            a.set_tag("RG", "synth")
            out.write(a)
    subprocess.check_call(["samtools", "index", str(bam_path)])
    return bam_path


def main():
    ref_seq = build_reference()
    fa = write_fasta(ref_seq)
    write_truth()
    reads = build_reads(ref_seq)
    bam = write_bam(reads)
    # Sanity: print coverage at one het site.
    print(f"wrote {fa} ({REF_LEN} bp), {bam} ({len(reads)} reads)", file=sys.stderr)
    print(f"truth: {len(TRUTH)} sites ({len(HET_SNPS)} het, {len(HOM_SNPS)} hom)", file=sys.stderr)


if __name__ == "__main__":
    main()
