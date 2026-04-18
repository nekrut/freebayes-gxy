#!/usr/bin/env python3
"""Build a synthetic FASTA + BAM fixture for the parity harness.

Produces:
  ref.fa, ref.fa.fai             — single 2kb contig `chrS`
  truth.tsv                      — ground-truth variants (chrom, pos, ref, alt, gt)
  sample.bam, sample.bam.bai     — ~20x coverage simulation covering
                                   the truth variants with Q40 bases
                                   and MAPQ 60

Truth set:
  SNPs (positions 0-based):
    - 9 het SNPs   at 100, 200, ..., 900  (A->G)
    - 3 hom SNPs   at 50, 550, 950        (A->T / C->G / G->A)
  Insertions (anchor = ref position of the base immediately preceding):
    - 2 het INS    at 1050 (insert "TT"), 1250 (insert "AAA")
    - 2 hom INS    at 1450 (insert "G"),  1650 (insert "CC")
  Deletions (anchor = 0-based ref position of first deleted base):
    - 2 het DEL    at 1150 (remove 1 bp), 1350 (remove 2 bp)
    - 2 hom DEL    at 1550 (remove 1 bp), 1750 (remove 3 bp)

All indels are placed in homopolymer-free regions (see the
ref-generator guard below) so upstream and the port don't disagree on
left-alignment — we are not testing that axis yet. The reference is
generated first, variants overwrite ref bases as needed, and the read
simulator applies CIGAR / sequence edits on a per-read basis.
"""

from __future__ import annotations

import random
import subprocess
import sys
from pathlib import Path

import pysam

random.seed(42)

CONTIG = "chrS"
REF_LEN = 2000
READ_LEN = 100
READ_STAGGER = 5  # one read started every 5 bp → ~20x coverage
OUT_DIR = Path(__file__).resolve().parent


# SNPs: (0-based pos, ref, alt, zygosity)
SNPS = [
    (50, "A", "T", "hom"),
    (100, "A", "G", "het"),
    (200, "A", "G", "het"),
    (300, "A", "G", "het"),
    (400, "A", "G", "het"),
    (500, "A", "G", "het"),
    (550, "C", "G", "hom"),
    (600, "A", "G", "het"),
    (700, "A", "G", "het"),
    (800, "A", "G", "het"),
    (900, "A", "G", "het"),
    (950, "G", "A", "hom"),
]

# Insertions: (0-based anchor pos, inserted seq, zygosity)
# Anchor = the ref position of the base immediately preceding the inserted bases.
# freebayes emits INS at 1-based `anchor + 1` with REF = anchor base, ALT =
# anchor base + inserted bases (anchor-base synthesis).
INSERTIONS = [
    (1050, "TT", "het"),
    (1250, "AAA", "het"),
    (1450, "G", "hom"),
    (1650, "CC", "hom"),
]

# Deletions: (0-based start pos of first deleted base, length, zygosity).
# freebayes emits DEL at 1-based `start`, REF = anchor + deleted bases, ALT =
# anchor base. We pick the anchor as `start - 1`.
DELETIONS = [
    (1150, 1, "het"),
    (1350, 2, "het"),
    (1550, 1, "hom"),
    (1750, 3, "hom"),
]


def build_reference() -> str:
    """Generate a reference that guarantees non-homopolymer neighbourhoods
    around every indel anchor, so left-alignment is unambiguous."""
    bases = ["A", "C", "G", "T"]
    seq = [random.choice(bases) for _ in range(REF_LEN)]
    # Pin SNP ref bases.
    for pos, ref, _alt, _ in SNPS:
        seq[pos] = ref

    # For each indel, force the 5 bp flanking both sides to alternate
    # bases so the inserted/deleted run never matches. Cheap and
    # avoids tripping upstream's left-align normalisation.
    def force_alternating(center: int, half: int = 5) -> None:
        rot = ["A", "C", "G", "T"]
        for i in range(center - half, center + half + 1):
            if 0 <= i < REF_LEN:
                seq[i] = rot[i % 4]

    for anchor, _ins, _ in INSERTIONS:
        force_alternating(anchor)
    for start, length, _ in DELETIONS:
        force_alternating(start + length // 2, half=5 + length)

    # Re-pin SNPs after the alternating overwrites in case they
    # clobbered a SNP site.
    for pos, ref, _alt, _ in SNPS:
        seq[pos] = ref
    return "".join(seq)


def write_fasta(ref_seq: str) -> Path:
    fa = OUT_DIR / "ref.fa"
    with open(fa, "w") as f:
        f.write(f">{CONTIG}\n{ref_seq}\n")
    subprocess.check_call(["samtools", "faidx", str(fa)])
    return fa


def write_truth(ref_seq: str) -> None:
    rows = []
    for pos, ref, alt, gt in SNPS:
        rows.append((pos + 1, ref, alt, gt, "snp"))
    for anchor, ins, gt in INSERTIONS:
        anchor_base = ref_seq[anchor]
        rows.append((anchor + 1, anchor_base, anchor_base + ins, gt, "ins"))
    for start, length, gt in DELETIONS:
        # Upstream's VCF convention: anchor = start - 1, REF includes anchor + deleted bases.
        anchor = start - 1
        anchor_base = ref_seq[anchor]
        deleted = ref_seq[start : start + length]
        rows.append((anchor + 1, anchor_base + deleted, anchor_base, gt, "del"))
    rows.sort(key=lambda r: r[0])
    with open(OUT_DIR / "truth.tsv", "w") as f:
        f.write("chrom\tpos\tref\talt\tgt\tkind\n")
        for pos, ref, alt, gt, kind in rows:
            f.write(f"{CONTIG}\t{pos}\t{ref}\t{alt}\t{gt}\t{kind}\n")


def apply_snp(base_seq: list[str], start: int, k: int) -> None:
    """Mutate the read sequence in place for SNPs that fall inside this read."""
    for pos, ref, alt, gt in SNPS:
        if start <= pos < start + READ_LEN:
            offset = pos - start
            if gt == "hom" or (start + k + pos) % 2 == 0:
                base_seq[offset] = alt


def build_reads(ref_seq: str):
    """Simulate reads. Each read is emitted as (qname, 0-based pos, seq, cigar)
    where `cigar` is a list of (op, len) tuples in pysam's numeric form
    (M=0, I=1, D=2)."""
    reads = []
    n_per_start = max(1, 20 * READ_STAGGER // READ_LEN)  # ~20x depth
    for start in range(0, REF_LEN - READ_LEN + 1, READ_STAGGER):
        for k in range(n_per_start):
            qname = f"r_{start}_{k}"
            base_seq = list(ref_seq[start : start + READ_LEN])
            apply_snp(base_seq, start, k)

            # Check for indels that overlap this read and decide whether
            # the read carries them (based on zygosity + deterministic
            # coin-flip on (start, k)).
            ins_events = []  # (ref_offset, ins_seq)
            del_events = []  # (ref_offset, del_len)

            for anchor, ins, gt in INSERTIONS:
                # Anchor must be fully inside the read window with some
                # slack for the inserted bases.
                if start <= anchor < start + READ_LEN:
                    carry = gt == "hom" or (start + k + anchor) % 2 == 0
                    if carry:
                        ins_events.append((anchor - start, ins))

            for del_start, length, gt in DELETIONS:
                # Need room to keep the whole deletion + at least one
                # anchor base on each side within the read.
                if start < del_start and del_start + length < start + READ_LEN:
                    carry = gt == "hom" or (start + k + del_start) % 2 == 0
                    if carry:
                        del_events.append((del_start - start, length))

            # Build the read bases + CIGAR by scanning left-to-right
            # through ref_seq[start:start+READ_LEN]. Insertions expand,
            # deletions skip.
            ins_events.sort()
            del_events.sort()
            read_bases: list[str] = []
            cigar: list[tuple[int, int]] = []
            ref_offset = 0

            def push_m(length: int) -> None:
                if length > 0:
                    if cigar and cigar[-1][0] == 0:
                        cigar[-1] = (0, cigar[-1][1] + length)
                    else:
                        cigar.append((0, length))

            while ref_offset < READ_LEN:
                # Insertion fires at the *end* of this ref_offset base
                # (after position ref_offset has been matched).
                ins = next((e for e in ins_events if e[0] == ref_offset), None)
                dele = next((e for e in del_events if e[0] == ref_offset + 1), None)

                # Consume the current ref base.
                read_bases.append(base_seq[ref_offset])
                push_m(1)
                ref_offset += 1

                if ins:
                    read_bases.extend(ins[1])
                    cigar.append((1, len(ins[1])))

                if dele:
                    length = dele[1]
                    cigar.append((2, length))
                    ref_offset += length  # skip those ref bases in the read

                if ref_offset >= READ_LEN:
                    break

            reads.append((qname, start, "".join(read_bases), cigar))
    return reads


def write_bam(reads) -> Path:
    bam_path = OUT_DIR / "sample.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": CONTIG, "LN": REF_LEN}],
        "RG": [{"ID": "synth", "SM": "SAMPLE", "PL": "ILLUMINA", "LB": "lib1"}],
    }
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out:
        for qname, pos, seq, cigar in sorted(reads, key=lambda r: r[1]):
            a = pysam.AlignedSegment()
            a.query_name = qname
            a.query_sequence = seq
            a.reference_id = 0
            a.reference_start = pos
            a.mapping_quality = 60
            a.cigar = cigar
            a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))  # Q40
            a.flag = 0
            a.set_tag("RG", "synth")
            out.write(a)
    subprocess.check_call(["samtools", "index", str(bam_path)])
    return bam_path


def main() -> None:
    ref_seq = build_reference()
    fa = write_fasta(ref_seq)
    write_truth(ref_seq)
    reads = build_reads(ref_seq)
    bam = write_bam(reads)
    print(
        f"wrote {fa} ({REF_LEN} bp), {bam} ({len(reads)} reads)",
        file=sys.stderr,
    )
    print(
        f"truth: {len(SNPS)} SNPs, {len(INSERTIONS)} INS, {len(DELETIONS)} DEL",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
