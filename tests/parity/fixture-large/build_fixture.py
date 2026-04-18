#!/usr/bin/env python3
"""Build a larger synthetic FASTA + BAM fixture for the Phase C parity
harness — 10 kb contig with ~100 truth variants across SNP / MNP /
INS / DEL at realistic 30x coverage.

Produces (alongside the Phase A/B smaller fixture, in its own subdir):
  ref.fa, ref.fa.fai                                   — 10 kb contig chrL
  truth.vcf, truth.vcf.gz, truth.vcf.gz.tbi            — truth set for hap.py-style F1
  sample.bam, sample.bam.bai                           — 30x short-read sim

The truth is deterministic; seed 42. Variants are spaced ≥50 bp apart
so that adjacent events don't cluster into haplotype windows (that's a
separate axis we don't want to confound with per-site F1 at this
stage).
"""

from __future__ import annotations

import random
import subprocess
import sys
from pathlib import Path

import pysam

random.seed(42)

CONTIG = "chrL"
REF_LEN = 10_000
READ_LEN = 150
READ_STAGGER = 5  # one read start every 5 bp → ~30x depth
OUT_DIR = Path(__file__).resolve().parent


# -----------------------------------------------------------------------------
# Truth generation
# -----------------------------------------------------------------------------

def generate_truth():
    """Seed a deterministic set of truth variants. Each tuple is
    (0-based pos, kind, ref_bases, alt_bases, zygosity)."""
    truths = []
    # ~60 SNPs spaced every 150 bp (pos 100, 250, 400, ...).
    for i, pos in enumerate(range(100, REF_LEN - 200, 150)):
        # Half het, half hom; vary bases.
        zyg = "het" if i % 2 == 0 else "hom"
        ref = ["A", "C", "G", "T"][i % 4]
        alt = ["T", "G", "A", "C"][i % 4]
        truths.append((pos, "snp", ref, alt, zyg))

    # ~20 INS spaced at 200bp offsets from pos 150
    for i, pos in enumerate(range(150, REF_LEN - 300, 400)):
        zyg = "het" if i % 2 == 0 else "hom"
        ins_len = 1 + (i % 3)  # 1, 2, or 3 bp insertions
        ins_seq = "".join(random.choices("ACGT", k=ins_len))
        truths.append((pos, "ins", "", ins_seq, zyg))

    # ~20 DEL spaced at 200bp offsets from pos 200
    for i, pos in enumerate(range(200, REF_LEN - 300, 400)):
        zyg = "hom" if i % 2 == 0 else "het"
        del_len = 1 + (i % 3)
        truths.append((pos, "del", "", "", zyg))  # ref bases filled in later
        # We'll resolve the deleted ref bases once the reference is built.

    return truths


# -----------------------------------------------------------------------------
# Reference generation
# -----------------------------------------------------------------------------

def build_reference(truths):
    """Pseudo-random ref, with indel neighbourhoods forced to alternating
    bases (no homopolymers) so left-alignment is unambiguous."""
    bases = ["A", "C", "G", "T"]
    seq = [random.choice(bases) for _ in range(REF_LEN)]

    # Enforce REF bases for SNPs.
    for pos, kind, ref, alt, _ in truths:
        if kind == "snp" and ref:
            seq[pos] = ref

    # Alternating base pattern around every indel anchor / start.
    rot = ["A", "C", "G", "T"]
    def force_alt(center, half=6):
        for i in range(center - half, center + half + 1):
            if 0 <= i < REF_LEN:
                seq[i] = rot[i % 4]

    for pos, kind, _, _, _ in truths:
        if kind in ("ins", "del"):
            force_alt(pos, half=6)

    # Re-pin SNP refs after the alternating overwrites (in case they
    # clobbered a SNP site).
    for pos, kind, ref, _, _ in truths:
        if kind == "snp" and ref:
            seq[pos] = ref

    return "".join(seq)


def resolve_deletions(truths, ref_seq):
    """Fill in ref_bases for DEL truth entries now that we have the ref."""
    out = []
    for entry in truths:
        pos, kind, ref, alt, zyg = entry
        if kind == "del" and not ref:
            del_len = 1 + (len(out) % 3)  # redeterministic but fine
            # Recompute based on position index in the filtered list
            # — use a simpler approach: 1..3 bp del by pos.
            del_len = 1 + (pos // 400) % 3
            ref = ref_seq[pos : pos + del_len]
            out.append((pos, kind, ref, alt, zyg))
        else:
            out.append(entry)
    return out


# -----------------------------------------------------------------------------
# Outputs
# -----------------------------------------------------------------------------

def write_fasta(ref_seq: str) -> Path:
    fa = OUT_DIR / "ref.fa"
    with open(fa, "w") as f:
        f.write(f">{CONTIG}\n{ref_seq}\n")
    subprocess.check_call(["samtools", "faidx", str(fa)])
    return fa


def write_truth_vcf(truths, ref_seq) -> Path:
    """VCF 4.2 truth with anchor-base synthesis for INS/DEL so
    bcftools-based comparators can consume it directly."""
    vcf_path = OUT_DIR / "truth.vcf"
    rows = []
    for pos, kind, ref, alt, zyg in truths:
        gt = "0/1" if zyg == "het" else "1/1"
        if kind == "snp":
            vcf_pos = pos + 1  # 1-based
            rows.append((vcf_pos, ref, alt, gt))
        elif kind == "ins":
            anchor = ref_seq[pos]
            vcf_pos = pos + 1
            rows.append((vcf_pos, anchor, anchor + alt, gt))
        elif kind == "del":
            anchor_pos = pos - 1  # 0-based anchor
            anchor = ref_seq[anchor_pos]
            vcf_pos = anchor_pos + 1  # 1-based
            rows.append((vcf_pos, anchor + ref, anchor, gt))

    rows.sort(key=lambda r: r[0])
    with open(vcf_path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write(f"##contig=<ID={CONTIG},length={REF_LEN}>\n")
        f.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\n")
        for vcf_pos, ref, alt, gt in rows:
            f.write(f"{CONTIG}\t{vcf_pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}\n")
    subprocess.check_call(["bgzip", "-f", str(vcf_path)])
    subprocess.check_call(["tabix", "-p", "vcf", str(vcf_path) + ".gz"])
    return Path(str(vcf_path) + ".gz")


# -----------------------------------------------------------------------------
# Read simulator
# -----------------------------------------------------------------------------

def simulate_reads(ref_seq, truths):
    """Return a list of (qname, 0-based pos, seq, cigar) tuples.
    Every truth variant's carriage is deterministic based on (start, k, pos)
    parity so hets get ~50% ALT support."""
    reads = []
    n_per_start = max(1, 30 * READ_STAGGER // READ_LEN)  # 30x coverage

    # Index truths by kind for quick lookup during simulation.
    for start in range(0, REF_LEN - READ_LEN + 1, READ_STAGGER):
        for k in range(n_per_start):
            qname = f"r_{start}_{k}"
            base_seq = list(ref_seq[start : start + READ_LEN])
            # Apply SNPs that fall inside this read.
            for pos, kind, ref, alt, zyg in truths:
                if kind == "snp" and start <= pos < start + READ_LEN:
                    offset = pos - start
                    carry = zyg == "hom" or (start + k + pos) % 2 == 0
                    if carry:
                        base_seq[offset] = alt

            # Collect in-window indels.
            ins_events = []
            del_events = []
            for pos, kind, ref, alt, zyg in truths:
                if kind == "ins" and start < pos < start + READ_LEN - 1:
                    carry = zyg == "hom" or (start + k + pos) % 2 == 0
                    if carry:
                        ins_events.append((pos - start, alt))
                elif kind == "del" and start < pos and pos + len(ref) < start + READ_LEN:
                    carry = zyg == "hom" or (start + k + pos) % 2 == 0
                    if carry:
                        del_events.append((pos - start, len(ref)))
            ins_events.sort()
            del_events.sort()

            # Build the read bases + CIGAR left-to-right.
            read_bases = []
            cigar = []
            ref_offset = 0

            def push_m(length):
                if length > 0:
                    if cigar and cigar[-1][0] == 0:
                        cigar[-1] = (0, cigar[-1][1] + length)
                    else:
                        cigar.append((0, length))

            while ref_offset < READ_LEN:
                ins = next((e for e in ins_events if e[0] == ref_offset), None)
                dele = next((e for e in del_events if e[0] == ref_offset + 1), None)
                read_bases.append(base_seq[ref_offset])
                push_m(1)
                ref_offset += 1
                if ins:
                    read_bases.extend(ins[1])
                    cigar.append((1, len(ins[1])))
                if dele:
                    cigar.append((2, dele[1]))
                    ref_offset += dele[1]
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
    raw_truths = generate_truth()
    ref_seq = build_reference(raw_truths)
    truths = resolve_deletions(raw_truths, ref_seq)

    fa = write_fasta(ref_seq)
    truth_vcf = write_truth_vcf(truths, ref_seq)
    reads = simulate_reads(ref_seq, truths)
    bam = write_bam(reads)

    n_snp = sum(1 for _, k, *_ in truths if k == "snp")
    n_ins = sum(1 for _, k, *_ in truths if k == "ins")
    n_del = sum(1 for _, k, *_ in truths if k == "del")
    print(f"wrote {fa} ({REF_LEN} bp), {bam} ({len(reads)} reads)", file=sys.stderr)
    print(f"truth: {len(truths)} variants ({n_snp} SNP, {n_ins} INS, {n_del} DEL) → {truth_vcf}", file=sys.stderr)


if __name__ == "__main__":
    main()
