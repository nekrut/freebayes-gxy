#!/usr/bin/env python3
"""Build a 100 kb benchmark fixture for M5 Phase B scaling
measurements. Same schema as ``fixture-large/build_fixture.py`` but
10x the contig length and proportionally more reads + variants so
parallel speedup is actually measurable.
"""

from __future__ import annotations

import random
import subprocess
import sys
from pathlib import Path

import pysam

random.seed(42)

CONTIG = "chrB"
REF_LEN = 100_000
READ_LEN = 150
READ_STAGGER = 5  # ~30× coverage
OUT_DIR = Path(__file__).resolve().parent


def generate_truth():
    truths = []
    # SNPs at 0-based positions ≡ 0 (mod 50) in [100, REF_LEN - 200).
    for i, pos in enumerate(range(100, REF_LEN - 200, 150)):
        zyg = "het" if i % 2 == 0 else "hom"
        ref = ["A", "C", "G", "T"][i % 4]
        alt = ["T", "G", "A", "C"][i % 4]
        truths.append((pos, "snp", ref, alt, zyg))
    # INS at ≡ 20 (mod 50).
    for i, pos in enumerate(range(120, REF_LEN - 300, 400)):
        zyg = "het" if i % 2 == 0 else "hom"
        ins_len = 1 + (i % 3)
        ins_seq = "".join(random.choices("ACGT", k=ins_len))
        truths.append((pos, "ins", "", ins_seq, zyg))
    # DEL at ≡ 35 (mod 50).
    for i, pos in enumerate(range(135, REF_LEN - 300, 400)):
        zyg = "hom" if i % 2 == 0 else "het"
        truths.append((pos, "del", "", "", zyg))
    return truths


def build_reference(truths):
    bases = ["A", "C", "G", "T"]
    seq = [random.choice(bases) for _ in range(REF_LEN)]
    for pos, kind, ref, _, _ in truths:
        if kind == "snp" and ref:
            seq[pos] = ref
    rot = ["A", "C", "G", "T"]
    def force_alt(center, half=6):
        for i in range(center - half, center + half + 1):
            if 0 <= i < REF_LEN:
                seq[i] = rot[i % 4]
    for pos, kind, _, _, _ in truths:
        if kind in ("ins", "del"):
            force_alt(pos, half=6)
    for pos, kind, ref, _, _ in truths:
        if kind == "snp" and ref:
            seq[pos] = ref
    return "".join(seq)


def resolve_deletions(truths, ref_seq):
    out = []
    for entry in truths:
        pos, kind, ref, alt, zyg = entry
        if kind == "del" and not ref:
            del_len = 1 + (pos // 400) % 3
            ref = ref_seq[pos : pos + del_len]
            out.append((pos, kind, ref, alt, zyg))
        else:
            out.append(entry)
    return out


def write_fasta(ref_seq):
    fa = OUT_DIR / "ref.fa"
    with open(fa, "w") as f:
        f.write(f">{CONTIG}\n{ref_seq}\n")
    subprocess.check_call(["samtools", "faidx", str(fa)])
    return fa


def write_truth_vcf(truths, ref_seq):
    vcf_path = OUT_DIR / "truth.vcf"
    rows = []
    for pos, kind, ref, alt, zyg in truths:
        gt = "0/1" if zyg == "het" else "1/1"
        if kind == "snp":
            rows.append((pos + 1, ref, alt, gt))
        elif kind == "ins":
            anchor = ref_seq[pos]
            rows.append((pos + 1, anchor, anchor + alt, gt))
        elif kind == "del":
            anchor_pos = pos - 1
            anchor = ref_seq[anchor_pos]
            rows.append((anchor_pos + 1, anchor + ref, anchor, gt))
    rows.sort(key=lambda r: r[0])
    with open(vcf_path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write(f"##contig=<ID={CONTIG},length={REF_LEN}>\n")
        f.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\n")
        for pos, ref, alt, gt in rows:
            f.write(f"{CONTIG}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}\n")
    subprocess.check_call(["bgzip", "-f", str(vcf_path)])
    subprocess.check_call(["tabix", "-p", "vcf", str(vcf_path) + ".gz"])
    return Path(str(vcf_path) + ".gz")


def simulate_reads(ref_seq, truths):
    reads = []
    n_per_start = max(1, 30 * READ_STAGGER // READ_LEN)
    for start in range(0, REF_LEN - READ_LEN + 1, READ_STAGGER):
        for k in range(n_per_start):
            qname = f"r_{start}_{k}"
            base_seq = list(ref_seq[start : start + READ_LEN])
            for pos, kind, ref, alt, zyg in truths:
                if kind == "snp" and start <= pos < start + READ_LEN:
                    offset = pos - start
                    carry = zyg == "hom" or (start + k + pos) % 2 == 0
                    if carry:
                        base_seq[offset] = alt
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


def write_bam(reads):
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
            a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
            a.flag = 0
            a.set_tag("RG", "synth")
            out.write(a)
    subprocess.check_call(["samtools", "index", str(bam_path)])
    return bam_path


def main():
    raw = generate_truth()
    ref_seq = build_reference(raw)
    truths = resolve_deletions(raw, ref_seq)
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
