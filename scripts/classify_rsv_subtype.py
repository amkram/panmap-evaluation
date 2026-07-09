#!/usr/bin/env python3
"""Assign each RSV panman leaf to subtype A or B by COMPETITIVE mapping: align the
leaf against a combined A+B reference and take the subtype whose reference accrues
more matching bases. This is robust where the whole-genome A-vs-B distance (~5%) is
comparable to within-subtype diversity, which defeats nearest-single-reference
classification (the tight within-B / far cross-A-B structure only emerges with the
competitive assignment). Writes: node<TAB>subtype.

Usage: classify_rsv_subtype.py <genomes.fa> <refA.fa> <refB.fa> <out.tsv>
"""
import os
import sys
import tempfile
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C

GEN, REFA, REFB, OUT = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
TMP = tempfile.mkdtemp(prefix="rsvsub_")


def _seq(p):
    return "".join(l.strip() for l in open(p) if not l.startswith(">"))


COMB = os.path.join(TMP, "AB.fa")
with open(COMB, "w") as o:                       # combined reference, contigs named A / B
    o.write(">A\n" + _seq(REFA) + "\n>B\n" + _seq(REFB) + "\n")
names = [l.split()[0] for l in open(GEN + ".fai")]


def classify(idx_n):
    i, n = idx_n
    fa = os.path.join(TMP, f"{i}.fa")
    C.faidx_seq(GEN, n, fa, "samtools")
    p = C.sh(["minimap2", "-cx", "asm20", COMB, fa])
    mA = mB = 0
    for line in p.stdout.splitlines():
        c = line.split("\t")
        if len(c) < 11:
            continue
        m = int(c[9])                            # PAF col 10: residue matches
        if c[5] == "A":
            mA += m
        elif c[5] == "B":
            mB += m
    os.remove(fa)
    return (n, "A" if mA >= mB else "B")


with Pool(min(24, os.cpu_count() or 4)) as pool:
    res = pool.map(classify, list(enumerate(names)))
with open(OUT, "w") as f:
    f.write("node\tsubtype\n")
    for n, s in res:
        f.write(f"{n}\t{s}\n")
sys.stderr.write(f"classify_rsv_subtype: A={sum(1 for _,s in res if s=='A')} "
                 f"B={sum(1 for _,s in res if s=='B')} -> {OUT}\n")
