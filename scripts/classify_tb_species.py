#!/usr/bin/env python3
"""Coarse MTBC species for each TB panman leaf from genome distance to two anchors:
  M. canettii        if d(H37Rv) > 10000   (smooth tubercle bacilli sit ~20-40k
                                            SNPs out; nothing else exceeds ~5.5k)
  M. bovis / animal  if d(M. bovis) < 2000 (the animal-adapted clade is
                                            <=~1100 from the bovis ref, the rest >=2600)
  M. tuberculosis    otherwise
H37Rv is an external reference; the M. bovis anchor (AF2122/97) is a panman leaf.
A coarse 3-way split, not a fine lineage call (L1-L7 / africanum vs bovis
need a SNP barcode). Writes: node<TAB>species.

Usage: classify_tb_species.py <genomes.fa> <h37rv.fa> <out.tsv>
"""
import os
import sys
import tempfile
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C

GEN, H37, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
BOVIS = "LT708304.1"          # M. bovis AF2122/97 reference (a leaf in the tb_400 panman)
EX = 250
TMP = tempfile.mkdtemp(prefix="tbsp_")
names = [l.split()[0] for l in open(GEN + ".fai")]

have_bovis = BOVIS in names
if have_bovis:
    C.faidx_seq(GEN, BOVIS, os.path.join(TMP, "bov.fa"), "samtools")
else:
    sys.stderr.write(f"classify_tb_species: WARNING bovis anchor {BOVIS} not a leaf; "
                     "animal clade will not be separated\n")


def prof(idx_n):
    i, n = idx_n
    fa = os.path.join(TMP, f"{i}.fa")
    C.faidx_seq(GEN, n, fa, "samtools")
    dh = C.genome_distance("minimap2", H37, fa, EX)
    db = C.genome_distance("minimap2", os.path.join(TMP, "bov.fa"), fa, EX) if have_bovis else 1e9
    os.remove(fa)
    return (n, dh, db)


with Pool(min(24, os.cpu_count() or 4)) as pool:
    res = pool.map(prof, list(enumerate(names)))
with open(OUT, "w") as f:
    f.write("node\tspecies\n")
    for n, dh, db in res:
        sp = ("M. canettii" if dh > 10000 else
              "M. bovis / animal" if db < 2000 else "M. tuberculosis")
        f.write(f"{n}\t{sp}\n")
from collections import Counter
sys.stderr.write("classify_tb_species: " + str(dict(Counter(
    ("M. canettii" if dh > 10000 else "M. bovis / animal" if db < 2000 else "M. tuberculosis")
    for _, dh, db in res))) + f" -> {OUT}\n")
