#!/usr/bin/env python3
"""Produce results/figure2_raw.tsv: the Fig 2 scores recomputed with the RAW edit
distance (indels weighted by base length, insertions included) instead of the
event-based distance. No placements are re-run: the placed node is re-derived from
each placement's --dump-all-scores (argmax logContainment, excluding self for real
LOO), exactly as the pipeline's place() does, then scored with genome_distance_raw.
Coverage-independent intermediates (desc/exp/truth) are reused from the 0.5x sibling.
The random column is likewise redrawn (pipeline seeding) and scored raw, keeping the
same taxon labels (RSV same/cross subtype, TB MTBC species).
"""
import os
import re
import sys
import glob
import yaml
import hashlib
import random as _random
from multiprocessing import Pool

sys.path.insert(0, "scripts")
import common as C

MM, SAM = "minimap2", "samtools"
CFG = yaml.safe_load(open("config.yaml"))
SEED, EX = CFG["seed"], CFG["exclude_bp"]
PANMAN = {sp: CFG["species"][sp]["panman"] for sp in CFG["species"]}
PU = os.path.join(os.path.dirname(CFG["panmap"]), "panmanUtils")
TMP = "/tmp/claude-1015/-scratch1-alex-poopdoop/e32e7b24-91a2-49f4-a2e6-8a1f8e802306/scratchpad/rawrescore"
os.makedirs(TMP, exist_ok=True)
DRAW = C.genome_distance_raw
_G = {}


def load_class(sp):
    f = {"rsv": "meta/rsv_subtype.tsv", "tb": "meta/tb_species.tsv"}.get(sp)
    m = {}
    if f and os.path.exists(f):
        for line in open(f).read().splitlines()[1:]:
            c = line.split("\t")
            if len(c) >= 2:
                m[c[0]] = c[1]
    return m


def hnode(fa):
    with open(fa) as f:
        return f.readline()[1:].strip()


def nodefa(node):
    o = os.path.join(TMP, _G["sp"] + "_" + str(abs(int(hashlib.md5(node.encode()).hexdigest()[:8], 16))) + ".fa")
    if not os.path.exists(o):
        C.faidx_seq(_G["genomes_fa"], node, o, SAM)
    return o


def pick_placed(scores_file, exclude):
    if not os.path.exists(scores_file):
        return None
    best, bv = None, None
    with open(scores_file) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        if "logContainment" not in hdr:
            return None
        mi = hdr.index("logContainment")
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) <= mi or c[0] == exclude:
                continue
            try:
                v = float(c[mi])
            except ValueError:
                continue
            if bv is None or v > bv:
                best, bv = c[0], v
    return best


def one(tsvpath):
    row = open(tsvpath).read().rstrip("\n").split("\t")
    if len(row) < 8:
        return None
    sp, kind, mut, cov, _score, _rand, wall, rss = row[:8]
    is_sim = kind == "sim"
    d = os.path.dirname(tsvpath)
    base = os.path.basename(tsvpath)[:-4]
    fcov = str(float(cov))                    # pipeline prefix uses float cov (1 -> 1.0)
    if is_sim:
        i, m = base.split("_")[0], base.split("_")[1]
        pre = f"{d}/{i}_{m}_{fcov}"; istem = f"{d}/{i}_{m}_0.5"
    else:
        ri = base.split("_")[0]
        pre = f"{d}/{ri}_{fcov}"; istem = f"{d}/{ri}_0.5"
    sample_fa = istem + (".desc.fa" if is_sim else ".truth.fa")
    exp_fa = istem + ".exp.fa"
    sf = istem + (".leaf.fa" if is_sim else ".truth.fa")
    if not all(os.path.exists(p) for p in (sample_fa, exp_fa, sf)):
        return None
    self_node = hnode(sf)
    ed = DRAW(MM, exp_fa, sample_fa, EX)
    # placement score (raw)
    placed = pick_placed(pre + ".scores.tsv", self_node if not is_sim else None)
    if placed is None:
        score = ""
    elif placed == self_node:
        score = 0                              # placed on truth leaf -> raw pd == ed
    else:
        score = max(0, DRAW(MM, nodefa(placed), sample_fa, EX) - ed)
    # random (raw), pipeline seeding + taxon label
    expected = self_node if is_sim else _G["par"].get(self_node)
    lv = [n for n in _G["leaves"] if n != (None if is_sim else self_node) and n != expected]
    rng = _random.Random(SEED + int(hashlib.md5(pre.encode()).hexdigest()[:8], 16))
    rn = rng.sample(lv, min(CFG["n_random"], len(lv)))[0]
    rd = max(0, DRAW(MM, nodefa(rn), sample_fa, EX) - ed)
    cm = _G["classmap"]
    if _G["sp"] == "rsv":
        ds, ss = cm.get(rn), cm.get(self_node)
        cls = "" if ds is None or ss is None else ("same" if ds == ss else "cross")
    else:
        cls = cm.get(rn, "")
    return f"{sp}\t{kind}\t{mut}\t{cov}\t{score}\t{rd}:{cls}\t{wall}\t{rss}\n"


def init(sp, leaves, par, classmap, genomes_fa):
    _G.update(sp=sp, leaves=leaves, par=par, classmap=classmap, genomes_fa=genomes_fa)


def main():
    out = open("results/figure2_raw.tsv", "w")
    out.write("species\tkind\tmut\tcoverage\tscore\trandom\twall_s\trss_mb\n")
    for sp in (sys.argv[1:] or list(CFG["species"])):
        par = C.newick_parents(PANMAN[sp], PU)
        leaves = sorted(C.leaves(par))
        cm = load_class(sp)
        gen = f"work/{sp}/genomes.fa"
        tsvs = []
        for sub, pat in [("fig2_sim", re.compile(r"^\d+_\d+_[\d.]+\.tsv$")),
                         ("fig2_real", re.compile(r"^\d+_[\d.]+\.tsv$"))]:
            tsvs += [t for t in sorted(glob.glob(f"work/{sp}/{sub}/*.tsv"))
                     if pat.match(os.path.basename(t))]
        with Pool(24, initializer=init, initargs=(sp, leaves, par, cm, gen)) as p:
            for r in p.map(one, tsvs):
                if r:
                    out.write(r)
        out.flush()
        sys.stderr.write(f"{sp}: {len(tsvs)} placements re-scored (raw)\n")
    out.close()


if __name__ == "__main__":
    main()
