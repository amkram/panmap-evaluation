#!/usr/bin/env python3
"""Apply the n_random=1 random baseline NOW without re-running placements: for each
per-sample Fig 2 tsv, redraw the random column exactly as the updated _score_row
would (same SEED + md5(prefix) seeding, one draw), as "dist:taxon", and write it
back. Then results/figure2.tsv is rebuilt by concatenation (fig2_table format).

Reuses the coverage-independent 0.5x intermediates (desc/exp) for scoring.
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
TMP = "/tmp/claude-1015/-scratch1-alex-poopdoop/e32e7b24-91a2-49f4-a2e6-8a1f8e802306/scratchpad/patchrnd"
os.makedirs(TMP, exist_ok=True)
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


def patch(tsvpath):
    """Redraw one placement's random column and write it back into tsvpath."""
    row = open(tsvpath).read().rstrip("\n").split("\t")
    if len(row) < 8:
        return None
    sp, kind, mut, cov, score, _oldrand, wall, rss = row[:8]
    base = os.path.basename(tsvpath)[:-4]
    is_sim = kind == "sim"
    # coverage-independent intermediates live at the 0.5x sibling
    if is_sim:
        i, m = base.split("_")[0], base.split("_")[1]
        istem = os.path.join(os.path.dirname(tsvpath), f"{i}_{m}_0.5")
    else:
        ri = base.split("_")[0]
        istem = os.path.join(os.path.dirname(tsvpath), f"{ri}_0.5")
    sample_fa = istem + (".desc.fa" if is_sim else ".truth.fa")
    exp_fa = istem + ".exp.fa"
    sf = istem + (".leaf.fa" if is_sim else ".truth.fa")
    if not all(os.path.exists(p) for p in (sample_fa, exp_fa, sf)):
        return ("MISS", tsvpath)
    self_node = hnode(sf)
    expected = self_node if is_sim else _G["par"].get(self_node)
    exclude = None if is_sim else self_node
    ed = C.genome_distance(MM, exp_fa, sample_fa, EX)
    lv = [n for n in _G["leaves"] if n != exclude and n != expected]
    pre = tsvpath[:-4]                       # matches _score_row's `pre` string
    rng = _random.Random(SEED + int(hashlib.md5(pre.encode()).hexdigest()[:8], 16))
    rn = rng.sample(lv, min(CFG["n_random"], len(lv)))[0]
    d = max(0, C.genome_distance(MM, nodefa(rn), sample_fa, EX) - ed)
    cm = _G["classmap"]
    if _G["sp"] == "rsv":                          # 'same'/'cross' subtype vs the sample
        ds, ss = cm.get(rn), cm.get(self_node)
        cls = "" if ds is None or ss is None else ("same" if ds == ss else "cross")
    else:
        cls = cm.get(rn, "")
    rstr = f"{d}:{cls}"
    with open(tsvpath, "w") as f:
        f.write(f"{sp}\t{kind}\t{mut}\t{cov}\t{score}\t{rstr}\t{wall}\t{rss}\n")
    return ("OK", tsvpath)


def init(sp, leaves, par, classmap, genomes_fa):
    _G.update(sp=sp, leaves=leaves, par=par, classmap=classmap, genomes_fa=genomes_fa)


def main():
    for sp in (sys.argv[1:] or list(CFG["species"])):
        par = C.newick_parents(PANMAN[sp], PU)
        leaves = sorted(C.leaves(par))
        classmap = load_class(sp)
        genomes_fa = f"work/{sp}/genomes.fa"
        tsvs = []
        for kind, sub, pat in [("sim", "fig2_sim", re.compile(r"^\d+_\d+_[\d.]+\.tsv$")),
                               ("real", "fig2_real", re.compile(r"^\d+_[\d.]+\.tsv$"))]:
            tsvs += [t for t in sorted(glob.glob(f"work/{sp}/{sub}/*.tsv"))
                     if pat.match(os.path.basename(t))]
        miss = 0
        with Pool(24, initializer=init, initargs=(sp, leaves, par, classmap, genomes_fa)) as p:
            for r in p.map(patch, tsvs):
                if r and r[0] == "MISS":
                    miss += 1
        sys.stderr.write(f"{sp}: patched {len(tsvs)-miss}/{len(tsvs)} ({miss} missing intermediates)\n")


if __name__ == "__main__":
    main()
