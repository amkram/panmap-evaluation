#!/usr/bin/env python3
"""Regenerate the Fig 2 random-placement baseline as ONE random node per placement
(n_random=1 semantics), recording each drawn node's classification so the plot can
color/shape the random points:
  RSV  -> subtype A / B          (meta/rsv_subtype.tsv)
  TB   -> MTBC species           (meta/tb_species.tsv)
  SARS -> (none; stays grey)

Random score = max(0, genome_distance(random_leaf, sample) - genome_distance(truth, sample)),
identical to a real placement's score but at a random node instead of panmap's pick.
Coverage-independent, so it reuses the surviving 0.5x intermediates.

Output: results/random_classed.tsv  (sp  kind  cov  dist  cls)
"""
import os
import re
import sys
import glob
import hashlib
import random as _random
from multiprocessing import Pool

sys.path.insert(0, "scripts")
import common as C

MM, SAM, EX = "minimap2", "samtools", 250
PANMAN = {"rsv": "../src/test/data/rsv_4K.panman",
          "sars": "../examples/data/sars_20000_twilight_dipper.panman",
          "tb": "/scratch1/alex/panstop/data/tb/tb_400.panman"}
TMP = "/tmp/claude-1015/-scratch1-alex-poopdoop/e32e7b24-91a2-49f4-a2e6-8a1f8e802306/scratchpad/randgen"
os.makedirs(TMP, exist_ok=True)

_G = {}                          # per-worker globals: leaves, classmap, genomes_fa


def load_class(sp):
    if sp == "rsv":
        return {r.split("\t")[0]: r.rstrip("\n").split("\t")[1]
                for r in open("meta/rsv_subtype.tsv").read().splitlines()[1:]}
    if sp == "tb":
        return {r.split("\t")[0]: r.rstrip("\n").split("\t")[1]
                for r in open("meta/tb_species.tsv").read().splitlines()[1:]}
    return {}


def header_node(fa):
    with open(fa) as f:
        return f.readline()[1:].strip()


def node_fa(node):
    o = os.path.join(TMP, _G["sp"] + "_g_" + str(abs(int(hashlib.md5(node.encode()).hexdigest()[:8], 16))) + ".fa")
    if not os.path.exists(o):
        C.faidx_seq(_G["genomes_fa"], node, o, SAM)
    return o


def one(task):
    """task = (kind, int_stem, cov). Draw 1 random leaf, score + classify."""
    kind, stem, cov = task
    is_sim = kind == "sim"
    sample_fa = stem + (".desc.fa" if is_sim else ".truth.fa")
    exp_fa = stem + ".exp.fa"
    sf = stem + (".leaf.fa" if is_sim else ".truth.fa")
    if not all(os.path.exists(p) for p in (sample_fa, exp_fa, sf)):
        return None
    self_node = header_node(sf)
    expected = self_node if is_sim else _G["par"].get(self_node)
    ed = C.genome_distance(MM, exp_fa, sample_fa, EX)
    lv = [n for n in _G["leaves"] if n != self_node and n != expected]
    key = f"{_G['sp']}/{kind}/{os.path.basename(stem)}/{cov}"
    rng = _random.Random(int(hashlib.md5(key.encode()).hexdigest()[:12], 16))
    rn = rng.choice(lv)
    d = max(0, C.genome_distance(MM, node_fa(rn), sample_fa, EX) - ed)
    cls = _G["classmap"].get(rn, "")
    return (_G["sp"], kind, cov, d, cls)


def init(sp, leaves, par, classmap, genomes_fa):
    _G.update(sp=sp, leaves=leaves, par=par, classmap=classmap, genomes_fa=genomes_fa)


def main():
    species = sys.argv[1:] or ["rsv", "sars", "tb"]
    out = open("results/random_classed.tsv", "a")
    for sp in species:
        par = C.newick_parents(PANMAN[sp], os.environ.get("PANMANUTILS", "panmanUtils"))
        leaves = sorted(C.leaves(par))
        classmap = load_class(sp)
        genomes_fa = f"work/{sp}/genomes.fa"
        # placements: every final per-sample tsv; map to its 0.5x intermediates
        tasks = []
        for kind, sub in [("sim", "fig2_sim"), ("real", "fig2_real")]:
            pat = (re.compile(r"^(\d+)_(\d+)_([\d.]+)\.tsv$") if kind == "sim"
                   else re.compile(r"^(\d+)_([\d.]+)\.tsv$"))
            for tsv in sorted(glob.glob(f"work/{sp}/{sub}/*.tsv")):
                m = pat.match(os.path.basename(tsv))
                if not m:
                    continue
                if kind == "sim":
                    i, mut, cov = m.groups()
                    int_stem = f"work/{sp}/{sub}/{i}_{mut}_0.5"   # coverage-independent intermediates
                else:
                    ri, cov = m.groups()
                    int_stem = f"work/{sp}/{sub}/{ri}_0.5"
                tasks.append((kind, int_stem, cov))
        with Pool(24, initializer=init, initargs=(sp, leaves, par, classmap, genomes_fa)) as p:
            for r in p.map(one, tasks):
                if r:
                    out.write("\t".join(str(x) for x in r) + "\n")
        out.flush()
        sys.stderr.write(f"{sp}: {len(tasks)} placements done\n")
    out.close()


if __name__ == "__main__":
    main()
