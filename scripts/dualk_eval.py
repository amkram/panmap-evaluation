#!/usr/bin/env python3
"""Dual-parameterization placement: combine dump-all-scores from index A (current
k/s/l) and index B (a second k/s/l) to resolve low-coverage ties.

For each sample: run panmap place with index_B, load A's cached scores + B's, and
select the placed node by several fusion rules, then score with genome_distance.

Fusion rules (all leaf-restricted; exclude self for real LOO):
  A_only     : argmax A.logContainment, first-in-file tie-break        (== current)
  B_only     : argmax B.logContainment
  tieBreakB  : argmax A.logContainment; break A ties by B.logContainment
  tieBreakA  : argmax B.logContainment; break B ties by A.logContainment
  rankSum    : min over nodes of rank_A + rank_B (by logContainment)
  sumContain : argmax (A.containment + B.containment)
  intersect  : first node in (A top-tie ∩ B top-tie); else A_only

Usage: dualk_eval.py <sp> <kind> <stem> [stem ...]  -> TSV: method sp kind mut cov score is_self
"""
import os
import sys

sys.path.insert(0, "scripts")
import common as C

MM, SAM, EX = "minimap2", "samtools", 250
PANMAN = {"rsv": "../src/test/data/rsv_4K.panman",
          "sars": "../examples/data/sars_20000_twilight_dipper.panman"}
IDX = os.environ.get("DUALK_IDX", "index_B")     # second parameterization's index name
RW = "/tmp/claude-1015/-scratch1-alex-poopdoop/e32e7b24-91a2-49f4-a2e6-8a1f8e802306/scratchpad/" + IDX
os.makedirs(RW, exist_ok=True)
_geno = {}


def node_fa(gen, node):
    if node not in _geno:
        o = os.path.join(RW, "g_" + str(abs(hash(node))) + ".fa")
        C.faidx_seq(gen, node, o, SAM)
        _geno[node] = o
    return _geno[node]


def load_scores(path):
    """return list of (node, containment, logContainment) in file order."""
    rows = []
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        ci, li = hdr.index("containment"), hdr.index("logContainment")
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) > li:
                rows.append((c[0], float(c[ci]), float(c[li])))
    return rows


def run_B(sp, stem, out):
    C.sh(["../build/bin/panmap", PANMAN[sp], stem + "_1.fq", stem + "_2.fq",
          "-i", f"work/{sp}/{IDX}.pmi", "--stop", "place", "-o", out,
          "--force-leaf", "--dump-all-scores", out + ".scores.tsv"])
    return out + ".scores.tsv"


def argmax_first(rows, keyidx, exclude):
    best, bv = None, None
    for r in rows:
        if r[0] == exclude:
            continue
        if bv is None or r[keyidx] > bv:
            best, bv = r[0], r[keyidx]
    return best, bv


def top_tie(rows, keyidx, exclude):
    _, bv = argmax_first(rows, keyidx, exclude)
    return [r[0] for r in rows if r[0] != exclude and r[keyidx] == bv], bv


def main():
    sp, kind = sys.argv[1], sys.argv[2]
    gen = f"work/{sp}/genomes.fa"
    is_sim = kind == "sim"
    for stem in sys.argv[3:]:
        rest = os.path.basename(stem)
        mut, cov = (rest.split("_")[1], rest.split("_")[2]) if is_sim else ("-1", rest.split("_")[1])
        sample_fa = stem + (".desc.fa" if is_sim else ".truth.fa")
        exp_fa = stem + ".exp.fa"
        sf = stem + (".leaf.fa" if is_sim else ".truth.fa")
        A_scores = stem + ".scores.tsv"
        if not all(os.path.exists(p) for p in (sample_fa, exp_fa, sf, A_scores)):
            continue
        self_node = open(sf).readline()[1:].strip()
        excl = None if is_sim else self_node
        try:
            A = load_scores(A_scores)
            B = load_scores(run_B(sp, stem, os.path.join(RW, f"{sp}_{kind}_{rest}")))
        except Exception as e:
            sys.stderr.write(f"FAIL {stem}: {e}\n")
            continue
        Bmap = {n: (c, l) for n, c, l in B}          # node -> (contain, logContain)
        ed = C.genome_distance(MM, exp_fa, sample_fa, EX)

        # rank maps by logContainment (0 = best); nodes absent in one index -> large rank
        def rankmap(rows):
            order = sorted([r for r in rows if r[0] != excl], key=lambda r: -r[2])
            return {r[0]: i for i, r in enumerate(order)}
        rA, rB = rankmap(A), rankmap(B)
        BIG = len(A) + len(B)

        picks = {}
        picks["A_only"], _ = argmax_first(A, 2, excl)
        picks["B_only"], _ = argmax_first(B, 2, excl)
        tieA, _ = top_tie(A, 2, excl)
        picks["tieBreakB"] = max(tieA, key=lambda n: Bmap.get(n, (0, 0))[1]) if tieA else None
        tieB, _ = top_tie(B, 2, excl)
        Amap = {n: l for n, c, l in A}
        picks["tieBreakA"] = max(tieB, key=lambda n: Amap.get(n, 0)) if tieB else None
        allnodes = set(rA) | set(rB)
        picks["rankSum"] = min(allnodes, key=lambda n: rA.get(n, BIG) + rB.get(n, BIG)) if allnodes else None
        cmax = max(((c + Bmap.get(n, (0, 0))[0]), n) for n, c, l in A if n != excl) if A else (0, None)
        picks["sumContain"] = cmax[1]
        inter = [n for n in tieA if n in set(tieB)]
        picks["intersect"] = inter[0] if inter else picks["A_only"]

        for method, node in picks.items():
            if not node:
                continue
            is_self = 1 if (not is_sim and node == self_node) else 0
            sc = 0 if (node == self_node and is_sim) else max(
                0, C.genome_distance(MM, node_fa(gen, node), sample_fa, EX) - ed)
            print(f"{method}\t{sp}\t{kind}\t{mut}\t{cov}\t{sc}\t{is_self}", flush=True)


if __name__ == "__main__":
    main()
