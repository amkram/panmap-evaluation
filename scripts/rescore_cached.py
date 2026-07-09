#!/usr/bin/env python3
"""Re-score the cached refine placement.tsv files per-sample (no panmap re-run),
emitting a sample id so methods can be compared on an identical subset.
Output rows: sample \t method \t sp \t kind \t mut \t cov \t score \t is_self
"""
import os
import sys
import glob

sys.path.insert(0, "scripts")
import common as C

MM, SAM, EX = "minimap2", "samtools", 250
RW = "/tmp/claude-1015/-scratch1-alex-poopdoop/e32e7b24-91a2-49f4-a2e6-8a1f8e802306/scratchpad/refine_work"
_geno = {}


def node_fa(gen, node):
    if node not in _geno:
        o = os.path.join(RW, "g_" + str(abs(hash(node))) + ".fa")
        C.faidx_seq(gen, node, o, SAM)
        _geno[node] = o
    return _geno[node]


def parse(pf):
    seed, refined = {}, {}
    for line in open(pf):
        c = line.rstrip("\n").split("\t")
        if len(c) < 3:
            continue
        name, nodes = c[0], c[2].split(",")
        (refined if name.startswith("refined_") else seed)[
            name[len("refined_"):] if name.startswith("refined_") else name] = (
            nodes[0] if name.startswith("refined_") else nodes)
    return seed, refined


def first_non_self(nodes, self_node):
    for n in nodes:
        if n != self_node:
            return n
    return nodes[0] if nodes else None


for pf in sorted(glob.glob(RW + "/*.placement.tsv")):
    b = os.path.basename(pf)[:-len(".placement.tsv")]
    parts = b.split("_")
    sp, kind, rest = parts[0], parts[1], "_".join(parts[2:])
    stem = f"work/{sp}/fig2_{kind}/{rest}"
    gen = f"work/{sp}/genomes.fa"
    is_sim = kind == "sim"
    if is_sim:
        bb = rest.split("_"); mut, cov = bb[1], bb[2]
    else:
        mut, cov = "-1", rest.split("_")[1]
    sample_fa = stem + (".desc.fa" if is_sim else ".truth.fa")
    exp_fa = stem + ".exp.fa"
    sf = stem + (".leaf.fa" if is_sim else ".truth.fa")
    if not (os.path.exists(sample_fa) and os.path.exists(exp_fa) and os.path.exists(sf)):
        continue
    self_node = open(sf).readline()[1:].strip()
    seed, refined = parse(pf)
    ed = C.genome_distance(MM, exp_fa, sample_fa, EX)
    excl = self_node if not is_sim else None
    picks = {
        "seed_logC":   first_non_self(seed.get("log_containment", []), excl),
        "seed_logCos": first_non_self(seed.get("log_cosine", []), excl),
        "ref_logC":    refined.get("log_containment"),
        "ref_logCos":  refined.get("log_cosine"),
        "ref_contain": refined.get("containment"),
        "ref_logRaw":  refined.get("log_raw"),
    }
    sample_id = f"{sp}/{kind}/{rest}"
    for method, node in picks.items():
        if not node:
            continue
        is_self = 1 if (not is_sim and node == self_node) else 0
        sc = 0 if (node == self_node and is_sim) else max(
            0, C.genome_distance(MM, node_fa(gen, node), sample_fa, EX) - ed)
        print(f"{sample_id}\t{method}\t{sp}\t{kind}\t{mut}\t{cov}\t{sc}\t{is_self}")
