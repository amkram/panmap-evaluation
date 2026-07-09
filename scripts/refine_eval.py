#!/usr/bin/env python3
"""Evaluate refine-mode (alignment-based) and other tie-break selections vs the
seed baseline. Runs `panmap --refine --stop place` per sample, reads placement.tsv,
and scores each selection's picked node with the harness genome_distance.

Selections compared (all leaf-restricted):
  seed_logC     : seed log_containment, first-non-self tied node   (== current Fig2)
  seed_logCos   : seed log_cosine best
  ref_logC      : refined_log_containment (alignment; logC seed tie-break)
  ref_logCos    : refined_log_cosine      (alignment; cosine seed tie-break)
  ref_contain   : refined_containment
  ref_logRaw    : refined_log_raw

Usage: refine_eval.py <sp> <kind> <stem> [stem ...]   -> TSV rows on stdout:
  method  sp  kind  mut  cov  score  is_self
"""
import os
import sys

sys.path.insert(0, "scripts")
import common as C

PANMAP = "../build/bin/panmap"
MM, SAM, EX = "minimap2", "samtools", 250
PANMAN = {"rsv": "../src/test/data/rsv_4K.panman",
          "sars": "../examples/data/sars_20000_twilight_dipper.panman",
          "tb": "/scratch1/alex/panstop/data/tb/tb_400.panman"}
TMP = "/tmp/claude-1015/-scratch1-alex-poopdoop/e32e7b24-91a2-49f4-a2e6-8a1f8e802306/scratchpad/refine_work"
os.makedirs(TMP, exist_ok=True)
_geno = {}


def node_fa(gen, node):
    if node not in _geno:
        o = os.path.join(TMP, str(abs(hash(node))) + ".fa")
        C.faidx_seq(gen, node, o, SAM)
        _geno[node] = o
    return _geno[node]


def run_panmap(sp, stem, out):
    r1, r2 = stem + "_1.fq", stem + "_2.fq"
    cmd = [PANMAP, PANMAN[sp], r1, r2, "-i", f"work/{sp}/index.pmi",
           "--stop", "place", "-o", out, "--force-leaf", "--refine",
           "--dump-all-scores", out + ".scores.tsv"]
    C.sh(cmd)
    return out + ".placement.tsv"


def parse(pf):
    seed, refined = {}, {}
    for line in open(pf):
        c = line.rstrip("\n").split("\t")
        if len(c) < 3:
            continue
        name, nodes = c[0], c[2].split(",")
        if name.startswith("refined_"):
            refined[name[len("refined_"):]] = nodes[0]
        else:
            seed[name] = nodes
    return seed, refined


def first_non_self(nodes, self_node):
    for n in nodes:
        if n != self_node:
            return n
    return nodes[0] if nodes else None


def main():
    sp, kind = sys.argv[1], sys.argv[2]
    gen = f"work/{sp}/genomes.fa"
    for stem in sys.argv[3:]:
        base = os.path.basename(stem).split("_")
        if kind == "sim":
            mut, cov = base[1], base[2]
        else:
            mut, cov = "-1", base[1]
        sample_fa = stem + (".desc.fa" if kind == "sim" else ".truth.fa")
        exp_fa = stem + ".exp.fa"
        self_node = None
        sf = stem + (".leaf.fa" if kind == "sim" else ".truth.fa")
        if os.path.exists(sf):
            self_node = open(sf).readline()[1:].strip()
        if not (os.path.exists(sample_fa) and os.path.exists(exp_fa)):
            continue
        out = os.path.join(TMP, f"{sp}_{kind}_" + os.path.basename(stem))
        try:
            pf = run_panmap(sp, stem, out)
            seed, refined = parse(pf)
        except Exception as e:
            sys.stderr.write(f"FAIL {stem}: {e}\n")
            continue
        ed = C.genome_distance(MM, exp_fa, sample_fa, EX)
        excl = self_node if kind == "real" else None
        picks = {
            "seed_logC":   first_non_self(seed.get("log_containment", []), excl),
            "seed_logCos": first_non_self(seed.get("log_cosine", []), excl),
            "ref_logC":    refined.get("log_containment"),
            "ref_logCos":  refined.get("log_cosine"),
            "ref_contain": refined.get("containment"),
            "ref_logRaw":  refined.get("log_raw"),
        }
        for method, node in picks.items():
            if not node:
                continue
            is_self = 1 if (kind == "real" and node == self_node) else 0
            if node == self_node and kind == "sim":
                sc = 0
            else:
                sc = max(0, C.genome_distance(MM, node_fa(gen, node), sample_fa, EX) - ed)
            print(f"{method}\t{sp}\t{kind}\t{mut}\t{cov}\t{sc}\t{is_self}", flush=True)


if __name__ == "__main__":
    main()
