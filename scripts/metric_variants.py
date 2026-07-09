#!/usr/bin/env python3
"""Rebuild Fig 2 score tables selecting the placed node by EACH ranking metric
(logRaw, logCosine, containment, logContainment) from the already-computed
--dump-all-scores TSVs. No panmap re-runs: the picked node's genome is an O(1)
faidx from work/{sp}/genomes.fa, scored with the harness's exact genome_distance.

The logContainment variant must reproduce results/figure2.tsv (validation).

Usage: metric_variants.py <sp> [sp2 ...]   (writes results/figure2_<metric>.<sp>.part)
"""
import os
import re
import sys
import glob

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
import common as C

MM, SAM = "minimap2", "samtools"
EXCLUDE_BP = 250                       # config.yaml exclude_bp
METRICS = ["logRaw", "logCosine", "containment", "logContainment", "logContainment+cosTB"]
CACHE = "/tmp/claude-1015/-scratch1-alex-poopdoop/e32e7b24-91a2-49f4-a2e6-8a1f8e802306/scratchpad/mv_genomes"
os.makedirs(CACHE, exist_ok=True)

_dist = {}                             # (node, sample_fa) -> distance
_geno = {}                             # node -> extracted fa path


def header_node(fa):
    with open(fa) as f:
        return f.readline()[1:].strip()


def node_fa(genomes_fa, node):
    if node not in _geno:
        safe = str(abs(hash(node)))
        out = os.path.join(CACHE, safe + ".fa")
        C.faidx_seq(genomes_fa, node, out, SAM)
        _geno[node] = out
    return _geno[node]


def dist(genomes_fa, node, sample_fa):
    key = (node, sample_fa)
    if key not in _dist:
        _dist[key] = C.genome_distance(MM, node_fa(genomes_fa, node), sample_fa, EXCLUDE_BP)
    return _dist[key]


def pick(scores_file, metric, exclude):
    """Replicate place(): argmax(metric), first-in-file on ties, drop `exclude`.
    metric 'logContainment+cosTB': argmax(logContainment), break exact ties by
    logCosine (keeps logContainment's unique picks, resolves only the ties)."""
    rows = []
    with open(scores_file) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(hdr)}
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < len(hdr) or c[0] == exclude:
                continue
            rows.append(c)
    if not rows:
        return None
    if metric == "logContainment+cosTB":
        lc, cos = idx["logContainment"], idx["logCosine"]
        mx = max(float(r[lc]) for r in rows)
        tied = [r for r in rows if float(r[lc]) == mx]
        return max(tied, key=lambda r: float(r[cos]))[0]
    mi = idx.get(metric)
    if mi is None:
        return None
    best_node, best_val = None, None
    for c in rows:
        try:
            v = float(c[mi])
        except ValueError:
            continue
        if best_val is None or v > best_val:          # strictly > keeps first on ties
            best_node, best_val = c[0], v
    return best_node


def process(sp):
    genomes_fa = f"work/{sp}/genomes.fa"
    out = {m: [] for m in METRICS}
    stats = {m: {} for m in METRICS}                  # (kind,cov) -> [scores]
    # per-sample score files only: sim = i_m_cov.tsv, real = ri_cov.tsv
    # (exclude sibling .scores.tsv / .placement.tsv that also end in .tsv)
    pat = {"fig2_sim": re.compile(r"^\d+_\d+_[\d.]+\.tsv$"),
           "fig2_real": re.compile(r"^\d+_[\d.]+\.tsv$")}
    for kind, sub, is_sim in [("sim", "fig2_sim", True), ("real", "fig2_real", False)]:
        for tsv in sorted(glob.glob(f"work/{sp}/{sub}/*.tsv")):
            if not pat[sub].match(os.path.basename(tsv)):
                continue
            stem = tsv[:-4]
            # per-sample final tsv rows look like: sp\tkind\tm\tcov\tscore\trand\twall\trss
            row = open(tsv).read().rstrip("\n").split("\t")
            if len(row) < 8:
                continue
            _, _, m, cov, score_orig, rand, wall, rss = row[:8]
            scores_file = stem + ".scores.tsv"
            sample_fa = stem + (".desc.fa" if is_sim else ".truth.fa")
            exp_fa = stem + ".exp.fa"
            base = f"{sp}\t{kind}\t{m}\t{cov}"
            tail = f"{rand}\t{wall}\t{rss}"
            # no scores / no sample -> keep original (covers no-calls)
            if (not os.path.exists(scores_file) or not os.path.exists(sample_fa)
                    or not os.path.exists(exp_fa)):
                for mm in METRICS:
                    out[mm].append(f"{base}\t{score_orig}\t{tail}\n")
                continue
            truth_node = header_node(stem + (".leaf.fa" if is_sim else ".truth.fa"))
            exclude = None if is_sim else truth_node
            ed = C.genome_distance(MM, exp_fa, sample_fa, EXCLUDE_BP)
            for metric in METRICS:
                node = pick(scores_file, metric, exclude)
                if node is None:
                    sc = ""
                elif is_sim and node == truth_node:
                    sc = 0                              # identical genome -> pd==ed
                else:
                    pd = dist(genomes_fa, node, sample_fa)
                    sc = max(0, pd - ed)
                out[metric].append(f"{base}\t{sc}\t{tail}\n")
                if sc != "":
                    stats[metric].setdefault((kind, cov), []).append(sc)
    for metric in METRICS:
        with open(f"results/figure2_{metric}.{sp}.part", "w") as f:
            f.writelines(out[metric])
    # brief per-metric summary at the noisy low-cov columns
    for metric in METRICS:
        s = stats[metric]
        cells = []
        for cov in ["0.5", "1.0"]:
            vals = s.get(("sim", cov), [])
            hi = sum(1 for v in vals if v > 10)
            cells.append(f"{cov}x:>10={hi}/{len(vals)} mean={sum(vals)/max(1,len(vals)):.1f}")
        print(f"[{sp}] {metric:15s} sim  " + "  ".join(cells))


if __name__ == "__main__":
    for sp in sys.argv[1:]:
        process(sp)
