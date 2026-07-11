#!/usr/bin/env python3
"""Supplementary figure: consensus quality, leaf-only vs internal-allowed placement.

Reads forceleaf_consensus.tsv (panmap consensus for every real leave-one-out sample x
coverage, under both placement modes) and draws, per species, two line plots vs
coverage with error bars, comparing the modes:
  - genome fraction  = aligned / interior                   (completeness; higher better)
  - per-base error   = (snps+del_bases+ins_bases) / aligned (accuracy; lower better)

Both modes share reads per sample, so each sample is paired. Line points are the mean
over samples with SEM error bars; each accuracy panel is annotated with a paired
two-sided sign test (the mean lines understate the effect because a few large
regressions offset many small improvements). Per-species panels (default) show the
claim holds on each dataset; --pool overlays all species into one pair of panels.

Usage: plot_forceleaf.py forceleaf_consensus.tsv out.pdf '<meta-json>' [pool]
"""
import csv
import sys
import json
from math import comb
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TSV, OUT = sys.argv[1], sys.argv[2]
META = json.loads(sys.argv[3]) if len(sys.argv) > 3 and sys.argv[3] else {}
POOL = len(sys.argv) > 4 and sys.argv[4] == "pool"
ORDER = META.get("order", ["rsv", "sars", "tb"])
LABELS = META.get("labels", {"rsv": "RSV", "sars": "SARS-CoV-2", "tb": "M. tuberculosis"})
COVS = [float(c) for c in META.get("cov", [0.5, 1, 10, 100])]
LEAF_C, INT_C = "#C0392B", "#1B9E77"   # leaf-only (red), internal-allowed (green)

gf = lambda r: 100.0 * float(r["aligned"]) / float(r["interior"])
pbe = lambda r: 100.0 * (float(r["snps"]) + float(r["del_bases"]) + float(r["ins_bases"])) / float(r["aligned"])

# unpaired series for the lines: [metric][sp][mode][cov] -> [values]
series = {"gf": defaultdict(lambda: defaultdict(lambda: defaultdict(list))),
          "pbe": defaultdict(lambda: defaultdict(lambda: defaultdict(list)))}
# paired by (sp,ri,cov) for the sign test
paired = defaultdict(dict)
for r in csv.DictReader(open(TSV), delimiter="\t"):
    if r["status"] != "ok" or r["aligned"] in ("", "nan") or float(r["aligned"]) <= 0:
        continue
    cov = float(r["coverage"])
    series["gf"][r["species"]][r["mode"]][cov].append(gf(r))
    series["pbe"][r["species"]][r["mode"]][cov].append(pbe(r))
    paired[(r["species"], r["ri"], r["coverage"])][r["mode"]] = r

species = [sp for sp in ORDER if sp in series["gf"]]


def mean_sem(v):
    a = np.array(v, float)
    return a.mean(), (a.std(ddof=1) / np.sqrt(len(a)) if len(a) > 1 else 0.0)


def sign_p(b, w):
    n = b + w
    return 1.0 if n == 0 else min(1.0, 2.0 * sum(comb(n, i) for i in range(min(b, w) + 1)) / 2 ** n)


def paired_test(sp_list, metric, cmp):
    """Return (better, worse, p) over all paired samples of the given species."""
    fn = gf if metric == "gf" else pbe
    b = w = 0
    for (s, ri, c), d in paired.items():
        if (sp_list and s not in sp_list) or "leaf" not in d or "internal" not in d:
            continue
        lv, iv = fn(d["leaf"]), fn(d["internal"])
        if abs(iv - lv) < 1e-9:
            continue
        if (cmp == ">" and iv > lv) or (cmp == "<" and iv < lv):
            b += 1
        else:
            w += 1
    return b, w, sign_p(b, w)


def lines(ax, metric, sp_list):
    means = []
    for mode, col, lab in (("leaf", LEAF_C, "Leaf-only (--force-leaf)"),
                           ("internal", INT_C, "Internal nodes allowed")):
        xs, ys, es = [], [], []
        for cov in COVS:
            pooled = [v for sp in sp_list for v in series[metric][sp][mode].get(cov, [])]
            if pooled:
                m, se = mean_sem(pooled)
                xs.append(cov); ys.append(m); es.append(se); means.append(m)
        ax.errorbar(xs, ys, yerr=es, marker="o", ms=5, color=col, label=lab, lw=1.9, capsize=3, zorder=3)
    ax.set_xscale("log")
    ax.set_xticks(COVS); ax.set_xticklabels([f"{c:g}×" for c in COVS])
    ax.set_xlabel("Read depth (coverage)")
    ax.grid(alpha=0.25); ax.margins(x=0.08)
    if metric == "gf" and means:                       # keep near-100% ceiling from auto-zooming
        ax.set_ylim(min(99.0, min(means) - 0.2), 100.15)
        ax.ticklabel_format(axis="y", useOffset=False, style="plain")
    return means


def annotate_pbe(ax, sp_list):
    b, w, p = paired_test(sp_list, "pbe", "<")
    ptxt = "p<1e-3" if p < 1e-3 else f"p={p:.2g}"
    verdict = "internal lower error" if b > w else ("tie" if b == w else "leaf lower error")
    ax.text(0.97, 0.95, f"{verdict}\n{b}/{b + w} paired, {ptxt}", transform=ax.transAxes,
            ha="right", va="top", fontsize=7.5,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7", alpha=0.85))


GF_YL = "Genome fraction (%)\n(completeness, higher better)"
PBE_YL = "Per-base error (%)\n(lower better)"

if POOL:
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.4))
    lines(ax[0], "gf", species); ax[0].set_ylabel(GF_YL); ax[0].set_title("Consensus completeness")
    lines(ax[1], "pbe", species); ax[1].set_ylabel(PBE_YL); ax[1].set_title("Consensus accuracy")
    annotate_pbe(ax[1], species)
    ax[0].legend(fontsize=8, loc="lower right")
else:
    fig, axes = plt.subplots(len(species), 2, figsize=(9.2, 3.4 * len(species)), squeeze=False)
    for i, sp in enumerate(species):
        n = max((len(series["gf"][sp][m].get(c, [])) for m in ("leaf", "internal") for c in COVS), default=0)
        lines(axes[i][0], "gf", [sp]); lines(axes[i][1], "pbe", [sp])
        axes[i][0].set_ylabel(GF_YL); axes[i][1].set_ylabel(PBE_YL)
        annotate_pbe(axes[i][1], [sp])
        axes[i][0].set_title(f"{LABELS.get(sp, sp)} — completeness  (n≤{n})", fontsize=10)
        axes[i][1].set_title(f"{LABELS.get(sp, sp)} — accuracy", fontsize=10)
    axes[0][0].legend(fontsize=8, loc="lower right")

fig.suptitle("Consensus quality: leaf-only vs internal-allowed placement", fontsize=12, y=1.002)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(OUT.rsplit(".", 1)[0] + "." + ext, dpi=200 if ext == "png" else None, bbox_inches="tight")

# console summary
print(f"{'species':10}{'metric':>5}{'better':>7}{'worse':>6}{'sign_p':>9}  (paired internal vs leaf)")
for sp in species + (["ALL"] if len(species) > 1 else []):
    spl = species if sp == "ALL" else [sp]
    for metric, cmp in (("gf", ">"), ("pbe", "<")):
        b, w, p = paired_test(spl, metric, cmp)
        print(f"{sp:10}{metric:>5}{b:>7}{w:>6}{p:>9.2g}")
print("wrote", OUT)
