#!/usr/bin/env python3
"""Composite of Fig 2 Row-1 (placement-accuracy panels), one ROW per ranking
metric x one COLUMN per species. Reads results/figure2_<metric>.<sp>.part
(produced by metric_variants.py). y-axis is symlog so the near-zero placement
scores and the high random-placement baseline are both legible; y is shared
across the metric rows within each species column for fair comparison.

Usage: plot_metric_grid.py out.pdf
"""
import csv
import sys
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

out = sys.argv[1] if len(sys.argv) > 1 else "results/figure2_metric_grid.pdf"
import os
METRICS = ["logRaw", "logCosine", "containment", "logContainment"]
# species columns: only those whose .part files exist (lets us render interim)
ORDER = [sp for sp in ["rsv", "sars", "tb"]
         if os.path.exists(f"results/figure2_logContainment.{sp}.part")]
LABELS = {"rsv": "RSV (4K)", "sars": "SARS-CoV-2 (20K)", "tb": "M. tb (400)"}
MUTS = {"rsv": [1e-4, 1e-3], "sars": [1e-4, 1e-3], "tb": [1e-5, 1e-4]}
COVS = [0.5, 1.0, 10.0, 100.0]
SIM_C = ["#9ECAE1", "#08519C"]
REAL_C = "#E6550D"
RAND_C = "#BBBBBB"

# data[metric][sp] -> (scores{(kind,mut)}{cov}:[...], rand[...])
data = defaultdict(dict)
for metric in METRICS:
    for sp in ORDER:
        scores = defaultdict(lambda: defaultdict(list))
        rand = []
        path = f"results/figure2_{metric}.{sp}.part"
        try:
            rows = [l.rstrip("\n").split("\t") for l in open(path)]
        except FileNotFoundError:
            data[metric][sp] = (scores, rand)
            continue
        for r in rows:
            if len(r) < 6:
                continue
            _, kind, mut, cov, score, rnd = r[0], r[1], int(r[2]), float(r[3]), r[4], r[5]
            if score not in ("", None):
                try:
                    scores[(kind, mut)][cov].append(float(score))
                except ValueError:
                    pass
            if rnd:
                rand += [float(x) for x in rnd.split(";")]
        data[metric][sp] = (scores, rand)

# shared y-limit per species column (max over all metrics of any plotted value)
ymax = {}
for sp in ORDER:
    mx = 1.0
    for metric in METRICS:
        sc, rnd = data[metric][sp]
        for cd in sc.values():
            for vals in cd.values():
                if vals:
                    mx = max(mx, max(vals))
        if rnd:
            mx = max(mx, max(rnd))
    ymax[sp] = mx


def violin(ax, x, vals, color, w=0.7):
    vals = [v for v in vals if v is not None]
    if not vals:
        return
    v = ax.violinplot([vals], positions=[x], widths=w, showextrema=False)
    for b in v["bodies"]:
        b.set_facecolor(color); b.set_alpha(0.35); b.set_edgecolor("none")
    jit = (np.arange(len(vals)) % 7 - 3) / 3.0 * 0.12      # deterministic jitter
    ax.scatter(x + jit, vals, s=5, color=color, alpha=0.55, edgecolors="none", zorder=3)
    ax.plot([x - 0.25, x + 0.25], [np.median(vals)] * 2, color="black", lw=1.4, zorder=4)


nrow, ncol = len(METRICS), len(ORDER)
fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 2.6 * nrow), squeeze=False)

for ri, metric in enumerate(METRICS):
    for ci, sp in enumerate(ORDER):
        ax = axes[ri][ci]
        sc, rnd = data[metric][sp]
        groups = [("sim", 0), ("sim", 1)] + ([("real", -1)] if ("real", -1) in sc else [])
        for k, cov in enumerate(COVS):
            base = k * (len(groups) + 0.6)
            for gi, (kind, mut) in enumerate(groups):
                col = REAL_C if kind == "real" else SIM_C[mut]
                violin(ax, base + gi, sc[(kind, mut)].get(cov, []), col, w=0.8)
        centers = [k * (len(groups) + 0.6) + (len(groups) - 1) / 2 for k in range(len(COVS))]
        rbase = len(COVS) * (len(groups) + 0.6) + 0.5
        violin(ax, rbase, rnd, RAND_C, w=1.2)
        ax.set_xticks(centers + [rbase])
        ax.set_xticklabels([f"{c:g}×" for c in COVS] + ["Rand"], fontsize=7)
        ax.set_yscale("symlog", linthresh=1)
        ax.set_ylim(-0.5, ymax[sp] * 1.3)
        ax.set_yticks([0, 1, 10, 100, 1000, 10000])
        ax.grid(axis="y", alpha=0.25)
        ax.tick_params(labelsize=7)
        if ri == 0:
            ax.set_title(LABELS[sp], fontsize=11, weight="bold")
        if ci == 0:
            ax.set_ylabel(f"{metric}\n\nscore", fontsize=9)
            ax.get_yaxis().get_label().set_weight("bold")

# legend (top-left cell)
handles = [plt.Line2D([], [], marker="o", ls="", color=SIM_C[0], label="Sim μ (low)"),
           plt.Line2D([], [], marker="o", ls="", color=SIM_C[1], label="Sim μ (high)"),
           plt.Line2D([], [], marker="o", ls="", color=REAL_C, label="Real"),
           plt.Line2D([], [], marker="o", ls="", color=RAND_C, label="Random"),
           plt.Line2D([], [], color="black", lw=1.4, label="Median")]
axes[0][0].legend(handles=handles, fontsize=6.5, loc="upper left", framealpha=0.9)

fig.suptitle("Fig 2 placement accuracy — one row per placement-ranking metric (symlog y)",
             fontsize=12, weight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.98])
fig.savefig(out, bbox_inches="tight")
fig.savefig(out.rsplit(".", 1)[0] + ".png", dpi=170, bbox_inches="tight")
print("wrote", out)
