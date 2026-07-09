#!/usr/bin/env python3
"""Reviewer figure: reference-selection accuracy vs candidate-genome batch size,
for panmap (pangenome, may pick internal/ancestral nodes) vs mash vs minimap2
(leaf genomes only). One panel per species; x = number of candidate genomes to
choose from (log scale); y = Fig-3 genotyping accuracy of the selected reference
against the sample's truth genome (median across samples, IQR band).
Usage: plot_reviewer_refsel.py <refsel.tsv> <out.pdf> '<labels_json>'"""
import csv
import json
import sys
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

tsv, out = sys.argv[1], sys.argv[2]
LABELS = json.loads(sys.argv[3]) if len(sys.argv) > 3 else {}

METHODS = [("panmap", "#6A3D9A", "Panmap (pangenome)"),
           ("minimap2", "#1B9E77", "minimap2 (closest leaf)"),
           ("mash", "#E08214", "Mash (closest leaf)")]

# acc[species][method][batch] -> list of accuracies
acc = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
sp_order = []
for r in csv.DictReader(open(tsv), delimiter="\t"):
    sp, m, b, a = r["species"], r["method"], int(r["batch"]), r["accuracy"]
    if sp not in sp_order:
        sp_order.append(sp)
    if a not in ("", "nan"):
        acc[sp][m][b].append(float(a))

ncol = len(sp_order)
fig, axes = plt.subplots(1, ncol, figsize=(5.0 * ncol, 4.2), squeeze=False)

for j, sp in enumerate(sp_order):
    ax = axes[0][j]
    for m, col, lab in METHODS:
        batches = sorted(acc[sp][m])
        if not batches:
            continue
        med = [np.median(acc[sp][m][b]) for b in batches]
        lo = [np.percentile(acc[sp][m][b], 25) for b in batches]
        hi = [np.percentile(acc[sp][m][b], 75) for b in batches]
        ax.plot(batches, med, "-o", color=col, label=lab, lw=2, ms=5, zorder=3)
        ax.fill_between(batches, lo, hi, color=col, alpha=0.15, zorder=1)
    ax.set_xscale("log", base=2)
    n = max((len(acc[sp][METHODS[0][0]][b]) for b in acc[sp][METHODS[0][0]]), default=0)
    ax.set_title(f"{LABELS.get(sp, sp)}  (n={n})", fontsize=11)
    ax.set_xlabel("Candidate genomes to select from")
    if j == 0:
        ax.set_ylabel("Selected-reference accuracy (%)")
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=8, loc="lower right")

fig.suptitle("Reference selection vs. database size: pangenome (panmap) vs. leaf-based (mash, minimap2)",
             fontsize=12, y=1.02)
fig.tight_layout()
fig.savefig(out, bbox_inches="tight")
fig.savefig(out.rsplit(".", 1)[0] + ".png", dpi=200, bbox_inches="tight")
print("wrote", out)
