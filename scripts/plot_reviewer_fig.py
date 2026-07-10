#!/usr/bin/env python3
"""Reviewer-response figure: leave-one-out genotyping, RSV + SARS side by side
(2 columns), genotyping accuracy on top and assembly runtime underneath. Uses the
internal+completeness-weight (cw) Fig 3 data. Three genotyping arms:
  - Panmap            : panmap selects the reference AND genotypes (native)
  - Panmap->BWA+iVar  : panmap selects the reference, BWA+iVar genotypes
  - BWA+iVar          : single standard reference + BWA+iVar
Accuracy = % of the held-out genome correctly reconstructed (250 bp flanks ignored);
each arm is a median line + IQR band vs coverage.
Usage: plot_reviewer_fig.py <fig3_cw.tsv> <out_prefix> ['<label_json>']
"""
import csv
import json
import sys
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

tsv, out = sys.argv[1], sys.argv[2]
meta = json.loads(sys.argv[3]) if len(sys.argv) > 3 else {}
LABELS = meta.get("labels", {"rsv": "RSV (4K panman)", "sars": "SARS-CoV-2 (20K panman)"})
COVS = meta.get("cov", [0.5, 1, 10, 100])
ORDER = meta.get("order", ["rsv", "sars"])

PANMAP_C, PANREF_C, REF_C = "#6A3D9A", "#1B9E77", "#C0392B"


def mclass(m):
    return "panmap" if m == "panmap" else "panmap_ref" if m.startswith("panmap_") else "ref"


SERIES = [("panmap", PANMAP_C, "Panmap (native genotyping)", "-"),
          ("panmap_ref", PANREF_C, "Panmap→BWA+iVar", "-."),
          ("ref", REF_C, "BWA+iVar (single ref)", "--")]

acc = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))   # [sp][cls][cov]
rt = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
for r in csv.DictReader(open(tsv), delimiter="\t"):
    sp, cls, cov = r["species"], mclass(r["method"]), float(r["coverage"])
    if r["accuracy"] not in ("", "nan"):
        acc[sp][cls][cov].append(float(r["accuracy"]))
    if r["wall_s"]:
        rt[sp][cls][cov].append(float(r["wall_s"]))

order = [sp for sp in ORDER if sp in acc]
fig, axes = plt.subplots(2, len(order), figsize=(5.4 * len(order), 8), squeeze=False)


def q(v, p):
    v = sorted(v)
    if not v:
        return np.nan
    idx = (len(v) - 1) * p
    lo = int(idx)
    return v[lo] if lo == idx else v[lo] + (v[lo + 1] - v[lo]) * (idx - lo)


for j, sp in enumerate(order):
    # top: accuracy (median line + IQR band)
    ax = axes[0][j]
    for cls, col, lab, ls in SERIES:
        xs = [c for c in COVS if acc[sp][cls].get(c)]
        if not xs:
            continue
        med = [np.median(acc[sp][cls][c]) for c in xs]
        lo = [q(acc[sp][cls][c], .25) for c in xs]
        hi = [q(acc[sp][cls][c], .75) for c in xs]
        ax.fill_between(xs, lo, hi, color=col, alpha=0.15, linewidth=0)
        ax.plot(xs, med, marker="o", color=col, label=lab, ls=ls, lw=2)
    ax.set_xscale("log"); ax.set_xticks(COVS); ax.set_xticklabels([f"{c}×" for c in COVS])
    ax.set_ylim(-3, 103)
    ax.set_title(LABELS.get(sp, sp), fontsize=11, fontweight="bold")
    if j == 0:
        ax.set_ylabel("Correctly genotyped (%)")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8, loc="lower right")
    # bottom: runtime
    ax2 = axes[1][j]
    for cls, col, lab, ls in SERIES:
        xs = [c for c in COVS if rt[sp][cls].get(c)]
        if not xs:
            continue
        med = [np.median(rt[sp][cls][c]) for c in xs]
        er = [np.std(rt[sp][cls][c]) for c in xs]
        ax2.errorbar(xs, med, yerr=er, marker="o", color=col, label=lab, ls=ls, lw=2, capsize=3)
    ax2.set_xscale("log"); ax2.set_xticks(COVS); ax2.set_xticklabels([f"{c}×" for c in COVS])
    ax2.set_xlabel("Read depth (coverage)")
    if j == 0:
        ax2.set_ylabel("Assembly runtime (s)")
    ax2.set_title(f"Runtime, {LABELS.get(sp, sp)}", fontsize=10)
    ax2.grid(alpha=0.25); ax2.legend(fontsize=8)

fig.tight_layout()
fig.savefig(out + ".pdf", bbox_inches="tight")
fig.savefig(out + ".png", dpi=200, bbox_inches="tight")
print("wrote", out + ".pdf/.png")
