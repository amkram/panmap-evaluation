#!/usr/bin/env python3
"""Alternative Fig-3 accuracy panels under three quality conventions, all showing
Panmap vs the field-standard pipeline (HaphPIPE / NCBI / Clockwork) per species x
coverage, from the cached assemblies (figure3_rescored.tsv). Emits:

  figure3_quast.*    QUAST-style decomposition (3 rows, linear y): genome fraction
                     (%), mismatches per 100 kbp, indel events per 100 kbp.
  figure3_baseerr.*  base error rate = (subs + del_bases + ins_bases) per 100 kbp
                     over alignment columns; charges every wrong/missing/spurious
                     base, so large indels and insertions are fully counted (symlog).
  figure3_baseacc.*  original base-based accuracy = matches/interior (%); excludes
                     insertions, uses truth length as denominator.

RSV's HaphPIPE arm is split into RSV-A / RSV-B (subtype-matched reference), as in the
main Fig 3; the per-sample subtype is recovered from results/figure3.tsv by matching
each sample's event-based accuracy. Genome fraction and per-100kbp mismatch/indel rates
match what QUAST reports (dnadiff/MUMmer agree).
Usage: plot_fig3_quality.py [rescored.tsv]
"""
import csv
import os
import sys
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RES = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    os.path.dirname(__file__), "..", "results", "figure3_rescored.tsv")
OUTDIR = os.path.dirname(os.path.abspath(RES))
PUB = os.path.join(OUTDIR, "figure3.tsv")            # for RSV subtype labels
SP = ["rsv", "sars", "tb"]
LAB = {"rsv": "RSV (4K)", "sars": "SARS-CoV-2 (20K)", "tb": "M. tb (400)"}
STDLAB = {"rsv": "HaphPIPE", "sars": "NCBI Pipeline", "tb": "Clockwork"}
COVS = [0.5, 1, 10, 100]
PANMAP_C, STD_C = "#6A3D9A", "#E08214"
STD_A_C, STD_B_C = "#E08214", "#FDBF6F"              # RSV-A (dark), RSV-B (light) -- matches Fig 3
np.random.seed(0)


def arms_for(sp):
    if sp == "rsv":
        return [("panmap", PANMAP_C, "Panmap"),
                ("standard_B", STD_B_C, "HaphPIPE (RSV-B, LR699737.1)"),
                ("standard_A", STD_A_C, "HaphPIPE (RSV-A, NC_038235.1)")]
    return [("panmap", PANMAP_C, "Panmap"), ("standard", STD_C, STDLAB[sp])]


def armof(m):
    if m == "panmap":
        return "panmap"
    if m.startswith("standard"):
        return "standard"
    return None                                      # baseline arms not shown here


# RSV subtype pool: per coverage, list of [event_accuracy, subtype] from published table.
sub_pool = defaultdict(list)
for r in csv.DictReader(open(PUB), delimiter="\t"):
    if r["species"] == "rsv" and r["method"].startswith("standard_") and r["accuracy"] not in ("", "nan"):
        sub_pool[float(r["coverage"])].append([float(r["accuracy"]), r["method"].split("_", 1)[1]])

# data[sp][arm][cov] = list of per-sample count dicts
data = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
for r in csv.DictReader(open(RES), delimiter="\t"):
    base_arm = armof(r["method"])
    if base_arm is None or r["interior"] in ("", "nan") or r["aligned"] in ("", "nan"):
        continue
    c = {k: float(r[k]) for k in ("snps", "indel_events", "del_bases",
                                  "ins_bases", "aligned", "interior")}
    if c["aligned"] <= 0 or c["interior"] <= 0:
        continue
    sp, cov, arm = r["species"], float(r["coverage"]), base_arm
    if sp == "rsv" and base_arm == "standard":       # recover A/B subtype by event-acc match
        ae = float(r["accuracy_event"])
        cand = sub_pool[cov]
        j = min(range(len(cand)), key=lambda k: abs(cand[k][0] - ae), default=None)
        sub = cand.pop(j)[1] if (j is not None and abs(cand[j][0] - ae) < 1e-3) else "A"
        arm = "standard_" + sub
    data[sp][arm][cov].append(c)


# ── per-sample metrics (c = count dict) ───────────────────────────────────────
def genome_fraction(c): return 100.0 * c["aligned"] / c["interior"]
def mism_100k(c):       return c["snps"] * 1e5 / c["aligned"]
def indels_100k(c):     return c["indel_events"] * 1e5 / c["aligned"]
def base_err_100k(c):   return (c["snps"] + c["del_bases"] + c["ins_bases"]) * 1e5 / (c["aligned"] + c["ins_bases"])
def acc_base(c):        return 100.0 * (c["aligned"] - c["snps"] - c["del_bases"]) / c["interior"]


def panel(ax, sp, fn, logy=False, legend=False):
    arms = arms_for(sp)
    pos = np.arange(len(COVS))
    offs = [-0.34 + 0.68 * (k + 0.5) / len(arms) for k in range(len(arms))]
    width = 0.68 / len(arms) * 0.85
    for (a, col, _), off in zip(arms, offs):
        vals = [[fn(c) for c in data[sp][a].get(cv, [])] for cv in COVS]
        xp = pos + off
        idx = [i for i in range(len(COVS)) if len(vals[i]) > 1]
        if idx:
            vp = ax.violinplot([vals[i] for i in idx], positions=[xp[i] for i in idx],
                               widths=width, showmedians=True, showextrema=False)
            for b in vp["bodies"]:
                b.set_facecolor(col); b.set_alpha(0.22)
                b.set_edgecolor(col); b.set_linewidth(0.6)
            vp["cmedians"].set_color("black"); vp["cmedians"].set_linewidth(1.2)
        for i in range(len(COVS)):
            v = vals[i]
            if not v:
                continue
            jit = np.random.uniform(-width * 0.35, width * 0.35, len(v))
            ax.scatter(xp[i] + jit, v, s=6, color=col, alpha=0.55, edgecolors="none", zorder=3)
    ax.set_xticks(pos); ax.set_xticklabels([f"{c}×" for c in COVS])
    if logy:
        ax.set_yscale("symlog", linthresh=1)
        ax.set_ylim(bottom=0)                        # rates are >=0; hide symlog's negative half
    else:
        ax.ticklabel_format(useOffset=False, axis="y")
    ax.grid(axis="y", alpha=0.25)
    if legend:
        ax.legend([plt.Rectangle((0, 0), 1, 1, fc=c, alpha=.55) for _, c, _ in arms],
                  [l for _, _, l in arms], fontsize=7, loc="best")


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUTDIR, f"{name}.{ext}"),
                    dpi=200 if ext == "png" else None, bbox_inches="tight")
    print("wrote", os.path.join(OUTDIR, name + ".{pdf,png}"))


# ── 1) QUAST-style decomposition: 3 metric rows x species (linear y) ──────────
rows = [("Genome fraction (%)", genome_fraction),
        ("Mismatches / 100 kbp", mism_100k),
        ("Indels / 100 kbp (events)", indels_100k)]
fig, axes = plt.subplots(len(rows), len(SP), figsize=(4.0 * len(SP), 2.9 * len(rows)),
                         squeeze=False)
for ri, (ylab, fn) in enumerate(rows):
    for ci, sp in enumerate(SP):
        ax = axes[ri][ci]
        panel(ax, sp, fn, logy=False, legend=(ri == 0))
        if ri == 0:                                  # genome fraction: fix SARS/TB floor at 99.5
            ax.set_title(LAB[sp], fontsize=10)
            if sp in ("sars", "tb"):
                ax.set_ylim(99.5, 100.02)
        else:                                        # rates: non-negative
            ax.set_ylim(bottom=0)
        if ci == 0:
            ax.set_ylabel(ylab)
        if ri == len(rows) - 1:
            ax.set_xlabel("Read depth (coverage)")
fig.suptitle("QUAST-style assembly quality: Panmap vs field-standard pipeline",
             fontsize=12, y=1.0)
fig.tight_layout()
save(fig, "figure3_quast")

# ── 2) base error rate (charges every erroneous base) ─────────────────────────
fig, axes = plt.subplots(1, len(SP), figsize=(4.0 * len(SP), 3.4), squeeze=False)
for ci, sp in enumerate(SP):
    ax = axes[0][ci]
    panel(ax, sp, base_err_100k, logy=True, legend=True)
    ax.set_title(LAB[sp], fontsize=10)
    if ci == 0:
        ax.set_ylabel("Base errors / 100 kbp\n(subs + indel bases)")
    ax.set_xlabel("Read depth (coverage)")
fig.suptitle("Base-based error rate: every wrong/missing/spurious base counted",
             fontsize=12, y=1.03)
fig.tight_layout()
save(fig, "figure3_baseerr")

# ── 3) original base-based accuracy (matches/interior) ────────────────────────
fig, axes = plt.subplots(1, len(SP), figsize=(4.0 * len(SP), 3.4), squeeze=False)
YL = {"rsv": (40, 101), "sars": (99.0, 100.1), "tb": (98.0, 100.1)}
for ci, sp in enumerate(SP):
    ax = axes[0][ci]
    panel(ax, sp, acc_base, legend=True)
    ax.set_ylim(*YL[sp])
    ax.set_title(LAB[sp], fontsize=10)
    if ci == 0:
        ax.set_ylabel("Correctly reconstructed (%)\n= matches / interior")
    ax.set_xlabel("Read depth (coverage)")
fig.suptitle("Base-based accuracy (matches / interior); insertions not charged",
             fontsize=12, y=1.03)
fig.tight_layout()
save(fig, "figure3_baseacc")
