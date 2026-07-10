#!/usr/bin/env python3
"""Figure 3 (revised): Panmap vs the field-standard pipeline (HaphPIPE / NCBI /
Clockwork), 4 rows x 3 species. The single blended "accuracy" row is replaced by
the two axes QUAST/dnadiff always report separately:

  Row 1  Genome fraction (%)       -- completeness; higher is better
  Row 2  Base errors / 100 kbp     -- per-base correctness incl. indel bases
                                      (subs + del_bases + ins_bases); lower is better
  Row 3  Runtime (s)
  Row 4  Peak memory (MB)

Rows 1-2 are computed from the cached assemblies (figure3_rescored.tsv, base-based
metric); rows 3-4 reuse the published runtime/memory (figure3.tsv). RSV's HaphPIPE
arm is split into RSV-A / RSV-B (subtype-matched reference), as in the original
Fig 3; the subtype per sample is the `subtype` column of figure3_rescored.tsv --
the competitive-mapping classification (classify_rsv_subtype.py) joined by node id,
an exact assignment rather than an accuracy-match heuristic.
Usage: plot_fig3_revised.py [rescored.tsv] [figure3.tsv] [out.pdf]
"""
import csv
import math
import os
import sys
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(__file__)
RES = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "results", "figure3_rescored.tsv")
PUB = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "..", "results", "figure3.tsv")
# Row-2 correctness metric (dimensionless, genome-size-independent): identity | qv | error
ROW2 = os.environ.get("ROW2METRIC", "identity")
OUT = sys.argv[3] if len(sys.argv) > 3 else os.path.join(HERE, "..", "results", f"figure3_revised_{ROW2}.pdf")

SP = ["rsv", "sars", "tb"]
LAB = {"rsv": "RSV (4K)", "sars": "SARS-CoV-2 (20K)", "tb": "M. tb (400)"}
SIZE = {"rsv": "4K", "sars": "20K", "tb": "400"}
# Field-standard pipeline label, with the single reference each uses (as in Fig 3).
STDLAB = {"rsv": "HaphPIPE", "sars": "NCBI Pipeline (NC_045512.2 ref)",
          "tb": "Clockwork (NC_000962.3 ref)"}
COVS = [0.5, 1, 10, 100]
PANMAP_C, STD_C = "#6A3D9A", "#E08214"
STD_A_C, STD_B_C = "#E08214", "#FDBF6F"          # RSV-A (dark), RSV-B (light)
np.random.seed(0)


def panmap_lab(sp):
    return f"Panmap ({SIZE[sp]}-sample pangenome)"


def arms_for(sp):
    if sp == "rsv":
        return [("panmap", PANMAP_C, panmap_lab(sp)),
                ("standard_B", STD_B_C, "HaphPIPE (RSV-B samples, LR699737.1 ref)"),
                ("standard_A", STD_A_C, "HaphPIPE (RSV-A samples, NC_038235.1 ref)")]
    return [("panmap", PANMAP_C, panmap_lab(sp)), ("standard", STD_C, STDLAB[sp])]


# ── rows 1-2: per-sample base counts from rescored table ──────────────────────
# RSV's HaphPIPE (standard) arm is split into RSV-A / RSV-B by the per-sample `subtype`
# column of figure3_rescored.tsv (competitive-mapping classification joined by node id),
# an exact assignment rather than matching accuracy values across tables.
data = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))    # [sp][arm][cov] -> [count dicts]
for r in csv.DictReader(open(RES), delimiter="\t"):
    m = r["method"]
    base_arm = "panmap" if m == "panmap" else ("standard" if m.startswith("standard") else None)
    if base_arm is None or r["interior"] in ("", "nan") or r["aligned"] in ("", "nan"):
        continue
    c = {k: float(r[k]) for k in ("snps", "indel_events", "del_bases", "ins_bases", "aligned", "interior")}
    if c["aligned"] <= 0 or c["interior"] <= 0:
        continue
    sp, cov, arm = r["species"], float(r["coverage"]), base_arm
    if sp == "rsv" and base_arm == "standard":
        sub = r.get("subtype", "").strip()
        arm = "standard_" + (sub if sub in ("A", "B") else "A")   # A fallback if unclassified
    data[sp][arm][cov].append(c)

# Only plot species actually present in the data (e.g. test mode = RSV+SARS, no TB),
# keeping the canonical order; fall back to all three if the table is unexpectedly empty.
SP = [sp for sp in SP if sp in data] or SP

# ── rows 3-4: runtime / peak memory from published table (standard* -> std) ────
rt = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
mem = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
for r in csv.DictReader(open(PUB), delimiter="\t"):
    m = r["method"]
    a = "panmap" if m == "panmap" else ("std" if m.startswith("standard") else None)
    if a is None:
        continue
    sp, cov = r["species"], float(r["coverage"])
    if r["wall_s"] not in ("", "nan"):
        rt[sp][a][cov].append(float(r["wall_s"]))
    if r.get("peak_mb") not in (None, "", "nan"):
        mem[sp][a][cov].append(float(r["peak_mb"]))


# ── per-sample metrics ────────────────────────────────────────────────────────
def genome_fraction(c): return 100.0 * c["aligned"] / c["interior"]


def _err_rate(c):       # per-base error fraction over alignment columns (subs + indel bases)
    return (c["snps"] + c["del_bases"] + c["ins_bases"]) / (c["aligned"] + c["ins_bases"])


def pct_identity(c):    return 100.0 * (1.0 - _err_rate(c))              # higher is better
def pct_error(c):       return 100.0 * _err_rate(c)                      # lower is better
def qv(c):              # Phred consensus quality; 0 observed errors -> 1-error floor
    r = _err_rate(c) or 1.0 / (c["aligned"] + c["ins_bases"])
    return -10.0 * math.log10(r)


# (fn, ylabel, lower_better) keyed by ROW2 choice
ROW2_METRIC = {
    "rate":     (_err_rate,    "Base error rate\n(erroneous / total bases)", True),
    "identity": (pct_identity, "Consensus identity (%)\n(subs + indel bases)", False),
    "qv":       (qv,           "Consensus QV (Phred)\n(subs + indel bases)", False),
    "error":    (pct_error,    "Base error (%)\n(subs + indel bases)", True),
}
row2_fn, row2_ylabel, row2_lower_better = ROW2_METRIC[ROW2]


def violin_panel(ax, sp, fn, legend=False):
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
                b.set_facecolor(col); b.set_alpha(0.22); b.set_edgecolor(col); b.set_linewidth(0.6)
            vp["cmedians"].set_color("black"); vp["cmedians"].set_linewidth(1.2)
        for i in range(len(COVS)):                    # overlay the actual points (jittered)
            v = vals[i]
            if not v:
                continue
            jit = np.random.uniform(-width * 0.35, width * 0.35, len(v))
            ax.scatter(xp[i] + jit, v, s=6, color=col, alpha=0.55, edgecolors="none", zorder=3)
    ax.set_xticks(pos); ax.set_xticklabels([f"{c}×" for c in COVS])
    ax.set_xlim(-0.5, len(COVS) - 0.5)
    ax.ticklabel_format(useOffset=False, axis="y")
    ax.grid(axis="y", alpha=0.25)
    if legend:
        ax.legend([plt.Rectangle((0, 0), 1, 1, fc=c, alpha=.55) for _, c, _ in arms],
                  [l for _, _, l in arms], fontsize=7, loc="lower left")


def line_panel(ax, sp, D, legend=False):
    for a, col, lab, ls in [("panmap", PANMAP_C, panmap_lab(sp), "-"), ("std", STD_C, STDLAB[sp], "--")]:
        xs = [c for c in COVS if D[sp][a].get(c)]
        if not xs:
            continue
        ys = [np.median(D[sp][a][c]) for c in xs]
        er = [np.std(D[sp][a][c]) for c in xs]
        ax.errorbar(xs, ys, yerr=er, marker="o", ms=4, color=col, label=lab, ls=ls, lw=1.8, capsize=3)
    ax.set_xscale("log"); ax.set_xticks(COVS); ax.set_xticklabels([f"{c}×" for c in COVS])
    ax.grid(alpha=0.25)
    if legend:
        ax.legend(fontsize=7)


fig, axes = plt.subplots(4, len(SP), figsize=(4.2 * len(SP), 11.6),
                         squeeze=False, gridspec_kw={"height_ratios": [1, 1, 0.48, 0.48]})
for ci, sp in enumerate(SP):
    # Row 1: genome fraction (completeness)
    ax = axes[0][ci]
    violin_panel(ax, sp, genome_fraction, legend=True)
    ax.set_title(LAB[sp], fontsize=11)
    gf_ymin = {"sars": 99.9, "tb": 99.5}              # RSV autoscales to show the HaphPIPE drop
    if sp in gf_ymin:
        # small headroom above 100 proportional to the panel's range (matches tb)
        ax.set_ylim(gf_ymin[sp], 100.0 + 0.04 * (100.0 - gf_ymin[sp]))
        ax.set_yticks([t for t in ax.get_yticks() if t <= 100.0 + 1e-9])
    if ci == 0:
        ax.set_ylabel("Genome fraction (%)")
    # Row 2: per-base correctness (identity / QV / error %)
    ax = axes[1][ci]
    violin_panel(ax, sp, row2_fn)
    if row2_lower_better:
        ax.set_ylim(bottom=0)
    if ci == 0:
        ax.set_ylabel(row2_ylabel)
    # Row 3: runtime
    ax = axes[2][ci]
    line_panel(ax, sp, rt, legend=True)
    if ci == 0:
        ax.set_ylabel("Runtime (s)")
    # Row 4: peak memory
    ax = axes[3][ci]
    line_panel(ax, sp, mem, legend=(ci == 0))
    if ci == 0:
        ax.set_ylabel("Peak memory (MB)")
    ax.set_xlabel("Read depth (coverage)")

fig.tight_layout()
for ext in ("pdf", "png"):
    p = OUT if OUT.endswith(ext) else OUT.rsplit(".", 1)[0] + "." + ext
    fig.savefig(p, dpi=200 if ext == "png" else None, bbox_inches="tight")
    print("wrote", p)
