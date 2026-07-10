#!/usr/bin/env python3
"""Figure 2: placement accuracy (A, per species) + runtime (B) + peak memory (C).
Usage: plot_fig2.py figure2.tsv out.pdf '<meta-json>'"""
import csv
import json
import os
import sys
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch, Rectangle

tsv, out, meta = sys.argv[1], sys.argv[2], json.loads(sys.argv[3])
# y-scale mode for the placement-score panels: linear (default) | symlog | split
mode = sys.argv[4] if len(sys.argv) > 4 else "linear"
order, labels, muts, covs = meta["order"], meta["labels"], meta["muts"], meta["cov"]
SIM_C = ["#9ECAE1", "#08519C"]          # two mutation rates (light->dark blue)
REAL_C = "#E6550D"                       # real (orange)
RAND_C = "#BBBBBB"                       # random (grey)
# Random-placement points, styled by the drawn node's taxon. Each tuple is
# (class, label, color, marker); most-divergent class last so it draws on top.
RAND_STYLES = {
    # RSV random bimodality is same-vs-cross subtype (within-B tight/low, cross A-B far/high)
    "rsv": [("same", "same subtype", "#B0B0B0", "o"), ("cross", "cross subtype", "#4D4D4D", "x")],
}   # TB/SARS random: plain grey, no legend
REST_LABEL = {}

rows = list(csv.DictReader(open(tsv), delimiter="\t"))
# scores[sp][(kind,mut)][cov] = [scores]; rand[sp] = [(dist,cls)]; perf[sp][cov]=(wall,rss)
scores = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
rand = defaultdict(list)
perf = defaultdict(lambda: defaultdict(lambda: ([], [])))
failed = defaultdict(int)        # placements that produced no node (low coverage)
for r in rows:
    sp, kind, mut, cov = r["species"], r["kind"], int(r["mut"]), float(r["coverage"])
    try:
        sc = float(r["score"])
    except ValueError:
        sc = float("nan")
    if sc != sc or sc >= 1e8:    # nan / placement-failure sentinel -> "no call"
        failed[(sp, cov)] += 1
    else:
        scores[sp][(kind, mut)][cov].append(sc)
    if r["wall_s"]:
        perf[sp][cov][0].append(float(r["wall_s"]))
        perf[sp][cov][1].append(float(r["rss_mb"]))
# Random baseline from the tsv 'random' column. Each token is "dist:taxon", one
# draw per placement; ":taxon" is optional so a plain number still parses.
for r in rows:
    if not r["random"]:
        continue
    for tok in r["random"].split(";"):
        p = tok.split(":", 1)
        try:
            rand[r["species"]].append((float(p[0]), p[1] if len(p) > 1 else ""))
        except ValueError:
            pass
if failed:
    print("placement no-calls (dropped):",
          ", ".join(f"{sp}@{cov}x={n}" for (sp, cov), n in sorted(failed.items())))

fig = plt.figure(figsize=(13, 11))
# row 0: accuracy + random ; row 1: zoom (coverage only, autoscaled) ; row 2: runtime/mem
gs = fig.add_gridspec(3, 6, height_ratios=[1.2, 0.9, 0.85], hspace=0.5, wspace=0.6)
axB = fig.add_subplot(gs[2, 0:3]); axC = fig.add_subplot(gs[2, 3:6])
first_ax = None                          # first score panel, for the shared legend


def violin(ax, x, vals, color, w=0.7):
    vals = [v for v in vals if v is not None]   # keep zeros (a 0 score is valid)
    if not vals:
        return
    ax.scatter(x + np.random.uniform(-0.28, 0.28, len(vals)) * w,
               vals, s=6, color=color, alpha=0.55, edgecolors="none", zorder=3)
    ax.plot([x - 0.25, x + 0.25], [np.median(vals)] * 2, color="black", lw=1.5, zorder=4)


def classed_random(ax, x, points, sp, w=1.2):
    """Random-placement column: one point per placement, coloured/shaped by the
    drawn node's taxon (RSV subtype / MTBC species). Grey violin shows the shape."""
    dists = [d for d, _ in points]
    if not dists:
        return
    v = ax.violinplot([dists], positions=[x], widths=w, showextrema=False)
    for b in v["bodies"]:
        b.set_facecolor(RAND_C); b.set_alpha(0.25); b.set_edgecolor("none")
    styles = RAND_STYLES.get(sp, [])
    styled = {cls for cls, *_ in styles}
    # default grey points: everything not in a highlighted class (incl. M. tb, SARS)
    rest = [d for d, c in points if c not in styled]
    if rest:
        ax.scatter(x + np.random.uniform(-0.30, 0.30, len(rest)), rest, s=3,
                   color=RAND_C, alpha=0.4, edgecolors="none", zorder=2)
    for cls, _lab, col, mk in styles:
        ys = [d for d, c in points if c == cls]
        if ys:                                     # no edgecolors="none" so 'x' renders
            ax.scatter(x + np.random.uniform(-0.30, 0.30, len(ys)), ys, s=6,
                       color=col, marker=mk, alpha=0.6, linewidths=0.8, zorder=3)
    ax.plot([x - 0.3, x + 0.3], [np.median(dists)] * 2, color="black", lw=1.5, zorder=4)


def rand_legend(ax, sp, loc="upper right"):
    styles = RAND_STYLES.get(sp)
    if not styles:
        return
    h = []
    if sp in REST_LABEL:
        h.append(plt.Line2D([], [], marker="o", ls="", color=RAND_C,
                            label=REST_LABEL[sp], markersize=4))
    h += [plt.Line2D([], [], marker=mk, ls="", color=col, label=lab, markersize=5)
          for _cls, lab, col, mk in styles]
    ax.add_artist(ax.legend(handles=h, fontsize=6, loc=loc, framealpha=0.9,
                            title="Random node", title_fontsize=6))


def plot_cov(ax, sp, groups):
    """Scatter the coverage-column strip points; return the per-coverage centers."""
    for ci, cov in enumerate(covs):
        base = ci * (len(groups) + 0.6)
        for gi, (kind, mut) in enumerate(groups):
            col = REAL_C if kind == "real" else SIM_C[mut]
            violin(ax, base + gi, scores[sp][(kind, mut)].get(cov, []), col, w=0.8)
    return [ci * (len(groups) + 0.6) + (len(groups) - 1) / 2 for ci in range(len(covs))]


zoom_links = []
for j, sp in enumerate(order):
    groups = [("sim", 0), ("sim", 1)] + ([("real", -1)] if ("real", -1) in scores[sp] else [])
    ncov = len(covs)
    x_left = -0.6
    x_right = (ncov - 1) * (len(groups) + 0.6) + (len(groups) - 1) + 0.6
    # ── row 0: coverage columns + random column ──
    ax = fig.add_subplot(gs[0, 2 * j:2 * j + 2])
    if j == 0:
        first_ax = ax
    centers = plot_cov(ax, sp, groups)
    rbase = ncov * (len(groups) + 0.6) + 0.6
    classed_random(ax, rbase, rand[sp], sp, w=1.2)
    rand_legend(ax, sp)
    ax.set_xticks(centers + [rbase])
    ax.set_xticklabels([f"{c}×" for c in covs] + ["Random\nplacement"], fontsize=8)
    ax.set_xlim(x_left - 0.4, rbase + 1.4)
    if mode == "symlog":
        ax.set_yscale("symlog", linthresh=1); ax.set_ylim(bottom=-0.5)
        ax.set_yticks([0, 1, 10, 100, 1000]); ax.set_yticklabels(["0", "1", "10", "100", "1000"])
    else:
        ax.set_ylim(-0.5, 6000 if sp == "tb" else None)   # cap TB (canettii runs off-scale)
    ax.set_title(labels[sp], fontsize=11, weight="bold")
    ax.grid(axis="y", alpha=0.25)
    if j == 0:
        ax.set_ylabel("Placement score\n(parsimony, clamped ≥ 0)")
    # ── row 1: zoom. Coverage columns only, symlog y. Low-score placements squashed
    #    in row 0 spread out; score 0 (correct) stays visible in the linear window
    #    below linthresh ──
    az = fig.add_subplot(gs[1, 2 * j:2 * j + 2])
    plot_cov(az, sp, groups)
    az.set_xticks(centers); az.set_xticklabels([f"{c}×" for c in covs], fontsize=8)
    az.set_xlim(x_left - 0.2, x_right + 0.2)                # tight to coverage -> wider columns
    pm = [v for (kind, mut) in groups for cov in covs for v in scores[sp][(kind, mut)].get(cov, [])]
    pmax = max(pm) if pm else 1
    az.set_yscale("symlog", linthresh=1)                   # log tail + linear 0–1 window (0 = correct)
    ytop = pmax * 1.6 + 1
    az.set_ylim(-0.6, ytop)
    ticks = [t for t in (0, 1, 10, 100, 1000, 10000) if t <= ytop]
    az.set_yticks(ticks); az.set_yticklabels([str(t) for t in ticks])
    az.grid(axis="y", alpha=0.25); az.tick_params(labelsize=8)
    az.set_xlabel("Coverage", fontsize=8)
    if j == 0:
        az.set_ylabel("Placement score\n(zoom, symlog, no random)", fontsize=9)
    zoom_links.append((ax, az, x_left, x_right, pmax))

# zoom connectors: dashed region box on the accuracy panel + splayed lines to the zoom
for ax, az, x_left, x_right, pmax in zoom_links:
    ax.add_patch(Rectangle((x_left, 0), x_right - x_left, pmax * 1.06, fill=False,
                           ec="0.45", lw=0.7, ls="--", zorder=6, clip_on=False))
    zx0, zx1 = az.get_xlim(); zy1 = az.get_ylim()[1]
    for xa, xz in [(x_left, zx0), (x_right, zx1)]:
        fig.add_artist(ConnectionPatch(xyA=(xa, 0), coordsA=ax.transData,
                                       xyB=(xz, zy1), coordsB=az.transData,
                                       color="0.45", lw=0.8, zorder=6))
handles = [plt.Line2D([], [], marker="o", ls="", color=SIM_C[0], label=f"Sim μ={muts[order[0]][0]:g}"),
           plt.Line2D([], [], marker="o", ls="", color=SIM_C[1], label=f"Sim μ={muts[order[0]][1]:g}"),
           plt.Line2D([], [], marker="o", ls="", color=REAL_C, label="Real"),
           plt.Line2D([], [], marker="o", ls="", color=RAND_C, label="Random"),
           plt.Line2D([], [], color="black", lw=1.5, label="Median")]
first_ax.legend(handles=handles, fontsize=7, loc="upper left", framealpha=0.9)

# B runtime, C memory: line plots vs coverage, per species (matches the manuscript)
SP_LINE = ["#377EB8", "#E41A1C", "#4DAF4A"]
for k, (ax, idx, ylab) in enumerate([(axB, 0, "Runtime (s)"), (axC, 1, "Peak RSS (MB)")]):
    for si, sp in enumerate(order):
        xs = [c for c in covs if perf[sp][c][idx]]
        ys = [np.median(perf[sp][c][idx]) for c in xs]
        er = [np.std(perf[sp][c][idx]) for c in xs]
        ax.errorbar(xs, ys, yerr=er, marker="o", color=SP_LINE[si], label=labels[sp], capsize=3)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xticks(covs); ax.set_xticklabels([f"{c}×" for c in covs])
    ax.set_xlabel("Coverage"); ax.set_ylabel(ylab); ax.grid(alpha=0.25)
    ax.legend(fontsize=8)

fig.savefig(out, bbox_inches="tight")
fig.savefig(out.rsplit(".", 1)[0] + ".png", dpi=200, bbox_inches="tight")
print("wrote", out)
