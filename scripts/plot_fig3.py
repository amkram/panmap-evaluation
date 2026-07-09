#!/usr/bin/env python3
"""Figure 3: leave-one-out consensus assembly. Two accuracy rows (paired arms)
plus a runtime row, one column per species:
  Row 1 (accuracy): Panmap (full, native genotyping)  vs  single-ref BWA+iVar + impute
  Row 2 (accuracy): Panmap->BWA+iVar (no impute)       vs  single-ref BWA+iVar (no impute)
  Row 3 (runtime):  the distinct pipelines
The single-reference BWA+iVar arms use, per sample, whichever reference the reads
align best to (for RSV, subtype A NC_038235.1 or B LR699737.1); all samples are
pooled into one baseline category. Accuracy (the `accuracy` column) = base-based
% of the held-out genome correctly reconstructed = matches/interior (250 bp flanks
ignored); the TSV also carries `accuracy_event` (legacy event-based) and raw error
counts. Usage: plot_fig3.py figure3.tsv out.pdf '<meta>'"""
import csv
import json
import re
import sys
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

tsv, out, meta = sys.argv[1], sys.argv[2], json.loads(sys.argv[3])
order, labels, baseline, covs = meta["order"], meta["labels"], meta["baseline"], meta["cov"]
# Optional per-species y-limits for the accuracy row (e.g. widen for base-based
# M. tb, which dips to ~98.9%). {sp: [lo, hi]}; falls back to the defaults below.
YLIMS = {sp: tuple(v) for sp, v in meta.get("ylims", {}).items()}
STYLE = sys.argv[4] if len(sys.argv) > 4 else "box"    # box | violin | scatter
LOGY = len(sys.argv) > 5 and sys.argv[5] == "logy"     # log-scale the runtime/memory rows
np.random.seed(0)                                       # reproducible scatter jitter

PANMAP_C, STD_C, PANREF_C, REF_C = "#6A3D9A", "#E08214", "#1B9E77", "#C0392B"
BLAB = {"bwa_ivar": "BWA+iVar", "clockwork": "Clockwork"}
# Field-standard Row-1 pipeline name, per species.
STDLAB = {"rsv": "HaphPIPE", "sars": "NCBI Pipeline (NC_045512.2 ref)",
          "tb": "Clockwork (NC_000962.3 ref)"}
STD_A_C, STD_B_C = "#E08214", "#FDBF6F"   # RSV HaphPIPE, subtype-matched (A / B)


def sp_name_size(sp):
    """Split labels[sp] like 'RSV (4K)' into ('RSV', '4K'): name is the size-free
    subplot title, size drives the 'N-sample pangenome' Panmap legend tag."""
    m = re.match(r"^(.*?)\s*\(([^)]*)\)\s*$", labels[sp])
    return (m.group(1), m.group(2)) if m else (labels[sp], "")


def arm(mth):
    """Classify a method into an arm. RSV's field-standard (HaphPIPE) rows carry a
    subtype tag (standard_A/standard_B) and stay split; other trailing _A/_B tags
    (per-sample baseline reference) pool into the single baseline category."""
    if mth.startswith("standard"):         # standard, standard_A, standard_B (kept distinct)
        return mth
    if mth.endswith("_A") or mth.endswith("_B"):
        mth = mth[:-2]
    if mth == "panmap":
        return "panmap"
    if mth.startswith("panmap_"):
        return "panmap_ref"
    return "ref"


acc = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))   # [sp][arm][cov]
rt = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
mem = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))   # peak RSS (MB)
for r in csv.DictReader(open(tsv), delimiter="\t"):
    sp, a, cov = r["species"], arm(r["method"]), float(r["coverage"])
    if r["accuracy"] not in ("", "nan"):
        acc[sp][a][cov].append(float(r["accuracy"]))
    if r["wall_s"]:
        rt[sp][a][cov].append(float(r["wall_s"]))
    if r.get("peak_mb") not in (None, "", "nan"):
        mem[sp][a][cov].append(float(r["peak_mb"]))

# Combine the (possibly subtype-split) standard arms into one "std" series for the
# runtime/memory rows, which compare Panmap vs the field-standard pipeline.
for D in (rt, mem):
    for sp in list(D):
        merged = defaultdict(list)
        for a, byc in list(D[sp].items()):
            if a.startswith("standard"):
                for c, vs in byc.items():
                    merged[c] += vs
        if merged:
            D[sp]["std"] = merged

order = [sp for sp in order if sp in acc]
ncol = len(order)
fig, axes = plt.subplots(3, ncol, figsize=(4.16 * ncol, 8.6), squeeze=False,
                         gridspec_kw={"height_ratios": [1, 0.5, 0.5]})


def boxrow(axrow, pairs_of, show_n=True, ylims=None, logy=False):
    """One accuracy row: N boxplot arms side by side per coverage, per species.
    pairs_of is a list (same arms for every species) or a callable
    sp -> [(arm, color, label), ...]. logy uses a symlog y-axis (arms hit exactly
    0 at low coverage, so plain log can't render them)."""
    getp = pairs_of if callable(pairs_of) else (lambda sp: pairs_of)
    for j, sp in enumerate(order):
        ax = axrow[j]
        pos = np.arange(len(covs))
        pairs = getp(sp)
        n = len(pairs)
        offs = [-0.35 + 0.70 * (k + 0.5) / n for k in range(n)]
        width = 0.70 / n * 0.85
        for (a, col, _), off in zip(pairs, offs):
            data = [acc[sp][a].get(c, []) for c in covs]
            xp = pos + off
            if STYLE == "box":
                bp = ax.boxplot(data, positions=xp, widths=width, patch_artist=True,
                                showfliers=True, flierprops=dict(marker=".", ms=3, alpha=.5))
                for b in bp["boxes"]:
                    b.set_facecolor(col); b.set_alpha(0.55)
                for med in bp["medians"]:
                    med.set_color("black")
            elif STYLE == "violin":
                idx = [i for i in range(len(covs)) if len(data[i]) > 1]
                if idx:
                    vp = ax.violinplot([data[i] for i in idx], positions=[xp[i] for i in idx],
                                       widths=width, showmedians=True, showextrema=False)
                    for b in vp["bodies"]:
                        b.set_facecolor(col); b.set_alpha(0.22)
                        b.set_edgecolor(col); b.set_linewidth(0.6)
                    vp["cmedians"].set_color("black"); vp["cmedians"].set_linewidth(1.2)
                for i in range(len(covs)):          # overlay the actual points (jittered)
                    vals = data[i]
                    if not vals:
                        continue
                    jit = np.random.uniform(-width * 0.35, width * 0.35, len(vals))
                    ax.scatter(xp[i] + jit, vals, s=6, color=col, alpha=0.55,
                               edgecolors="none", zorder=3)
            else:                                    # scatter: jittered dots + median bar
                for i in range(len(covs)):
                    vals = data[i]
                    if not vals:
                        continue
                    jit = np.random.uniform(-width * 0.45, width * 0.45, len(vals))
                    ax.scatter(xp[i] + jit, vals, s=9, color=col, alpha=0.5, edgecolors="none")
                    ax.plot([xp[i] - width * 0.45, xp[i] + width * 0.45],
                            [np.median(vals)] * 2, color="black", lw=1.2)
        ax.set_xticks(pos); ax.set_xticklabels([f"{c}×" for c in covs])
        if logy:
            ax.set_yscale("symlog", linthresh=1)
            ax.set_ylim(0, 130)
            ax.set_yticks([0, 1, 10, 100]); ax.set_yticklabels(["0", "1", "10", "100"])
        else:
            ax.set_ylim(*(ylims or {}).get(sp, (-3, 103)))
        nn = max((len(acc[sp][pairs[0][0]].get(c, [])) for c in covs), default=0)
        ax.set_title(sp_name_size(sp)[0] + (f"  (n={nn})" if show_n else ""), fontsize=10)
        if j == 0:
            ax.set_ylabel("Correctly genotyped (%)")
        ax.grid(axis="y", alpha=0.25)
        legl = [(l.get(sp, next(iter(l.values()))) if isinstance(l, dict) else l)
                for _, _, l in pairs]
        ax.legend([plt.Rectangle((0, 0), 1, 1, fc=c, alpha=.55) for _, c, _ in pairs],
                  legl, fontsize=7, loc="lower right")


blab = BLAB.get(baseline[order[0]], "BWA+iVar") if order else "BWA+iVar"
# Row 1: Panmap (full) vs the field-standard reference-based pipeline (HaphPIPE for
# RSV, NCBI SC2VC for SARS, Clockwork for TB) ;  Row 2: Panmap->ref vs single-ref.
# RSV row-1 arms cluster at ~98-100% with divergent-A low outliers -> start at 40;
# SARS row-1 arms both sit at ~99.8-100% -> zoom to 99.0.
def row0_pairs(sp):
    size = sp_name_size(sp)[1]
    p = [("panmap", PANMAP_C, "Panmap" + (f" ({size}-sample pangenome)" if size else ""))]
    if sp == "rsv":                        # HaphPIPE split by subtype (subtype-matched ref)
        p += [("standard_B", STD_B_C, "HaphPIPE (RSV-B samples, LR699737.1 ref)"),
              ("standard_A", STD_A_C, "HaphPIPE (RSV-A samples, NC_038235.1 ref)")]
    else:
        p += [("standard", STD_C, STDLAB[sp])]
    return p


# Main figure: Row 1 (accuracy vs field-standard pipeline) + runtime row.
boxrow(axes[0], row0_pairs,
       ylims=YLIMS or {"rsv": (40, 101), "sars": (99.0, 100.1), "tb": (99.0, 100.1)})

# Half-height rows: Panmap vs the single-reference pipeline, per coverage.
def linerow(axrow, vals, ylabel, title, xlabel=True, legend=True, ymax=None, logy=False):
    for j, sp in enumerate(order):
        ax = axrow[j]
        size = sp_name_size(sp)[1]
        series = [("panmap", PANMAP_C, "Panmap" + (f" ({size}-sample pangenome)" if size else ""), "-"),
                  ("std", STD_C, STDLAB[sp], "--")]
        for a, col, lab, ls in series:
            xs = [c for c in covs if vals[sp][a].get(c)]
            if not xs:
                continue
            ys = [np.median(vals[sp][a][c]) for c in xs]
            er = [np.std(vals[sp][a][c]) for c in xs]
            if logy:                                       # keep the lower whisker positive
                er = [min(e, y * 0.999) for e, y in zip(er, ys)]
            ax.errorbar(xs, ys, yerr=er, marker="o", ms=4, color=col, label=lab, ls=ls, lw=1.8, capsize=3)
        ax.set_xscale("log"); ax.set_xticks(covs); ax.set_xticklabels([f"{c}×" for c in covs])
        if logy:
            ax.set_yscale("log")
        if xlabel:
            ax.set_xlabel("Read depth (coverage)")
        if j == 0:
            ax.set_ylabel(ylabel)
        if ymax is not None:
            ax.set_ylim(0, ymax)
        ax.set_title(f"{title}, {sp_name_size(sp)[0]}", fontsize=9)
        ax.grid(alpha=0.25)
        if legend:
            ax.legend(fontsize=7)


linerow(axes[1], rt, "Runtime (s)", "Runtime", xlabel=False, logy=LOGY)
linerow(axes[2], mem, "Peak memory (MB)", "Memory", logy=LOGY)

fig.tight_layout()
fig.savefig(out, bbox_inches="tight")
fig.savefig(out.rsplit(".", 1)[0] + ".png", dpi=200, bbox_inches="tight")
print("wrote", out)

# Supplementary figure (figure_S): the read-pileup arms, panmap-selected reference
# vs single reference (both BWA+iVar, no impute) -- the coverage-limited comparison.
figS, axesS = plt.subplots(1, ncol, figsize=(5.2 * ncol, 4.2), squeeze=False)
pan_ref_lab = {sp: f"Panmap→{blab}" + (f" ({sz}-sample pangenome)" if (sz := sp_name_size(sp)[1]) else "")
               for sp in order}
boxrow(axesS[0], [("panmap_ref", PANREF_C, pan_ref_lab),
                  ("ref", REF_C, f"{blab} (single ref)")])
figS.suptitle("Read-pileup genotyping: panmap-selected vs single reference "
              "(BWA+iVar, no impute)", fontsize=12, y=1.02)
figS.tight_layout()
outS = out.replace("figure3", "figure_S")
figS.savefig(outS, bbox_inches="tight")
figS.savefig(outS.rsplit(".", 1)[0] + ".png", dpi=200, bbox_inches="tight")
print("wrote", outS)
