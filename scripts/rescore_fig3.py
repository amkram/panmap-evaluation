#!/usr/bin/env python3
"""Re-score cached Fig-3 consensus assemblies under both accuracy conventions
without re-running the assemble rule. For every held-out sample in work/<sp>/fig3/,
align each method's consensus to the truth genome and report base-based vs legacy
event-based accuracy plus raw error counts. Writes results/figure3_rescored.tsv
and prints an event-vs-base summary.

Scoring identical consensus files isolates the metric choice from run-to-run
assembly variation.

Usage: python3 scripts/rescore_fig3.py [minimap2] [mask_bp]
"""
import glob
import os
import statistics
import sys
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(__file__))
import common as C

MINIMAP2 = sys.argv[1] if len(sys.argv) > 1 else "minimap2"
MASK_BP = int(sys.argv[2]) if len(sys.argv) > 2 else 250
ROOT = os.path.join(os.path.dirname(__file__), "..")
# RSV A/B subtype per node, from the competitive-mapping classifier (classify_rsv_subtype.py).
# Joined to each consensus by the node id in its truth.fa header, this lets plot_fig3_revised
# split the HaphPIPE arm into RSV-A/RSV-B by an exact node join rather than an accuracy match.
SUBMAP_PATH = sys.argv[3] if len(sys.argv) > 3 else os.path.join(ROOT, "meta", "rsv_subtype.tsv")


def load_submap(path):
    m = {}
    if os.path.exists(path):
        for line in open(path).read().splitlines()[1:]:
            c = line.split("\t")
            if len(c) >= 2:
                m[c[0]] = c[1]
    return m


def node_of(truth_fa):
    """The node id in a *.truth.fa header (samtools faidx writes '>node')."""
    with open(truth_fa) as f:
        h = f.readline()
    return h[1:].split()[0].strip() if h.startswith(">") else ""


SUBMAP = load_submap(SUBMAP_PATH)

# method label -> consensus-file suffix, per species (see Snakefile assemble rule).
STD_SUFFIX = {"rsv": "hp.consensus.fa", "sars": "ncbi.consensus.fa", "tb": "cw.consensus.fa"}
def methods(sp):
    return {
        "panmap":           "pm.consensus.fa",
        "standard":         STD_SUFFIX[sp],
        "bwa_ivar":         "bwa.fa",
        "panmap_bwa_ivar":  "pbwa.fa",
    }

HEADER = ("species\tmethod\tcoverage\taccuracy\taccuracy_event\t"
          "snps\tindel_events\tdel_bases\tins_bases\taligned\tinterior\tsubtype")


def jobs():
    for sp in ("rsv", "sars", "tb"):
        d = os.path.join(ROOT, "work", sp, "fig3")
        for truth in sorted(glob.glob(os.path.join(d, "*.truth.fa"))):
            pre = truth[:-len(".truth.fa")]
            ri_cov = os.path.basename(pre)
            cov = ri_cov.rsplit("_", 1)[1]
            for mth, suf in methods(sp).items():
                cons = f"{pre}.{suf}"
                if os.path.exists(cons) and os.path.getsize(cons) > 0:
                    yield (sp, mth, cov, cons, truth)


def score(job):
    sp, mth, cov, cons, truth = job
    s = C.assembly_scores(MINIMAP2, cons, truth, MASK_BP)
    sub = SUBMAP.get(node_of(truth), "") if sp == "rsv" else ""
    return (sp, mth, cov, sub, s)


def main():
    js = list(jobs())
    print(f"re-scoring {len(js)} cached consensus assemblies "
          f"(minimap2={MINIMAP2}, mask_bp={MASK_BP})...", file=sys.stderr)
    with ThreadPoolExecutor(max_workers=8) as ex:
        rows = list(ex.map(score, js))

    out = os.path.join(ROOT, "results", "figure3_rescored.tsv")
    with open(out, "w") as o:
        o.write(HEADER + "\n")
        for sp, mth, cov, sub, s in rows:
            o.write(f"{sp}\t{mth}\t{cov}\t" + "\t".join(str(x) for x in s) + f"\t{sub}\n")
    print(f"wrote {out}", file=sys.stderr)

    # summary: per species x method x coverage, event vs base
    agg = {}
    for sp, mth, cov, sub, s in rows:
        agg.setdefault((sp, mth, cov), []).append(s)
    covs = ["0.5", "1", "1.0", "10", "10.0", "100", "100.0"]
    covkey = lambda c: float(c)
    print("\n=== Fig 3 accuracy: event-based vs base-based (median over samples) ===")
    print(f"{'species':7} {'method':17} {'cov':>6} {'n':>4} "
          f"{'acc_event':>10} {'acc_base':>10} {'Δ(ev-base)':>11} "
          f"{'med_del_bp':>10} {'med_ins_bp':>10}")
    for sp in ("rsv", "sars", "tb"):
        for mth in ("panmap", "standard", "bwa_ivar", "panmap_bwa_ivar"):
            keys = sorted((k for k in agg if k[0] == sp and k[1] == mth),
                          key=lambda k: covkey(k[2]))
            for k in keys:
                ss = agg[k]
                me = statistics.median(x.acc_event for x in ss)
                mb = statistics.median(x.acc_base for x in ss)
                dd = statistics.median(x.acc_event - x.acc_base for x in ss)
                mdel = statistics.median(x.del_bases for x in ss)
                mins = statistics.median(x.ins_bases for x in ss)
                print(f"{sp:7} {mth:17} {k[2]:>6} {len(ss):>4} "
                      f"{me:>10.3f} {mb:>10.3f} {dd:>11.3f} "
                      f"{mdel:>10.0f} {mins:>10.0f}")
        print()


if __name__ == "__main__":
    main()
