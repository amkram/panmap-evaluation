#!/usr/bin/env python3
"""Supplementary figure driver: panmap consensus quality, leaf-only vs internal-allowed.

Runs only the panmap arm of the Fig-3 consensus pipeline (the sole arm whose reference
is chosen by placement, hence the only one force_leaf affects), for every real
leave-one-out sample x coverage, under both placement modes:
  leaf     = --force-leaf   (restrict placement to observed leaves)
  internal = default        (internal/ancestral nodes allowed in the candidate set)

Both modes share the same enriched, subsampled reads, so each sample is a paired
observation. Reuses common.py's exact place/assemble/assembly_scores, and the shared
per-species index (LOO is enforced at placement via exclude_self, and assemble scores
the placed node, so a per-sample index rebuild would be identical). The standard /
baseline arms are mode-independent and deliberately skipped.

Emits results/forceleaf_consensus.tsv:
  species  ri  node  coverage  mode  status  best  aligned  interior  snps  del_bases  ins_bases  acc_base
from which plot_forceleaf.py derives genome fraction (aligned/interior) and per-base
error ((snps+del_bases+ins_bases)/aligned).
"""
import csv
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from shutil import which

import yaml

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
import common as C

CFG = yaml.safe_load(open("config.yaml"))
PANMAP = sys.argv[1] if len(sys.argv) > 1 else CFG["panmap"]
OUT = sys.argv[2] if len(sys.argv) > 2 else "results/forceleaf_consensus.tsv"
SPECIES = (sys.argv[3].split(",") if len(sys.argv) > 3 and sys.argv[3] else list(CFG["species"]))
MAX_REAL = int(sys.argv[4]) if len(sys.argv) > 4 else 0
JOBS = int(sys.argv[5]) if len(sys.argv) > 5 else 48

PANMANUTILS = os.path.join(os.path.dirname(PANMAP), "panmanUtils")
SP = CFG["species"]
COV = [str(c) for c in CFG["coverages"]]
SEED = CFG["seed"]
MASK = CFG["mask_bp"]
MM = True                       # panmap genotype prior (mutation_spectrum default)
CW = 0.0                        # completeness_weight: 0 for both modes (only force_leaf differs)
GT_MIN_COV = float(CFG.get("gt_min_cov_frac", 0.80))
MT = 1


def bin_(tool):
    p = which(tool)
    if not p:
        raise FileNotFoundError(tool)
    return p


def real_samples(sp):
    """Locked real set after ground-truth QC + duplicate-group dedup (mirrors Snakefile)."""
    gen = f"work/{sp}/genomes.fa"
    C.ensure_genomes_fa(PANMANUTILS, SP[sp]["panman"], gen, bin_("samtools"))
    max_n = float(CFG.get("gt_max_n_frac", 0.01))
    min_len = float(CFG.get("gt_min_len_frac", 0.90)) * SP[sp]["genome_size"]
    dgf = os.path.join(os.path.dirname(SP[sp]["samples_tsv"]), "duplicate_groups.tsv")
    dup = ({r["node_id"]: r["duplicate_group"] for r in
            csv.DictReader(open(dgf), delimiter="\t")} if os.path.exists(dgf) else {})
    kept, seen = [], set()
    with open(SP[sp]["samples_tsv"]) as f:
        rows = [(r["node"], r["run"]) for r in csv.DictReader(f, delimiter="\t")]
    for node, run in rows:
        seq = C.faidx_get(gen, node, bin_("samtools"))
        if not seq or seq.upper().count("N") / len(seq) > max_n or len(seq) < min_len:
            continue
        g = dup.get(node)
        if g is not None:
            if g in seen:
                continue
            seen.add(g)
        kept.append((node, run))
    return kept[:MAX_REAL] if MAX_REAL else kept


def dup_members(sp):
    dgf = os.path.join(os.path.dirname(SP[sp]["samples_tsv"]), "duplicate_groups.tsv")
    n2m = {}
    if os.path.exists(dgf):
        g2n = {}
        for r in csv.DictReader(open(dgf), delimiter="\t"):
            g2n.setdefault(r["duplicate_group"], set()).add(r["node_id"])
        for members in g2n.values():
            for n in members:
                n2m[n] = members
    return n2m


COLS = ["species", "ri", "node", "coverage", "mode", "status", "best",
        "aligned", "interior", "snps", "del_bases", "ins_bases", "acc_base"]


def run_unit(u):
    """One (species, ri, cov): enrich+subsample once, then place+assemble in both modes."""
    sp, ri, node, run, cov, excl, idx = u
    cfg, pan = SP[sp], SP[sp]["panman"]
    d = f"work/{sp}/fig_forceleaf/{ri}_{cov}"
    os.makedirs(os.path.dirname(d), exist_ok=True)
    pre = d
    rows = []

    def row(mode, status, best="", s=None):
        r = dict(species=sp, ri=ri, node=node, coverage=cov, mode=mode, status=status, best=best,
                 aligned="", interior="", snps="", del_bases="", ins_bases="", acc_base="")
        if s is not None:
            r.update(aligned=s.aligned, interior=s.interior, snps=s.snps,
                     del_bases=s.del_bases, ins_bases=s.ins_bases, acc_base=s.acc_base)
        return r

    truth = C.get_seq(PANMANUTILS, PANMAP, pan, node, pre + ".truth.fa",
                      f"work/{sp}/genomes.fa", bin_("samtools")) and pre + ".truth.fa"
    rr1, rr2 = C.ena_fastqs(run, f"work/{sp}/reads/{run}")
    if not rr1:
        return [row(m, "nocall_reads") for m in ("leaf", "internal")]
    er1, er2 = C.ontarget_reads(bin_("minimap2"), bin_("samtools"), bin_("seqtk"),
                                truth, rr1, rr2, f"work/{sp}/reads/{run}", MT)
    rl = C.mean_read_len(er1)
    if 2 * C.count_reads(er1) * rl / cfg["genome_size"] < GT_MIN_COV * float(cov):
        return [row(m, "nocall_reads") for m in ("leaf", "internal")]
    r1, r2 = C.subsample(bin_("seqtk"), er1, er2, float(cov), cfg["genome_size"], rl,
                         SEED, pre + "_1.fq", pre + "_2.fq")
    for mode, fl in (("leaf", True), ("internal", False)):
        mp = pre + "." + mode
        best, _, _ = C.place(PANMAP, pan, r1, r2, idx, mp + ".pl", MT,
                             exclude_self=excl, force_leaf=fl, completeness_weight=CW)
        if not best:
            rows.append(row(mode, "nocall_place"))
            continue
        cons, _, _ = C.assemble(PANMAP, pan, r1, r2, idx, best, mp + ".pm", MT, mutation_spectrum=MM)
        if not cons:
            rows.append(row(mode, "nocall_cons", best))
            continue
        s = C.assembly_scores(bin_("minimap2"), cons, pre + ".truth.fa", MASK)
        rows.append(row(mode, "ok", best, s))
    return rows


def main():
    # Build a fresh index per species with the CURRENT panmap. The on-disk shared
    # work/{sp}/index.pmi can be a stale format version (TB's was v0 vs the v4 this
    # binary expects), which makes every placement error out to a no-call. LOO is
    # enforced at placement (exclude_self), so one index per species is correct.
    idxs = {}
    for sp in SPECIES:
        idx = f"work/{sp}/fig_forceleaf/index.pmi"
        os.makedirs(os.path.dirname(idx), exist_ok=True)
        print(f"building fresh {sp} index (k,s,l={SP[sp]['ksl']}) ...", flush=True)
        C.build_index(PANMAP, SP[sp]["panman"], *SP[sp]["ksl"], idx, min(16, JOBS))
        idxs[sp] = idx
    units = []
    for sp in SPECIES:
        samples = real_samples(sp)
        n2m = dup_members(sp)
        for ri, (node, run) in enumerate(samples):
            excl = n2m.get(node, {node})
            for cov in COV:
                units.append((sp, ri, node, run, cov, excl, idxs[sp]))
    print(f"{len(units)} work units ({', '.join(f'{sp}:{len(real_samples(sp))}' for sp in SPECIES)} real samples x {len(COV)} cov), jobs={JOBS}", flush=True)
    os.makedirs(os.path.dirname(OUT) or ".", exist_ok=True)
    done = 0
    with open(OUT, "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=COLS, delimiter="\t")
        w.writeheader()
        with ProcessPoolExecutor(max_workers=JOBS) as ex:
            futs = {ex.submit(run_unit, u): u for u in units}
            for fut in as_completed(futs):
                u = futs[fut]
                try:
                    for r in fut.result():
                        w.writerow(r)
                except Exception as e:
                    sys.stderr.write(f"FAIL {u[0]} ri={u[1]} cov={u[4]}: {e}\n")
                done += 1
                if done % 20 == 0:
                    fo.flush()
                    print(f"  {done}/{len(units)} units done", flush=True)
    print(f"wrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
