#!/usr/bin/env python3
"""Reference-selection accuracy of panmap vs mash vs minimap2 as the candidate-genome
pool grows.

Samples are the Fig-3 QC-passed real samples (work/<sp>/qc_pass.tsv). Leave-one-out:
all test genomes are removed from the candidate pool before subsampling, so no method
can select a test sample's own genome. For each batch size x we draw the same seeded
x-genome subset and hand it to every method:
  - mash, minimap2 : the x leaf genomes as FASTA (can only pick a leaf)
  - panmap         : a --subnet PanMAN of those x leaves, which also exposes their
                     induced internal/ancestral nodes (panmap may pick those)
Each method selects one reference from the sample's reads; we score it against the
sample's truth genome with the Fig-3 genotyping-accuracy metric (completeness x
correctness, 250 bp flanks masked).

Usage: reviewer_refsel.py <out_tsv> [species_csv] [batches_csv] [depth] [workers]
"""
import concurrent.futures as cf
import csv
import os
import subprocess
import sys

import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C

CFG = yaml.safe_load(open(os.path.join(os.path.dirname(__file__), "..", "config.yaml")))
OUT = sys.argv[1]
SPECIES = (sys.argv[2].split(",") if len(sys.argv) > 2 else ["rsv", "sars", "tb"])
BATCHES = [int(x) for x in (sys.argv[3].split(",") if len(sys.argv) > 3
                            else "2,4,8,16,32,64,128,256".split(","))]
DEPTH = float(sys.argv[4]) if len(sys.argv) > 4 else 30.0      # reads subsample depth
WORKERS = int(sys.argv[5]) if len(sys.argv) > 5 else 6
RL, SEED = 150, 42

# Run from evaluation/: config paths are relative to that cwd.
PANMAP = os.path.abspath(CFG["panmap"])                    # ../build/bin/panmap
# Patched panmanUtils (fixes 6 --subnet serialization bugs); kept outside the git
# worktree so a repo cleanup can't delete it. Falls back to the in-tree build.
_PU_PATCHED = "/tmp/claude-1015/-scratch1-alex-poopdoop/638cc60a-73b8-44fb-97a7-8108b1b218c5/scratchpad/panmanUtils.patched"
PU = _PU_PATCHED if os.path.exists(_PU_PATCHED) else os.path.abspath("../build/bin/panmanUtils")
CFG_PANMAN = {sp: os.path.abspath(CFG["species"][sp]["panman"]) for sp in CFG["species"]}
BIN = "/home/alex/micromamba/envs/panmap-eval/bin"
MASH = "/home/alex/micromamba/envs/refsel/bin/mash"
B = lambda x: os.path.join(BIN, x)
MM = B("minimap2"); SAM = B("samtools"); SEQTK = B("seqtk")


def prep_reads(sp, node, run, truth_fa, depth):
    """On-target-enriched reads for a sample, subsampled to `depth`x. Cached."""
    d = f"work/{sp}/reads/{run}"
    r1, r2 = C.ena_fastqs(run, d)
    if not r1:
        return None, None
    er1, er2 = C.ontarget_reads(MM, SAM, SEQTK, truth_fa, r1, r2, d, threads=4)
    od = f"work/{sp}/refsel"; os.makedirs(od, exist_ok=True)
    s1, s2 = f"{od}/{run}_{depth:g}_1.fq", f"{od}/{run}_{depth:g}_2.fq"
    if not (os.path.exists(s1) and os.path.getsize(s1)):
        C.subsample(SEQTK, er1, er2, depth, CFG["species"][sp]["genome_size"], RL,
                    SEED, s1, s2)
    return s1, s2


def build_candidates(sp, pool, x, seed, gen):
    """Draw x genomes from pool (seeded) and build the shared candidate set: a --subnet
    panman (+ index) for panmap, a concatenated FASTA + mash sketch for minimap2/mash.
    Returns dict of paths, or None on failure."""
    import random
    rng = random.Random(seed)
    picks = rng.sample(pool, x)
    d = f"work/{sp}/refsel/cand_{x}_{seed}"; os.makedirs(d, exist_ok=True)
    # leaf FASTA (concatenated) for mash + minimap2
    concat = f"{d}/genomes.fa"
    if not (os.path.exists(concat) and os.path.getsize(concat)):
        with open(concat, "w") as o:
            for n in picks:
                seq = C.faidx_get(gen, n, SAM)
                if seq:
                    o.write(f">{n}\n{seq}\n")
    msh = f"{d}/genomes.msh"
    if not os.path.exists(msh):
        C.sh([MASH, "sketch", "-o", f"{d}/genomes", concat])
    # subnet panman (P_x leaves + their ancestral internal nodes) + panmap index.
    # panmanUtils --subnet input is "<treeId> <leaf1> <leaf2> ..." on one line and
    # writes to ./panman/<name>.panman relative to the process cwd, so run it in `d`.
    subpan = f"{d}/panman/sub.panman"
    if not os.path.exists(subpan):
        open(f"{d}/nodes.txt", "w").write("0 " + " ".join(picks) + "\n")
        subprocess.run([PU, "-I", CFG_PANMAN[sp], "--subnet", "-i", "nodes.txt", "-o", "sub"],
                       cwd=d, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if not (os.path.exists(subpan) and os.path.getsize(subpan)):
        return None
    idx = f"{d}/sub.pmi"
    if not os.path.exists(idx):
        k, s, l = CFG["species"][sp]["ksl"]
        C.build_index(PANMAP, subpan, k, s, l, idx, 4)
    return {"picks": picks, "concat": concat, "msh": msh, "subpan": subpan,
            "idx": idx, "dir": d}


def select_mash(cand, s1, pre):
    """Closest candidate genome by Mash distance (reads sketch vs genome sketches)."""
    rmsh = pre + ".reads"
    C.sh([MASH, "sketch", "-r", "-m", "2", "-o", rmsh, s1])
    out = C.sh([MASH, "dist", cand["msh"], rmsh + ".msh"]).stdout
    best, bd = None, 1e9
    for ln in out.splitlines():
        f = ln.split("\t")
        if len(f) >= 3:
            try:
                dv = float(f[2])
            except ValueError:
                continue
            if dv < bd:
                bd, best = dv, f[0]
    return best


def select_minimap2(cand, s1, s2, pre):
    """Candidate genome that captures the most mapped reads (minimap2 -ax sr)."""
    bam = pre + ".mm.bam"
    p = subprocess.Popen([MM, "-ax", "sr", "-t", "4", cand["concat"], s1, s2],
                         stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    vf = subprocess.Popen([SAM, "view", "-b", "-F", "0x904", "-"], stdin=p.stdout,
                          stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    p.stdout.close()
    with open(bam, "wb") as bf:
        srt = subprocess.Popen([SAM, "sort", "-o", "-", "-"], stdin=vf.stdout,
                               stdout=bf, stderr=subprocess.DEVNULL)
        vf.stdout.close(); srt.wait()
    p.wait(); vf.wait()
    C.sh([SAM, "index", bam])
    best, bc = None, -1
    for ln in C.sh([SAM, "idxstats", bam]).stdout.splitlines():
        f = ln.split("\t")
        if len(f) >= 3 and f[0] != "*":
            c = int(f[2])
            if c > bc:
                bc, best = c, f[0]
    return best


def acc_of(sp, selected_seq_fa, truth_fa):
    return C.genotyping_accuracy(MM, selected_seq_fa, truth_fa,
                                 CFG["species"][sp]["genome_size"], CFG["mask_bp"])


def one_sample(sp, node, run, x, seed, cand, gen):
    """All three selections for one (sample, batch). Returns list of result rows."""
    od = f"work/{sp}/refsel"; pre = f"{od}/{run}_{x}_{seed}"
    truth = f"{od}/{run}.truth.fa"
    if not (os.path.exists(truth) and os.path.getsize(truth)):
        seq = C.faidx_get(gen, node, SAM)
        if not seq:
            return []
        open(truth, "w").write(f">{node}\n{seq}\n")
    s1, s2 = prep_reads(sp, node, run, truth, DEPTH)
    if not s1:
        return []
    rows = []
    # panmap: place on the subnet panman, allow internal nodes
    best, _, _ = C.place(PANMAP, cand["subpan"], s1, s2, cand["idx"], pre + ".pm",
                         4, exclude_self=node, force_leaf=False)
    sel = {"panmap": best,
           "mash": select_mash(cand, s1, pre),
           "minimap2": select_minimap2(cand, s1, s2, pre)}
    for meth, snode in sel.items():
        a = ""
        if snode:
            sfa = pre + f".{meth}.sel.fa"
            if meth == "panmap":
                # selected node may be an internal/ancestral node -> dump from the subnet
                got = C.dump_seq(PANMAP, cand["subpan"], snode, sfa)
            else:
                got = C.faidx_seq(gen, snode, sfa, SAM)   # a leaf in the full genomes.fa
            if got:
                a = acc_of(sp, sfa, truth)
        rows.append((sp, meth, x, seed, run, snode or "", a))
    return rows


def main():
    out_rows = [("species", "method", "batch", "seed", "run", "selected", "accuracy")]
    for sp in SPECIES:
        qcp = f"work/{sp}/qc_pass.tsv"
        if not os.path.exists(qcp):
            print(f"{sp}: no qc_pass, skipping", flush=True); continue
        samples = [(r["node"], r["run"]) for r in csv.DictReader(open(qcp), delimiter="\t")]
        gen = f"work/{sp}/genomes.fa"
        C.ensure_genomes_fa(PU, CFG_PANMAN[sp], gen, SAM)
        test_nodes = {n for n, _ in samples}
        pool = [l.split("\t")[0] for l in open(gen + ".fai")
                if l.split("\t")[0] not in test_nodes]           # leave-one-out pool
        print(f"{sp}: {len(samples)} samples, pool={len(pool)}", flush=True)
        for x in BATCHES:
            if x > len(pool):
                continue
            cand = build_candidates(sp, pool, x, SEED, gen)
            if not cand:
                print(f"  {sp} x={x}: subnet build failed", flush=True); continue
            with cf.ThreadPoolExecutor(max_workers=WORKERS) as ex:
                futs = [ex.submit(one_sample, sp, n, r, x, SEED, cand, gen)
                        for n, r in samples]
                for fu in cf.as_completed(futs):
                    out_rows.extend(fu.result())
            done = sum(1 for r in out_rows if r[2] == x and r[0] == sp) // 3
            print(f"  {sp} x={x}: {done}/{len(samples)} samples done", flush=True)
            with open(OUT, "w") as f:                            # incremental write
                for r in out_rows:
                    f.write("\t".join(str(v) for v in r) + "\n")
    print("wrote", OUT, flush=True)


if __name__ == "__main__":
    main()
