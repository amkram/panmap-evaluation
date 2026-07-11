# Reproduce Figure 2 (phylogenetic placement) and Figure 3 (consensus assembly).
#
# Leave-one-out uses the current panmap build (no --remove-node): the held-out
# leaf is dropped from the candidate set after scoring (placement) and never used
# as the reference (assembly) -- equivalent to pruning the leaf.
#
#   snakemake -j8 --configfile config.yaml            # everything
#   snakemake -j8 results/figure2.pdf                 # Fig 2 only
#   snakemake -n                                      # dry run / DAG check
import csv
import math
import os
import sys

from snakemake.utils import min_version
min_version("8.0")          # older snakemake mangles f-strings in its parser

sys.path.insert(0, os.path.join(workflow.basedir, "scripts"))
import common as C

configfile: "config.yaml"

# panmap binary. Built separately (not by this workflow) -- slurm_run.sbatch runs
# scripts/build_panmap.sh to build it from the latest main before the run, or build
# it yourself and point `panmap:` at it. See the README.
PANMAP = config["panmap"]
PANMANUTILS = os.path.join(os.path.dirname(PANMAP), "panmanUtils")
# Thread count for every benchmarked step (panmap index/place/assemble and the
# baselines). Held at 1 for a thread-fair Fig-3 runtime/memory comparison: the
# standard pipelines are all single-threaded (HaphPIPE --ncpu 1, NCBI -threads 1,
# Clockwork --cpus 1, BWA -t 1), so panmap is too. With every tool at 1 thread,
# each timed rule reserves exactly `MT` Snakemake core(s), so `snakemake --cores N`
# runs up to N jobs at once with no CPU oversubscription: each job owns a physical
# core and its wall-clock is contention-free even at full-node parallelism (see
# slurm_run.sbatch).
MT = int(config.get("method_threads", 1))

# Peak-RSS budget per job (MB). With `--resources mem_mb=<~node RAM>` it caps
# concurrency so total residency never exceeds RAM; swapping would corrupt both the
# wall-clock and the peak-RSS the figures report. TB's panman index dominates
# residency (~20 GB); the standard-pipeline JVMs (GATK/Picard, HaphPIPE) want a few
# GB. On a 1 TB node this allows ~27 TB jobs, or hundreds of RSV/SARS jobs
# concurrently: memory-bound for TB, core-bound (--cores) for the small genomes.
_JOB_MEM = {"index":    {"tb": 34000, "_": 4000},
            "place":    {"tb": 34000, "_": 3000},
            "assemble": {"tb": 34000, "_": 8000}}
def job_mem(rule, sp):
    d = _JOB_MEM[rule]
    return d.get(sp, d["_"])
SP = config["species"]
COV = [str(c) for c in config["coverages"]]
RL = config["read_length"]
SEED = config["seed"]
# Optional scope knobs (defaults reproduce the full manuscript run):
#   --config only_species=rsv          restrict to a subset (comma-separated)
#   --config sim_replicates=6          fewer simulated leaves per species
#   --config max_real=3                cap real (SRA) samples per species
#   --config test_mode=true            fast smoke run: 20% of RSV+SARS, all plots
# Test mode (or subset_frac on its own) takes a fraction of both the simulated leaves
# and the real samples, per species, so a run exercises every rule + plot at ~1/5 cost.
TEST = str(config.get("test_mode", "")).strip().lower() in ("true", "1", "yes", "on")
SUBSET = float(config.get("subset_frac", 0.2 if TEST else 0) or 0)   # 0 = full data
_only = str(config.get("only_species", "")).strip()
_default_species = [s for s in ("rsv", "sars") if s in SP] if TEST else list(SP)
SPECIES = [s for s in _only.split(",") if s] if _only else _default_species
MAX_REAL = int(config.get("max_real", 0))      # 0 = all real samples
MUT = [0, 1]                                   # two per-species mutation rates
# sim_leaf still draws from the full replicate list, so the subset is a reproducible
# prefix of the full run's leaves.
REPS = range(max(1, math.ceil(SUBSET * config["sim_replicates"])) if SUBSET > 0
             else config["sim_replicates"])


def _cap_real(lst):
    """Cap the real-sample list: hard MAX_REAL if set, else a SUBSET fraction (ceil,
    >=1) per species for test/smoke runs, else the whole list."""
    if MAX_REAL:
        return lst[:MAX_REAL]
    if SUBSET > 0:
        return lst[:max(1, math.ceil(SUBSET * len(lst)))]
    return lst
RESULT = "results"
WORK = "work"

# Placement-mode knobs -> separate, tagged output trees so variants coexist and
# share the (mode-independent) index/genomes/reads:
#   --config force_leaf=false               allow ancestral (internal) placement
#   --config completeness_weight=1.0        re-rank internal nodes by completeness
def _as_bool(v, default):
    if isinstance(v, bool):
        return v
    if v is None:
        return default
    return str(v).strip().lower() not in ("false", "0", "no", "off", "")
FORCE_LEAF = _as_bool(config.get("force_leaf", True), True)
CW = float(config.get("completeness_weight", 0) or 0)
MM = _as_bool(config.get("mutation_spectrum", True), True)   # panmap genotype prior
_tag = []
if not FORCE_LEAF:
    _tag.append("internal")
if CW > 0:
    _tag.append("cw" + f"{CW:g}".replace(".", "p"))
TAG = ("_" + "_".join(_tag)) if _tag else ""       # placement tag (Fig 2 + Fig 3)
MM_TAG = "" if MM else "_nomm"                      # genotyping tag (Fig 3 only)
F2SIM, F2REAL = f"fig2_sim{TAG}", f"fig2_real{TAG}"
FIG2 = f"figure2{TAG}"
F3 = f"fig3{TAG}{MM_TAG}"                           # Fig 3 depends on placement and genotyping
FIG3 = f"figure3{TAG}{MM_TAG}"


_REAL_SAMPLES = {}
_DUP_MEMBERS = {}
def real_samples(sp):
    """[(node, run)] from the species samples TSV, after ground-truth QC on the
    PanMAN-leaf truth genome:
      (1) drop truth with > gt_max_n_frac Ns or length < gt_min_len_frac*genome_size,
      (2) dedup by duplicate_groups.tsv (one sample per group).
    Read-support (3) is enforced as a runtime no-call in place_real/assemble.
    Cached; capped at MAX_REAL when set (>0) for scoped / smoke runs."""
    if sp not in _REAL_SAMPLES:
        with open(SP[sp]["samples_tsv"]) as f:
            rows = [(r["node"], r["run"]) for r in csv.DictReader(f, delimiter="\t")]
        gen = f"work/{sp}/genomes.fa"
        C.ensure_genomes_fa(PANMANUTILS, SP[sp]["panman"], gen, "samtools")
        max_n = float(config.get("gt_max_n_frac", 0.01))
        min_len = float(config.get("gt_min_len_frac", 0.90)) * SP[sp]["genome_size"]
        dgf = os.path.join(os.path.dirname(SP[sp]["samples_tsv"]), "duplicate_groups.tsv")
        dup = ({r["node_id"]: r["duplicate_group"] for r in
                csv.DictReader(open(dgf), delimiter="\t")} if os.path.exists(dgf) else {})
        kept, seen = [], set()
        for node, run in rows:
            seq = C.faidx_get(gen, node, "samtools")
            if not seq or seq.upper().count("N") / len(seq) > max_n or len(seq) < min_len:
                continue                                          # (1) completeness/ambiguity
            g = dup.get(node)
            if g is not None:
                if g in seen:
                    continue                                      # (2) dedup
                seen.add(g)
            kept.append((node, run))
        _REAL_SAMPLES[sp] = kept
    return _cap_real(_REAL_SAMPLES[sp])


def dup_group_members(sp, node):
    """All node_ids in `node`'s duplicate_group (>=1, incl. node itself). LOO must
    exclude every byte-identical twin from the candidate set, not just the held-out
    node, else panmap can place on an identical copy (dist ~= 0) and inflate
    placement/assembly accuracy. Falls back to {node} if node is in no group.
    Cached per species."""
    if sp not in _DUP_MEMBERS:
        dgf = os.path.join(os.path.dirname(SP[sp]["samples_tsv"]), "duplicate_groups.tsv")
        n2m = {}
        if os.path.exists(dgf):
            g2n = {}
            for r in csv.DictReader(open(dgf), delimiter="\t"):
                g2n.setdefault(r["duplicate_group"], set()).add(r["node_id"])
            for members in g2n.values():
                for n in members:
                    n2m[n] = members
        _DUP_MEMBERS[sp] = n2m
    return _DUP_MEMBERS[sp].get(node, {node})


def sim_leaf(sp, i):
    """Deterministic i-th of N randomly selected leaves for species sp."""
    import random
    par = C.newick_parents(SP[sp]["panman"], PANMANUTILS)
    lv = sorted(C.leaves(par))
    return random.Random(SEED).sample(lv, config["sim_replicates"])[i], par


def fig2_real_inputs():
    return [f"work/{sp}/{F2REAL}/{ri}_{cov}.tsv"
            for sp in SPECIES for ri in range(len(real_samples(sp))) for cov in COV]


def fig3_inputs():
    return [f"work/{sp}/{F3}/{ri}_{cov}.tsv"
            for sp in SPECIES for ri in range(len(real_samples(sp))) for cov in COV]


def fig2_sim_inputs():
    return expand("work/{sp}/" + F2SIM + "/{i}_{m}_{cov}.tsv",
                  sp=SPECIES, i=REPS, m=MUT, cov=COV)


wildcard_constraints:
    sp="|".join(SPECIES), cov="|".join(COV), i=r"\d+", m=r"\d+", ri=r"\d+"


rule all:
    input: f"results/{FIG2}.pdf", f"results/{FIG3}.pdf", f"results/{FIG3}_revised_rate.pdf"


# ── shared index ──────────────────────────────────────────────────────────────
rule index:
    output: idx="work/{sp}/index.pmi", t="work/{sp}/index.time"
    threads: MT
    resources: mem_mb=lambda wc: job_mem("index", wc.sp)
    run:
        k, s, l = SP[wildcards.sp]["ksl"]
        wall, rss = C.build_index(PANMAP, SP[wildcards.sp]["panman"], k, s, l,
                                  output.idx, MT)
        open(output.t, "w").write(f"{wall}\t{rss}\n")


# ── random-node taxon maps (Fig 2A colouring) ─────────────────────────────────
# _score_row reads these at placement time to label each random draw. place_sim/
# place_real depend on the species' map so a clean build generates it before any
# placement; ancient() keeps that dependency from retriggering cached placements.
_CLASS_MAP = {"rsv": "meta/rsv_subtype.tsv", "tb": "meta/tb_species.tsv"}


def class_map(wildcards):
    m = _CLASS_MAP.get(wildcards.sp)
    return [ancient(m)] if m else []


rule rsv_subtype_map:            # competitive A/B mapping to the two reference genomes
    output: "meta/rsv_subtype.tsv"
    threads: 24                  # classify_rsv_subtype.py fans out over a Pool(24); reserve to match
    run:
        gen = "work/rsv/genomes.fa"
        C.ensure_genomes_fa(PANMANUTILS, SP["rsv"]["panman"], gen, SP_BIN("samtools"))
        C.sh(["python3", os.path.join(workflow.basedir, "scripts/classify_rsv_subtype.py"),
              gen, SP["rsv"]["ref_a"], SP["rsv"]["ref_b"], output[0]])


rule tb_species_map:             # coarse MTBC species by distance to H37Rv + M. bovis anchors
    output: "meta/tb_species.tsv"
    threads: 24                  # classify_tb_species.py fans out over a Pool(24); reserve to match
    run:
        gen = "work/tb/genomes.fa"
        C.ensure_genomes_fa(PANMANUTILS, SP["tb"]["panman"], gen, SP_BIN("samtools"))
        C.sh(["python3", os.path.join(workflow.basedir, "scripts/classify_tb_species.py"),
              gen, SP["tb"]["ref"], output[0]])


# ── FIGURE 2 : placement ──────────────────────────────────────────────────────
# Simulated: mutate leaf -> descendant -> reads -> place on full panman; score vs
# the true parent (the mutated leaf). Real: leave-one-out, score vs parent(leaf).
rule place_sim:
    input: idx="work/{sp}/index.pmi", maps=class_map
    output: "work/{sp}/" + F2SIM + "/{i}_{m}_{cov}.tsv"
    threads: MT
    resources: mem_mb=lambda wc: job_mem("place", wc.sp)
    run:
        sp, i, m, cov = wildcards.sp, int(wildcards.i), int(wildcards.m), float(wildcards.cov)
        cfg, pan = SP[sp], SP[wildcards.sp]["panman"]
        leaf, par = sim_leaf(sp, i)
        d = os.path.dirname(output[0]); os.makedirs(d, exist_ok=True)
        pre = f"{d}/{i}_{m}_{cov}"
        leaf_fa = C.get_seq(PANMANUTILS, PANMAP, pan, leaf, pre + ".leaf.fa",
                            f"work/{sp}/genomes.fa", SP_BIN("samtools")) and pre + ".leaf.fa"
        desc = pre + ".desc.fa"
        C.mutate_genome(pre + ".leaf.fa", cfg["mut_rates"][m], config["indel_fraction"],
                        SEED + i, desc, ts_tv=config.get("ts_tv", 2.0))
        r1, r2 = C.sim_reads(SP_BIN("wgsim"), desc, cov, cfg["genome_size"], RL,
                             config["seq_error"], SEED + i, pre + "_1.fq", pre + "_2.fq")
        best, wall, rss = C.place(PANMAP, pan, r1, r2, input.idx, pre, MT,
                                  force_leaf=FORCE_LEAF, completeness_weight=CW)
        row = _score_row(sp, pan, par, best, leaf, desc, cov, "sim", m, wall, rss, pre)
        open(output[0], "w").write(row)


rule place_real:
    input: idx="work/{sp}/index.pmi", maps=class_map
    output: "work/{sp}/" + F2REAL + "/{ri}_{cov}.tsv"
    threads: MT
    resources: mem_mb=lambda wc: job_mem("place", wc.sp)
    run:
        sp, ri, cov = wildcards.sp, int(wildcards.ri), float(wildcards.cov)
        cfg, pan = SP[sp], SP[wildcards.sp]["panman"]
        node, run = real_samples(sp)[ri]
        par = C.newick_parents(pan, PANMANUTILS)
        d = os.path.dirname(output[0]); os.makedirs(d, exist_ok=True)
        pre = f"{d}/{ri}_{cov}"
        rr1, rr2 = C.ena_fastqs(run, f"work/{sp}/reads/{run}")
        if not rr1:                              # no paired reads on ENA -> no-call
            open(output[0], "w").write(f"{sp}\treal\t-1\t{cov}\t\t\t\t\n")
            return
        # Enrich for on-target reads (map to the sample's own genome) so `cov` is
        # depth of that genome, not of the mostly-off-target raw library.
        truth = C.get_seq(PANMANUTILS, PANMAP, pan, node, pre + ".truth.fa",
                          f"work/{sp}/genomes.fa", SP_BIN("samtools")) and pre + ".truth.fa"
        er1, er2 = C.ontarget_reads(SP_BIN("minimap2"), SP_BIN("samtools"), SP_BIN("seqtk"),
                                    pre + ".truth.fa", rr1, rr2, f"work/{sp}/reads/{run}", MT)
        rl_real = C.mean_read_len(er1)   # actual read length; real runs vary ~44-250 bp
        avail_cov = 2 * C.count_reads(er1) * rl_real / cfg["genome_size"]
        if avail_cov < config.get("gt_min_cov_frac", 0.80) * cov:   # (3) read support -> no-call
            open(output[0], "w").write(f"{sp}\treal\t-1\t{cov}\t\t\t\t\n")
            return
        r1, r2 = C.subsample(SP_BIN("seqtk"), er1, er2, cov, cfg["genome_size"], rl_real,
                             SEED, pre + "_1.fq", pre + "_2.fq")
        best, wall, rss = C.place(PANMAP, pan, r1, r2, input.idx, pre, MT,
                                  exclude_self=dup_group_members(sp, node), force_leaf=FORCE_LEAF,
                                  completeness_weight=CW)
        row = _score_row(sp, pan, par, best, node, pre + ".truth.fa", cov, "real",
                         -1, wall, rss, pre, exclude=node)
        open(output[0], "w").write(row)


rule fig2_table:
    input: fig2_sim_inputs() + fig2_real_inputs()
    output: f"results/{FIG2}.tsv"
    run:
        hdr = "species\tkind\tmut\tcoverage\tscore\trandom\twall_s\trss_mb\n"
        with open(output[0], "w") as o:
            o.write(hdr)
            for f in input:
                o.write(open(f).read())


rule plot_fig2:
    input: f"results/{FIG2}.tsv"
    output: f"results/{FIG2}.pdf", f"results/{FIG2}.png"
    run:
        labels = {sp: SP[sp]["label"] for sp in SPECIES}
        muts = {sp: SP[sp]["mut_rates"] for sp in SPECIES}
        C.sh(["python3", os.path.join(workflow.basedir, "scripts/plot_fig2.py"),
              input[0], output[0], json_dumps({"order": SPECIES, "labels": labels,
              "muts": muts, "cov": config["coverages"]})])


# ── FIGURE 2 (raw): edit-distance variant (indels weighted by base length) ─────
# On-demand target: snakemake results/figure2_raw.pdf . Re-derives the placed node
# from each dump-all-scores and re-scores with genome_distance_raw; no re-placement.
rule fig2_raw_table:
    input: fig2_sim_inputs() + fig2_real_inputs() + \
           [ancient(_CLASS_MAP[sp]) for sp in SPECIES if sp in _CLASS_MAP]
    output: "results/figure2_raw.tsv"
    run:
        C.sh(["python3", os.path.join(workflow.basedir, "scripts/rescore_raw.py")] + list(SPECIES))


rule plot_fig2_raw:
    input: "results/figure2_raw.tsv"
    output: "results/figure2_raw.pdf", "results/figure2_raw.png"
    run:
        labels = {sp: SP[sp]["label"] for sp in SPECIES}
        muts = {sp: SP[sp]["mut_rates"] for sp in SPECIES}
        C.sh(["python3", os.path.join(workflow.basedir, "scripts/plot_fig2.py"),
              input[0], output[0], json_dumps({"order": SPECIES, "labels": labels,
              "muts": muts, "cov": config["coverages"]}), "symlog"])


# ── FIGURE 3 : consensus assembly (real data, leave-one-out) ──────────────────
# Ground truth = the leaf's GenBank sequence (dumped from the panman).
rule assemble:
    input: idx="work/{sp}/index.pmi"
    output: "work/{sp}/" + F3 + "/{ri}_{cov}.tsv"
    threads: MT
    resources: mem_mb=lambda wc: job_mem("assemble", wc.sp)   # TB panman index ~20GB; SARS/RSV JVMs a few GB
    run:
        sp, ri, cov = wildcards.sp, int(wildcards.ri), float(wildcards.cov)
        cfg, pan = SP[sp], SP[wildcards.sp]["panman"]
        node, run = real_samples(sp)[ri]
        d = os.path.dirname(output[0]); os.makedirs(d, exist_ok=True)
        pre = f"{d}/{ri}_{cov}"
        truth = C.get_seq(PANMANUTILS, PANMAP, pan, node, pre + ".truth.fa",
                          f"work/{sp}/genomes.fa", SP_BIN("samtools")) and pre + ".truth.fa"
        rr1, rr2 = C.ena_fastqs(run, f"work/{sp}/reads/{run}")
        if not rr1:                              # no paired reads on ENA -> no-call
            with open(output[0], "w") as o:
                for m in ("panmap", "standard", cfg["baseline"], "panmap_" + cfg["baseline"]):
                    o.write(f"{sp}\t{m}\t{cov}" + "\t" * 10 + "\n")
            return
        # Enrich for on-target reads (map to the sample's own genome) so `cov` is
        # depth of that genome, not of the mostly-off-target raw library.
        er1, er2 = C.ontarget_reads(SP_BIN("minimap2"), SP_BIN("samtools"), SP_BIN("seqtk"),
                                    truth, rr1, rr2, f"work/{sp}/reads/{run}", MT)
        rl_real = C.mean_read_len(er1)   # actual read length; real runs vary ~44-250 bp
        avail_cov = 2 * C.count_reads(er1) * rl_real / cfg["genome_size"]
        if avail_cov < config.get("gt_min_cov_frac", 0.80) * cov:   # (3) read support -> no-call
            with open(output[0], "w") as o:
                for m in ("panmap", "standard", cfg["baseline"], "panmap_" + cfg["baseline"]):
                    o.write(f"{sp}\t{m}\t{cov}" + "\t" * 10 + "\n")
            return
        r1, r2 = C.subsample(SP_BIN("seqtk"), er1, er2, cov, cfg["genome_size"], rl_real,
                             SEED, pre + "_1.fq", pre + "_2.fq")
        # panmap (leave-one-out): build the index per sample (counted in runtime),
        # place excluding self -> best node -> consensus. best is None when
        # placement produces no node (too few reads) -> no-call. No --min-read-support
        # flag is passed, so the build's adaptive default applies: singleton seeds are
        # filtered (min-read-support 2) only when estimated coverage > 3x, else all
        # seeds are kept (min-read-support 1) -- matters most at the low coverages here.
        pmi = pre + ".pmi"
        iw, irss = C.build_index(PANMAP, pan, *cfg["ksl"], pmi, MT)
        best, pw, prss = C.place(PANMAP, pan, r1, r2, pmi, pre + ".pl", MT,
                                 exclude_self=dup_group_members(sp, node), force_leaf=FORCE_LEAF,
                                 completeness_weight=CW)
        # Score a consensus against the truth genome -> both accuracy conventions
        # (base-based headline + legacy event-based) plus raw base/event counts.
        scr = lambda fa: (C.assembly_scores(SP_BIN("minimap2"), fa, pre + ".truth.fa",
                                            config["mask_bp"]) if fa else None)
        arss = pbrss = float("nan")
        if best:
            cons, aw, arss = C.assemble(PANMAP, pan, r1, r2, pmi, best, pre + ".pm", MT,
                                        mutation_spectrum=MM)
            s_pm = scr(cons)
            # panmap selects the reference (its placed leaf), then genotype that
            # reference with BWA+iVar (no impute): isolates reference selection from
            # the genotyping method. This arm feeds figure_S, not the main Fig 3.
            selref = C.get_seq(PANMANUTILS, PANMAP, pan, best, pre + ".selref.fa",
                               f"work/{sp}/genomes.fa", SP_BIN("samtools")) and pre + ".selref.fa"
            cons_pb, pbw, pbrss = C.bwa_ivar(pre + ".selref.fa", r1, r2, pre + ".pbwa",
                                             MT, SP_BIN("minimap2"),
                                             bindir=os.path.dirname(SP_BIN("samtools")))
            s_pb = scr(cons_pb)
        else:
            cons, aw, s_pm = None, 0.0, None
            cons_pb, pbw, s_pb = None, 0.0, None
        # Row 2 baseline: single standard reference + BWA+iVar (no impute). A
        # baseline that produces no consensus is 0% reconstructed (ZERO_SCORES),
        # a real coverage-limited data point rather than a dropped no-call.
        cons_b, bw, ref_b, subtype, brss = _baseline(sp, cfg, r1, r2, pre)
        s_b = scr(cons_b) if cons_b else C.ZERO_SCORES
        # Row 1 standard pipeline: the field-standard reference-based assembly for
        # the species (HaphPIPE for RSV, NCBI SC2VC for SARS, Clockwork for TB),
        # run exactly as published, then bcftools consensus.
        cons_std, stdw, std_sub, stdrss = _standard_pipeline(sp, cfg, r1, r2, pre)
        s_std = scr(cons_std)
        pk = lambda *xs: max([x for x in xs if x == x], default="")   # nan-safe peak MB
        # 8 score fields (accuracy=base, accuracy_event, snps, indel_events,
        # del_bases, ins_bases, aligned, interior); None -> all blank (no-call).
        def fmt(s):
            if s is None:
                return "\t".join([""] * 8)
            return "\t".join(str(x) for x in s)
        with open(output[0], "w") as o:
            o.write(f"{sp}\tpanmap\t{cov}\t{fmt(s_pm)}\t{iw + pw + aw}\t{pk(irss, prss, arss)}\n")
            o.write(f"{sp}\tstandard{('_' + std_sub) if std_sub else ''}\t{cov}\t{fmt(s_std)}\t{stdw}\t{pk(stdrss)}\n")
            o.write(f"{sp}\t{cfg['baseline']}\t{cov}\t{fmt(s_b)}\t{bw}\t{pk(brss)}\n")
            o.write(f"{sp}\tpanmap_{cfg['baseline']}\t{cov}\t{fmt(s_pb)}\t{iw + pw + pbw}\t{pk(irss, prss, pbrss)}\n")


rule fig3_table:
    input: fig3_inputs()
    output: f"results/{FIG3}.tsv"
    run:
        with open(output[0], "w") as o:
            # accuracy = base-based (% of the genome correctly reconstructed);
            # accuracy_event = legacy event-based metric; the count columns
            # (snps, indel_events, del_bases, ins_bases, aligned, interior) back
            # the supplementary event- vs base-based error table.
            o.write("species\tmethod\tcoverage\taccuracy\taccuracy_event\t"
                    "snps\tindel_events\tdel_bases\tins_bases\taligned\tinterior\t"
                    "wall_s\tpeak_mb\n")
            for f in input:
                o.write(open(f).read())


rule plot_fig3:
    input: f"results/{FIG3}.tsv"
    output: f"results/{FIG3}.pdf", f"results/{FIG3}.png"
    run:
        meta = {"order": SPECIES, "labels": {sp: SP[sp]["label"] for sp in SPECIES},
                "baseline": {sp: SP[sp]["baseline"] for sp in SPECIES},
                "cov": config["coverages"]}
        C.sh(["python3", os.path.join(workflow.basedir, "scripts/plot_fig3.py"),
              input[0], output[0], json_dumps(meta), "violin"])


# Re-score the cached Fig-3 consensus assemblies under the base-based counts
# (common.assembly_scores; minimap2 asm20 cs). Cheap: re-aligns already-written
# consensuses, no re-assembly, so it also propagates the common._cs_counts
# union-of-covered-intervals fix (overlapping alignment blocks no longer double-
# count, so genome fraction can never exceed 100%). Depends on the assemble outputs
# so it runs only once every held-out consensus exists.
rule rescore_fig3:
    # submap: RSV A/B competitive-mapping classification, so rescore_fig3 can tag each
    # RSV consensus with its subtype (joined by the truth.fa node id) for the revised
    # figure's HaphPIPE A/B split. Built first when RSV is in scope.
    input: res=fig3_inputs(),
           submap=(["meta/rsv_subtype.tsv"] if "rsv" in SPECIES else [])
    output: "results/figure3_rescored.tsv"
    run:
        cmd = ["python3", os.path.join(workflow.basedir, "scripts/rescore_fig3.py"),
               SP_BIN("minimap2"), str(config["mask_bp"])]
        if input.submap:
            cmd.append(input.submap[0])
        C.sh(cmd)


# Revised Fig 3 (final): the single blended accuracy row is split into the two axes
# QUAST/dnadiff report separately, genome fraction (completeness) and base error rate
# (per-base correctness incl. indel bases), over the panmap vs field-standard arms,
# plus runtime/memory. Reads the rescored base counts for rows 1-2 and the published
# table for runtime/memory (rows 3-4). ROW2METRIC=rate selects the base error rate for
# row 2 (vs the script default of consensus identity).
rule plot_fig3_revised:
    input: res="results/figure3_rescored.tsv", pub=f"results/{FIG3}.tsv"
    output: f"results/{FIG3}_revised_rate.pdf", f"results/{FIG3}_revised_rate.png"
    run:
        C.sh(["python3", os.path.join(workflow.basedir, "scripts/plot_fig3_revised.py"),
              input.res, input.pub, output[0]],
             env={**os.environ, "ROW2METRIC": "rate"})


# ── helpers used inside run-blocks ────────────────────────────────────────────
def json_dumps(x):
    import json
    return json.dumps(x)


def SP_BIN(tool):
    """Resolve a helper binary (bwa/samtools/ivar/minimap2/seqtk/wgsim) on PATH."""
    from shutil import which
    p = which(tool)
    if not p:
        raise FileNotFoundError(f"{tool} not found on PATH")
    return p


def _score_row(sp, pan, par, placed, expected_leaf, sample_fa, cov, kind, m,
               wall, rss, pre, exclude=None):
    """Fig 2 row: parsimony(placed,sample) - parsimony(expected,sample), clamped."""
    mm = SP_BIN("minimap2")
    gen, sam = f"work/{sp}/genomes.fa", SP_BIN("samtools")
    expected = par.get(expected_leaf) if kind == "real" else expected_leaf
    placed_fa = C.get_seq(PANMANUTILS, PANMAP, pan, placed, pre + ".placed.fa", gen, sam) and pre + ".placed.fa" if placed else None
    exp_fa = C.get_seq(PANMANUTILS, PANMAP, pan, expected, pre + ".exp.fa", gen, sam) and pre + ".exp.fa" if expected else None
    ex = config["exclude_bp"]
    ed = C.genome_distance(mm, exp_fa, sample_fa, ex) if exp_fa else 0
    # placed_fa is None when placement produced no node (e.g. too few reads at low
    # coverage) -> emit an empty score ("no call"), never a sentinel distance.
    pd = C.genome_distance(mm, placed_fa, sample_fa, ex) if placed_fa else None
    score = max(0, pd - ed) if pd is not None else ""
    # random-placement baseline: one random node per placement (n_random=1). Each
    # draw is emitted as "dist:taxon" so Fig 2A can colour it (RSV subtype / MTBC
    # species, from meta/); taxon is empty for species without a class map (SARS).
    import random, hashlib
    # LOO exclusion for the random pool must drop every exact duplicate (byte-identical
    # twin) of the held-out node and of the true parent, not just the single labels, so
    # a random draw can't land on an identical copy (dist ~= 0) and deflate the random
    # baseline. Mirrors the dup-group exclusion used for panmap's candidate set.
    excl_rand = set()
    if exclude is not None:
        excl_rand |= set(dup_group_members(sp, exclude))
    if expected is not None:
        excl_rand |= set(dup_group_members(sp, expected))
    lv = [n for n in C.leaves(par) if n not in excl_rand]
    rng = random.Random(SEED + int(hashlib.md5(pre.encode()).hexdigest()[:8], 16))
    rand = []
    for rn in rng.sample(lv, min(config["n_random"], len(lv))):
        rf = C.get_seq(PANMANUTILS, PANMAP, pan, rn, pre + ".rnd.fa", gen, sam) and pre + ".rnd.fa"
        d = max(0, C.genome_distance(mm, rf, sample_fa, ex) - ed)
        rand.append(f"{d}:{_rand_class(sp, rn, expected_leaf)}")
    rstr = ";".join(rand)
    return f"{sp}\t{kind}\t{m}\t{cov}\t{score}\t{rstr}\t{wall}\t{rss}\n"


_RAND_CLASS = {}
def _rand_class(sp, drawn, sample):
    """Random-baseline taxon label. RSV: 'same'/'cross' subtype of the drawn node
    vs the placement's sample (competitive-mapping subtypes explain the bimodal
    random distribution). TB: MTBC species of the drawn node. Empty otherwise."""
    if sp not in _RAND_CLASS:
        f = {"rsv": "meta/rsv_subtype.tsv", "tb": "meta/tb_species.tsv"}.get(sp)
        m = {}
        if f and os.path.exists(f):
            for line in open(f).read().splitlines()[1:]:
                c = line.split("\t")
                if len(c) >= 2:
                    m[c[0]] = c[1]
        _RAND_CLASS[sp] = m
    m = _RAND_CLASS[sp]
    if sp == "rsv":
        ds, ss = m.get(drawn), m.get(sample)
        return "" if ds is None or ss is None else ("same" if ds == ss else "cross")
    return m.get(drawn, "")


def _baseline(sp, cfg, r1, r2, pre):
    """Fig 3 baseline consensus + runtime + reference used. RSV/SARS: BWA+iVar
    (no impute); TB: Clockwork. Returns (consensus_fa, wall_s, ref_fa)."""
    bd = os.path.dirname(SP_BIN("samtools"))
    if cfg["baseline"] == "bwa_ivar":
        ref, subtype = cfg.get("ref"), None
        if sp == "rsv":                                  # subtype-matched reference
            ref, subtype = _rsv_ref(cfg, r1, r2, pre)
        fa, w, rss = C.bwa_ivar(ref, r1, r2, pre + ".bwa", MT,
                                SP_BIN("minimap2"), bindir=bd)
        return fa, w, ref, subtype, rss
    else:                                                # TB: Clockwork (Docker)
        fa, w, rss = _clockwork(cfg, r1, r2, pre)
        return fa, w, cfg.get("ref"), None, rss


def _standard_pipeline(sp, cfg, r1, r2, pre):
    """Fig 3 Row 1 field-standard reference-based assembly, run exactly as
    published, then bcftools consensus. Returns (consensus_fa, wall, subtype);
    subtype is 'A'/'B' for RSV (subtype-matched HaphPIPE reference) else None."""
    if sp == "rsv":                                      # HaphPIPE, subtype-matched ref (A: NC_038235.1, B: LR699737.1)
        ref, subtype = _rsv_ref(cfg, r1, r2, pre)
        fa, w, rss = C.haphpipe_consensus(config["haphpipe_bin"], ref, r1, r2, pre,
                                          bcftools=SP_BIN("bcftools"))
        return fa, w, subtype, rss
    if sp == "sars":                                     # NCBI SC2VC Illumina (NC_045512.2)
        fa, w, rss = C.ncbi_sars_consensus(config["ncbi_sars_script"], config["ncbi_sars_ref"],
                                           config["ncbi_sars_hisat2"], r1, r2, pre,
                                           config["ncbi_custom_filter"])
        return fa, w, None, rss
    fa, w, rss = C.clockwork_consensus(cfg["clockwork_image"], cfg["clockwork_ref_dir"],  # TB: Clockwork (H37Rv)
                                       cfg["ref"], r1, r2, pre, bcftools=SP_BIN("bcftools"))
    return fa, w, None, rss


def _rsv_ref(cfg, r1, r2, pre):
    """Pick RSV-A vs RSV-B reference by which maps more reads. Returns
    (ref_fa, subtype) where subtype is 'A' (ref_a) or 'B' (ref_b)."""
    mm = SP_BIN("minimap2")
    best, bestn, subtype = cfg["ref_a"], -1, "A"
    for ref, st in ((cfg["ref_a"], "A"), (cfg["ref_b"], "B")):
        n = sum(1 for ln in C.sh([mm, "-x", "sr", ref, r1, r2]).stdout.splitlines()
                if "tp:A:P" in ln)
        if n > bestn:
            best, bestn, subtype = ref, n, st
    return best, subtype


def _clockwork(cfg, r1, r2, pre):
    return C.clockwork_consensus(cfg["clockwork_image"], cfg["clockwork_ref_dir"],
                                 cfg["ref"], r1, r2, pre, bcftools=SP_BIN("bcftools"))


# ── Supplementary figure: leaf-only vs internal-allowed placement (consensus) ──
# force_leaf only affects panmap's placement, so this isolates the panmap arm: for
# every real leave-one-out sample x coverage, build the consensus under both
# --force-leaf and internal-allowed placement (shared reads -> each sample paired) and
# score genome fraction + per-base error. The driver builds its own per-species index
# and parallelises across samples internally, so it never re-runs the mode-independent
# standard/baseline arms. Run: snakemake results/figureSX_forceleaf.pdf --cores N
rule figure_forceleaf:
    output: pdf="results/figureSX_forceleaf.pdf", png="results/figureSX_forceleaf.png",
            tsv="results/forceleaf_consensus.tsv"
    run:
        jobs = min(getattr(workflow, "cores", 8) or 8, 64)
        C.sh(["python3", os.path.join(workflow.basedir, "scripts/forceleaf_consensus.py"),
              PANMAP, output.tsv, ",".join(SPECIES), "0", str(jobs)])
        meta = {"order": SPECIES, "labels": {sp: SP[sp]["label"] for sp in SPECIES},
                "cov": config["coverages"]}
        C.sh(["python3", os.path.join(workflow.basedir, "scripts/plot_forceleaf.py"),
              output.tsv, output.pdf, json_dumps(meta)])
