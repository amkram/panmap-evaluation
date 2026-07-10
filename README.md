# Panmap evaluation — Figures 2 & 3

Reproduces the manuscript's **Figure 2** (phylogenetic placement: accuracy, runtime,
memory) and **Figure 3** (consensus assembly accuracy + runtime) for RSV, SARS-CoV-2,
and M. tuberculosis.

## Run it

```bash
snakemake -j8                          # both figures -> results/
snakemake -j8 results/figure2.pdf      # Figure 2 only
snakemake -n results/figure2.pdf       # dry run / DAG check
```

`config.yaml` loads automatically. One snakemake ≥ 9 gotcha: `--config` and
`--configfile` are greedy and will eat a target that follows them, so put the target
**first** — `snakemake -n results/figure2.pdf --config only_species=rsv`, not the
other way around.

Scoped runs (defaults reproduce the full manuscript run):

```bash
snakemake -j8 results/figure2.pdf --config only_species=rsv sim_replicates=6 max_real=3
```

- `only_species=rsv[,sars,tb]` — restrict to a subset
- `sim_replicates=N` — simulated leaves per species (default 70)
- `max_real=N` — cap real SRA samples per species (0 = all)

## panmap is built for you

By default the eval builds `panmap` + `panmanUtils` itself from the latest **main** of
the panmap repo, so every run uses a current, correct binary. `scripts/build_panmap.sh`
checks out main into an isolated git worktree (`panmap_build_dir`, default
`../../.panmap-build`) — your own checkout and branch are never touched — and rebuilds
only when main has moved. The first build compiles the C++ dependencies and takes a
while; after that it's incremental.

You need `git`, `cmake`, a C++ compiler, and `make` on `PATH` in addition to the eval
tools. To use your own binary instead, set in `config.yaml`:

```yaml
panmap_auto_build: false
panmap: /path/to/panmap
```

The eval scores placement by the `logContainment` column of `--dump-all-scores` — the
exact metric panmap places by. That column landed in panmap#80; `build_panmap.sh`
refuses any ref that predates it, and `place()` errors rather than silently scoring a
different column.

## Requirements

- **snakemake ≥ 8** (older versions mangle f-strings; the Snakefile enforces this).
- On `PATH`: `minimap2 bwa samtools ivar seqtk wgsim bcftools bedtools`; `docker` for
  the TB Clockwork baseline only; plus the build tools above.
- Python: `numpy matplotlib dendropy pyyaml`.

The `panmap-eval` conda env provides the tool + Python stack:

```bash
micromamba create -n panmap-eval -c conda-forge -c bioconda \
  python=3.11 "snakemake-minimal>=8" numpy matplotlib pandas dendropy pyyaml \
  seqtk ivar wgsim minimap2 bwa samtools bcftools bedtools
```

## Data

Point `config.yaml` at, per species: a current-format panman (readable by the current
panmap build), a samples TSV, and the single-reference baseline FASTA(s). The samples
TSV is one row per real sample:

```
node	run
USA/.../|OQ782918.1|2023-03-23	SRR24110076
```

`node` is the sample's leaf label in the panman (its GenBank genome is the ground
truth); `run` is a paired-Illumina SRA/ENA accession.

## What it measures

**Figure 2 — placement.** *Simulated:* for each of 70 random leaves × 2 mutation rates
× 4 coverages, mutate the leaf (HKY substitutions + indels), simulate 150 bp reads,
place on the full panman, and score `parsimony(placed) − parsimony(true parent)`,
clamped ≥ 0; a random-placement baseline draws one random node per placement. *Real:*
leave-one-out placement of subsampled SRA reads, scored against the parent of the
held-out leaf. Panels B/C are wall-time and peak RSS vs coverage.

**Figure 3 — assembly**, leave-one-out, three methods per sample:

- **panmap** — places (excluding the held-out leaf), then genotypes a consensus against
  its placed reference (native behaviour).
- **panmap→BWA+iVar** — panmap picks the reference (its placed leaf), but BWA-MEM + iVar
  call the consensus. Isolates reference *selection* from the genotyping method.
- **BWA+iVar** (baseline) — one standard reference (subtype-matched for RSV, Wuhan-1 for
  SARS) + iVar; or **Clockwork** for TB.

Both BWA+iVar variants leave iVar's N-gaps at uncovered positions (no reference
imputation), so at low depth they score most of the genome as unreconstructed — that's
the coverage-limited baseline panmap's pangenome imputation is measured against.
Accuracy is the % of the leaf's GenBank genome correctly reconstructed (150 bp ends
excluded).

**On-target enrichment.** Raw SRA runs are mostly off-target (RSV can be <1% of a
library), so subsampling the whole library makes the coverage axis meaningless. Before
subsampling, both figures map the raw reads to the sample's *own* GenBank assembly and
keep only the pairs that map (cached under `work/<sp>/reads/<run>/`), so "0.5×–100×" is
depth of that genome. panmap and the baseline get the identical enriched reads; the
baseline still assembles against its own standard reference, so no ground truth leaks
into the reconstruction.

**Leave-one-out** uses the current panmap build (no `--remove-node`): the held-out leaf,
and every byte-identical twin in its duplicate group, are dropped from the candidate set
after scoring and never used as the reference — equivalent to pruning the leaf.

## Timing is contention-free at full parallelism

The runtime/memory rows compare panmap against the standard pipelines, so no job's
wall-clock or peak-RSS may be perturbed by a neighbour. Rather than serialize with
`-j1`, the pipeline runs many jobs at once but keeps each honest:

- **`method_threads: 1`** — every benchmarked tool (panmap, BWA+iVar, HaphPIPE, the NCBI
  pipeline, Clockwork) runs single-threaded, so each Snakemake job owns exactly one core.
- **per-job `mem_mb` budgets** (`Snakefile: job_mem`; TB ≈ 34 GB, others 3–8 GB) paired
  with `--resources mem_mb=<~node RAM>` cap concurrency to physical RAM, so nothing
  swaps — swapping would corrupt both the wall-clock and the peak-RSS.

Net: RSV/SARS run core-bound (hundreds at once), TB runs memory-bound (~27 at once), and
the accuracy panels are deterministic regardless of `-j`.

## On a cluster

`slurm_run.sbatch` runs the whole thing on one exclusive node with that contention-free
guarantee (`--exclusive --hint=nomultithread` for physical cores, `--cores` and
`mem_mb` set from the node):

```bash
sbatch -p <partition> -A <account> slurm_run.sbatch                     # rule all
sbatch -p <partition> -A <account> slurm_run.sbatch results/figure2.pdf # a subtarget
```

The full run is ~2400 jobs (1680 `place_sim` + 356 `place_real` + 356 `assemble` +
tables/plots); `--rerun-incomplete` resumes a killed allocation. A few things are
hardcoded to the dev host and will no-call (not crash, thanks to `--keep-going`) if they
don't resolve on the compute node — check before a long run:

- the data root in `config.yaml` must be visible cluster-wide;
- the SARS standard pipeline (`scripts/ncbi_sars_pipeline.py`) hardcodes absolute paths
  to Trimmomatic/Picard/HISAT2/GATK;
- the RSV standard pipeline needs `config.yaml: haphpipe_bin`;
- TB's Clockwork baseline runs via Docker, which most SLURM nodes lack — port it to
  Apptainer, or run without it and TB's "standard" row becomes no-calls.

## Layout

```
config.yaml            species data paths + parameters
Snakefile              index -> place/assemble -> tables -> plots (auto-builds panmap)
slurm_run.sbatch       single-node cluster runner
scripts/
  common.py            panmap commands, scoring, BWA+iVar, sim, ENA reads (shared lib)
  build_panmap.sh      build panmap from latest main into an isolated worktree
  plot_fig2.py         plot_fig3.py  plot_fig3_revised.py     figure rendering
  rescore_fig3.py  rescore_raw.py                             alternative score conventions
  classify_{rsv_subtype,tb_species}.py                        RSV A/B, MTBC species labels
  ncbi_sars_pipeline.py  ncbi_custom_vcf_filter.py            SARS standard pipeline
  screen_qc.py  qc_precompute.py  find_tb_illumina.py         real-sample search + QC
  plot_reviewer_*.py  reviewer_refsel.py  plot_fig3_quality.py  reviewer-response figures
results/figure{2,3}.{tsv,pdf,png}
```
