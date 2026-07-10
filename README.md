# Panmap evaluation (Figures 2 & 3)

Reproduces Figure 2 (phylogenetic placement) and Figure 3 (consensus assembly) —
accuracy, runtime, and memory — for RSV, SARS-CoV-2, and M. tuberculosis.

## Run

Build panmap first (below), then:

```bash
snakemake -j8                        # everything -> results/
snakemake -j8 results/figure2.pdf    # one target
snakemake -n                         # dry run
```

`config.yaml` loads automatically. Paper figures are `results/figure2.pdf` and
`results/figure3_revised_rate.pdf`. Scope a quick run with `--config` (put targets
before it — snakemake treats a target after `--config` as a config entry):

```bash
snakemake -j8 --config only_species=rsv sim_replicates=6 max_real=3
```

`only_species` (subset), `sim_replicates` (default 70), `max_real` (0 = all).
`test_mode=true` is a fast smoke run — 20% of RSV+SARS (no TB), every rule, all three
plots (~370 jobs).

## Building panmap

The workflow doesn't build panmap, so no timed job ever waits on a compile.
`scripts/build_panmap.sh` checks out the latest `main` into an isolated worktree
(`panmap_build_dir`, default `../../.panmap-build`) and builds panmap + panmanUtils
there, leaving your own checkout alone. It rebuilds only when main moves, and
refuses any ref missing the `logContainment` dump column (panmap#80) — the metric
`place()` scores by. `slurm_run.sbatch` runs it automatically before snakemake.

```bash
scripts/build_panmap.sh .. ../../.panmap-build main   # config panmap: points at the output
```

Needs `git`, `cmake`, a C++ compiler, and `make`. First build is slow, then
incremental. To use your own binary, point `config.yaml: panmap` at it and skip this.

## Requirements

snakemake ≥ 8, plus on PATH: `minimap2 bwa samtools ivar seqtk wgsim bcftools
bedtools` (and `docker` for TB's Clockwork). Python: numpy, matplotlib, dendropy,
pyyaml.

```bash
micromamba create -n panmap-eval -c conda-forge -c bioconda \
  python=3.11 "snakemake-minimal>=8" numpy matplotlib dendropy pyyaml \
  minimap2 bwa samtools ivar seqtk wgsim bcftools bedtools
```

## Data

Everything a run reads lives under `data/<sp>/`. The small authoritative pieces —
`samples.tsv` (the locked real-sample set) and `duplicate_groups.tsv` — are committed;
the large panmans, reference FASTAs, and TB's `clockwork_ref/` are git-ignored and
fetched separately. `samples.tsv` is one row per real sample — `node` is the leaf's
label in the panman (its GenBank genome is the ground truth), `run` a paired Illumina
SRA/ENA accession:

```
node	run
USA/.../|OQ782918.1|2023-03-23	SRR24110076
```

It's already deduplicated to one sample per byte-identical group and is the single
source of truth (no separate QC cache).

### Fetch the large inputs

Not in git — on a fresh checkout, pull them from the host that has them, into `data/`
(run from the repo root):

```bash
DEV=alex@silverbullet.ucsc.edu

# refs + TB panman/clockwork_ref (~90 MB) from the pantry -> data/<sp>/
rsync -av --files-from=- "$DEV:/scratch1/alex/panstop/data/" data/ <<'EOF'
rsv/NC_038235.1.fa
rsv/LR699737.1.fa
sars/wuhan1.fa
tb/tb_400.panman
tb/NC_000962.3.fa
tb/clockwork_ref
EOF
# the RSV + SARS panmans live in the panmap source tree, not the pantry
rsync -av "$DEV:/scratch1/alex/poopdoop/src/test/data/rsv_4K.panman" data/rsv/
rsync -av "$DEV:/scratch1/alex/poopdoop/examples/data/sars_20000_twilight_dipper.panman" data/sars/
```

Reads aren't transferred — the pipeline pulls them from ENA on demand and caches them
under `work/<sp>/reads/`.

## What it measures

**Figure 2 (placement).** Simulated: 70 random leaves × 2 mutation rates × 4
coverages, each mutated (HKY substitutions + indels), sequenced to 150 bp reads,
and placed on the full panman. Score is `parsimony(placed) − parsimony(true parent)`
clamped at 0, against a one-draw random-node baseline. Real: leave-one-out placement
of subsampled SRA reads, scored against the held-out leaf's parent. Panels B/C are
wall time and peak RSS vs coverage.

**Figure 3 (assembly), leave-one-out.** Four arms per sample:

- **panmap** — place, then genotype a consensus against the placed reference.
- **standard** — the field-standard pipeline for the species (HaphPIPE for RSV, NCBI
  SC2VC for SARS, Clockwork for TB), run as published.
- **BWA+iVar** — one fixed reference (subtype-matched RSV, Wuhan-1 SARS, H37Rv TB) +
  iVar.
- **panmap→BWA+iVar** — panmap picks the reference, iVar genotypes it.

The revised figure plots panmap against the standard pipeline on genome fraction,
per-base error rate, runtime, and memory. The BWA+iVar arms leave iVar's N-gaps
unfilled, so at low depth they reconstruct little — the coverage-limited baseline
panmap's pangenome imputation is measured against. Accuracy is over the leaf's
GenBank genome, 150 bp of each end masked.

**On-target enrichment.** Raw SRA runs are mostly off-target (RSV can be under 1% of
a library), so before subsampling, reads are mapped to the sample's own GenBank
assembly and only mapped pairs are kept (cached under `work/<sp>/reads/<run>/`).
"0.5×–100×" is then depth of that genome, and every method gets the same reads. The
baseline still assembles against its own reference, so no ground truth leaks in.

**Leave-one-out** drops the held-out leaf and every byte-identical twin in its
duplicate group from the candidate set (no `--remove-node`; equivalent to pruning).

## Timing

Runtime and memory rows compare tools, so no job may perturb another's measurement.
Every tool runs single-threaded (`method_threads: 1`), so one Snakemake job uses one
core; per-job `mem_mb` budgets (`job_mem`) with `--resources mem_mb` cap concurrency
to physical RAM so nothing swaps. Small genomes run core-bound (hundreds at once), TB
memory-bound (~27). Accuracy is deterministic regardless of `-j`.

## Cluster

`slurm_run.sbatch` builds panmap and runs everything on one exclusive node
(`--exclusive --hint=nomultithread` for physical cores; `--cores` and `mem_mb` from
the node):

```bash
sbatch -p <partition> -A <account> slurm_run.sbatch           # full run
sbatch -p <partition> -A <account> slurm_run.sbatch --test     # 20% RSV+SARS smoke run
```

~2400 jobs (`--test` is ~370); `--rerun-incomplete` resumes a killed allocation. A few things are
host-specific and will no-call (not crash, via `--keep-going`) if missing on the node:

- the `config.yaml` data root must be visible cluster-wide;
- `scripts/ncbi_sars_pipeline.py` hardcodes Trimmomatic/Picard/HISAT2/GATK paths;
- RSV's HaphPIPE needs `config.yaml: haphpipe_bin`;
- TB's Clockwork runs under Docker, which most clusters lack — port it to Apptainer,
  or run without it and TB's standard row becomes no-calls.

## Layout

```
config.yaml        paths + parameters
Snakefile          index -> place/assemble -> tables -> plots
slurm_run.sbatch   builds panmap, then runs on one node
scripts/           common.py (shared lib), build_panmap.sh, plot_*, rescore_*,
                   classify_*, ncbi_sars_pipeline.py, screen_qc.py, ...
results/           figure2.pdf, figure3_revised_rate.pdf, *.tsv
```
