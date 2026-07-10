# Running the panmap evaluation on SLURM

Reproduces Figures 2 (placement) and 3 (assembly) from scratch on a single
exclusive node, with **maximum parallelism but contention-free timing/memory**.

```bash
sbatch -p <partition> -A <account> slurm_run.sbatch          # full run (rule all)
sbatch -p <partition> -A <account> slurm_run.sbatch results/figure2.pdf   # a subtarget
```

The final figure is `results/figure3_revised_rate.{pdf,png}` (plus `figure2.pdf`).

## Scale

~2400 jobs: 1680 `place_sim` + 356 `place_real` + 356 `assemble` + tables/plots.
Every benchmarked tool runs **single-threaded**, so throughput comes from running
many jobs at once, not from per-job threads.

## Why the timings stay honest at full parallelism

The runtime/memory rows of Fig 3 compare panmap against the standard pipelines, so
no job's wall-clock or peak-RSS may be perturbed by its neighbours. Four things
guarantee that (see the header of `slurm_run.sbatch` and the `MT` / `job_mem`
comments in the `Snakefile`):

1. **1 thread per tool** (`config method_threads: 1`) — panmap index/place/assemble,
   BWA+iVar, HaphPIPE, the NCBI pipeline, Clockwork, and the (untimed) ontarget
   enrichment. Each Snakemake job therefore uses exactly one core for its lifetime.
2. **`--exclusive --hint=nomultithread`** — we own the node and are given *physical*
   cores only (no SMT sibling sharing a core and skewing timing).
3. **`--cores = physical_cores − headroom`** with each rule reserving 1 core → no CPU
   oversubscription; every concurrent job owns a core.
4. **`--resources mem_mb ≈ 90% RAM`** + per-job `mem_mb` budgets (`Snakefile: job_mem`,
   TB ≈ 34 GB, others 3–8 GB) → total residency never exceeds RAM, so nothing swaps
   (swap would corrupt both wall-clock and peak-RSS). Peak-RSS is per-process, so with
   no swap it is inherently contention-free.

Net: RSV/SARS jobs run core-bound (hundreds at once); TB runs memory-bound (~27 at
once). Residual shared-L3 / memory-bandwidth effects are minor for these CPU-bound,
seconds-to-minutes jobs; for zero residual contention, pin with `taskset`/`numactl`
or drop `--cores` to one job per socket.

## Pre-flight checklist (edit for your site)

These are **hardcoded to the dev host** and must resolve on the compute node, or the
affected arm silently no-calls (the run continues via `--keep-going`):

- [ ] **Conda env** exposing `snakemake>=8`, `minimap2 samtools bwa ivar seqtk wgsim
      bcftools bedtools panmanUtils`, and python `numpy matplotlib dendropy pyyaml`.
      Set `PANMAP_EVAL_ENV` (default `panmap-eval`).
- [ ] **`panmap` binary** — `config.yaml: panmap` (currently `../build/bin/panmap`).
- [ ] **Data root** `config.yaml: /scratch1/alex/panstop/data/...` must be visible
      cluster-wide (Ceph primary storage, not a login-only mount).
- [ ] **SARS standard pipeline** — `scripts/ncbi_sars_pipeline.py` hardcodes absolute
      paths: Trimmomatic/Picard in `.../envs/haphpipe/bin`, HISAT2 in `.../envs/telexpr/bin`,
      `GATK=/usr/local/bin/gatk`. Verify these exist on the node.
- [ ] **RSV standard pipeline** — `config.yaml: haphpipe_bin`.
- [ ] **TB standard pipeline (Clockwork) uses Docker** (`common.py: clockwork_consensus`
      → `docker run`). Most SLURM nodes have **no Docker daemon**. Options: (a) ensure
      Docker is available on the node; (b) port `clockwork_consensus` to Apptainer/
      Singularity; or (c) run without it — TB "standard" rows become no-calls and the
      other three arms still populate. This node currently *does* have `/usr/bin/docker`.

## Notes

- **TB ground truth = the 10 samples** in `work/tb/qc_pass.tsv` (the `real_samples`
  source). Stale 7-sample `work/tb/{fig3,fig2_real}` partials were removed so the
  rebuild aligns `ri` → sample with the current set. The `work/*/reads/` ontarget
  cache and indexes were kept.
- **Resume**: `--rerun-incomplete` in the sbatch recovers a killed allocation; just
  re-submit. Snakemake reuses every finished job.
- **Timing** rows only reflect this run once every job has re-run; a partial run mixes
  vintages. Inspect `logs/panmap-eval.<jobid>.out` for any `--keep-going` failures.
