#!/bin/bash
# ============================================================================
# Panmap evaluation (Figures 2 & 3) — single-node, contention-free full run.
#
#   sbatch -p <partition> -A <account> slurm_run.sbatch            # everything
#   sbatch ... slurm_run.sbatch results/figure2.pdf               # a subtarget
#
# Design goals:
#   * Reproduce the manuscript figures from scratch after the intermediates were
#     cleaned (work/ is rebuilt; the ontarget read cache is reused if present).
#   * MAXIMUM parallelism across experiments, but every *timed* job must measure
#     its wall-clock / peak-RSS with NO influence from co-running jobs.
#
# How the "fast but contention-free" guarantee is met:
#   1. Every benchmarked tool is single-threaded (config method_threads=1): panmap
#      index/place/assemble, BWA+iVar, HaphPIPE, the NCBI pipeline, Clockwork, and
#      even the (untimed) ontarget enrichment. So each Snakemake job uses exactly
#      ONE core for its whole lifetime.
#   2. --exclusive + --hint=nomultithread => we own the whole node and are handed
#      PHYSICAL cores only (no SMT sibling sharing a core, which would skew timing).
#   3. snakemake --cores = (physical cores - headroom) => at most that many 1-core
#      jobs run at once => NO CPU oversubscription. Each concurrent job owns a core,
#      so its wall-clock is unaffected by the others.
#   4. --resources mem_mb = ~90% of node RAM, and each rule declares a per-job RSS
#      budget (Snakefile job_mem: TB ~34 GB, others 3-8 GB). Snakemake caps concurrency
#      so total residency never exceeds RAM => nothing swaps. Swapping would corrupt
#      BOTH wall-clock and peak-RSS, so this is what keeps the memory axis honest.
#      (Memory is per-process, so with no swap, peak-RSS is inherently contention-free.)
#   Net effect: small-genome (RSV/SARS) jobs run core-bound (hundreds at once); TB
#   runs memory-bound (~27 at once). No job ever shares a core or induces swap.
#
# Residual contention (shared L3 / memory bandwidth across distinct cores) is minor
# for these CPU-bound, multi-second-to-minutes jobs; if you need to eliminate it
# entirely, drop --cores to one-job-per-socket, or add taskset/numactl pinning.
# ============================================================================
#SBATCH --job-name=panmap-eval
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=0                       # all node RAM (paired with --resources mem_mb below)
#SBATCH --hint=nomultithread          # physical cores only -> no SMT co-runner timing skew
#SBATCH --time=48:00:00               # TB single-threaded (incl. Clockwork) is slow; be generous
#SBATCH --output=logs/panmap-eval.%j.out
#SBATCH --error=logs/panmap-eval.%j.out
# Set your scheduler specifics here or on the command line:
##SBATCH --partition=CHANGE_ME
##SBATCH --account=CHANGE_ME

set -euo pipefail
mkdir -p logs
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"

# --- environment: must expose snakemake>=8 AND every eval tool on PATH ----------
# (minimap2, samtools, bwa, ivar, seqtk, wgsim, bcftools, bedtools, panmanUtils,
#  plus python w/ numpy/matplotlib/dendropy/pyyaml). Adjust to your site.
ENV_NAME="${PANMAP_EVAL_ENV:-panmap-eval}"
if command -v micromamba >/dev/null 2>&1; then
    eval "$(micromamba shell hook -s bash)"; micromamba activate "$ENV_NAME"
elif command -v conda >/dev/null 2>&1; then
    source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate "$ENV_NAME"
fi

# --- concurrency: physical cores (minus a little headroom) and ~90% of RAM -------
CORES_AVAIL="$(nproc)"                                   # = physical cores (nomultithread)
JOBS="${SNAKE_JOBS:-$(( CORES_AVAIL > 8 ? CORES_AVAIL - 4 : CORES_AVAIL ))}"
MEM_TOTAL_MB="$(free -m | awk '/^Mem:/{print $2}')"
MEM_MB="${SNAKE_MEM_MB:-$(( MEM_TOTAL_MB * 90 / 100 ))}"

# Keep temp on node-local NVMe (fast, and off the shared FS) so I/O can't contend.
export TMPDIR="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}"

echo "node=$(hostname)  physical_cores=$CORES_AVAIL  snakemake_jobs=$JOBS  mem_mb=$MEM_MB  TMPDIR=$TMPDIR"

# panmap is auto-built from latest main by the Snakefile (config: panmap_auto_build),
# which needs a C++ toolchain + git on PATH on top of the eval tools. Check now so a
# missing compiler fails here, not deep in the first job.
for t in git cmake c++ make; do
    command -v "$t" >/dev/null || { echo "missing build tool '$t' (needed to auto-build panmap; or set panmap_auto_build: false and point panmap: at a prebuilt binary)"; exit 1; }
done

# --rerun-incomplete: recover cleanly if a prior allocation was killed mid-job.
# --keep-going: one failed job (e.g. a Clockwork/Docker TB-standard sample) must not
#               abort the other ~2000; inspect gaps afterward.
snakemake \
    --cores "$JOBS" \
    --resources mem_mb="$MEM_MB" \
    --configfile config.yaml \
    --rerun-incomplete \
    --keep-going \
    --printshellcmds \
    "$@"

echo "snakemake exit=$?"
