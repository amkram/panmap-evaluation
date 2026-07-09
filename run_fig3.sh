#!/bin/bash
# Fig 3 (RSV + SARS) with published Row-1 pipelines (HaphPIPE / NCBI SC2VC) on
# the QC'd real ground-truth sets (22 RSV + 62 SARS). Needs the panmap-eval env
# bin on PATH for seqtk/bwa/minimap2/ivar/samtools (bcftools/bedtools -> /usr/bin).
set -o pipefail
cd /scratch1/alex/poopdoop/evaluation
export PATH=/home/alex/micromamba/envs/panmap-eval/bin:$PATH
exec snakemake results/figure3.pdf -j25 --keep-going --rerun-triggers mtime \
     --config only_species=rsv,sars
