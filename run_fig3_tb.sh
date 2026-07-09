#!/bin/bash
# Fig 3 with TB added (Clockwork Row-1). RSV+SARS tsvs are cached (mtime triggers),
# so only the missing TB jobs compute; then the 3-species figure regenerates.
# Higher -j: Clockwork is mostly single-threaded, so oversubscribe job slots to
# keep several Cortex/minos runs in flight (128 cores, ~750 GB RAM free).
set -o pipefail
cd /scratch1/alex/poopdoop/evaluation
export PATH=/home/alex/micromamba/envs/panmap-eval/bin:$PATH
exec snakemake results/figure3.pdf -j48 --keep-going --rerun-triggers mtime \
     --config only_species=rsv,sars,tb
