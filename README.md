# Panmap evaluation — Figures 2 & 3

Reproduces the manuscript's **Figure 2** (phylogenetic placement: accuracy +
runtime + memory) and **Figure 3** (consensus assembly accuracy + runtime) for
RSV, SARS-CoV-2, and M. tuberculosis.

`config.yaml` is loaded automatically. **Put targets/`--config` LAST** — snakemake ≥ 9
treats `--configfile`/`--config` as greedy, so `--configfile config.yaml results/…` swallows
the target. Just use:

```
snakemake -j8                                 # both figures -> results/
snakemake -j8 results/figure2.pdf             # Fig 2 only
snakemake -n results/figure2.pdf              # dry run / DAG check
```

### Scoped / smoke runs
Optional `--config` knobs (defaults reproduce the full manuscript run) — keep them last:
```
snakemake -j8 --config only_species=rsv sim_replicates=6 max_real=3
```
- `only_species=rsv[,sars,tb]` — restrict to a subset of species
- `sim_replicates=N`           — fewer simulated leaves per species (default 70)
- `max_real=N`                 — cap real SRA samples per species (0 = all)

## Requirements
- **snakemake ≥ 8** (older versions mangle f-strings; the Snakefile enforces this).
- `panmap` + `panmanUtils` (current build; path in `config.yaml`, default `../build/bin`).
- On `PATH`: `minimap2`, `bwa`, `samtools`, `ivar`, `seqtk`, `wgsim`; `docker` (TB Clockwork baseline only).
- Python: `numpy`, `matplotlib`, `dendropy`.

All of the above are provided by the `panmap-eval` conda env:
```
micromamba create -n panmap-eval -c conda-forge -c bioconda \
  python=3.11 "snakemake-minimal>=8" numpy matplotlib pandas dendropy \
  seqtk ivar wgsim minimap2 bwa samtools
micromamba activate panmap-eval
```
Build `panmap`/`panmanUtils` from the repo root first:
`cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DUSE_SYSTEM_LIBS=OFF && cmake --build build -j`.

## Data (fill in `config.yaml` paths)
Per species: a **current-format** panman (readable by the current panmap build;
old published panmans need a matching panmanUtils), a samples TSV, and the
single-reference baseline FASTA(s).

`samples.tsv` (real data, leave-one-out) — tab-separated, one row per sample:
```
node	run
USA/.../|OQ782918.1|2023-03-23	SRR24110076
```
`node` is the sample's leaf label in the panman (its GenBank genome = ground
truth); `run` is the SRA/ENA accession (paired Illumina).

## What it does
**Fig 2** — *Simulated*: for each of 70 random leaves × 2 mutation rates × 4
coverages, build a descendant genome (`wgsim`-style SNP/indel), simulate 150 bp
reads, place on the full panman, score = parsimony(placed) − parsimony(true
parent), clamped ≥ 0; plus a random-placement baseline. *Real*: leave-one-out
placement of subsampled SRA reads, scored against the parent of the held-out
leaf. Panels B/C report wall-time and peak RSS vs coverage.

**Fig 3** — leave-one-out assembly, three methods per sample:
- **Panmap** — panmap places (excluding the held-out leaf), then genotypes a
  consensus against its placed reference (native behavior).
- **Panmap→BWA+iVar** — panmap selects the reference (its placed leaf), but the
  consensus is called by BWA-MEM + iVar against that reference. Isolates the
  contribution of reference *selection* from the genotyping method.
- **BWA+iVar** (baseline) — single standard reference (subtype-matched for RSV,
  Wuhan-1 for SARS) + iVar; or **Clockwork** (TB, Docker).

Both BWA+iVar variants leave iVar's **N-gaps** at uncovered / low-depth positions
(no reference imputation — the coverage-limited baseline), so at low depth they
score most of the genome as unreconstructed, exposing panmap's advantage from
imputing across the pangenome. (`common.bwa_ivar(..., impute=True)` can fill N's
from the reference if a like-for-like comparison is wanted instead.) Accuracy =
% of the leaf's GenBank genome correctly reconstructed (150 bp ends excluded);
runtime is reported per method.

**On-target enrichment (real data).** Raw SRA runs are mostly off-target
(host/background — RSV can be <1% of a library), so subsampling from the whole
library makes the coverage axis meaningless and starves the read-support-based
baseline. Before subsampling, both figures map the raw reads to the *sample's own
GenBank assembly* and keep only the pairs that map (`common.ontarget_reads`,
mapped once per sample and cached under `work/<sp>/reads/<run>/`). "0.5×–100×"
therefore means depth of that genome. Both panmap and the baseline receive the
identical enriched read set; the baseline still assembles against its own
standard reference, so there is no ground-truth leakage into the reconstruction.

Leave-one-out uses the current panmap build (no `--remove-node`): the held-out
leaf is dropped from the candidate set after scoring and never used as the
reference — equivalent to pruning the leaf. The placed node is chosen by
`logContainment` — exactly the metric panmap's own `[ok] place … LogC` decision
uses (added as a column to `--dump-all-scores`), so it matches panmap's real
placement node-for-node.

**Runtimes** (Fig 2 B/C, Fig 3 bottom) are measured per job against a prebuilt
index. They are only meaningful with no CPU contention — run the timing pass
with `-j1` (accuracy panels are deterministic and unaffected by `-j`). Index
build time is recorded separately in `work/<sp>/index.time`.

## Layout
```
config.yaml          species data paths + parameters
Snakefile            rules: index -> place/assemble -> tables -> plots
scripts/common.py    panmap commands, scoring, BWA+iVar, sim, ENA reads
scripts/plot_fig2.py  scripts/plot_fig3.py
results/figure{2,3}.{tsv,pdf,png}
```
