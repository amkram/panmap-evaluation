#!/usr/bin/env python3
"""Ground-truth QC precompute for the real-data assembly benchmark, applying the
published sample-selection criteria and writing work/<sp>/qc_pass.tsv (the runs
kept for Fig 3 real data). Criteria, per candidate sample:
  (1) Assembly completeness/ambiguity:
        viral (RSV, SARS): truth genome has <= 5 N characters
        M. tuberculosis  : complete assembly (>= 0.98*genome_size, 0 N) + paired Illumina
  (2) Estimated coverage depth >= 500x  (total raw SRA bases / genome_size)
  (3) Concordance: align the sample's reads to its assembly and discard if
        bcftools call reports > 5 variants (reads must derive from that genome)
Reuses the eval's ENA download cache. Concordance mapping uses a subsample
(<=100x) since the variant count (<=5) is depth-stable.
Usage: qc_precompute.py <sp> <samples_tsv> <genomes_fa> <genome_size> <viral 0|1> <reads_root> <out_tsv> <bindir>
"""
import csv
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C

sp, samples_tsv, genomes_fa, genome_size, viral, reads_root, out_tsv, bindir = sys.argv[1:9]
genome_size = int(genome_size)
viral = viral == "1"
# Thresholds (env-overridable). Defaults are the published-strict criteria; the
# PanMAN-leaf-as-truth leave-one-out design makes reads differ from the leaf
# consensus by ~0.1-0.2% (10-25 variants over a viral genome), so the concordance
# bound is relaxed to the value that curated the manuscript's real-data set.
MIN_COV = int(os.environ.get("QC_MIN_COV", 500))
MAX_N_VIRAL = int(os.environ.get("QC_MAX_N", 5))
MAX_VARIANTS = int(os.environ.get("QC_MAX_VAR", 5))
CONC_COV = 100                      # subsample depth for the concordance check
RL = 150
B = lambda x: os.path.join(bindir, x) if bindir else x


def seqtk_bases(fq):
    """(n_reads, total_bases) via `seqtk size`."""
    out = C.sh([B("seqtk"), "size", fq]).stdout.split()
    return (int(out[0]), int(out[1])) if len(out) >= 2 else (0, 0)


def concordant(truth_fa, r1, r2, pre):
    """bwa reads -> truth assembly, bcftools call; return number of variants."""
    C.sh([B("bwa"), "index", truth_fa])
    bam = pre + ".conc.bam"
    p = subprocess.Popen([B("bwa"), "mem", "-t", "4", truth_fa, r1, r2],
                         stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    with open(bam, "wb") as bf:
        subprocess.run([B("samtools"), "sort", "-o", "-", "-"], stdin=p.stdout,
                       stdout=bf, stderr=subprocess.DEVNULL)
    p.wait()
    subprocess.run([B("samtools"), "index", bam], stderr=subprocess.DEVNULL)
    mp = subprocess.Popen([B("bcftools"), "mpileup", "-f", truth_fa, bam],
                          stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    call = subprocess.run([B("bcftools"), "call", "-mv"], stdin=mp.stdout,
                          capture_output=True, text=True)
    mp.wait()
    return sum(1 for ln in call.stdout.splitlines() if ln and not ln.startswith("#"))


rows = [(r["node"], r["run"]) for r in csv.DictReader(open(samples_tsv), delimiter="\t")]
C.ensure_genomes_fa(os.environ.get("PANMANUTILS", "panmanUtils"),
                    os.environ["QC_PANMAN"], genomes_fa, B("samtools"))
kept = []
for node, run in rows:
    reason = "pass"
    seq = C.faidx_get(genomes_fa, node, B("samtools"))
    if not seq:
        reason = "no-tip"
    else:
        ncount = seq.upper().count("N")
        if viral and ncount > MAX_N_VIRAL:
            reason = f"N={ncount}>5"
        elif not viral and (ncount > 0 or len(seq) < 0.98 * genome_size):
            reason = f"incomplete(N={ncount},len={len(seq)})"
    if reason == "pass":
        d = os.path.join(reads_root, run)
        r1, r2 = C.ena_fastqs(run, d)
        if not r1:
            reason = "no-paired-reads"
        else:
            (n1, b1), (n2, b2) = seqtk_bases(r1), seqtk_bases(r2)
            cov = (b1 + b2) / genome_size
            if cov < MIN_COV:
                reason = f"cov={cov:.0f}<500"
            else:
                truth = os.path.join(d, "qc_truth.fa")
                with open(truth, "w") as f:
                    f.write(f">{node}\n{seq}\n")
                s1, s2 = C.subsample(B("seqtk"), r1, r2, CONC_COV, genome_size, RL,
                                     42, os.path.join(d, "qc_s1.fq"), os.path.join(d, "qc_s2.fq"))
                nvar = concordant(truth, s1, s2, os.path.join(d, "qc"))
                if nvar > MAX_VARIANTS:
                    reason = f"discordant({nvar}var)"
    print(f"  {sp} {run} {node}: {reason}", flush=True)
    if reason == "pass":
        kept.append((node, run))

os.makedirs(os.path.dirname(out_tsv) or ".", exist_ok=True)
with open(out_tsv, "w") as f:
    f.write("node\trun\n")
    for node, run in kept:
        f.write(f"{node}\t{run}\n")
print(f"{sp}: {len(kept)}/{len(rows)} pass -> {out_tsv}", flush=True)
