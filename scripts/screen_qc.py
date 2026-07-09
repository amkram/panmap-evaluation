#!/usr/bin/env python3
"""Broader real-sample search over PanMAN leaves + published ground-truth QC.
For each candidate leaf node (a GenBank accession) in a random order, until
`target` samples pass:
  0) resolve the node accession -> SRA run via NCBI (nuccore -> biosample -> sra)
  1) completeness/ambiguity: viral truth <= 5 N; M.tb complete (>=0.98*len, 0 N) + paired
  2) on-target depth >= 500x  (map reads to the assembly; estimate on-target depth)
  3) concordance: <= 5 variants from bcftools call of reads vs the assembly
Writes out_tsv (node<TAB>run) incrementally; resumable (skips nodes already logged).
Usage: screen_qc.py <sp> <genomes_fa> <genome_size> <viral 0|1> <nodes_file> <target> <reads_root> <out_tsv> <log_tsv> <edirect_bin> <bindir>
"""
import csv
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C

(sp, genomes_fa, genome_size, viral, nodes_file, target, reads_root,
 out_tsv, log_tsv, edirect_bin, bindir) = sys.argv[1:12]
genome_size = int(genome_size); viral = viral == "1"; target = int(target)
# QC thresholds (overridable via env for loosened runs)
MIN_COV = int(os.environ.get("QC_MIN_COV", 500))
MAX_N = int(os.environ.get("QC_MAX_N", 5))
MAX_VAR = int(os.environ.get("QC_MAX_VAR", 5))
RL, SUB = 150, 1_500_000
B = lambda x: os.path.join(bindir, x) if bindir else x
E = lambda x: os.path.join(edirect_bin, x) if edirect_bin else x


import urllib.request
import xml.etree.ElementTree as ET
import time
_EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def _get(url):
    for _ in range(4):
        try:
            with urllib.request.urlopen(url, timeout=45) as r:
                return r.read()
        except Exception:
            time.sleep(1.5)
    return b""


def _linked(dbfrom, db, uid):
    """Linked UIDs from elink (only <LinkSetDb> links, not the input IdList)."""
    try:
        root = ET.fromstring(_get(f"{_EU}/elink.fcgi?dbfrom={dbfrom}&db={db}&id={uid}"))
        return [e.text for ls in root.iter("LinkSetDb") for e in ls.iter("Id")]
    except Exception:
        return []


def resolve_run(acc):
    """GenBank accession -> paired-Illumina SRA run (nuccore->biosample->sra)."""
    try:
        uid = ET.fromstring(_get(f"{_EU}/esearch.fcgi?db=nuccore&term={acc}")).find(".//IdList/Id").text
    except Exception:
        return None
    time.sleep(0.34)
    for bs in _linked("nuccore", "biosample", uid):
        for sr in _linked("biosample", "sra", bs):
            info = _get(f"{_EU}/efetch.fcgi?db=sra&id={sr}&rettype=runinfo&retmode=text").decode(errors="ignore")
            for row in csv.DictReader(info.splitlines()):
                if row.get("Run", "").startswith(("SRR", "ERR", "DRR")) and \
                   row.get("LibraryLayout", "") == "PAIRED" and "ILLUMINA" in row.get("Platform", "").upper():
                    return row["Run"]
            time.sleep(0.34)
    return None


def map_qc(truth_fa, r1, r2, total_pairs, pre):
    """Map (subsampled) reads to the assembly; return (est_on_target_depth, n_variants)."""
    C.sh([B("bwa"), "index", truth_fa])
    s1, s2 = r1, r2
    if total_pairs > SUB:
        s1, s2 = pre + ".s1.fq", pre + ".s2.fq"
        for src, dst in ((r1, s1), (r2, s2)):
            with open(dst, "w") as f:
                subprocess.run([B("seqtk"), "sample", "-s42", src, str(SUB)],
                               stdout=f, stderr=subprocess.DEVNULL)
    sub_pairs = min(total_pairs, SUB)
    bam = pre + ".bam"
    p = subprocess.Popen([B("bwa"), "mem", "-t", "4", truth_fa, s1, s2],
                         stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    with open(bam, "wb") as bf:
        subprocess.run([B("samtools"), "sort", "-o", "-", "-"], stdin=p.stdout,
                       stdout=bf, stderr=subprocess.DEVNULL)
    p.wait(); subprocess.run([B("samtools"), "index", bam], stderr=subprocess.DEVNULL)
    fs = C.sh([B("samtools"), "flagstat", bam]).stdout
    mapped = next((int(l.split()[0]) for l in fs.splitlines() if " mapped (" in l), 0)
    frac = mapped / (2 * sub_pairs) if sub_pairs else 0
    est_depth = frac * total_pairs * 2 * RL / genome_size
    mp = subprocess.Popen([B("bcftools"), "mpileup", "-f", truth_fa, bam],
                          stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    call = subprocess.run([B("bcftools"), "call", "-mv"], stdin=mp.stdout,
                          capture_output=True, text=True); mp.wait()
    nvar = sum(1 for l in call.stdout.splitlines() if l and not l.startswith("#"))
    return est_depth, nvar


done = set()
if os.path.exists(log_tsv):
    for r in csv.DictReader(open(log_tsv), delimiter="\t"):
        done.add(r["node"])
kept = []
if os.path.exists(out_tsv):
    for r in csv.DictReader(open(out_tsv), delimiter="\t"):
        kept.append((r["node"], r["run"]))
nodes = [l.strip() for l in open(nodes_file) if l.strip()]
new_log = not os.path.exists(log_tsv) or os.path.getsize(log_tsv) == 0
logf = open(log_tsv, "a")
if new_log:
    logf.write("node\trun\tresult\n")

for node in nodes:
    if len(kept) >= target:
        break
    if node in done:
        continue
    acc = node.split("|")[1] if "|" in node else node          # SARS labels: name|ACC|date
    run, res = None, "pass"
    seq = C.faidx_get(genomes_fa, node, B("samtools"))
    if not seq:
        res = "no-tip"
    else:
        nc = seq.upper().count("N")
        if viral and nc > MAX_N:
            res = f"N={nc}"
        elif not viral and (nc > 0 or len(seq) < 0.98 * genome_size):
            res = f"incomplete"
    if res == "pass":
        run = resolve_run(acc)
        if not run:
            res = "no-sra"
    if res == "pass":
        d = os.path.join(reads_root, run); os.makedirs(d, exist_ok=True)
        r1, r2 = C.ena_fastqs(run, d)
        if not r1:
            res = "no-paired"
        else:
            total = int(C.sh([B("seqtk"), "size", r1]).stdout.split()[0] or 0)
            truth = os.path.join(d, "truth.fa")
            open(truth, "w").write(f">{node}\n{seq}\n")
            depth, nvar = map_qc(truth, r1, r2, total, os.path.join(d, "qc"))
            if depth < MIN_COV:
                res = f"depth={depth:.0f}"
            elif nvar > MAX_VAR:
                res = f"discordant={nvar}"
    logf.write(f"{node}\t{run or ''}\t{res}\n"); logf.flush()
    print(f"  [{len(kept)}/{target}] {node} {acc} run={run or '-'} -> {res}", flush=True)
    if res == "pass":
        kept.append((node, run))
        with open(out_tsv, "w") as f:
            f.write("node\trun\n")
            for n, rr in kept:
                f.write(f"{n}\t{rr}\n")

print(f"{sp}: {len(kept)}/{target} passed (screened up to {len(done)+1} nodes)", flush=True)
