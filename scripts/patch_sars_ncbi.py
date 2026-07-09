#!/usr/bin/env python3
"""Recompute the SARS Row-1 'standard' (NCBI SC2VC) arm for Fig 3, using the cached
subsampled reads/truth, and patch it into the per-sample fig3 tsvs. Needed because
work/sars/ncbi_ref (HISAT2 index) was wiped, leaving the arm empty. Reuses cached
work/sars/fig3/{ri}_{covf}_{1,2}.fq + {ri}_{covf}.truth.fa; no re-enrichment."""
import concurrent.futures as cf
import glob
import os
import sys

sys.path.insert(0, "scripts")
import common as C

NCBI_SCRIPT = "scripts/ncbi_sars_pipeline.py"
NCBI_REF = "work/sars/ncbi_ref/ref.fa"
NCBI_HISAT2 = "work/sars/ncbi_ref/ref"
CUSTOM = "scripts/ncbi_custom_vcf_filter.py"
GSIZE, MASK = 29903, 250
COVS = [("0.5", "0.5"), ("1", "1.0"), ("10", "10.0"), ("100", "100.0")]   # (tsv wildcard, float)
NSAMPLES = 57
WORKERS = int(sys.argv[1]) if len(sys.argv) > 1 else 8


def one(ri, cov_wc, cov_f):
    d = "work/sars/fig3"
    r1, r2 = f"{d}/{ri}_{cov_f}_1.fq", f"{d}/{ri}_{cov_f}_2.fq"
    truth = f"{d}/{ri}_{cov_f}.truth.fa"
    tsv = f"{d}/{ri}_{cov_wc}.tsv"
    if not os.path.exists(tsv):
        return (ri, cov_wc, "no-tsv")
    acc = ""
    if os.path.exists(r1) and os.path.exists(r2) and os.path.exists(truth) and os.path.getsize(r1):
        cons, _ = C.ncbi_sars_consensus(NCBI_SCRIPT, NCBI_REF, NCBI_HISAT2, r1, r2,
                                        f"{d}/{ri}_{cov_wc}.pncbi", CUSTOM)
        if cons:
            acc = C.genotyping_accuracy("minimap2", cons, truth, GSIZE, MASK)
    # patch the 'standard' line's accuracy (field 4) in place
    lines = open(tsv).read().splitlines()
    out = []
    for ln in lines:
        f = ln.split("\t")
        if len(f) >= 4 and f[1] == "standard":
            f[3] = str(acc)
            ln = "\t".join(f)
        out.append(ln)
    with open(tsv, "w") as o:
        o.write("\n".join(out) + "\n")
    return (ri, cov_wc, acc)


def main():
    tasks = [(ri, wc, f) for ri in range(NSAMPLES) for (wc, f) in COVS]
    done = 0
    with cf.ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futs = [ex.submit(one, *t) for t in tasks]
        for fu in cf.as_completed(futs):
            ri, wc, acc = fu.result()
            done += 1
            if done % 25 == 0 or acc in ("no-tsv", ""):
                print(f"  [{done}/{len(tasks)}] ri={ri} cov={wc}: acc={acc}", flush=True)
    print("patched SARS standard arm", flush=True)


if __name__ == "__main__":
    main()
