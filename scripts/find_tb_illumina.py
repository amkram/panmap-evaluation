#!/usr/bin/env python3
"""Scan TB panman leaves for those backed by deep (>=100x on 4.4Mb, >=100bp) Illumina
SRA runs, so Fig 3 LOO has samples that pass the 100x on-target coverage gate.
leaf accession -> GenBank DBLINK BioSample -> ENA read_run (best Illumina run).
Writes a ranked TSV: node, biosample, run, platform, base_count, read_len, raw_depth_x
"""
import sys, time, re, urllib.request, urllib.parse
from concurrent.futures import ThreadPoolExecutor
GS = 4411532
def get(url):
    for _ in range(3):
        try:
            return urllib.request.urlopen(url, timeout=30).read().decode("utf-8", "replace")
        except Exception:
            time.sleep(1)
    return ""
accs = [l.strip() for l in open(sys.argv[1]) if l.strip()]
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
# 1) batched efetch: accession -> biosample (seq_stop=1 keeps it to the header)
acc2bs = {}
for i in range(0, len(accs), 20):
    batch = accs[i:i+20]
    txt = get(EUTILS + "?" + urllib.parse.urlencode(
        {"db": "nuccore", "id": ",".join(batch), "rettype": "gb", "retmode": "text", "seq_stop": "1"}))
    # split records on VERSION lines; map each record's accession.version to its BioSample
    for rec in re.split(r"\n(?=LOCUS )", txt):
        m = re.search(r"VERSION\s+(\S+)", rec)
        b = re.search(r"BioSample:\s*(SAM[END][A-Z]?\d+)", rec)
        if m and b:
            acc2bs[m.group(1)] = b.group(1)
    time.sleep(0.34)
    sys.stderr.write(f"efetch {min(i+20,len(accs))}/{len(accs)} (biosamples so far {len(acc2bs)})\n")
# 2) biosample -> best Illumina run via ENA (threaded; ENA tolerates concurrency)
def best_illumina(item):
    acc, bs = item
    txt = get("https://www.ebi.ac.uk/ena/portal/api/filereport?" + urllib.parse.urlencode(
        {"accession": bs, "result": "read_run",
         "fields": "run_accession,instrument_platform,library_strategy,base_count,read_count,library_layout",
         "format": "tsv"}))
    best = None
    for ln in txt.splitlines()[1:]:
        c = ln.split("\t")
        if len(c) < 6 or c[1] != "ILLUMINA":
            continue
        try:
            base, reads = int(c[3]), int(c[4])
        except ValueError:
            continue
        rl = base / reads if reads else 0
        if c[5] == "PAIRED" and rl >= 100 and (best is None or base > best[4]):
            best = (acc, bs, c[0], c[1], base, int(rl), base / GS)
    return best
rows = []
with ThreadPoolExecutor(max_workers=12) as ex:
    for r in ex.map(best_illumina, list(acc2bs.items())):
        if r:
            rows.append(r)
rows.sort(key=lambda r: -r[4])
with open("meta/tb_illumina_candidates.tsv", "w") as f:
    f.write("node\tbiosample\trun\tplatform\tbase_count\tread_len\traw_depth_x\n")
    for acc, bs, run, plat, base, rl, dep in rows:
        f.write(f"{acc}\t{bs}\t{run}\t{plat}\t{base}\t{rl}\t{dep:.0f}\n")
deep = [r for r in rows if r[6] >= 100]
sys.stderr.write(f"\nDONE: {len(acc2bs)}/{len(accs)} leaves linked to biosample; "
                 f"{len(rows)} have paired Illumina >=100bp; {len(deep)} are >=100x raw depth\n")
