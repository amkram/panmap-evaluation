#!/usr/bin/env python3
"""Build two plot-ready Fig-3 tables that differ ONLY in the accuracy column:
  results/figure3_event.tsv  -- legacy event-based accuracy (== published)
  results/figure3_base.tsv   -- base-based accuracy (matches/interior)
Runtime/memory columns are identical in both (the metric change doesn't touch
them). Base values come from figure3_rescored.tsv, matched onto each published
row by its exact acc_event value (verified to reproduce the published accuracy),
so the per-(species,method,coverage) distribution is preserved. RSV's subtype-
split standard_A/standard_B labels are kept from the published table."""
import csv
import os
import sys
from collections import defaultdict

ROOT = os.path.join(os.path.dirname(__file__), "..")
PUB = os.path.join(ROOT, "results", "figure3.tsv")
RES = os.path.join(ROOT, "results", "figure3_rescored.tsv")


def mapmeth(m):
    return "standard" if m.startswith("standard") else m


# rescore -> {(sp, mapped_method, cov): [(acc_event, acc_base), ...]}
pool = defaultdict(list)
for r in csv.DictReader(open(RES), delimiter="\t"):
    if r["accuracy_event"] in ("", "nan"):
        continue
    pool[(r["species"], r["method"], float(r["coverage"]))].append(
        (float(r["accuracy_event"]), float(r["accuracy"])))

pub_rows = list(csv.DictReader(open(PUB), delimiter="\t"))
HEADER = "species\tmethod\tcoverage\taccuracy\twall_s\tpeak_mb\n"
unmatched = 0
ev = open(os.path.join(ROOT, "results", "figure3_event.tsv"), "w")
ba = open(os.path.join(ROOT, "results", "figure3_base.tsv"), "w")
ev.write(HEADER); ba.write(HEADER)
for r in pub_rows:
    sp, m, cov = r["species"], r["method"], float(r["coverage"])
    base = r["accuracy"]
    if r["accuracy"] not in ("", "nan"):
        cand = pool.get((sp, mapmeth(m), cov), [])
        ae = float(r["accuracy"])
        j = min(range(len(cand)), key=lambda k: abs(cand[k][0] - ae), default=None)
        if j is not None and abs(cand[j][0] - ae) < 1e-3:
            base = f"{cand.pop(j)[1]:.6f}"
        else:
            unmatched += 1
            base = r["accuracy"]                       # fall back to event value
    row_ev = f"{sp}\t{m}\t{r['coverage']}\t{r['accuracy']}\t{r['wall_s']}\t{r['peak_mb']}\n"
    row_ba = f"{sp}\t{m}\t{r['coverage']}\t{base}\t{r['wall_s']}\t{r['peak_mb']}\n"
    ev.write(row_ev); ba.write(row_ba)
ev.close(); ba.close()
print(f"wrote results/figure3_event.tsv and results/figure3_base.tsv "
      f"({len(pub_rows)} rows, {unmatched} unmatched -> kept event value)",
      file=sys.stderr)
