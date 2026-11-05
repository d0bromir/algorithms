#!/usr/bin/env python3
"""vc_table.py DIR GENOME EVAL NAME=LABEL...: markdown tables of hap.py's
results (PASS calls) for each aligner run: errors (FN + FP), recall,
precision and F1 for SNVs and indels, in the whole evaluated region and in
strata, with the mapping cost (CPU seconds as user + sys, peak memory) of
the 30x run from DIR/NAME.time."""
import csv, sys

d, genome, ev = sys.argv[1:4]
runs = [a.split("=", 1) for a in sys.argv[4:]]
strata = [("*", "all"), ("segdups", "segmental duplications"), ("lowmap_or_segdup", "low mappability or segdups"),
          ("alldifficult", "all difficult"), ("notinalldifficult", "not difficult"), ("CMRG_genes", "CMRG genes"),
          ("chr20", "chr20")]


def load(name):
    rows = {}
    with open(f"{d}/{name}.{ev}.extended.csv") as f:
        for r in csv.DictReader(f):
            if r["Filter"] == "PASS" and r["Subtype"] == "*" and r["Subset"] in dict(strata):
                rows[(r["Subset"], r["Type"])] = r
    return rows


def cost(name):
    try:
        w, u, s, m = open(f"{d}/{name}.time").read().split()[:4]
        return f"{float(u) + float(s):,.0f} s, {int(m) / 1024 ** 2:.1f} GiB"
    except (OSError, ValueError):
        return "not measured"


data = {n: load(n) for n, _ in runs}
for sub, title in strata:
    if not any((sub, "SNP") in data[n] for n, _ in runs):
        continue
    print(f"\n{genome}, {ev}, {title}:\n")
    print("| aligner | SNV errors | SNV recall | SNV precision | indel errors | indel recall | indel precision | mapping CPU, memory |")
    print("|---|---|---|---|---|---|---|---|")
    for n, label in runs:
        cells = []
        for t in ("SNP", "INDEL"):
            r = data[n].get((sub, t))
            if r is None or not r["TRUTH.TOTAL"]:
                cells += ["", "", ""]
                continue
            err = int(r["TRUTH.FN"]) + int(r["QUERY.FP"])
            cells += [f"{err:,}", f"{100 * float(r['METRIC.Recall']):.3f} %", f"{100 * float(r['METRIC.Precision']):.3f} %"]
        print(f"| {label} | " + " | ".join(cells) + f" | {cost(n)} |")
