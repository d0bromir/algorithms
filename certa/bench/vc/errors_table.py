#!/usr/bin/env python3
"""errors_table.py DIR EVAL NAME=LABEL...: hap.py PASS errors (FN, FP) per run for
SNVs and indels in the whole region, segdups and all difficult regions."""
import csv
import sys

d, ev = sys.argv[1:3]
runs = [a.split("=", 1) for a in sys.argv[3:]]
strata = [("*", "all"), ("segdups", "segdups"), ("alldifficult", "all difficult"), ("notinalldifficult", "not difficult")]


def load(name):
    rows = {}
    with open(f"{d}/{name}.{ev}.extended.csv") as f:
        for r in csv.DictReader(f):
            if r["Filter"] == "PASS" and r["Subtype"] == "*":
                rows[(r["Subset"], r["Type"])] = r
    return rows


data = {n: load(n) for n, _ in runs}
for sub, title in strata:
    if not any((sub, "SNP") in data[n] for n, _ in runs):
        continue
    print(f"\n{ev}, {title}:\n")
    print("| run | SNV FN | SNV FP | indel FN | indel FP | errors |")
    print("|---|---|---|---|---|---|")
    for n, label in runs:
        cells, tot = [], 0
        for t in ("SNP", "INDEL"):
            r = data[n].get((sub, t))
            if r is None:
                cells += ["", ""]
                continue
            fn, fp = int(r["TRUTH.FN"]), int(r["QUERY.FP"])
            tot += fn + fp
            cells += [f"{fn:,}", f"{fp:,}"]
        print(f"| {label} | " + " | ".join(cells) + f" | {tot:,} |")
