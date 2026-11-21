#!/usr/bin/env python3
"""Couple carrier screening on real couples, checked against their child (README section 7).

    couple_real_eval.py DIR --couple wes_couple.py --carrier-genes carrier_genes.tsv

DIR holds, per sample, SAMPLE.tsv (clinical report findings table, report >= 0.13) and SAMPLE.vcf (its annotated
VCF). The three GIAB families are HG002 (HG003 x HG004), HG005 (HG006 x HG007) and HG001 (NA12891 x NA12892).
For every pair of the six parents the couple tool is run; for the real couples each parent's carrier alleles are
looked up in the child's calls: a real germline carrier allele reaches the child with probability 1/2.
"""
import argparse
import collections
import csv
import itertools
import os
import subprocess
import sys

FAMILIES = {("HG003", "HG004"): "HG002", ("HG006", "HG007"): "HG005", ("NA12891", "NA12892"): "HG001"}
SEX = {"HG003": "male", "HG004": "female", "HG006": "male", "HG007": "female", "NA12891": "male", "NA12892": "female"}


def calls(vcf):
    out = {}
    with open(vcf) as f:
        for line in f:
            if line[0] == "#":
                continue
            t = line.split("\t", 10)
            gt = t[9].split(":")[0] if len(t) > 9 else "./."
            for i, alt in enumerate(t[4].split(","), start=1):
                n = sum(1 for x in gt.replace("|", "/").split("/") if x == str(i))
                if n:
                    out[(t[0], t[1], t[3], alt)] = n
    return out


def carrier_alleles(tsv_path, couple_tsv, partner):
    """(chrom, pos, ref, alt) of the partner's carrier alleles in the couple result."""
    names = set()
    for r in csv.DictReader(open(couple_tsv), delimiter="\t"):
        if partner in r["partner"].split("+"):
            names.add(r["gene"])
    keys = []
    for r in csv.DictReader(open(tsv_path), delimiter="\t"):
        if r.get("gene") in names and r.get("pos") and r.get("ref"):
            keys.append((r["gene"], (r["chrom"], r["pos"], r["ref"], r["alt"])))
    return keys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dir")
    ap.add_argument("--couple", required=True)
    ap.add_argument("--carrier-genes", required=True)
    a = ap.parse_args()
    parents = [p for pair in FAMILIES for p in pair]
    print("pair\trelation\tboth_carriers\tto_confirm\tx_linked_female\tone_carrier\ttreatable_both")
    for x, y in itertools.combinations(parents, 2):
        out = os.path.join(a.dir, f"couple_{x}_{y}")
        subprocess.run([sys.executable, a.couple, "--a", os.path.join(a.dir, f"{x}.tsv"), "--b",
                        os.path.join(a.dir, f"{y}.tsv"), "--a-sex", SEX[x], "--b-sex", SEX[y],
                        "--carrier-genes", a.carrier_genes, "--html", out + ".html",
                        "--tsv", out + ".tsv"], check=True, capture_output=True)
        rows = list(csv.DictReader(open(out + ".tsv"), delimiter="\t"))
        c = collections.Counter(r["category"] for r in rows)
        xlf = sum(1 for r in rows if r["category"] == "x_linked" and r["risk"] == "1 in 2 sons")
        tb = sum(1 for r in rows if r["category"] == "both_carriers" and r["treatable"].startswith("yes"))
        rel = "couple" if (x, y) in FAMILIES else "unrelated"
        print(f"{x} x {y}\t{rel}\t{c['both_carriers']}\t{c['both_carriers_confirm']}\t{xlf}\t{c['one_carrier']}\t{tb}"
              + ("\t" + "; ".join(f"{r['gene']} ({r['alleles']})" for r in rows if r["category"].startswith("both"))
                 if c["both_carriers"] or c["both_carriers_confirm"] else ""))
    print("\nTransmission to the child of each parent's carrier alleles (expected 50 %):")
    tot = collections.Counter()
    for (x, y), child in FAMILIES.items():
        ch = calls(os.path.join(a.dir, f"{child}.vcf"))
        out = os.path.join(a.dir, f"couple_{x}_{y}.tsv")
        for p, who in ((x, "A"), (y, "B")):
            keys = carrier_alleles(os.path.join(a.dir, f"{p}.tsv"), out, who)
            got = sum(1 for _g, k in keys if k in ch)
            tot["alleles"] += len(keys)
            tot["transmitted"] += got
            print(f"  {p} -> {child}: {got} of {len(keys)} carrier alleles transmitted"
                  + (f" ({100 * got / len(keys):.0f} %)" if keys else ""))
    if tot["alleles"]:
        n, k = tot["alleles"], tot["transmitted"]
        p = k / n
        se = (0.25 / n) ** 0.5
        print(f"  all: {k} of {n} ({100 * p:.0f} %; 95 % range under 50 %: {100 * (0.5 - 1.96 * se):.0f}-"
              f"{100 * (0.5 + 1.96 * se):.0f} %)")


if __name__ == "__main__":
    main()
