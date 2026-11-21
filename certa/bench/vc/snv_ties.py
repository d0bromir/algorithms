#!/usr/bin/env python3
"""snv_ties.py SAMTOOLS REF SITES.tsv CERTA.cram: CERTA's false-allele reads (MAPQ >= 5) at CERTA-only false SNV
calls, by tier and XB, the number of equally best single-end loci (CERTA writes no XS for most tiers)."""
import re
import subprocess
import sys
from collections import Counter

st, ref, sites_path, certa = sys.argv[1:5]
OP = re.compile(r"(\d+)([MIDNSHP=X])")


def base_at(c, pos):
    r, q = int(c[3]), 0
    for ln, op in OP.findall(c[5]):
        ln = int(ln)
        if op in "M=X":
            if r <= pos < r + ln:
                return c[9][q + pos - r]
            r += ln
            q += ln
        elif op in "IS":
            q += ln
        elif op in "DN":
            if r <= pos < r + ln:
                return None
            r += ln
    return None


tab = Counter()
for chrom, pos, kind, rb, ab, tgt, qgt in (l.rstrip("\n").split("\t") for l in open(sites_path) if l.strip()):
    pos = int(pos)
    excess = rb if (kind == "geno" and tgt == "homalt") else ab.split(",")[0]
    out = subprocess.run([st, "view", "-T", ref, "-F", "0x904", certa, f"{chrom}:{pos}-{pos}"], capture_output=True,
                         text=True, check=True).stdout.splitlines()
    for rec in out:
        c = rec.split("\t", 11)
        if base_at(c, pos) != excess or int(c[4]) < 5:
            continue
        t = re.search(r"\tXT:Z:(\w+)", rec)
        a = re.search(r"\tAS:i:(-?\d+)", rec)
        x = re.search(r"\tXS:i:(-?\d+)", rec)
        tier = t.group(1) if t else "-"
        b = re.search(r"\tXB:i:(-?\d+)", rec)
        nb = int(b.group(1)) if b else -1
        tab[(tier, "no XB" if nb < 0 else ("one best locus" if nb <= 1 else f"{min(nb, 8)}{'+' if nb >= 8 else ''} equal best loci"))] += 1
for (t, k), v in sorted(tab.items(), key=lambda kv: (kv[0][0], -kv[1])):
    print(f"   {t:4} {k:22} {v}")
