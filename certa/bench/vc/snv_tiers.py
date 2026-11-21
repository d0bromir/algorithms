#!/usr/bin/env python3
"""snv_tiers.py SAMTOOLS REF SITES.tsv CERTA.cram MINIBWA.cram: at CERTA-only false SNV calls, the reads that carry
the false allele (the called alt base; the ref base where the truth is homalt) with MAPQ >= 5 in CERTA
(DeepVariant's default floor), by CERTA's tier and MAPQ band, against what minibwa does with the same read: at the
site with MAPQ >= 5, at the site below 5, elsewhere in the targets, or absent (outside the targets or unmapped)."""
import re
import subprocess
import sys
from collections import Counter, defaultdict

st, ref, sites_path, certa, minibwa = sys.argv[1:6]
OP = re.compile(r"(\d+)([MIDNSHP=X])")


def view(cram, region=None):
    cmd = [st, "view", "-T", ref, "-F", "0x904", cram] + ([region] if region else [])
    return subprocess.run(cmd, capture_output=True, text=True, check=True).stdout.splitlines()


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


def tier(rec):
    m = re.search(r"\tXT:Z:(\w+)", rec)
    return m.group(1) if m else "-"


sites = [l.rstrip("\n").split("\t") for l in open(sites_path) if l.strip()]
cases = []  # (site, name, mate, CERTA tier, CERTA MAPQ, minibwa at site: MAPQ or None)
for chrom, pos, kind, rb, ab, tgt, qgt in sites:
    pos = int(pos)
    excess = rb if (kind == "geno" and tgt == "homalt") else ab.split(",")[0]
    mb = {}
    for rec in view(minibwa, f"{chrom}:{pos}-{pos}"):
        c = rec.split("\t", 11)
        if base_at(c, pos) is not None:
            mb[(c[0], int(c[1]) & 192)] = int(c[4])
    for rec in view(certa, f"{chrom}:{pos}-{pos}"):
        c = rec.split("\t", 11)
        if base_at(c, pos) == excess and int(c[4]) >= 5:
            k = (c[0], int(c[1]) & 192)
            cases.append((f"{chrom}:{pos}", k[0], k[1], tier(rec), int(c[4]), mb.get(k)))
names = {x[1] for x in cases if x[5] is None}
where = {}
proc = subprocess.Popen([st, "view", "-T", ref, "-F", "0x904", minibwa], stdout=subprocess.PIPE, text=True)
for rec in proc.stdout:
    n = rec.split("\t", 1)[0]
    if n in names:
        c = rec.split("\t", 6)
        where[(n, int(c[1]) & 192)] = int(c[4])
proc.wait()


def band(q):
    return "5-19" if q < 20 else ("20-49" if q < 50 else "50-60")


table = Counter()
for site, n, mate, t, q, m in cases:
    if m is not None:
        mb_state = "site, MAPQ >= 5" if m >= 5 else "site, MAPQ < 5"
    elif (n, mate) in where:
        mb_state = "elsewhere"
    else:
        mb_state = "absent"
    table[(t, band(q), mb_state)] += 1
print(f"{len(sites)} sites, {len(cases)} false-allele reads with CERTA MAPQ >= 5")
by_t = defaultdict(int)
for (t, b, s), v in table.items():
    by_t[t] += v
for t in sorted(by_t, key=lambda x: -by_t[x]):
    print(f"\n{t}: {by_t[t]}")
    for (tt, b, s), v in sorted(table.items(), key=lambda x: -x[1]):
        if tt == t:
            print(f"   CERTA MAPQ {b:6}  minibwa {s:16} {v}")
