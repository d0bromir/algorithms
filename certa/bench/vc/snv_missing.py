#!/usr/bin/env python3
"""snv_missing.py SAMTOOLS REF SITES.tsv CERTA.cram MINIBWA.cram: at CERTA-only false SNV calls, the reads minibwa
has at the site with MAPQ >= 5 that CERTA's pileup lacks: where CERTA puts them (elsewhere in the targets, by tier
and distance, or absent: outside the targets or unmapped), and the sites' depth in both."""
import re
import subprocess
import sys
from collections import Counter

st, ref, sites_path, certa, minibwa = sys.argv[1:6]
OP = re.compile(r"(\d+)([MIDNSHP=X])")


def view(cram, region=None):
    cmd = [st, "view", "-T", ref, "-F", "0x904", cram] + ([region] if region else [])
    return subprocess.run(cmd, capture_output=True, text=True, check=True).stdout.splitlines()


def covers(c, pos):
    r = int(c[3])
    for ln, op in OP.findall(c[5]):
        ln = int(ln)
        if op in "M=XDN":
            if r <= pos < r + ln:
                return op in "M=X"
            r += ln
    return False


sites = [l.rstrip("\n").split("\t") for l in open(sites_path) if l.strip()]
missing, depth = [], []
for chrom, pos, *_ in sites:
    pos = int(pos)
    cs = set()
    for rec in view(certa, f"{chrom}:{pos}-{pos}"):
        c = rec.split("\t", 11)
        if covers(c, pos) and int(c[4]) >= 5:
            cs.add((c[0], int(c[1]) & 192))
    nm = 0
    for rec in view(minibwa, f"{chrom}:{pos}-{pos}"):
        c = rec.split("\t", 11)
        if covers(c, pos) and int(c[4]) >= 5:
            nm += 1
            k = (c[0], int(c[1]) & 192)
            if k not in cs:
                missing.append((chrom, pos, k))
    depth.append((len(cs), nm))
names = {k[0] for _, _, k in missing}
where = {}
proc = subprocess.Popen([st, "view", "-T", ref, "-F", "0x904", certa], stdout=subprocess.PIPE, text=True)
for rec in proc.stdout:
    n = rec.split("\t", 1)[0]
    if n in names:
        c = rec.split("\t", 11)
        t = re.search(r"\tXT:Z:(\w+)", rec)
        where[(n, int(c[1]) & 192)] = (c[2], int(c[3]), int(c[4]), t.group(1) if t else "-")
proc.wait()
tab = Counter()
for chrom, pos, k in missing:
    w = where.get(k)
    if w is None:
        tab[("absent from CERTA's targets", "-", "-")] += 1
        continue
    d = "other contig" if w[0] != chrom else ("< 1 kb" if abs(w[1] - pos) < 1000 else ">= 1 kb")
    tab[(d, w[3], "MAPQ >= 5" if w[2] >= 5 else "MAPQ < 5")] += 1
dc = sum(a for a, _ in depth)
dm = sum(b for _, b in depth)
print(f"{len(sites)} sites; reads with MAPQ >= 5 over them: CERTA {dc}, minibwa {dm}; minibwa's missing from CERTA: {len(missing)}")
print(f"sites where CERTA has < 70 % of minibwa's depth (MAPQ >= 5): {sum(a < 0.7 * b for a, b in depth)}")
for (d, t, q), v in tab.most_common(25):
    print(f"   {d:28} {t:4} {q:10} {v}")
