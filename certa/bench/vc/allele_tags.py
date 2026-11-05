#!/usr/bin/env python3
"""allele_tags.py SITES < SAM records over the sites: per allele (ALT/REF at
the site), the tier (XT), XD presence and MAPQ classes of the reads."""
import re, sys
from collections import Counter, defaultdict
sites = {}
for line in open(sys.argv[1]):
    c, p, r, a = line.split()[:4]
    sites.setdefault(c, []).append((int(p), r.upper(), a.upper()))
cig = re.compile(r"(\d+)([MIDNSHP=X])")
stat = defaultdict(Counter)
seen = set()
for line in sys.stdin:
    f = line.rstrip("\n").split("\t")
    if int(f[1]) & 0x904:
        continue
    c, pos, mq, cg, seq = f[2], int(f[3]), int(f[4]), f[5], f[9]
    tags = dict((t[:2], t[5:]) for t in f[11:])
    for p, r, a in sites.get(c, []):
        if (f[0], f[1], p) in seen:
            continue
        rp, qp, base = pos, 0, None
        for n, op in cig.findall(cg):
            n = int(n)
            if op in "M=X":
                if rp <= p < rp + n:
                    base = seq[qp + p - rp]
                    break
                rp += n; qp += n
            elif op in "IS":
                qp += n
            elif op in "DN":
                if rp <= p < rp + n:
                    break
                rp += n
        if base is None:
            continue
        al = "alt" if base == a else "ref" if base == r else None
        if al is None:
            continue
        seen.add((f[0], f[1], p))
        s = stat[al]
        s["n"] += 1
        s["tier " + tags.get("XT", "?")] += 1
        s["XD set" if tags.get("XD", "-1") != "-1" else "XD none"] += 1
        s["MAPQ >= 40" if mq >= 40 else "MAPQ 20-39" if mq >= 20 else "MAPQ < 20"] += 1
for al in ("ref", "alt"):
    s = stat[al]
    n = s.pop("n")
    print(al, f"{n:,} reads:", ", ".join(f"{k} {100 * v / n:.1f} %" for k, v in sorted(s.items())))
