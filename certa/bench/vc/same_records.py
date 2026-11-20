#!/usr/bin/env python3
"""same_records.py A.sam B.sam: A and B hold the same records in the same
order except for MAPQ and the XQ tag (a MAPQ change in CERTA). Counts the
records whose MAPQ (and XQ) differ, by tier and direction; exits 1 at the
first record that differs otherwise."""
import re
import sys
from collections import Counter

XQ = re.compile(r"\tXQ:i:-?\d+")
TIER = re.compile(r"\tXT:Z:(\w+)")
n = 0
mapq, xq = Counter(), Counter()
with open(sys.argv[1]) as fa, open(sys.argv[2]) as fb:
    for la, lb in zip(fa, fb):
        n += 1
        a, b = la.split("\t", 5), lb.split("\t", 5)
        if a[:4] != b[:4] or XQ.sub("", a[5]) != XQ.sub("", b[5]):
            sys.exit(f"record {n} differs beyond MAPQ and XQ:\n{la}{lb}")
        t = TIER.search(la)
        t = t.group(1) if t else "-"
        qa, qb = int(a[4]), int(b[4])
        if qa != qb:
            mapq[(t, "up" if qb > qa else "down")] += 1
        if XQ.search(la) and XQ.search(la).group() != (XQ.search(lb).group() if XQ.search(lb) else ""):
            xq[t] += 1
    if fa.readline() or fb.readline():
        sys.exit(f"record counts differ after {n}")
print(f"same_records: {n:,} records; MAPQ differs: "
      + (", ".join(f"{t} {d} {c:,}" for (t, d), c in sorted(mapq.items())) or "none")
      + "; XQ differs: " + (", ".join(f"{t} {c:,}" for t, c in sorted(xq.items())) or "none"))
