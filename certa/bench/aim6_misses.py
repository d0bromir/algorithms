"""Aim 6, descriptive: what the misses look like in minimap2's alignments.

  aim6_misses.py run_tau16.tsv minimap2.paf

Per miss class (aim6_summary.py): split alignments, the largest gap and
the identity of matched bases in the primary path, MAPQ; and the contigs
the misses map to.
"""
import collections
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from aim6_summary import classify, load, q  # noqa: E402

rows = load(sys.argv[1])
cls = {r["read"]: ("fit" if r["fits"] else classify(r)) for r in rows}
info = {}
nlines = collections.Counter()
nprim = collections.Counter()
with open(sys.argv[2]) as f:
    for line in f:
        t = line.rstrip("\n").split("\t")
        name = t[0]
        nlines[name] += 1
        tags = dict((x[:2], x[5:]) for x in t[12:])
        if tags.get("tp") == "P":
            nprim[name] += 1
        if name in info or tags.get("tp") != "P":
            continue
        cg = tags.get("cg", "")
        ops = re.findall(r"(\d+)([=XID])", cg)
        maxi = max([int(n) for n, o in ops if o == "I"], default=0)
        maxd = max([int(n) for n, o in ops if o == "D"], default=0)
        eq = sum(int(n) for n, o in ops if o == "=")
        x = sum(int(n) for n, o in ops if o == "X")
        info[name] = (t[5], int(t[7]), maxi, maxd, eq / max(1, eq + x), int(t[11]))

print("per class: reads, with >1 primary-type line (split), median largest insertion / deletion, "
      "share with a gap >= 500, median identity of matched bases, median MAPQ")
for c in ("fit", "clipped", "high loss", "repeat-limited", "unmapped"):
    names = [n for n, k in cls.items() if k == c]
    have = [info[n] for n in names if n in info]
    split = sum(nprim[n] > 1 for n in names)
    if not have:
        print(f"  {c}: {len(names)} reads")
        continue
    big = sum(max(h[2], h[3]) >= 500 for h in have)
    print(f"  {c}: {len(names)} reads; split {split} ({100 * split / len(names):.1f} %); "
          f"largest ins {q([h[2] for h in have], 0.5)} / del {q([h[3] for h in have], 0.5)}; "
          f"gap >= 500: {big} ({100 * big / len(have):.1f} %); identity {q([h[4] for h in have], 0.5):.4f}; "
          f"MAPQ {q([h[5] for h in have], 0.5)}")
contigs = collections.Counter()
for n, k in cls.items():
    if k in ("clipped", "high loss", "repeat-limited") and n in info:
        contigs[info[n][0]] += 1
allc = collections.Counter(info[n][0] for n in info)
print("misses by contig (share of that contig's mapped reads), top 12:")
print("  " + ", ".join(f"{c} {k} ({100 * k / allc[c]:.1f} %)" for c, k in contigs.most_common(12)))
