#!/usr/bin/env python3
"""t2t_tiers.py FEAT.tsv.zst SEGDUPS.bed.gz LOWMAP_SEGDUPS.bed.gz: where T2T-CHM13 disagrees with CERTA's GRCh38
placement, by CERTA's tier and region. FEAT: per primary record of a --t2t --feature-tags run (t2t_feat.sh:
qname, mate, contig, pos, MAPQ, XT, XC). A disagreement: T2T's placement lifts elsewhere in GRCh38 or does not lift
(XC class 2 or 3) with T2T MAPQ >= 20; confident when the record's MAPQ is >= 20 too. Regions (GIAB v3.6, renamed
to the reference's contigs): segdups, else low mappability (alllowmapandsegdupregions), else the rest. If T2T
disagreements concentrate in a few tiers and regions, a --t2t run could map only those records on CHM13."""
import bisect
import gzip
import subprocess
import sys
from collections import defaultdict

feat_path, seg_path, low_path = sys.argv[1:4]


def load(path):
    iv = defaultdict(list)
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith(("#", "track")):
                continue
            c, s, e = line.split("\t", 3)[:3]
            iv[c].append((int(s), int(e)))
    out = {}
    for c, lst in iv.items():
        lst.sort()
        starts, ends = [], []
        for s, e in lst:
            if starts and s <= ends[-1]:
                ends[-1] = max(ends[-1], e)
            else:
                starts.append(s)
                ends.append(e)
        out[c] = (starts, ends)
    return out


def inside(iv, c, pos0):
    t = iv.get(c)
    if t is None:
        return False
    i = bisect.bisect_right(t[0], pos0) - 1
    return i >= 0 and pos0 < t[1][i]


seg, low = load(seg_path), load(low_path)
# (tier, region) -> [records, MAPQ >= 20, disagreements, confident disagreements]
n = defaultdict(lambda: [0, 0, 0, 0])
z = subprocess.Popen(["zstd", "-dc", feat_path], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
for line in z.stdout:
    q, mate, contig, pos, mq, xt, xc = line.rstrip("\n").split("\t")
    if contig == "*":
        continue
    p0 = int(pos) - 1
    region = "segdups" if inside(seg, contig, p0) else ("low map" if inside(low, contig, p0) else "rest")
    mq = int(mq)
    dis = False
    if xc != "-":
        f = xc.split(",")
        dis = f[0] in ("2", "3") and int(f[2]) >= 20
    v = n[(xt, region)]
    v[0] += 1
    v[1] += mq >= 20
    v[2] += dis
    v[3] += dis and mq >= 20
z.wait()
tot = [sum(v[i] for v in n.values()) for i in range(4)]
print(f"{tot[0]:,} records, {tot[1]:,} with MAPQ >= 20; T2T disagrees on {tot[2]:,} ({tot[3]:,} of them MAPQ >= 20)")
print(f"{'tier':4} {'region':8} {'records':>11} {'% rec':>6} {'disagree':>9} {'% dis':>6} {'conf dis':>9} {'% conf':>6}"
      f" {'rate conf':>9}")
for k, v in sorted(n.items(), key=lambda x: -x[1][3]):
    print(f"{k[0]:4} {k[1]:8} {v[0]:11,} {100 * v[0] / tot[0]:6.2f} {v[2]:9,} {100 * v[2] / max(tot[2], 1):6.2f}"
          f" {v[3]:9,} {100 * v[3] / max(tot[3], 1):6.2f} {100 * v[3] / max(v[1], 1):8.3f}%")
