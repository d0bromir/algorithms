#!/usr/bin/env python3
"""A run's confident near misses (bench/VALIDATION_RESULTS.md, "Matched
sensitivity"): mates with MAPQ >= 40 at the truth within 150 bp but not
within 20, by CIGAR.

  samtools view X.bam | truth_nearmiss.py truth.npz X.w20.npz X.w150.npz

Each near miss is counted once: with a soft clip of >= 10 bases (offset
from the truth less the longest clip), else an indel of >= 10, else neither
(offset). The offset is that of truth_eval.py: unclipped leftmost base
against the truth's.
"""
import re
import sys

import numpy as np

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
LEAD = re.compile(r"^(\d+)[SH]")
tp = np.load(sys.argv[1])["tp"]
a, b = np.load(sys.argv[2]), np.load(sys.argv[3])
near = b["correct"] & ~a["correct"] & (b["mapq"] >= 40)
kinds = {"soft clip ≥ 10": [], "indel ≥ 10": [], "neither": []}
for line in sys.stdin:
    f = line.split("\t", 6)
    flag = int(f[1])
    if flag & (4 | 256 | 2048):
        continue
    i, m = int(f[0][1:]), 1 if flag & 128 else 0
    if not near[i, m]:
        continue
    ops = CIGAR.findall(f[5])
    lc = LEAD.match(f[5])
    d = abs(int(f[3]) - (int(lc.group(1)) if lc else 0) - int(tp[i, m]))
    clip = max([int(k) for k, o in ops if o in "SH"] or [0])
    indel = max([int(k) for k, o in ops if o in "ID"] or [0])
    if clip >= 10:
        kinds["soft clip ≥ 10"].append(d - clip)
    else:
        kinds["indel ≥ 10" if indel >= 10 else "neither"].append(d)
print(f"{int(near.sum()):,} confident near misses:")
for k, v in kinds.items():
    q = ", ".join(f"{x:g}" for x in np.percentile(v, [25, 50, 75])) if v else "–"
    print(f"- {k}: {len(v):,}; {'offset less the clip' if k.startswith('soft') else 'offset'}, quartiles {q}")
