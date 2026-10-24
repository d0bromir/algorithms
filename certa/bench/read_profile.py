#!/usr/bin/env python3
"""Profile recent-Illumina read specifics in a FASTQ (first N reads).

Per file: the quality alphabet (binning), mean quality per cycle (first,
middle and last cycles), reads with a 3' poly-G tail (>= 10 / >= 20 G, two-
colour no-signal), reads with the Illumina adapter prefix AGATCGGAAGAGC
(read-through on short inserts) and where it starts, reads with N, lengths.
  read_profile.py N FASTQ[.gz]...
"""
import gzip
import sys
from collections import Counter

ADAPTER = "AGATCGGAAGAGC"


def records(path, n):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as f:
        for i in range(n):
            h = f.readline()
            if not h:
                return
            s = f.readline().rstrip("\n")
            f.readline()
            q = f.readline().rstrip("\n")
            yield s, q


def profile(path, n):
    reads = 0
    qalpha = Counter()
    qsum, qcnt = Counter(), Counter()
    lens = Counter()
    g10 = g20 = withn = adapt = 0
    apos = Counter()
    lowq_tail = 0
    for s, q in records(path, n):
        reads += 1
        lens[len(s)] += 1
        qs = [ord(c) - 33 for c in q]
        qalpha.update(qs)
        for i, v in enumerate(qs):
            qsum[i] += v
            qcnt[i] += 1
        t = len(s) - len(s.rstrip("G"))
        g10 += t >= 10
        g20 += t >= 20
        withn += "N" in s
        a = s.find(ADAPTER)
        if a >= 0:
            adapt += 1
            apos[min(a // 25 * 25, 150)] += 1
        lowq_tail += len(qs) >= 10 and sum(v < 20 for v in qs[-10:]) >= 5
    L = max(qcnt)
    mean = [qsum[i] / qcnt[i] for i in range(L + 1)]
    print(f"== {path}: {reads} reads")
    print(f"quality values used: {sorted(qalpha)} (share of bases: "
          + ", ".join(f"Q{k} {100 * v / sum(qalpha.values()):.2f} %" for k, v in sorted(qalpha.items())) + ")")
    cyc = lambda a, b: sum(mean[a:b]) / (b - a)
    print(f"mean quality: cycles 1-10 {cyc(0, 10):.1f}, 70-80 {cyc(70, 80):.1f}, last 10 {cyc(L - 9, L + 1):.1f}")
    print(f"reads with >= 5 of the last 10 bases below Q20: {lowq_tail} ({100 * lowq_tail / reads:.3f} %)")
    print(f"3' poly-G >= 10: {g10} ({100 * g10 / reads:.3f} %), >= 20: {g20} ({100 * g20 / reads:.3f} %)")
    print(f"adapter {ADAPTER}: {adapt} ({100 * adapt / reads:.3f} %), starting at (25-bp bins): "
          + ", ".join(f"{k}: {v}" for k, v in sorted(apos.items())))
    print(f"reads with N: {withn} ({100 * withn / reads:.3f} %)")
    print("lengths (>= 0.1 %): " + ", ".join(f"{k}: {v}" for k, v in sorted(lens.items()) if v >= reads / 1000))


if __name__ == "__main__":
    n = int(sys.argv[1])
    for p in sys.argv[2:]:
        profile(p, n)
