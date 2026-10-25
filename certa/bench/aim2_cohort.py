"""Follow the reads the default pipeline (q22 + SL) leaves uncertified through
the later tiers: q22 large budget (pass 3), q15 (Aim 1), tier R (Aim 2).

  aim2_cohort.py unc_def.fq unc_best.fq resid.fq tierR.tsv aligner.sam

Cohorts by the aligner (BWA-MEM2): repeat-like = MAPQ < 20; near-perfect =
aligner score >= 120 (loses <= 30 points, the ceiling of 6 exact parts).
"""
import collections
import re
import sys

CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def score(cig, nm):
    ops = [(int(n), o) for n, o in CIG.findall(cig)]
    m = sum(n for n, o in ops if o in "M=X")
    gaps = [n for n, o in ops if o in "ID"]
    mism = max(0, nm - sum(gaps))
    clips = sum(1 for n, o in ops if o in "SH")
    return (m - mism) - 4 * mism - sum(6 + g for g in gaps) - 5 * clips


def names(fq):
    out = set()
    with open(fq) as f:
        for i, line in enumerate(f):
            if i % 4 == 0:
                out.add(line[1:].split()[0])
    return out


unc_def, unc_best, resid = names(sys.argv[1]), names(sys.argv[2]), names(sys.argv[3])
tier = {}
for line in open(sys.argv[4]):
    if line.startswith("read\t"):
        continue
    t = line.rstrip("\n").split("\t")
    tier[t[0]] = (int(t[1]), int(t[2]), int(t[3]), int(t[4]))
al = {}
for line in open(sys.argv[5]):
    if line[0] == "@":
        continue
    t = line.split("\t")
    if t[0] not in unc_def or int(t[1]) & 0x900:
        continue
    if int(t[1]) & 4:
        al[t[0]] = None
        continue
    al[t[0]] = (int(t[4]), score(t[5], int(re.search(r"NM:i:(\d+)", line).group(1))))


def status(n):
    if n not in unc_best:
        return "exact: q22 large budget"
    if n not in resid:
        return "exact: q15 (Aim 1)"
    e, lo, hi, copies = tier.get(n, (0, -1, 0, 0))
    if e:
        return "exact: tier R" + (" (proven tie)" if copies >= 2 else "")
    if lo != -1 and copies >= 2 and hi - lo <= 5:
        return "proven ambiguous: tier R"
    return "interval only"


cohorts = {
    "repeat-like (aligner MAPQ < 20)": lambda a: a is not None and a[0] < 20,
    "repeat-like and near-perfect (MAPQ < 20, score >= 120)": lambda a: a is not None and a[0] < 20 and a[1] >= 120,
    "repeat-like and lossy (MAPQ < 20, score < 120)": lambda a: a is not None and a[0] < 20 and a[1] < 120,
    "unique-like (MAPQ >= 20)": lambda a: a is not None and a[0] >= 20,
}
for label, pred in cohorts.items():
    members = [n for n in unc_def if pred(al.get(n))]
    c = collections.Counter(status(n) for n in members)
    N = len(members)
    good = sum(v for k, v in c.items() if k.startswith("exact") or k.startswith("proven"))
    print(f"{label}: {N} reads; exact or proven ambiguous {good} ({100 * good / max(N, 1):.1f}%)")
    for k in sorted(c):
        print(f"   {k:40s} {c[k]:7d}  {100 * c[k] / N:6.2f}%")
