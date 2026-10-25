"""Aim 2 on its pre-registered cohort.

The cohort: reads the default pipeline (q22 + SL) leaves uncertified while
searching fewer parts than the read has (|S| < P: repeat-limited). They are
followed through the later tiers: q22 large budget (pass 3), q15 (Aim 1),
tier R (Aim 2).

  aim2_prereg.py tagged_unc.fq unc_best.fq resid.fq tierR.tsv aligner.sam
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


tagged, unc_best, resid, tsv, sam = sys.argv[1:6]
cohort = set()
with open(tagged) as f:
    for i, line in enumerate(f):
        if i % 4:
            continue
        t = line[1:].split()
        kv = dict(x.split("=") for x in t[1:])
        if int(kv["used"]) < int(kv["parts"]):
            cohort.add(t[0])
best, res = names(unc_best), names(resid)
tier = {}
for line in open(tsv):
    if line.startswith("read\t"):
        continue
    t = line.rstrip("\n").split("\t")
    tier[t[0]] = (int(t[1]), int(t[2]), int(t[3]), int(t[4]))
al = {}
for line in open(sam):
    if line[0] == "@":
        continue
    t = line.split("\t")
    if t[0] not in cohort or int(t[1]) & 0x900:
        continue
    al[t[0]] = None if int(t[1]) & 4 else score(t[5], int(re.search(r"NM:i:(\d+)", line).group(1)))


def status(n):
    if n not in best:
        return "exact: q22 large budget"
    if n not in res:
        return "exact: q15 (Aim 1)"
    e, lo, hi, copies = tier.get(n, (0, -1, 0, 0))
    if e:
        return "exact: tier R" + (" (proven tie)" if copies >= 2 else "")
    if lo != -1 and copies >= 2 and hi - lo <= 5:
        return "proven ambiguous: tier R"
    return "interval only"


for label, pred in [("all repeat-limited (|S| < P)", lambda s: True),
                    ("  of which aligner score >= 120 (near-perfect)", lambda s: s is not None and s >= 120),
                    ("  of which aligner score < 120 or unmapped", lambda s: s is None or s < 120)]:
    members = [n for n in cohort if pred(al.get(n))]
    c = collections.Counter(status(n) for n in members)
    N = len(members)
    good = sum(v for k, v in c.items() if not k.startswith("interval"))
    print(f"{label}: {N} reads; exact or proven ambiguous {good} ({100 * good / max(N, 1):.1f}%)")
    for k in sorted(c):
        print(f"      {k:36s} {c[k]:7d}  {100 * c[k] / N:6.2f}%")
