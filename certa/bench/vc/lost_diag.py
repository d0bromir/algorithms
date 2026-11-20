#!/usr/bin/env python3
"""lost_diag.py SAMTOOLS REF MINIBWA.sam CERTA.sam: why CERTA loses reads that
minibwa places with MAPQ >= 20 (lost_reads.py's lost reads: CERTA > 1 kb
away, on another contig, or unmapped). MINIBWA.sam and CERTA.sam hold both
mates of the lost reads (CERTA's with --feature-tags --witness).

minibwa's alignment is re-scored under CERTA's objective (+1/-4, gap 6+k,
5 per clipped end; differences recounted against the reference) and set
against CERTA's record: its score (AS), the interval [XL, XU] that holds the
best score over the reference, and its tier R parts (XW) on that strand.
  search:  minibwa's score is above CERTA's XL and <= XU (CERTA's interval is
           honest, but its search did not reach the locus);
  bound:   minibwa's score is above XU (the bound is wrong, or the
           objectives differ);
  score:   minibwa's score is <= CERTA's (CERTA's placement scores as well
           under its objective; the pair decides, or minibwa's does not
           score the clip/mismatch trade-off the same).
Search misses are split by whether some tier R part occurs exactly at
minibwa's locus, and by the longest exact run there (>= 22: a q-mer seed
exists; < 22: the index cannot seed the locus). The other mate's placement
by CERTA is set against minibwa's locus (within 2 kb: a mate window would
reach it)."""
import re
import subprocess
import sys
from collections import Counter

st, ref, mb_path, ce_path = sys.argv[1:5]
OPS = re.compile(r"(\d+)([MIDNSHP=X])")


def recs(path):
    out = {}
    for line in open(path):
        if line[0] == "@":
            continue
        f = line.rstrip("\n").split("\t")
        if int(f[1]) & 0x900:
            continue
        out[(f[0], (int(f[1]) >> 6) & 3)] = f
    return out


def tag(f, name, default=None):
    for t in f[11:]:
        if t.startswith(name + ":"):
            v = t.split(":", 2)[2]
            return int(v) if t[3] == "i" else v
    return default


mb, ce = recs(mb_path), recs(ce_path)
lost = []
for k, m in mb.items():
    flag = int(m[1])
    if flag & 4 or int(m[4]) < 20:
        continue
    c = ce.get(k)
    if c is None:
        continue
    cf = int(c[1])
    if not cf & 4 and c[2] == m[2] and abs(int(c[3]) - int(m[3])) <= 1000:
        continue
    lost.append(k)
print(f"{len(lost):,} lost reads (minibwa MAPQ >= 20; CERTA > 1 kb away, another contig or unmapped)", file=sys.stderr)

# the reference under minibwa's alignments, one samtools call
regions = {}
for k in lost:
    m = mb[k]
    span = sum(int(n) for n, o in OPS.findall(m[5]) if o in "MD=XN")
    regions[k] = f"{m[2]}:{m[3]}-{int(m[3]) + span - 1}"
with open("/tmp/lost_regions.txt", "w") as f:
    f.write("\n".join(sorted(set(regions.values()))) + "\n")
fa = subprocess.run([st, "faidx", ref, "-r", "/tmp/lost_regions.txt"], capture_output=True, text=True).stdout
seqs, name, buf = {}, None, []
for line in fa.split("\n"):
    if line.startswith(">"):
        if name:
            seqs[name] = "".join(buf).upper()
        name, buf = line[1:], []
    elif line:
        buf.append(line)
if name:
    seqs[name] = "".join(buf).upper()

tab = Counter()
runs_hist = Counter()
examples = {}
for k in lost:
    m, c = mb[k], ce[k]
    flag, cig, seq = int(m[1]), m[5], m[9]
    rs = seqs[regions[k]]
    bad, qi, ri, mism, gaps, gapb = [], 0, 0, 0, 0, 0
    ops = OPS.findall(cig)
    for n, o in ops:
        n = int(n)
        if o in "M=X":
            for j in range(n):
                if seq[qi + j] != rs[ri + j]:
                    bad.append(qi + j); mism += 1
            qi += n; ri += n
        elif o == "I":
            bad.extend(range(qi, qi + n)); qi += n; gaps += 1; gapb += n
        elif o == "S":
            bad.extend(range(qi, qi + n)); qi += n
        elif o in "DN":
            bad.append(qi); ri += n; gaps += 1; gapb += n
    aligned = sum(int(n) for n, o in ops if o in "M=X")
    ends = (ops[0][1] in "SH") + (len(ops) > 1 and ops[-1][1] in "SH")
    mscore = aligned - 5 * mism - 6 * gaps - gapb - 5 * ends
    L = len(seq)
    last, longest = -1, 0
    for x in sorted(set(bad)) + [L]:
        longest = max(longest, x - last - 1); last = x
    cf = int(c[1])
    xt = tag(c, "XT", "-")
    xl, xu = tag(c, "XL", -999), tag(c, "XU", -999)
    cas = tag(c, "AS", -999) if not cf & 4 else -999
    best_found = max(cas, xl)
    if mscore > xu:
        cause = "bound: minibwa scores above XU"
    elif mscore > best_found:
        cause = "search: CERTA did not reach the locus"
    else:
        cause = "score: CERTA's placement scores as well"
    strand = 1 if flag & 16 else 0
    xw = tag(c, "XW", "")
    exact_part = "-"
    if cause.startswith("search"):
        parts = xw.split(";")[1 + strand].split(",") if xw.count(";") >= 2 else []
        if not parts:
            exact_part = "no parts"
        else:
            hit = any(not any(int(p.split("-")[0]) <= x < int(p.split("-")[1]) for x in bad) for p in parts if "-" in p)
            exact_part = "a part exact there" if hit else "no part exact there"
        runs_hist["run >= 22" if longest >= 22 else "run < 22"] += 1
    mk = (k[0], 3 - k[1])
    cm = ce.get(mk)
    mate = "-"
    if cm is not None:
        if int(cm[1]) & 4:
            mate = "mate unmapped"
        else:
            near = cm[2] == m[2] and abs(int(cm[3]) - int(m[3])) <= 2000
            mq = int(cm[4])
            mate = ("mate within 2 kb" if near else "mate elsewhere") + (", MAPQ >= 40" if mq >= 40 else ", MAPQ < 40")
    state = "unmapped" if cf & 4 else xt
    nmb = "NM 0-4" if mism + gapb < 5 else "NM 5-9" if mism + gapb < 10 else "NM >= 10"
    key = (cause, exact_part, "run >= 22" if longest >= 22 else "run < 22", state, mate)
    tab[key] += 1
    tab2 = tab.setdefault(("_nm", cause, nmb), 0)
    tab[("_nm", cause, nmb)] = tab2 + 1
    if key not in examples:
        examples[key] = (f"{k[0]}/{k[1]} minibwa {m[2]}:{m[3]} {cig} score {mscore} | CERTA {c[2]}:{c[3]} {c[5]} {xt} "
                         f"AS {cas} XL {xl} XU {xu} | mate {cm[2] + ':' + cm[3] + ' ' + cm[5] if cm else '-'}")
tot = len(lost)
print("| cause | parts at minibwa's locus | longest exact run | CERTA's record | CERTA's mate | reads | share |")
print("|---|---|---|---|---|---|---|")
for key, n in sorted(((k, v) for k, v in tab.items() if k[0] != "_nm"), key=lambda x: -x[1]):
    print("| " + " | ".join(key) + f" | {n:,} | {n / tot:.1%} |")
print("\nby cause and minibwa's differences:")
for key, n in sorted((k, v) for k, v in tab.items() if k[0] == "_nm"):
    print(f"  {key[1]} {key[2]}: {n:,}")
print("\nexamples:")
for key, n in sorted(((k, v) for k, v in tab.items() if k[0] != "_nm"), key=lambda x: -x[1])[:12]:
    print(f"  [{' / '.join(key)}] {examples[key]}")
