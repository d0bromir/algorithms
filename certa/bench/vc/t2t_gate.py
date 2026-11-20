#!/usr/bin/env python3
"""t2t_gate.py SAMTOOLS BEFORE.bam AFTER.bam FEAT.tsv.zst OUT_PREFIX: the T2T-gated
MAPQ for DeepVariant. BEFORE and AFTER are the same records (same order) with
the MAPQ before and after PR mates' point MAPQ from the best other pair
found; FEAT (t2t_feat.sh) has per (qname, mate) the placement and MAPQ of the
same build's --t2t --feature-tags run and its XC (T2T class, MAPQ there).
here: XC class 1 (T2T's placement lifts within 150 bp) with T2T MAPQ >= 40;
away: class 2 (lifts elsewhere) with T2T MAPQ >= 40. Writes three BAMs of
BEFORE's records with the MAPQ of each rule:
  gate:  0 if away, else AFTER's for PR mates where here, else BEFORE's;
  raise: AFTER's for PR mates where here, else BEFORE's;
  veto:  0 if away, else BEFORE's."""
import subprocess
import sys
from collections import Counter

st, bam_b, bam_a, feat_path, prefix = sys.argv[1:6]
feat = {}
z = subprocess.Popen(["zstd", "-dc", feat_path], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
for line in z.stdout:
    q, mate, contig, pos, mq, xt, xc = line.rstrip("\n").split("\t")
    if xc == "-":
        cls, tmq = -1, -1
    else:
        p = xc.split(",")
        cls, tmq = int(p[0]), int(p[2])
    feat[(q, int(mate))] = (contig, int(pos), int(mq), cls, tmq)
z.wait()
print(f"{len(feat):,} records in {feat_path}", file=sys.stderr)


def view(path, header=False):
    p = subprocess.Popen([st, "view", *(["-h"] if header else []), path], stdout=subprocess.PIPE, text=True,
                         bufsize=1 << 20)
    return p.stdout


rules = ("gate", "raise", "veto")
outs = {r: subprocess.Popen([st, "view", "-b", "-o", f"{prefix}{r}.bam", "-"], stdin=subprocess.PIPE, text=True,
                            bufsize=1 << 20) for r in rules}
fb, fa = view(bam_b, True), view(bam_a)
line_b = next(fb)
while line_b[0] == "@":
    for o in outs.values():
        o.stdin.write(line_b)
    line_b = next(fb)
n = Counter()
tab = Counter()  # (tier, class, here/away/-, direction of AFTER vs BEFORE) -> records
first = True
for line_a in fa:
    if not first:
        line_b = next(fb)
    first = False
    b = line_b.split("\t", 11)
    a = line_a.split("\t", 11)
    if b[0] != a[0] or b[1] != a[1] or b[3] != a[3]:
        sys.exit(f"records differ: {b[:4]} {a[:4]}")
    n["records"] += 1
    qb, qa = int(b[4]), int(a[4])
    tier = "-"
    i = a[11].find("XT:Z:")
    if i >= 0:
        tier = a[11][i + 5:].split("\t", 1)[0].rstrip("\n")
    f = feat.get((a[0], (int(a[1]) >> 6) & 3))
    if f is None:
        n["no T2T record"] += 1
        cls, tmq = -1, -1
    else:
        if f[0] != a[2] or f[1] != int(a[3]):
            n["placement differs from the --t2t run"] += 1
        if f[2] != qa:
            n["MAPQ differs from the --t2t run"] += 1
        cls, tmq = f[3], f[4]
    here = cls == 1 and tmq >= 40
    away = cls == 2 and tmq >= 40
    pr_here = tier == "PR" and here
    q = {"gate": 0 if away else (qa if pr_here else qb), "raise": qa if pr_here else qb, "veto": 0 if away else qb}
    d = "same" if qa == qb else "up" if qa > qb else "down"
    tab[(tier if tier == "PR" else "other", "here" if here else "away" if away else "-", d)] += 1
    for r in rules:
        if q[r] != qb:
            n[f"{r}: MAPQ changed"] += 1
            if q[r] < 5 <= qb:
                n[f"{r}: MAPQ to below 5"] += 1
            if q[r] >= 5 > qb:
                n[f"{r}: MAPQ from below 5"] += 1
        o = outs[r].stdin
        o.write(line_b if q[r] == qb else "\t".join(b[:4] + [str(q[r])] + b[5:]))
for o in outs.values():
    o.stdin.close()
    o.wait()
for k, v in n.items():
    print(f"{k}\t{v:,}")
print("tier\tT2T\tAFTER vs BEFORE\trecords")
for k, v in sorted(tab.items()):
    print("\t".join(k) + f"\t{v:,}")
