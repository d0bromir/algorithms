#!/usr/bin/env python3
"""lost_reads.py MINIBWA.sam CERTA.sam: reads (primary records) that minibwa
places with MAPQ >= 20 and CERTA places more than 1 kb away, on another
contig, or not at all; by minibwa's NM, with CERTA's tier and clipping."""
import re
import sys
from collections import Counter

CLIP = re.compile(r"(\d+)S")


def records(path):
    for line in open(path):
        if line[0] == "@":
            continue
        f = line.split("\t", 11)
        flag = int(f[1])
        if flag & 0x900:
            continue
        tags = f[11] if len(f) > 11 else ""
        nm = re.search(r"NM:i:(\d+)", tags)
        xt = re.search(r"XT:Z:(\w+)", tags)
        yield (f[0], (flag >> 6) & 3), (flag, f[2], int(f[3]), int(f[4]), f[5],
                                        int(nm.group(1)) if nm else -1, xt.group(1) if xt else "-", abs(int(f[8])))


mb = {k: v for k, v in records(sys.argv[1]) if not v[0] & 4 and v[3] >= 20}
by_nm, lost_nm, how, why = Counter(), Counter(), Counter(), Counter()
seen = 0
cer = dict(records(sys.argv[2]))
rescue = Counter()
examples = []
for k, c in cer.items():
    m = mb.get(k)
    if m is None:
        continue
    seen += 1
    b = min(m[5], 10)
    by_nm[b] += 1
    lost = bool(c[0] & 4) or c[1] != m[1] or abs(c[2] - m[2]) > 1000
    if lost:
        o = cer.get((k[0], 3 - k[1]))
        near = o is not None and not o[0] & 4 and o[3] >= 20 and o[1] == m[1] and abs(o[2] - m[2]) <= 2000
        rescue["mate placed by CERTA with MAPQ >= 20 within 2 kb" if near else "no such mate"] += 1
        if near and len(examples) < 40 and hash(k[0]) % 50 == 0:
            examples.append((k, m, c, o))
        mclip = sum(int(x) for x in CLIP.findall(m[4]))
        short = 0 < m[7] < 151
        why[("insert < 151 (adapter), " if short else "insert >= 151, ") + ("minibwa clips >= 10" if mclip >= 10 else "minibwa clips < 10")] += 1
    if c[0] & 4:
        lost_nm[b] += 1
        how["unmapped"] += 1
    elif c[1] != m[1] or abs(c[2] - m[2]) > 1000:
        lost_nm[b] += 1
        clip = sum(int(x) for x in CLIP.findall(c[4]))
        how[f"elsewhere, tier {c[6]}, {'clipped >= 30' if clip >= 30 else 'clipped < 30'}, MAPQ {'>= 20' if c[3] >= 20 else '< 20'}"] += 1
print(f"{seen:,} reads that minibwa places with MAPQ >= 20 (in the targets)")
print("| minibwa's NM | reads | CERTA elsewhere or unmapped | per million |")
print("|---|---|---|---|")
for b in sorted(by_nm):
    print(f"| {b if b < 10 else '>= 10'} | {by_nm[b]:,} | {lost_nm[b]:,} | {1e6 * lost_nm[b] / by_nm[b]:,.0f} |")
tot = sum(lost_nm.values())
print(f"| all | {seen:,} | {tot:,} | {1e6 * tot / max(seen, 1):,.0f} |")
print("the lost reads, by minibwa's record:")
for k, v in why.most_common():
    print(f"  {v:,} {k}")
print("rescue by the mate's window:")
for k, v in rescue.most_common():
    print(f"  {v:,} {k}")
print("examples (rescuable): minibwa's record | CERTA's | CERTA's mate")
for k, m, c, o in examples:
    print(f"  {k[0][-16:]}/{k[1]}: {m[1]}:{m[2]} {m[4]} NM{m[5]} MQ{m[3]} TLEN{m[7]} | "
          f"{'unmapped' if c[0] & 4 else f'{c[1]}:{c[2]} {c[4]} {c[6]} NM{c[5]} MQ{c[3]}'} | {o[1]}:{o[2]} {o[4]} {o[6]} MQ{o[3]}")
print("and CERTA's:")
for k, v in how.most_common():
    print(f"  {v:,} {k}")
