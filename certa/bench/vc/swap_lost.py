#!/usr/bin/env python3
"""swap_lost.py [--mapq] MINIBWA.sam CERTA.sam > OUT.sam: CERTA's records,
except that each read lost_reads.py counts as lost (minibwa places it with
MAPQ >= 20, CERTA more than 1 kb away, on another contig, or not at all)
gets minibwa's record instead. With --mapq, every other CERTA record within
1 kb of minibwa's placement of the read (any minibwa MAPQ) also gets
minibwa's MAPQ, its alignment staying CERTA's. Diagnostics: DeepVariant on
OUT tells how much of the gap to minibwa those reads, and then the MAPQ,
explain."""
import sys

args = [a for a in sys.argv[1:] if a != "--mapq"]
mapq_too = len(args) < len(sys.argv) - 1


def key(f):
    return f[0], (int(f[1]) >> 6) & 3


mb = {}
for line in open(args[0]):
    if line[0] == "@":
        continue
    f = line.split("\t", 6)
    flag = int(f[1])
    if flag & 0x904:
        continue
    mb[key(f)] = (f[2], int(f[3]), int(f[4]), line)
out = sys.stdout
swapped = requal = 0
for line in open(args[1]):
    if line[0] == "@":
        out.write(line)
        continue
    f = line.split("\t", 6)
    m = None if int(f[1]) & 0x900 else mb.get(key(f))
    if m is not None:
        away = int(f[1]) & 4 or f[2] != m[0] or abs(int(f[3]) - m[1]) > 1000
        if away and m[2] >= 20:
            out.write(m[3])
            swapped += 1
            continue
        if mapq_too and not away and int(f[4]) != m[2]:
            f[4] = str(m[2])
            out.write("\t".join(f))
            requal += 1
            continue
    out.write(line)
print(f"swap_lost: {swapped:,} records swapped, {requal:,} given minibwa's MAPQ", file=sys.stderr)
