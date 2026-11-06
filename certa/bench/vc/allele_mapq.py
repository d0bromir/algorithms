#!/usr/bin/env python3
"""allele_mapq.py SITES < mpileup --output-MAPQ: at SNV sites (CHR POS REF ALT),
mean MAPQ and share below 5 of the reads showing ALT and of those showing REF."""
import sys
alt_of = {}
for line in open(sys.argv[1]):
    c, p, r, a = line.split()[:4]
    alt_of[(c, int(p))] = a.upper()
sums = {"alt": [0, 0, 0], "ref": [0, 0, 0]}  # reads, MAPQ sum, MAPQ < 5
for line in sys.stdin:
    f = line.rstrip("\n").split("\t")
    c, p, ref, bases, mq = f[0], int(f[1]), f[2].upper(), f[4], f[6]
    a = alt_of.get((c, p))
    if a is None:
        continue
    i, k, calls = 0, 0, []
    while i < len(bases):  # strip read starts/ends and indels from the pileup string
        ch = bases[i]
        if ch == "^":
            i += 2
            continue
        if ch == "$":
            i += 1
            continue
        if ch in "+-":
            j = i + 1
            while j < len(bases) and bases[j].isdigit():
                j += 1
            i = j + int(bases[i + 1:j])
            continue
        calls.append(ch)
        i += 1
    for ch, q in zip(calls, mq):
        q = ord(q) - 33
        b = ref if ch in ".," else ch.upper()
        key = "alt" if b == a else "ref" if b == ref else None
        if key:
            s = sums[key]
            s[0] += 1
            s[1] += q
            s[2] += q < 5
for key, (n, t, low) in sums.items():
    print(f"{key}: {n:,} reads, mean MAPQ {t / max(n, 1):.1f}, MAPQ < 5 {100 * low / max(n, 1):.1f} %")
