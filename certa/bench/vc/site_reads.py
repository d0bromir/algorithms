#!/usr/bin/env python3
"""site_reads.py SAMTOOLS REF SITES.tsv NAME=CRAM...: per indel site, the reads of each run that
cover it: how many, how many with MAPQ >= 20, how many carry an indel within 5 bp of the site in
their CIGAR (and of those, how many with MAPQ >= 20), how many are soft-clipped within 10 bp of it."""
import re
import subprocess
import sys

st, ref, sites = sys.argv[1:4]
runs = [a.split("=", 1) for a in sys.argv[4:]]
OP = re.compile(r"(\d+)([MIDNSHP=X])")


def summary(cram, chrom, pos):
    out = subprocess.run([st, "view", "-T", ref, "-F", "0x904", cram, f"{chrom}:{pos}-{pos}"],
                         capture_output=True, text=True, check=True).stdout
    n = hq = ind = ind_hq = clip = 0
    for line in out.splitlines():
        c = line.split("\t")
        mq, start, cig = int(c[4]), int(c[3]), c[5]
        n += 1
        hq += mq >= 20
        r = start
        has_ind = has_clip = False
        ops = OP.findall(cig)
        for k, (ln, op) in enumerate(ops):
            ln = int(ln)
            if op in "M=XDN":
                if op == "D" and abs(r - pos) <= 5:
                    has_ind = True
                r += ln
            elif op == "I" and abs(r - pos) <= 5:
                has_ind = True
            elif op == "S" and abs(r - pos) <= 10:
                has_clip = True
        ind += has_ind
        ind_hq += has_ind and mq >= 20
        clip += has_clip
    return n, hq, ind, ind_hq, clip


print("site\tkind\tref\talt\t" + "\t".join(f"{name}: reads, MAPQ>=20, indel, indel MAPQ>=20, clipped" for name, _ in runs))
for line in open(sites):
    chrom, pos, kind, r, a = line.rstrip("\n").split("\t")
    cells = [",".join(str(v) for v in summary(cram, chrom, int(pos))) for _, cram in runs]
    print(f"{chrom}:{pos}\t{kind}\t{r[:10]}\t{a[:12]}\t" + "\t".join(cells))
