#!/usr/bin/env python3
"""B2 Aganezov collapsed-loci sinks: the genuine missing copies.

The chain-uncovered gaps (b2_extract) are mostly satellites; the missing *copies*
are CHM13 segmental-duplication sequence absent from GRCh38. This intersects the
CHM13-only gaps with the CHM13 v2.0 SEDEF segdup annotation (chm13v2.0_SD.bed),
excludes Cen/Sat, and extracts the resulting missing-copy regions from hs1.fa as
sink contigs. These are the Aganezov-style collapsed-locus copies (euchromatic
SDs GRCh38 lacks), not satellites.

  python3 b2_aganezov.py HS1.fa GAPS.bed CHM13_SD.bed CENSAT.bed OUT_PREFIX [--min 1000]
"""
import sys
from collections import defaultdict


def read_intervals(path):
    iv = defaultdict(list)
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) < 3 or p[0].startswith(("#", "track", "browser")):
                continue
            iv[p[0]].append((int(p[1]), int(p[2])))
    for c in iv:
        iv[c] = merge(iv[c])
    return iv


def merge(iv):
    iv.sort()
    out = []
    for a, b in iv:
        if out and a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return out


def overlaps(iv, c, s, e):
    for a, b in iv.get(c, []):
        if a >= e:
            break
        if b > s:
            return True
    return False


def read_fasta(path):
    seqs = {}
    name = None
    chunks = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def main():
    hs1, gaps_bed, sd_bed, censat_bed, out = sys.argv[1:6]
    minlen = 1000
    if "--min" in sys.argv:
        minlen = int(sys.argv[sys.argv.index("--min") + 1])
    sd = read_intervals(sd_bed)
    censat = read_intervals(censat_bed)
    print(f"CHM13 segdup merged intervals: {sum(len(v) for v in sd.values())}", flush=True)
    print("reading hs1.fa into RAM...", flush=True)
    seqs = read_fasta(hs1)
    pieces = []
    n_gap = 0
    with open(gaps_bed) as f:
        for line in f:
            c, s, e, g = line.split()
            s, e = int(s), int(e)
            n_gap += 1
            for ss, se in sd.get(c, []):
                if ss >= e:
                    break
                if se > s:
                    a, b = max(s, ss), min(e, se)
                    if b - a >= minlen and not overlaps(censat, c, a, b):
                        pieces.append((c, a, b))
    # merge adjacent/overlapping pieces per contig
    byc = defaultdict(list)
    for c, a, b in pieces:
        byc[c].append((a, b))
    kept = []
    for c in sorted(byc):
        for a, b in merge(byc[c]):
            kept.append((c, a, b))
    total = sum(b - a for c, a, b in kept)
    with open(out + ".fa", "w") as fa, open(out + ".bed", "w") as bed:
        for i, (c, a, b) in enumerate(kept):
            sub = seqs[c][a:b].upper()
            fa.write(f">sink_{i} {c}:{a}-{b} len={b-a}\n")
            for j in range(0, len(sub), 80):
                fa.write(sub[j:j + 80] + "\n")
            bed.write(f"{c}\t{a}\t{b}\tsink_{i}\n")
    print(f"\nGaps: {n_gap}; missing-copy sink regions (segdup ∩ gap, ¬Cen/Sat, >={minlen}bp): "
          f"{len(kept)}, {total:,} bp ({total/1e6:.2f} Mb)")
    print(f"FASTA: {out}.fa  BED: {out}.bed")
    print("Largest 12 missing-copy sinks:")
    for c, a, b in sorted(kept, key=lambda r: -(r[2] - r[1]))[:12]:
        print(f"  {c}:{a}-{b}  {b-a:,} bp")


if __name__ == "__main__":
    main()
