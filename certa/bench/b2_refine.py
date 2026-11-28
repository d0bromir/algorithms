#!/usr/bin/env python3
"""B2 refinement: exclude satellite sinks using the T2T Cen/Sat v2.1 annotation.

The heuristic-filtered sinks (b2_build.py) stole 916 real reads because they
included pericentromeric/telomeric satellite fragments. This excludes any sink
candidate overlapping a Cen/Sat region (centromeric satellite, higher-order
repeat, monomeric, rDNA, ct), plus the size/acrocentric/chrY heuristic, leaving
euchromatic missing-copy candidates. Extracts the refined sinks from hs1.fa.

  python3 b2_refine.py HS1.fa GAPS.bed CENSAT.bed OUT_PREFIX
"""
import sys
from collections import defaultdict


ACRO = {"chr13", "chr14", "chr15", "chr21", "chr22"}


def read_intervals(path):
    iv = defaultdict(list)
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) < 3 or p[0].startswith(("#", "track", "browser")):
                continue
            iv[p[0]].append((int(p[1]), int(p[2])))
    for c in iv:
        iv[c].sort()
    return iv


def overlaps(iv, c, s, e):
    for a, b in iv.get(c, []):
        if a >= e:
            break
        if b > s:
            return True
    return False


def heuristic_sat(c, s, e):
    g = e - s
    return g >= 500000 or (c in ACRO and s < 20000000) or (c == "chrY" and s > 26000000)


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
    hs1, gaps_bed, censat_bed, out = sys.argv[1:5]
    censat = read_intervals(censat_bed)
    print(f"Cen/Sat regions: {sum(len(v) for v in censat.values())}", flush=True)
    print("reading hs1.fa into RAM...", flush=True)
    seqs = read_fasta(hs1)
    kept = []
    n_heur = n_censat = 0
    with open(gaps_bed) as f:
        for line in f:
            c, s, e, g = line.split()
            s, e = int(s), int(e)
            if heuristic_sat(c, s, e):
                n_heur += 1
                continue
            if overlaps(censat, c, s, e):
                n_censat += 1
                continue
            kept.append((c, s, e))
    kept.sort(key=lambda r: (r[0], r[1]))
    total = sum(e - s for c, s, e in kept)
    with open(out + ".fa", "w") as fa, open(out + ".bed", "w") as bed:
        for i, (c, s, e) in enumerate(kept):
            sub = seqs[c][s:e].upper()
            fa.write(f">sink_{i} {c}:{s}-{e} len={e-s}\n")
            for j in range(0, len(sub), 80):
                fa.write(sub[j:j + 80] + "\n")
            bed.write(f"{c}\t{s}\t{e}\tsink_{i}\n")
    print(f"\nExcluded: {n_heur} heuristic-satellite, {n_censat} Cen/Sat-overlapping gaps")
    print(f"Refined sinks: {len(kept)} contigs, {total:,} bp ({total/1e6:.2f} Mb)")
    print(f"FASTA: {out}.fa  BED: {out}.bed")
    import collections
    hist = collections.Counter()
    for c, s, e in kept:
        g = e - s
        hist["5-20k" if g < 20000 else "20-100k" if g < 100000 else "100-500k"] += 1
    print("Size distribution:", dict(hist))
    print("Largest 10 refined sinks:")
    for c, s, e in sorted(kept, key=lambda r: -(r[2] - r[1]))[:10]:
        print(f"  {c}:{s}-{e}  {e-s:,} bp")


if __name__ == "__main__":
    main()
