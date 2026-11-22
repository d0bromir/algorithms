#!/usr/bin/env python3
"""B2-build step 1: filter CHM13-only gaps to sink candidates and extract them.

Satellite exclusion (heuristic; the precise T2T Cen/Sat v2.1 annotation is a
follow-up download): drop gaps >= 500 kb (centromeric satellite arrays), gaps in
the acrocentric short arms (chr13/14/15/21/22, start < 20 Mb), and chrY
heterochromatin (> 26 Mb). The remaining moderate gaps are missing-copy
candidates (a segdup-flanking filter with the Aganezov/Cen-Sat annotations would
refine further). Extracts their sequence from hs1.fa into a sink FASTA (ordinary
contigs, named sink_<n>) + a BED + stats. Reads hs1.fa into RAM; outputs are small
(sink content ~tens of Mb). No index built here (B2-build step 2, disk-heavy).

  python3 b2_build.py HS1.fa GAPS.bed OUT_PREFIX
"""
import sys


ACRO = {"chr13", "chr14", "chr15", "chr21", "chr22"}


def is_satellite(c, s, e):
    g = e - s
    if g >= 500000:
        return True
    if c in ACRO and s < 20000000:
        return True
    if c == "chrY" and s > 26000000:
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
    hs1, gaps_bed, out = sys.argv[1], sys.argv[2], sys.argv[3]
    print("reading hs1.fa into RAM...", flush=True)
    seqs = read_fasta(hs1)
    print(f"  {len(seqs)} contigs, {sum(len(s) for s in seqs.values()):,} bp", flush=True)
    kept = []
    n_sat = 0
    with open(gaps_bed) as f:
        for line in f:
            c, s, e, g = line.split()
            s, e = int(s), int(e)
            if is_satellite(c, s, e):
                n_sat += 1
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
    print(f"\nExcluded {n_sat} satellite gaps; kept {len(kept)} sink candidates, {total:,} bp ({total/1e6:.1f} Mb)")
    print(f"Sink FASTA: {out}.fa  BED: {out}.bed")
    print("Largest 10 sinks:")
    for c, s, e in sorted(kept, key=lambda r: -(r[2] - r[1]))[:10]:
        print(f"  {c}:{s}-{e}  {e-s:,} bp")
    # N content (assembly gaps within sinks) and size distribution
    import collections
    hist = collections.Counter()
    for c, s, e in kept:
        g = e - s
        hist["5-20k" if g < 20000 else "20-100k" if g < 100000 else "100-500k"] += 1
    print("Sink size distribution:", dict(hist))


if __name__ == "__main__":
    main()
