#!/usr/bin/env python3
"""B2-extract: characterize the CHM13-only sequence (sink candidates).

Parses the hs1ToHg38 lift-over chain, computes the hs1 (CHM13) intervals it
covers, and reports the uncovered gaps -- sequence present in CHM13 but absent
from GRCh38. These are the sink candidates for B2 (one-pass sink contigs). Large
gaps are typically centromeric satellite arrays (to exclude); moderate gaps
flanked by segmental duplications are the missing-copy candidates.

  python3 b2_extract.py CHAIN.gz OUT_BED [--min-gap 5000]

Writes gaps >= min-gap to OUT_BED (chrom start end size) and prints per-contig
coverage, gap totals, and a size histogram. Read-only on the chain; the BED is
small. No reference/index built here (that is B2-build, after this quantifies
the sink content and disk cost).
"""
import argparse
import gzip
from collections import defaultdict


def parse_chain(path):
    # The hs1ToHg38 chain is oriented hg38(query) -> hs1(target): header
    # tName/tSize are CHM13 (hs1), qName/qSize are GRCh38. CHM13-only sequence
    # is the uncovered TARGET (hs1) intervals, so track tName/tStart and dt.
    covered = defaultdict(list)
    sizes = {}
    with gzip.open(path, "rt") as f:
        tname = None
        tpos = 0
        tsize = 0
        tminus = False
        for line in f:
            line = line.strip()
            if not line:
                continue
            p = line.split()
            if p[0] == "chain":
                # chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
                tname = p[2]
                tsize = int(p[3])
                tminus = p[4] == "-"
                tpos = int(p[5])
                sizes[tname] = tsize
            else:
                size = int(p[0])
                dt = int(p[1]) if len(p) > 1 else 0
                if tminus:
                    covered[tname].append((tsize - tpos - size, tsize - tpos))
                else:
                    covered[tname].append((tpos, tpos + size))
                tpos += size + dt
    return covered, sizes


def merge(iv):
    iv.sort()
    out = []
    for a, b in iv:
        if out and a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("chain")
    ap.add_argument("out_bed")
    ap.add_argument("--min-gap", type=int, default=5000)
    a = ap.parse_args()
    covered, sizes = parse_chain(a.chain)
    hist = defaultdict(int)
    total_gap = 0
    total_kept = 0
    n_kept = 0
    rows = []
    print(f"{'contig':8} {'size':>12} {'covered':>12} {'gap_bp':>12} {'gaps>=min':>10} {'kept_bp':>12}")
    for c in sorted(sizes, key=lambda x: (len(x), x)):
        cov = merge(covered.get(c, []))
        gaps = []
        prev = 0
        for s, e in cov:
            if s > prev:
                gaps.append((prev, s))
            prev = max(prev, e)
        if prev < sizes[c]:
            gaps.append((prev, sizes[c]))
        gap_bp = sum(e - s for s, e in gaps)
        kept = [(s, e) for s, e in gaps if e - s >= a.min_gap]
        kept_bp = sum(e - s for s, e in kept)
        total_gap += gap_bp
        total_kept += kept_bp
        n_kept += len(kept)
        for s, e in gaps:
            g = e - s
            bucket = "<1k" if g < 1000 else "1-5k" if g < 5000 else "5-50k" if g < 50000 else "50-500k" if g < 500000 else ">=500k"
            hist[bucket] += 1
        for s, e in kept:
            rows.append((c, s, e, e - s))
        print(f"{c:8} {sizes[c]:12d} {sizes[c]-gap_bp:12d} {gap_bp:12d} {len(kept):10d} {kept_bp:12d}")
    with open(a.out_bed, "w") as f:
        for c, s, e, g in rows:
            f.write(f"{c}\t{s}\t{e}\t{g}\n")
    print(f"\nTotal CHM13-only gap bp: {total_gap:,} ({total_gap/1e6:.1f} Mb)")
    print(f"Gaps >= {a.min_gap} bp (sink candidates): {n_kept}, {total_kept:,} bp ({total_kept/1e6:.1f} Mb)")
    print("Gap size histogram (all gaps):")
    for b in ["<1k", "1-5k", "5-50k", "50-500k", ">=500k"]:
        print(f"  {b:8}: {hist[b]}")
    print(f"\nLargest 15 sink candidates:")
    for c, s, e, g in sorted(rows, key=lambda r: -r[3])[:15]:
        print(f"  {c}:{s}-{e}  {g:,} bp")
    print(f"\nWrote {len(rows)} gaps to {a.out_bed}")


if __name__ == "__main__":
    main()
