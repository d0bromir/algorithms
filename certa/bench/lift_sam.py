#!/usr/bin/env python3
"""Lifts an aligner's primary SAM records to another assembly through a UCSC
chain file, for truth_eval.py (bench/VALIDATION_RESULTS.md, "Matched
sensitivity").

  bwa-mem2 mem hs1.fa R1.fq R2.fq | lift_sam.py hs1ToHg38.over.chain.gz \\
      | truth_eval.py score truth.npz GRCh38.fai out.npz
  vg giraffe ... -o SAM --ref-name GRCh38 | lift_sam.py --names | truth_eval.py score ...

Each record keeps its name, flag, MAPQ and span; its position becomes the
lifted unclipped leftmost base, so truth_eval scores it as on the target.
Target contigs get galaxy's Ensembl names (truth_lift.ensembl).
The read's first base on the source (its unclipped leftmost, or rightmost
when the chain is reversed) is lifted; if it falls in a chain gap, the other
end is lifted and the position offset by the span. A record neither end of
which lifts is written unmapped and counted on stderr. Secondary and
supplementary records are dropped.
"""
import bisect
import gzip
import os
import re
import sys

sys.path.insert(0, os.path.dirname(__file__))
from truth_lift import ensembl  # noqa: E402

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def load_chain(path):
    """Per source contig: sorted block starts, and (start, end, target, strand, size, target start)."""
    blocks = {}
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("chain"):
                h = line.split()
                t, tpos, q, qsize, qstrand, qpos = h[2], int(h[5]), ensembl(h[7]), int(h[8]), h[9], int(h[10])
                cur = blocks.setdefault(t, [])
            elif line.strip():
                p = line.split()
                n = int(p[0])
                cur.append((tpos, tpos + n, q, qstrand, qsize, qpos))
                tpos += n
                qpos += n
                if len(p) == 3:
                    tpos += int(p[1])
                    qpos += int(p[2])
    index = {}
    for t, b in blocks.items():
        b.sort()
        index[t] = ([x[0] for x in b], b)
    return index


def lift(index, contig, pos):
    """(target, strand, forward position) of a 0-based source position, or None."""
    e = index.get(contig)
    if e is None:
        return None
    i = bisect.bisect_right(e[0], pos) - 1
    if i < 0:
        return None
    s, end, q, strand, qsize, qs = e[1][i]
    if pos >= end:
        return None
    p = qs + (pos - s)
    return (q, strand, p if strand == "+" else qsize - 1 - p)


def main():
    out = sys.stdout
    if sys.argv[1] == "--names":  # no lift: PanSN names (GRCh38#0#chr1) to Ensembl ones
        for line in sys.stdin:
            if line[0] != "@":
                t = line.split("\t", 3)
                t[2] = ensembl(t[2].rsplit("#", 1)[-1])
                out.write("\t".join(t))
        return
    index = load_chain(sys.argv[1])
    lifted = unlifted = 0
    for line in sys.stdin:
        if line[0] == "@":
            continue
        t = line.split("\t", 11)
        flag = int(t[1])
        if flag & (256 | 2048):
            continue
        if flag & 4:
            out.write("\t".join(t[:11]).rstrip("\n") + "\n")
            continue
        ops = CIGAR.findall(t[5])
        lead = int(ops[0][0]) if ops and ops[0][1] in "SH" else 0
        span = sum(int(k) for k, o in ops if o in "MDN=XSH")
        left = int(t[3]) - 1 - lead
        right = left + span - 1
        a, b = lift(index, t[2], left), lift(index, t[2], right)
        if a and a[1] == "+":
            q, new, flip = a[0], a[2], False
        elif b and b[1] == "-":
            q, new, flip = b[0], b[2], True
        elif b and b[1] == "+":
            q, new, flip = b[0], b[2] - (span - 1), False
        elif a and a[1] == "-":
            q, new, flip = a[0], a[2] - (span - 1), True
        else:
            unlifted += 1
            out.write(f"{t[0]}\t{(flag | 4) & ~16}\t*\t0\t0\t*\t*\t0\t0\t*\t*\n")
            continue
        lifted += 1
        out.write(f"{t[0]}\t{flag ^ 16 if flip else flag}\t{q}\t{new + 1}\t{t[4]}\t{span}M\t*\t0\t0\t*\t*\n")
    print(f"[lift_sam] lifted {lifted}, not lifted {unlifted}", file=sys.stderr)


if __name__ == "__main__":
    main()
