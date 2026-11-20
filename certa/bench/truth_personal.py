#!/usr/bin/env python3
"""The limit of any placement: the reads mapped to the sample's own genome
(bench/VALIDATION_RESULTS.md, "How far from the limit").

  ALIGNER HAP.idx R1.fq R2.fq | truth_personal.py place blocks.npz N HAP.npz
  truth_personal.py score truth.npz MAT.npz PAT.npz [--classes NAME=CLS.npz ...]

`place` keeps each mate's primary record on one haplotype genome (reads
named <letter><pair index>, N pairs) and projects it to GRCh38 through the
dipcall blocks, as truth_lift.py projects the truth: contig, unclipped
leftmost base, strand, MAPQ, AS, and whether it projects at all.

`score` takes, per pair, the haplotype whose two records score more (AS
summed; the maternal one on a tie) and counts the mates that project within
150 bp of their truth: what an aligner that knew the sample's genome would
place. A mate away from its truth there comes from sequence that repeats in
the sample itself, or that the projection cannot follow. The pair is
"settled" when both haplotypes project it to the same loci, or one scores
at least 10 more, and its MAPQ there is >= 40.
--classes: truth_oracle.py's classes of a GRCh38 run; for each class, the
mates the sample's own genome places at the truth.
"""
import argparse
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from truth_lift import CLS_TRUTH, Lift  # noqa: E402

LOCUS = 150
CIG = re.compile(r"(\d+)([MIDNSHP=X])")
AS = re.compile(r"\tAS:i:(-?\d+)")


def cmd_place(a):
    lift = Lift(a.blocks)
    n = a.n
    h = np.full((n, 2), -1, np.int32)
    s = np.zeros((n, 2), np.int64)
    e = np.zeros((n, 2), np.int64)
    st = np.zeros((n, 2), np.int8)
    lead = np.zeros((n, 2), np.int32)
    trail = np.zeros((n, 2), np.int32)
    mapq = np.zeros((n, 2), np.uint8)
    sc = np.zeros((n, 2), np.int16)
    for line in sys.stdin:
        if line[0] == "@":
            continue
        f = line.split("\t", 11)
        flag = int(f[1])
        if flag & (4 | 256 | 2048):
            continue
        hi = lift.hidx.get(f[2])
        if hi is None:
            continue
        i, m = int(f[0][1:]), 1 if flag & 128 else 0
        ops = CIG.findall(f[5])
        span = sum(int(k) for k, o in ops if o in "MD=XN")
        h[i, m] = hi
        s[i, m] = int(f[3]) - 1
        e[i, m] = s[i, m] + span
        st[i, m] = 1 if flag & 16 else 0
        lead[i, m] = int(ops[0][0]) if ops and ops[0][1] in "SH" else 0
        trail[i, m] = int(ops[-1][0]) if len(ops) > 1 and ops[-1][1] in "SH" else 0
        mapq[i, m] = min(int(f[4]), 255)
        x = AS.search(line)
        sc[i, m] = int(x.group(1)) if x else 0
    on = h >= 0
    valid, rc, left, rs = lift.project(np.where(on, h, 0), s, np.where(on, e, s + 1), st)
    valid &= on
    # The clip at the end that is leftmost on GRCh38: the leading one unless
    # the haplotype contig aligns to GRCh38 in reverse.
    flipped = rs != st
    left = left - np.where(flipped, trail, lead)
    np.savez(a.out, mapped=on, placed=valid, rc=rc.astype(np.int16), left=left.astype(np.int32),
             strand=rs.astype(np.int8), mapq=mapq, score=sc)
    print(f"[personal place] mapped {int(on.sum()):,}, projected {int(valid.sum()):,}", file=sys.stderr)


def cmd_score(a):
    t = np.load(a.truth)
    truth = t["cls"] == CLS_TRUTH
    nt = int(truth.sum())
    hp = [dict(np.load(p)) for p in (a.mat, a.pat)]
    for z in hp:
        z["left"] = z["left"].astype(np.int64)
    tot = [np.where(z["mapped"], z["score"].astype(np.int64), 0).sum(axis=1) for z in hp]
    pick = (tot[1] > tot[0])[:, None]  # the paternal genome scores more
    pk = lambda k: np.where(pick, hp[1][k], hp[0][k])  # noqa: E731
    ot = lambda k: np.where(pick, hp[0][k], hp[1][k])  # noqa: E731
    placed, rc, left, strand, mapq = pk("placed"), pk("rc"), pk("left"), pk("strand"), pk("mapq")
    tl = t["tp"].astype(np.int64) - 1
    ok = placed & truth & (rc == t["tc"]) & (strand == t["ts"]) & (np.abs(left - tl) <= LOCUS)
    same = (placed & ot("placed") & (rc == ot("rc")) & (strand == ot("strand"))
            & (np.abs(left - ot("left")) <= LOCUS))
    margin = np.abs(tot[1] - tot[0])[:, None]
    settled = (mapq >= 40) & ((same & (ot("mapq") >= 40)) | (~same & (margin >= 10)))
    mis = truth & ~ok
    print(f"Mates with a truth: {nt:,}. Mapped to the sample's own genome and projected: at the truth "
          f"{int(ok.sum()):,} ({100 * ok.sum() / nt:.3f} %); away from it {int((mis & placed).sum()):,} "
          f"({1e6 * (mis & placed).sum() / nt:,.0f} per million), not projected {int((mis & ~placed).sum()):,} "
          f"({1e6 * (mis & ~placed).sum() / nt:,.0f} per million).")
    st = truth & settled
    print(f"Settled pairs (MAPQ >= 40, the haplotypes agree or one scores 10 more): {int(st.sum()):,} mates "
          f"({100 * st.sum() / nt:.3f} %); away from the truth {int((st & ~ok).sum()):,} "
          f"({1e6 * (st & ~ok).sum() / max(1, st.sum()):,.1f} per million of them).")
    for spec in a.classes:
        name, p = spec.split("=", 1)
        cls = np.load(p)["cls"]
        print(f"\n{name}: of each class of its mates away from the truth, those the sample's own genome places "
              "at the truth (and settled):")
        for c, cn in ((1, "missed"), (2, "tie"), (3, "reference")):
            x = cls == c
            print(f"  {cn}: {int(x.sum()):,} -> {int((x & ok).sum()):,} ({100 * (x & ok).sum() / max(1, x.sum()):.1f} %), "
                  f"settled {int((x & ok & settled).sum()):,} ({100 * (x & ok & settled).sum() / max(1, x.sum()):.1f} %)")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    p = sp.add_parser("place")
    p.add_argument("blocks"), p.add_argument("n", type=int), p.add_argument("out")
    p = sp.add_parser("score")
    p.add_argument("truth"), p.add_argument("mat"), p.add_argument("pat")
    p.add_argument("--classes", action="append", default=[])
    a = ap.parse_args()
    {"place": cmd_place, "score": cmd_score}[a.cmd](a)


if __name__ == "__main__":
    main()
