#!/usr/bin/env python3
"""A1: the offline ceiling of a gap-graded MAPQ (research guideline, Pillar A).

Re-ranks a CERTA run's mates by min(its MAPQ, a gap-graded MAPQ fitted from a
competitor aligner's AS-XS gap) and reports misplaced per million at matched N
(the N most confident mates with a truth, boundary MAPQ pro rata, as
truth_mapq_pair.py / truth_headroom.py). No mapper code is touched: the gap
grades confidence only, the certificate is unchanged. The competitor gap is
BWA-MEM2's AS-XS by default (its XS is a *local* score, so it counts clipped
partial paralogs that CERTA's end-to-end second-best misses). Fit on HG002's
even pairs and report the odd pairs held out (--cv / --half), per the
calibrate-then-freeze rule.

  samtools view BWA.bam | truth_ceiling.py feat TRUTH.npz CERTA.npz OUT.npz
      (or stream the aligner, writing no BAM:
       bwa-mem2 mem -t 64 IDX R1.fq R2.fq | truth_ceiling.py feat ...)
  truth_ceiling.py fit [--cv] OUT.npz MODEL.npz
  truth_ceiling.py eval [--half] [--shares S,S,...] MODEL.npz OUT.npz

`feat` keeps per mate: CERTA's placed / at-locus (ok) / MAPQ, whether it has a
GRCh38 truth, and the competitor's gap = AS-XS (0..GAPMAX; GAPMAX when the
competitor mapped with no suboptimal hit, i.e. no competitor; -1 when the
competitor did not map the mate, so it adds no penalty). `fit` pools the gap
into bins and runs weighted PAVA (isotonic, non-increasing error in gap),
writing a gap->MAPQ lookup (Phred of the pooled misplaced share, +0.5/+1, at
most 60, as truth_mapq_fit.py). `eval` reports the baseline MAPQ, the re-ranked
min(MAPQ, gap MAPQ), and the gap alone, at matched N.
"""
import argparse
import re
import sys

import numpy as np

sys.path.insert(0, __import__("os").path.dirname(__file__))
from truth_lift import CLS_TRUTH  # noqa: E402

GAPMAX = 1000  # gap clip and the "mapped, no competitor" sentinel
AS = re.compile(r"\tAS:i:(-?\d+)")
XS = re.compile(r"\tXS:i:(-?\d+)")


def read_sam_gap(stream, n):
    """Per-mate competitor gap (n pairs x 2), -1 where it has no alignment score."""
    gap = np.full((n, 2), -1, np.int32)
    for line in stream:
        if line[0] == "@":
            continue
        t = line.split("\t", 11)
        flag = int(t[1])
        if flag & (256 | 2048 | 4):
            continue
        i = int(t[0][1:])
        m = 1 if flag & 128 else 0
        if len(t) > 11:
            rest = "\t" + t[11]
            a = AS.search(rest)
            if a:
                x = XS.search(rest)
                g = int(a.group(1)) - int(x.group(1)) if x else GAPMAX
                gap[i, m] = min(max(g, 0), GAPMAX)
    return gap


def cmd_feat(a):
    t = np.load(a.truth)
    n = t["cls"].shape[0]
    truth = t["cls"] == CLS_TRUTH
    r = np.load(a.certa)
    gap = read_sam_gap(sys.stdin, n)
    np.savez(a.out, truth=truth, placed=r["placed"], ok=r["correct"],
             mapq=r["mapq"].astype(np.int64), gap=gap)
    print(f"[feat] {int(truth.sum()):,} mates with a truth; competitor scored "
          f"{int((gap >= 0).sum()):,}; no competitor (gap={GAPMAX}) "
          f"{int((gap == GAPMAX).sum()):,}", file=sys.stderr)


def phred(e, n):
    return min(60, int(round(-10 * np.log10((e + 0.5) / (n + 1)))))


def pava_decreasing(counts, errors):
    """Weighted pool-adjacent-violators: pool bins until errors/counts is
    non-increasing. Returns per-bin pooled (counts, errors) of equal length."""
    bc, be, bn = [], [], []
    for c, e in zip(counts, errors):
        bc.append(int(c))
        be.append(int(e))
        bn.append(1)
        while len(bc) > 1 and be[-2] * bc[-1] < be[-1] * bc[-2]:
            bc[-2] += bc[-1]
            be[-2] += be[-1]
            bn[-2] += bn[-1]
            bc.pop()
            be.pop()
            bn.pop()
    out_c = np.empty(len(counts), np.int64)
    out_e = np.empty(len(counts), np.int64)
    i = 0
    for c, e, k in zip(bc, be, bn):
        out_c[i:i + k] = c
        out_e[i:i + k] = e
        i += k
    return out_c, out_e


def build_qmap(counts, errors):
    """gap->MAPQ lookup (GAPMAX+1) from per-gap counts/errors: isotonic
    non-increasing error rate, then Phred of the pooled share."""
    counts = np.asarray(counts, np.int64)
    errors = np.asarray(errors, np.int64)
    valid = counts > 0
    if not valid.any():
        return np.full(GAPMAX + 1, 60, np.uint8)
    gv = np.flatnonzero(valid)
    cp, ep = pava_decreasing(counts[valid], errors[valid])
    idx = np.searchsorted(gv, np.arange(GAPMAX + 1), side="right") - 1
    idx = np.clip(idx, 0, len(gv) - 1)
    q = np.array([phred(ep[j], cp[j]) for j in idx], np.uint8)
    return q


def cmd_fit(a):
    z = np.load(a.out)
    n = z["truth"].shape[0]
    sel = z["placed"] & z["truth"] & (z["gap"] >= 0)
    if a.cv:
        sel &= (np.arange(n) % 2 == 0)[:, None]
    g = z["gap"][sel]
    bad = ~z["ok"][sel]
    counts = np.bincount(g, minlength=GAPMAX + 1)
    errors = np.bincount(g[bad], minlength=GAPMAX + 1)
    qmap = build_qmap(counts, errors)
    np.savez(a.model, qmap=qmap, counts=counts, errors=errors,
             n_fit=int(sel.sum()), n_err=int(bad.sum()), cv=bool(a.cv))
    print(f"[fit] {int(sel.sum()):,} mates ({int(bad.sum()):,} misplaced), "
          f"{'even pairs' if a.cv else 'all pairs'}; gap MAPQ at "
          + ", ".join(f"{g}:{qmap[g]}" for g in (0, 1, 2, 5, 10, 20, 50, 100, GAPMAX)),
          file=sys.stderr)


def at_n(mq, bad, n):
    """Misplaced among the n most confident (the boundary MAPQ pro rata)."""
    order = np.argsort(-mq, kind="stable")
    m, b = mq[order], bad[order]
    if n >= len(m):
        return float(b.sum())
    q = m[n - 1]
    above = m > q
    at = m == q
    return float(b[above].sum()) + (n - int(above.sum())) * float(b[at].sum()) / float(at.sum())


def cmd_eval(a):
    z = np.load(a.out)
    qmap = np.load(a.model)["qmap"]
    n = z["truth"].shape[0]
    sel = z["placed"] & z["truth"]
    if a.half:
        sel &= (np.arange(n) % 2 == 1)[:, None]
    nt = int(sel.sum())
    mq = z["mapq"][sel]
    bad = ~z["ok"][sel]
    gap = z["gap"][sel]
    gq = np.where(gap < 0, 60, qmap[np.clip(gap, 0, GAPMAX)])
    rerank = np.minimum(mq, gq)
    shares = [float(s) for s in a.shares.split(",")]
    runs = {"baseline MAPQ": mq, "min(MAPQ, gap MAPQ)": rerank, "gap MAPQ alone": gq}
    print("| run | " + " | ".join(f"{100 * s:.2f} %" for s in shares) + " | MAPQ >= 40 |")
    print("|---|" + "---|" * (len(shares) + 1))
    for name, q in runs.items():
        row = [at_n(q, bad, int(round(s * nt))) / int(round(s * nt)) * 1e6 for s in shares]
        print(f"| {name} | " + " | ".join(f"{v:,.1f}" for v in row)
              + f" | {100 * np.mean(q >= 40):.2f} % |")
    print(f"\n{nt:,} mates with a truth ({'odd pairs held out' if a.half else 'all pairs'}); "
          "misplaced per million among each run's N most confident.")
    print("\n**At each run's own MAPQ >= 40**: confident %, misplaced per million confident.\n")
    print("| run | confident % | misplaced per M |")
    print("|---|---|---|")
    for name, q in runs.items():
        conf = q >= 40
        nc = int(conf.sum())
        print(f"| {name} | {100 * nc / nt:.2f} | {1e6 * int((conf & bad).sum()) / max(1, nc):,.1f} |")
    ce = bad & (mq >= 40)
    print(f"\nConfident (MAPQ >= 40) misplaced: {int(ce.sum()):,}; of them competitor gap < 5: "
          f"{int((ce & (gap >= 0) & (gap < 5)).sum()):,}; gap < 50: "
          f"{int((ce & (gap >= 0) & (gap < 50)).sum()):,}; no competitor score: "
          f"{int((ce & (gap < 0)).sum()):,}.")
    q = np.sort(rerank)[::-1]
    for s in shares:
        N = int(round(s * nt))
        if N < len(q):
            qb = q[N - 1]
            print(f"  N={100 * s:.2f}%: boundary MAPQ {qb}, tie group {int((rerank == qb).sum()):,}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    p = sp.add_parser("feat")
    p.add_argument("truth"), p.add_argument("certa"), p.add_argument("out")
    p = sp.add_parser("fit")
    p.add_argument("--cv", action="store_true")
    p.add_argument("out"), p.add_argument("model")
    p = sp.add_parser("eval")
    p.add_argument("--half", action="store_true")
    p.add_argument("--shares", default="0.9975,0.9929,0.98,0.95")
    p.add_argument("model"), p.add_argument("out")
    a = ap.parse_args()
    {"feat": cmd_feat, "fit": cmd_fit, "eval": cmd_eval}[a.cmd](a)


if __name__ == "__main__":
    main()
