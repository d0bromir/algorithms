#!/usr/bin/env python3
"""How --iv-rescue local changes placements against the truth
(bench/VALIDATION_RESULTS.md, "A local mate rescue"). The inputs are
truth_t2t.py place files of two runs of the same reads.

  truth_rescue.py matched truth.npz SHARES NAME=RUN.npz ...
  truth_rescue.py moved truth.npz BEFORE.npz AFTER.npz

`matched`: misplaced mates (at the locus, 150 bp) among each run's N most
confident, N = share x the mates with a truth: per million to one decimal,
and the count (truth_mapq_pair.py eval's at_n).
`moved`: the mates placed differently, those with a truth at the locus
before and after; the confident mates with a truth (MAPQ >= 20, >= 40) and
their misplaced; and the mates with a truth whose MAPQ >= 40 status
changes, by whether they moved and are at the locus after.
"""
import argparse
import os
import sys
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from truth_lift import CLS_TRUTH  # noqa: E402
from truth_mapq_pair import at_n  # noqa: E402


def at_locus(r, t):
    return r["placed"] & (r["rc"] == t["tc"]) & (r["strand"] == t["ts"]) & \
        (np.abs(r["left"].astype(np.int64) - (t["tp"].astype(np.int64) - 1)) <= 150)


def cmd_matched(a):
    t = np.load(a.truth)
    truth = t["cls"] == CLS_TRUTH
    nt = int(truth.sum())
    shares = [float(s) for s in a.shares.split(",")]
    print("| run | " + " | ".join(f"{100 * s:.2f} %" for s in shares) + " |")
    print("|---|" + "---|" * len(shares))
    for spec in a.runs:
        name, p = spec.split("=", 1)
        r = np.load(p)
        sel = r["placed"] & truth
        mq, bad = r["mapq"].astype(np.int64)[sel], ~at_locus(r, t)[sel]
        cells = []
        for s in shares:
            n = int(round(s * nt))
            b = at_n(mq, bad, n)
            cells.append(f"{b / n * 1e6:,.1f} ({b:,.0f})")
        print(f"| {name} | " + " | ".join(cells) + " |")
    print(f"\n{nt:,} mates with a truth; misplaced per million among each run's N most confident (count).")


def cmd_moved(a):
    t, x, y = np.load(a.truth), np.load(a.before), np.load(a.after)
    truth = t["cls"] == CLS_TRUTH
    ox, oy = at_locus(x, t), at_locus(y, t)
    d = (x["placed"] != y["placed"]) | (x["rc"] != y["rc"]) | (x["strand"] != y["strand"]) | (x["left"] != y["left"])
    dt = d & truth
    print(f"mates placed differently: {int(d.sum()):,}; with a truth {int(dt.sum()):,}: at the locus before "
          f"{int((dt & ox).sum()):,}, after {int((dt & oy).sum()):,} (wrong to right {int((dt & ~ox & oy).sum()):,}, "
          f"right to wrong {int((dt & ox & ~oy).sum()):,}); unplaced before {int((dt & ~x['placed']).sum()):,}, "
          f"after {int((dt & ~y['placed']).sum()):,}")
    for q in (20, 40):
        cx = truth & x["placed"] & (x["mapq"] >= q)
        cy = truth & y["placed"] & (y["mapq"] >= q)
        print(f"MAPQ >= {q}, with a truth: {int(cx.sum()):,} -> {int(cy.sum()):,} mates; misplaced "
              f"{int((cx & ~ox).sum()):,} -> {int((cy & ~oy).sum()):,}")
    cx = x["placed"] & (x["mapq"] >= 40)
    cy = y["placed"] & (y["mapq"] >= 40)
    c = Counter()
    for p, m in np.argwhere(truth & (cx != cy)):
        c[("becomes" if cy[p, m] else "stops being", "moved" if d[p, m] else "same place",
           "at the locus" if oy[p, m] else "elsewhere")] += 1
    print("MAPQ >= 40 status changes (mates with a truth):")
    for k, v in sorted(c.items(), key=lambda kv: -kv[1]):
        print(f"  {v:,} {k[0]} confident, {k[1]}, {k[2]}")


def main():
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    p = sp.add_parser("matched")
    p.add_argument("truth"), p.add_argument("shares"), p.add_argument("runs", nargs="+")
    p = sp.add_parser("moved")
    p.add_argument("truth"), p.add_argument("before"), p.add_argument("after")
    a = ap.parse_args()
    {"matched": cmd_matched, "moved": cmd_moved}[a.cmd](a)


if __name__ == "__main__":
    main()
