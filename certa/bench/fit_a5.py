#!/usr/bin/env python3
"""A5: score the phantom parameter sweep on HG002 and freeze the best by
calibration error (sum |misplaced - expected| over MAPQ bands, expected =
Sum 10^(-MAPQ/10)). Also reports matched-N errors and MAPQ>=40 errors. Writes
a5_params.json (frozen delta/pi/pi_unique) and a5_fit.md.

  python3 fit_a5.py TRUTH.npz DIR_WITH_fit_d*_p*.npz
"""
import glob
import json
import os
import re
import sys

import numpy as np


def at_n(mq, bad, n):
    order = np.argsort(-mq, kind="stable")
    m, b = mq[order], bad[order]
    if n >= len(m):
        return float(b.sum())
    q = m[n - 1]
    above = m > q
    at = m == q
    return float(b[above].sum()) + (n - int(above.sum())) * float(b[at].sum()) / float(at.sum())


def score(truth_cls, z):
    sel = z["placed"] & (truth_cls == 0)
    nt = int(sel.sum())
    mq = z["mapq"][sel].astype(np.int64)
    bad = ~z["correct"][sel]
    cal = 0.0
    for lo in range(0, 60, 10):
        hi = 60 if lo == 50 else lo + 9
        band = sel & (z["mapq"] >= lo) & (z["mapq"] <= hi)
        mis = int((band & ~z["correct"]).sum())
        exp = float(np.power(10.0, -z["mapq"][band].astype(np.float64) / 10.0).sum())
        cal += abs(mis - exp)
    n98 = at_n(mq, bad, int(round(0.98 * nt))) / int(round(0.98 * nt)) * 1e6
    n95 = at_n(mq, bad, int(round(0.95 * nt))) / int(round(0.95 * nt)) * 1e6
    e40 = int((bad & (mq >= 40)).sum())
    return cal, n98, n95, e40, nt


def main():
    truth = np.load(sys.argv[1])["cls"]
    d = sys.argv[2]
    rows = []
    for f in sorted(glob.glob(os.path.join(d, "fit_d*_p*.npz"))):
        m = re.search(r"fit_d([0-9.]+)_p([0-9.]+)\.npz", os.path.basename(f))
        if not m:
            continue
        delta, pi = float(m.group(1)), float(m.group(2))
        cal, n98, n95, e40, nt = score(truth, np.load(f))
        rows.append((delta, pi, cal, n98, n95, e40, nt))
    if not rows:
        sys.exit("no fit_d*_p*.npz found")
    rows.sort(key=lambda r: r[2])  # best calibration first
    best = rows[0]
    lines = [
        "# A5 fit: phantom parameter sweep on HG002 (4M single-end subsample)",
        "",
        "Ranked by calibration error (Sum|misplaced-expected| over MAPQ bands);",
        "pi_unique=0 throughout (locus prior, A4). N=98%/95% are misplaced per",
        "million among the N most confident mates with a truth.",
        "",
        "| delta | pi | calib err | N=98% /M | N=95% /M | MAPQ>=40 err | mates |",
        "|---|---|---|---|---|---|---|",
    ]
    for delta, pi, cal, n98, n95, e40, nt in rows:
        lines.append(f"| {delta} | {pi} | {cal:,.0f} | {n98:,.1f} | {n95:,.1f} | {e40:,} | {nt:,} |")
    lines += ["", f"**Frozen: delta={best[0]}, pi={best[1]}, pi_unique=0** (calib err {best[2]:,.0f})."]
    with open(os.path.join(d, "a5_fit.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    with open(os.path.join(d, "a5_params.json"), "w") as fh:
        json.dump({"delta": best[0], "pi": best[1], "pi_unique": 0.0}, fh)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
