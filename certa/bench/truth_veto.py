#!/usr/bin/env python3
"""The T2T-CHM13 veto of truth_t2t.py, applied to one run's placements
(bench/VALIDATION_RESULTS.md, "A MAPQ from T2T-CHM13"): a placement gets
MAPQ 0 when the T2T placement has MAPQ >= 40 and lifts more than 150 bp away,
or to another contig or strand. Placements are unchanged.

  truth_veto.py RUN.npz T2T.npz OUT.npz

RUN.npz and T2T.npz are truth_t2t.py place files (T2T.npz with --chain).
"""
import sys

import numpy as np


def main(a):
    r = dict(np.load(a[0]))
    t2 = np.load(a[1])
    left = r["left"].astype(np.int64)
    same = r["placed"] & t2["placed"] & (r["rc"] == t2["rc"]) & (r["strand"] == t2["strand"]) & \
        (np.abs(left - t2["left"].astype(np.int64)) <= 150)
    away = r["placed"] & (t2["mapq"] >= 40) & t2["placed"] & ~same
    r["mapq"] = np.where(away, 0, r["mapq"])
    np.savez(a[2], **r)
    print(f"[truth_veto] {int(away.sum()):,} placements set to MAPQ 0", file=sys.stderr)


if __name__ == "__main__":
    main(sys.argv[1:])
