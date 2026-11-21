#!/usr/bin/env python3
"""low_levels.py: misplaced share (150 bp locus rule) per MAPQ value 0..29 for minibwa and BWA-MEM2, and per XP
(the pair model's Phred) 0..12 for CERTA's --mapq pair runs (placed by XP), on tests A and B and on HG005."""
import sys
import numpy as np
sys.path.insert(0, "/home/mpiuser/certa-lim/bench")
from truth_lift import CLS_TRUTH  # noqa: E402

def ok_of(r, t):
    return r["placed"] & (r["rc"] == t["tc"]) & (r["strand"] == t["ts"]) & \
        (np.abs(r["left"].astype(np.int64) - (t["tp"].astype(np.int64) - 1)) <= 150)

T = "/home/mpiuser/certa-truth"
runs = {
    "A": {"minibwa": "/home/mpiuser/certa-trace/runs/t2t/A/minibwa.npz", "bwa-mem2": "/home/mpiuser/certa-trace/runs/t2t/A/bwa-mem2.npz",
          "certa XP": "/home/mpiuser/certa-lim/work/A.fpair.npz"},
    "B": {"minibwa": "/home/mpiuser/certa-trace/runs/t2t/B/minibwa.npz", "bwa-mem2": "/home/mpiuser/certa-trace/runs/t2t/B/bwa-mem2.npz",
          "certa XP": "/home/mpiuser/certa-lim/work/B.fpair.npz"},
    "H": {"minibwa": "/home/mpiuser/certa-trace/runs/hg005/place/minibwa.npz", "bwa-mem2": "/home/mpiuser/certa-trace/runs/hg005/place/bwa-mem2.npz",
          "certa XP": "/home/mpiuser/certa-lim/refit/out/H.pair.xp.npz"},
}
for X, rr in runs.items():
    t = np.load(f"{T}/{X}/truth.npz")
    truth = t["cls"] == CLS_TRUTH
    print(f"== {X}")
    for name, path in rr.items():
        r = np.load(path)
        sel = r["placed"] & truth
        mq = r["mapq"][sel].astype(np.int64)
        bad = ~ok_of(r, t)[sel]
        top = 13 if name == "certa XP" else 30
        cells = []
        for q in range(top):
            n = int((mq == q).sum())
            cells.append(f"{q}:{100 * bad[mq == q].mean():.0f}%/{n}" if n else f"{q}:-")
        print(f"  {name:9} " + "  ".join(cells))
