#!/usr/bin/env python3
"""truth_t2t_gate.py truth.npz SHARES BEFORE.npz AFTER.npz T2TF.npz: MAPQ
rules for PR mates that read T2T-CHM13 (bench/VALIDATION_RESULTS.md, "The
MAPQ of mates placed by the pair"). BEFORE and AFTER are truth_t2t.py place
files of the same placements with the MAPQ before and after the change (PR
mates' MAPQ from the best other pair found); T2TF is truth_mapq_pair.py feat
of a --t2t --feature-tags run of the same reads (tier, and where the T2T
placement lifts: XC). Misplaced mates per million among each rule's N most
confident, N = share x the mates with a truth (truth_rescue.py matched's
scoring)."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from truth_lift import CLS_TRUTH  # noqa: E402
from truth_mapq_pair import TIERS, at_n  # noqa: E402
from truth_rescue import at_locus  # noqa: E402

t = np.load(sys.argv[1])
shares = [float(s) for s in sys.argv[2].split(",")]
b, a, f = (np.load(p) for p in sys.argv[3:6])
truth = t["cls"] == CLS_TRUTH
ok = at_locus(b, t)
bad = truth & b["placed"] & ~ok
qb = np.where(b["placed"], b["mapq"], -1).astype(np.int64)
qa = np.where(a["placed"], a["mapq"], -1).astype(np.int64)
pr = f["tier"] == TIERS.index("PR")
here = (f["t_cls"] == 1) & (f["t_mq"] >= 40)
away = (f["t_cls"] == 2) & (f["t_mq"] >= 40)
rules = {
    "before": qb,
    "after (pair found)": qa,
    "after where T2T agrees (PR)": np.where(pr & here, qa, qb),
    "+ 0 where T2T is elsewhere": np.where(away, np.minimum(qb, 0), np.where(pr & here, qa, qb)),
    "before, 0 where T2T is elsewhere": np.where(away, np.minimum(qb, 0), qb),
}
n_truth = int(truth.sum())
print("| rule | " + " | ".join(f"{s:.2%}" for s in shares) + " |")
print("|---|" + "---|" * len(shares))
for name, q in rules.items():
    cells = []
    for s in shares:
        n = int(round(s * n_truth))
        nb = at_n(q[truth], bad[truth], n)
        cells.append(f"{nb / n * 1e6:,.1f}")
    print(f"| {name} | " + " | ".join(cells) + " |")
changed = truth & (qa != qb)
print(f"\n{int(changed.sum()):,} mates with a truth change MAPQ; PR {int((changed & pr).sum()):,}, "
      f"of which T2T agrees {int((changed & pr & here).sum()):,}, elsewhere {int((changed & pr & away).sum()):,}")
for name, m in (("T2T agrees", changed & pr & here), ("T2T elsewhere", changed & pr & away),
                ("neither", changed & pr & ~here & ~away)):
    print(f"  {name}: misplaced {int((m & bad).sum()):,} of {int(m.sum()):,}")
