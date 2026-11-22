#!/usr/bin/env python3
"""B1 ceiling (offline): of CERTA's confident GRCh38 misplacements, how many
does a CHM13 lift-and-commit place at the truth (RESCUE), abstain on, or still
misplace? This quantifies what reference content (Pillar B) adds over the
phantom abstention (Pillar A): the lift can PLACE missing-copy reads correctly,
not just abstain.

Uses existing scored npz (no mapper change, no new data):
  certa.npz  CERTA on GRCh38 (default mode)
  t2t.npz    BWA-MEM2 on T2T-CHM13 (hs1), lifted to GRCh38 (headroom/t2t.sh)
  truth.npz  the assembly-projected truth

  python3 b1_ceiling.py CERTA_TRUTH_DIR HEADROOM_DIR
"""
import sys

import numpy as np


def breakdown(name, mask, t2):
    n = int(mask.sum())
    if n == 0:
        print(f"{name}: 0")
        return
    resc = int((mask & t2["correct"]).sum())
    abst = int((mask & ~t2["placed"]).sum())
    wrong = int((mask & t2["placed"] & ~t2["correct"]).sum())
    print(f"{name}: {n:,}")
    print(f"  CHM13 lift at truth (RESCUE):      {resc:,} ({100 * resc / n:.1f}%)")
    print(f"  CHM13 lift abstains (not lifted):  {abst:,} ({100 * abst / n:.1f}%)")
    print(f"  CHM13 lift still misplaced:        {wrong:,} ({100 * wrong / n:.1f}%)")


def main():
    W, H = sys.argv[1], sys.argv[2]
    t = np.load(f"{W}/truth.npz")
    truth = t["cls"] == 0
    nolin = t["cls"] == 1  # missing-copy class: no GRCh38 truth
    st = t["st"]
    segdup = ((st >> 4) & 1) == 1
    c = np.load(f"{W}/certa.npz")
    t2 = np.load(f"{H}/t2t.npz")

    print("=== CERTA confident (MAPQ>=40) GRCh38 misplacements ===")
    ce = truth & c["placed"] & ~c["correct"] & (c["mapq"] >= 40)
    breakdown("all confident errors", ce, t2)
    breakdown("  in segdup", ce & segdup, t2)
    breakdown("  outside segdup", ce & ~segdup, t2)

    print("\n=== Missing-copy class (no linear GRCh38 truth) ===")
    # Any confident GRCh38 placement of these is an error; abstention is correct.
    nl_conf = nolin & c["placed"] & (c["mapq"] >= 40)
    print(f"mates with no GRCh38 truth: {int(nolin.sum()):,}")
    print(f"  CERTA places confidently (error): {int(nl_conf.sum()):,}")
    print(f"  of those, CHM13 lift abstains (correct): {int((nl_conf & ~t2['placed']).sum()):,}")
    print(f"  CHM13 lift placed (any): {int((nolin & t2['placed']).sum()):,}")

    print("\n=== Overall placement (mates with a GRCh38 truth) ===")
    for nm, z in (("CERTA GRCh38", c), ("CHM13 lift", t2)):
        pl = truth & z["placed"]
        corr = truth & z["correct"]
        print(f"{nm}: placed {int(pl.sum()):,}, at truth {int(corr.sum()):,} "
              f"({100 * corr.sum() / truth.sum():.3f}%), confident errors "
              f"{int((pl & ~z['correct'] & (z['mapq'] >= 40)).sum()):,}")


if __name__ == "__main__":
    main()
