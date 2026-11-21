#!/usr/bin/env python3
"""t2t_skip.py X.f4.npz...: which read pairs a lean --t2t could leave without a T2T-CHM13 lookup. A pair is
skipped when both mates are of tier S0 or S1 with a gap (the certified distance from the best alignment to any
other locus: XD, else XO, minus XE, clipped to 10; 11 when CERTA records neither, as the pair model reads it) of
at least G. Per G: the mapped mates in skipped pairs, and among them the T2T disagreements (XC class 2 or 3 with
T2T MAPQ >= 20) and the misplaced mates (a truth, more than 150 bp from it), with and without such a
disagreement, as shares of all mapped mates of each kind."""
import sys

import numpy as np

for path in sys.argv[1:]:
    z = np.load(path)
    mapped = (z["flag"] & 4) == 0
    tier = z["tier"].astype(np.int32)
    second = np.where(z["XD"] >= 0, z["XD"], z["XO"])
    gap = np.where(second >= 0, np.clip(second - z["XE"], 0, 10), 11)
    s01 = mapped & ((tier == 0) | (tier == 1))
    dis = mapped & np.isin(z["t_cls"], (2, 3)) & (z["t_mq"] >= 20)
    mis = mapped & z["truth"] & z["placed"] & ~z["ok"]
    n_m, n_d, n_x, n_xd = mapped.sum(), dis.sum(), mis.sum(), (mis & dis).sum()
    print(f"\n{path}: {n_m:,} mapped mates; {n_d:,} T2T disagreements; {n_x:,} misplaced, {n_xd:,} of them with a"
          f" T2T disagreement; S0/S1 {s01.sum() / n_m:.1%} of the mapped")
    print("| skip pairs of S0/S1 mates with gap >= | mates skipped | disagreements lost | misplaced skipped"
          " | misplaced with a disagreement lost |")
    print("|---|---|---|---|---|")
    for g in (0, 3, 6, 11):
        s = s01 & (gap >= g)
        sk = mapped & (s[:, 0] & s[:, 1])[:, None]
        print(f"| {g} | {sk.sum() / n_m:.2%} | {(sk & dis).sum() / max(n_d, 1):.2%} ({(sk & dis).sum():,})"
              f" | {(sk & mis).sum() / max(n_x, 1):.2%} ({(sk & mis).sum():,})"
              f" | {(sk & mis & dis).sum() / max(n_xd, 1):.2%} ({(sk & mis & dis).sum():,}) |")
    print("\ngap of S0/S1 mates: " + ", ".join(f"{g}: {(s01 & (gap == g)).sum() / max(s01.sum(), 1):.1%}"
                                               for g in range(12)))
