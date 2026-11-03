#!/usr/bin/env python3
"""MAPQ of the mates CERTA does not certify (bench/VALIDATION_RESULTS.md,
"MAPQ of the uncertified mates"), from truth_complete.sh's npz files.

  truth_caps.py TRUTH_DIR RUN_DIR > caps.md

TRUTH_DIR holds truth.npz and phase 1's competitor runs; RUN_DIR the npz
files of truth_complete.sh. A MAPQ rule changes MAPQ only, never a
placement, so each variant is scored by rewriting the npz arrays:
- fallback mates (placed by minibwa: in the complete pipeline, absent from
  CERTA's own SAM): capped at 39, or divided by 4;
- interval mode: XQ (its run with --iv-mapq cons) for tiers IV and MR,
  the point estimate for the rest; and its run with --iv-mapq cert (XQ for
  tiers R, IV and MR), checked against the other two runs;
- --bio with interval mode placing its fallback mates: composed from the
  two runs, as interval mode places a pair from the pair alone.
"""
import sys

import numpy as np

TIERS = ("S0", "S1", "SR", "S2", "SL", "PR", "R", "IV", "MR", "BM")
OTHERS = ("minibwa", "bwa-mem2", "strobealign", "minimap2")
W, D = sys.argv[1], sys.argv[2]
t = np.load(f"{W}/truth.npz")
truth = t["cls"] == 0
easy = truth & (((t["st"] >> 1) & 1) == 1)
nt = int(truth.sum())


def load(p):
    r = dict(np.load(p))
    r["mapq"] = r["mapq"].astype(np.int64)
    return r


runs = {"certa --bio + minibwa": load(f"{D}/bio.complete.npz"), "certa + minibwa": load(f"{D}/def.complete.npz"),
        "certa interval mode": load(f"{D}/iv.npz"), "certa --bio": load(f"{D}/bio.npz"), "certa": load(f"{D}/def.npz")}
runs.update({a: load(f"{W}/{a}.npz") for a in OTHERS})
ivc = load(f"{D}/ivc.npz")
cert = load(f"{D}/ivcert.npz")


def with_mapq(r, sel, mq):
    return dict(r, mapq=np.where(sel, mq, r["mapq"]))


def row(name, r, base=None):
    """base: the report run whose placements r shares, left out of 'another run'."""
    others = np.zeros(truth.shape, bool)
    for o, ro in runs.items():
        if o != base:
            others |= ro["correct"]
    pl, ok, mq = r["placed"] & truth, r["correct"], r["mapq"]
    bad = pl & ~ok
    cells = [f"{100 * ok.sum() / nt:.3f}", f"{100 * (pl & (mq >= 40)).sum() / nt:.2f}"]
    for q in (40, 60):
        n = (pl & (mq >= q)).sum()
        cells.append(f"{1e6 * (bad & (mq >= q)).sum() / n:,.0f} / {1e6 * (bad & others & (mq >= q)).sum() / n:,.1f}")
    cells.append(f"{int((easy & r['placed'] & ~ok & (mq >= 20)).sum()):,}")
    print(f"| {name} | " + " | ".join(cells) + " |")


print("**minibwa's MAPQ on the fallback mates** (with a truth): mates, misplaced, "
      "expected from MAPQ, and the Phred of the misplaced fraction.\n")
print("| MAPQ | `certa` + minibwa | `certa --bio` + minibwa |")
print("|---|---|---|")
fbs = {}
for m, c, a in (("def", runs["certa + minibwa"], runs["certa"]), ("bio", runs["certa --bio + minibwa"], runs["certa --bio"])):
    fbs[m] = c["placed"] & ~a["placed"]
for lo, hi in ((0, 0), (1, 9), (10, 19), (20, 29), (30, 39), (40, 49), (50, 59), (60, 60)):
    cells = []
    for m, c in (("def", runs["certa + minibwa"]), ("bio", runs["certa --bio + minibwa"])):
        b = fbs[m] & truth & (c["mapq"] >= lo) & (c["mapq"] <= hi)
        n, w = int(b.sum()), int((b & ~c["correct"]).sum())
        exp = float((10.0 ** (-c["mapq"][b] / 10.0)).sum())
        ph = f"{-10 * np.log10(w / n):.1f}" if w else "–"
        cells.append(f"{n:,}: {w:,} / {exp:,.1f}, Phred {ph}")
    print(f"| {lo}–{hi} | " + " | ".join(cells) + " |")

print("\n**The same, divided by 4** (MAPQ / 4, rounding down).\n")
print("| MAPQ / 4 | `certa` + minibwa | `certa --bio` + minibwa |")
print("|---|---|---|")
for lo, hi in ((0, 0), (1, 4), (5, 9), (10, 14), (15, 19)):
    cells = []
    for m, c in (("def", runs["certa + minibwa"]), ("bio", runs["certa --bio + minibwa"])):
        q = c["mapq"] // 4
        b = fbs[m] & truth & (q >= lo) & (q <= hi)
        n, w = int(b.sum()), int((b & ~c["correct"]).sum())
        exp = float((10.0 ** (-q[b] / 10.0)).sum())
        cells.append(f"{n:,}: {w:,} / {exp:,.1f}" + (f", {w / exp:.1f}×" if exp else ""))
    print(f"| {lo}–{hi} | " + " | ".join(cells) + " |")

print("\n**Interval mode by tier** (mates with a truth): placed, misplaced; with MAPQ ≥ 40; "
      "with XQ (`--iv-mapq cons`) ≥ 20.\n")
print("| tier | placed: misplaced | MAPQ ≥ 40 | XQ ≥ 20 |")
print("|---|---|---|---|")
iv = runs["certa interval mode"]
for k, nm in enumerate(TIERS):
    s = iv["placed"] & truth & (iv["tier"] == k)
    if s.any():
        bad = s & ~iv["correct"]
        cells = [f"{int(s.sum()):,}: {int(bad.sum()):,}"]
        cells += [f"{int((s & (q >= 40 if j == 0 else q >= 20)).sum()):,}: {int((bad & (q >= 40 if j == 0 else q >= 20)).sum()):,}"
                  for j, q in enumerate((iv["mapq"], ivc["mapq"]))]
        print(f"| {nm} | " + " | ".join(cells) + " |")

unc = np.isin(iv["tier"], [TIERS.index(x) for x in ("R", "IV", "MR")])
same = all((cert[k] == iv[k]).all() for k in ("placed", "correct", "tier"))
same &= bool((cert["mapq"] == np.where(unc, ivc["mapq"], iv["mapq"])).all())
print(f"\n`--iv-mapq cert` against the runs with the default MAPQ and `--iv-mapq cons`: the same placements, "
      f"and XQ for tiers R, IV and MR, the point estimate otherwise: {'yes' if same else 'NO'}.\n")

print("\n**MAPQ rules** (mates with a truth). Confident: MAPQ ≥ 40. Avoidable: misplaced while "
      "another run of the report placed the mate at the truth.\n")
print("| run | at truth % | confident % | misplaced / avoidable per M, MAPQ ≥ 40 | ≥ 60 | easy regions, misplaced at MAPQ ≥ 20 |")
print("|---|---|---|---|---|---|")
for a in OTHERS:
    row(a, runs[a], a)
for name, m in (("certa + minibwa", "def"), ("certa --bio + minibwa", "bio")):
    c, fb = runs[name], fbs[m]
    row(f"`{name}`", c, name)
    row(f"`{name}`, fallback MAPQ ≤ 39", with_mapq(c, fb, np.minimum(c["mapq"], 39)), name)
    row(f"`{name}`, fallback MAPQ / 4", with_mapq(c, fb, c["mapq"] // 4), name)
iname = "certa interval mode"
row(f"`{iname}`", iv, iname)
row(f"`{iname}`, `--iv-mapq cons`", ivc, iname)
row(f"`{iname}`, XQ for IV, MR", with_mapq(iv, np.isin(iv["tier"], [TIERS.index("IV"), TIERS.index("MR")]),
                                          ivc["mapq"]), iname)
row(f"`{iname}`, `--iv-mapq cert`", cert, iname)
b = runs["certa --bio"]
fbm = ~b["placed"]
comp = {k: np.where(fbm, iv[k], b[k]) for k in ("placed", "correct", "mapq")}
row("`certa --bio` + interval mode (composed)", comp)
row("`certa --bio` + interval mode, XQ", with_mapq(comp, fbm, ivc["mapq"]))
row("`certa --bio` + interval mode, MAPQ / 4", with_mapq(comp, fbm, iv["mapq"] // 4))
for name in ("certa --bio", "certa"):
    row(f"*`{name}` alone*", runs[name], name)
