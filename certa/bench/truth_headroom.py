#!/usr/bin/env python3
"""Aligners at matched sensitivity, and how much a reference richer than
GRCh38 gains (bench/VALIDATION_RESULTS.md, "Matched sensitivity"), from
truth_eval.py's npz files.

  truth_headroom.py TRUTH_DIR RUN_DIR [NAME=OUT.npz ...] > headroom.md

TRUTH_DIR holds truth.npz and phase 1's runs of the other aligners; RUN_DIR
the npz files of truth_complete.sh; NAME=OUT.npz are runs against another
reference, scored in GRCh38 coordinates (lift_sam.py, or a pangenome mapper's
projection onto GRCh38). All are scored with the same truth_eval.py --window.
Reported:
- the mates that no GRCh38 run places at the truth, and how many each new
  run places there, by MAPQ (with NAME=OUT.npz only);
- every run at a matched number of confident mates: its N most confident
  mates with a truth (ties at the boundary MAPQ counted pro rata), so a run
  is not credited for calling fewer mates confident. Avoidable: misplaced
  while some other GRCh38 run places the mate at the truth;
- every run at its own MAPQ >= 40.
The fallback mates of the minibwa pipelines get MAPQ / 4, as
drop_certified_mates.awk -v div=4 does.
"""
import sys

import numpy as np

W, D = sys.argv[1], sys.argv[2]
t = np.load(f"{W}/truth.npz")
truth = t["cls"] == 0
easy = truth & (((t["st"] >> 1) & 1) == 1)
nt = int(truth.sum())


def load(p):
    r = dict(np.load(p))
    r["mapq"] = r["mapq"].astype(np.int64)
    return r


CERT, BIO = "certa interval mode, --iv-mapq cert", "certa --bio + minibwa (MAPQ / 4)"
lin = {CERT: load(f"{D}/ivcert.npz"),
       "certa + minibwa (MAPQ / 4)": load(f"{D}/def.complete.npz"),
       BIO: load(f"{D}/bio.complete.npz")}
alone = {f: load(f"{D}/{f}.npz") for f in ("iv", "bio", "def")}
for name, f in (("certa + minibwa (MAPQ / 4)", "def"), (BIO, "bio")):
    r = lin[name]
    fb = r["placed"] & ~alone[f]["placed"]
    r["mapq"] = np.where(fb, r["mapq"] // 4, r["mapq"])
for a in ("minibwa", "bwa-mem2", "strobealign", "minimap2"):
    lin[a] = load(f"{W}/{a}.npz")
extra = {}
for spec in sys.argv[3:]:
    name, path = spec.split("=", 1)
    extra[name] = load(path)
runs = dict(lin, **extra)
grch38 = dict(lin, **{f"{f} alone": r for f, r in alone.items()})


def others_ok(name):
    """Mates some other GRCh38 run places at the truth."""
    ok = np.zeros(truth.shape, bool)
    for o, r in grch38.items():
        if o != name:
            ok |= r["correct"]
    return ok


nobody = truth & ~others_ok(None)
print(f"Mates with a truth: {nt:,}. Placed at the truth by no GRCh38 run: {int(nobody.sum()):,} "
      f"({100 * nobody.sum() / nt:.3f} %).\n")
if extra:
    print("**The mates no GRCh38 run places at the truth**: placed there by each new run, and with MAPQ ≥ 20 / 40.\n")
    print("| run | at the truth | MAPQ ≥ 20 | MAPQ ≥ 40 | misplaced with MAPQ ≥ 40 |")
    print("|---|---|---|---|---|")
    for name, r in extra.items():
        ok = nobody & r["correct"]
        bad = nobody & r["placed"] & ~r["correct"]
        print(f"| {name} | {int(ok.sum()):,} ({100 * ok.sum() / max(1, nobody.sum()):.1f} %) | "
              f"{int((ok & (r['mapq'] >= 20)).sum()):,} | {int((ok & (r['mapq'] >= 40)).sum()):,} | "
              f"{int((bad & (r['mapq'] >= 40)).sum()):,} |")
    print()


def at_n(r, n, sel):
    """Selected mates among the n most confident mates with a truth (boundary level pro rata)."""
    pl = r["placed"] & truth
    mq, s = r["mapq"][pl], sel[pl]
    order = np.argsort(-mq, kind="stable")
    mq, s = mq[order], s[order]
    if n >= len(mq):
        return float(s.sum())
    q = mq[n - 1]
    above = mq > q
    lvl = mq == q
    return float(s[above].sum() + (n - int(above.sum())) * (s[lvl].mean() if lvl.any() else 0.0))


def confident(r, q=40):
    return int((r["placed"] & truth & (r["mapq"] >= q)).sum())


refs = {f"{CERT} ≥ 40": confident(lin[CERT]), f"{BIO} ≥ 40": confident(lin[BIO]),
        "minibwa ≥ 40": confident(lin["minibwa"]), "minibwa ≥ 60": confident(lin["minibwa"], 60)}
refs.update({f"{name} ≥ 40": confident(r) for name, r in extra.items() if confident(r)})
print("**Matched confident mates**: misplaced / avoidable per million among each run's N most confident mates "
      "with a truth, N = the confident mates (MAPQ ≥ 40, or 60) of the run named.\n")
print("| run | at truth % | " + " | ".join(f"N = {k} ({100 * n / nt:.2f} %)" for k, n in refs.items()) + " |")
print("|---|---|" + "---|" * len(refs))
for name, r in runs.items():
    mis = truth & r["placed"] & ~r["correct"]
    avoid = mis & others_ok(name)
    cells = [f"{1e6 * at_n(r, n, mis) / n:,.0f} / {1e6 * at_n(r, n, avoid) / n:,.1f}" for n in refs.values()]
    print(f"| {name} | {100 * (r['correct'] & truth).sum() / nt:.3f} | " + " | ".join(cells) + " |")

print("\n**At each run's own MAPQ ≥ 40**: confident %, misplaced per million confident, and misplaced in easy "
      "regions at MAPQ ≥ 20.\n")
print("| run | confident % | misplaced per M | easy, MAPQ ≥ 20 |")
print("|---|---|---|---|")
for name, r in runs.items():
    n = confident(r)
    mis = truth & r["placed"] & ~r["correct"] & (r["mapq"] >= 40)
    e = easy & r["placed"] & ~r["correct"] & (r["mapq"] >= 20)
    print(f"| {name} | {100 * n / nt:.2f} | {1e6 * mis.sum() / max(1, n):,.0f} | {int(e.sum()):,} |")
