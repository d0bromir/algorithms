#!/usr/bin/env python3
"""Scores aligners' GRCh38 placements against truth (bench/VALIDATION_PLAN.md).

  samtools view ALN.bam | truth_eval.py score truth.npz REF.fai OUT.npz
  truth_eval.py report truth.npz NAME=OUT.npz ... [--bio NAME] [--default NAME]
                [--scores NAME=bio_score.tsv ...] > report.md

`score` reads one aligner's primary records (read names <letter><pair index>,
as truth_lift.py / truth_label.py number them) and keeps, per mate: placed,
at truth (right contig and strand, unclipped leftmost base within 20 bp;
chrY PAR placements moved to chrX first), MAPQ, CERTA's tier, and for --bio
records YA and the certified optimum p1 (YP, or YS when the mate was moved
to another locus: loci[0] is the best one). `report` prints Markdown tables:
placement, misplaced reads by MAPQ, strata, calibration, --bio against the
default mode, reads with no linear truth, and the soundness check (bio_score
penalties of other aligners' alignments against --bio's p1).
"""
import argparse
import re
import sys

import numpy as np

sys.path.insert(0, __import__("os").path.dirname(__file__))
from truth_lift import CLS_NOLIN, CLS_TRUTH, STRATA, Ref  # noqa: E402

WINDOW = 20
TIERS = ("S0", "S1", "SR", "S2", "SL", "PR", "R", "IV", "MR", "BM")
TIER = {t: i for i, t in enumerate(TIERS)}
LEAD_CLIP = re.compile(r"^(\d+)[SH]")
TAG = re.compile(r"\t(XT:Z|YA:i|YP:i|YS:i):([^\t\n]+)")


def read_sam(stream, n, contigs, tags=False):
    """Per-mate arrays (n pairs x 2) of an aligner's primary records."""
    out = dict(placed=np.zeros((n, 2), bool), rc=np.full((n, 2), -1, np.int64), left=np.zeros((n, 2), np.int64),
               strand=np.zeros((n, 2), np.int8), mapq=np.zeros((n, 2), np.uint8), span=np.zeros((n, 2), np.int32))
    if tags:
        out.update(tier=np.full((n, 2), -1, np.int8), ya=np.full((n, 2), -1, np.int16),
                   p1=np.full((n, 2), -1, np.int32))
    records = 0
    for line in stream:
        if line[0] == "@":
            continue
        t = line.split("\t", 11)
        flag = int(t[1])
        if flag & (256 | 2048):
            continue
        i = int(t[0][1:])
        m = 1 if flag & 128 else 0
        records += 1
        if flag & 4:
            continue
        c = contigs.get(t[2])
        if c is None:
            continue
        cig = t[5]
        lc = LEAD_CLIP.match(cig)
        pos = int(t[3]) - 1 - (int(lc.group(1)) if lc else 0)  # unclipped leftmost, 0-based
        out["placed"][i, m] = True
        out["rc"][i, m] = c
        out["left"][i, m] = pos
        out["strand"][i, m] = 1 if flag & 16 else 0
        out["mapq"][i, m] = min(int(t[4]), 255)
        out["span"][i, m] = sum(int(k) for k, o in re.findall(r"(\d+)([MDN=XSH])", cig))
        if tags and len(t) > 11:
            tg = dict(TAG.findall("\t" + t[11]))
            out["tier"][i, m] = TIER.get(tg.get("XT:Z", ""), -1)
            if "YA:i" in tg:
                out["ya"][i, m] = int(tg["YA:i"])
            if "YP:i" in tg:
                yp, ys = int(tg["YP:i"]), int(tg.get("YS:i", "-1"))
                out["p1"][i, m] = min(yp, ys) if ys >= 0 else yp
    print(f"[read_sam] {records} primary records", file=sys.stderr)
    return out


def cmd_score(a):
    t = dict(np.load(a.truth))
    n = t["cls"].shape[0]
    ref = Ref(a.fai)
    r = read_sam(sys.stdin, n, ref.index, tags=True)
    ref.par_to_x(r["rc"], r["left"])
    correct = r["placed"] & (t["cls"] == CLS_TRUTH) & (r["rc"] == t["tc"]) & (r["strand"] == t["ts"])
    correct &= np.abs(r["left"] + 1 - t["tp"]) <= WINDOW
    np.savez(a.out, placed=r["placed"], correct=correct, mapq=r["mapq"], tier=r["tier"], ya=r["ya"], p1=r["p1"])


def pct(k, n):
    return f"{100.0 * k / n:.3f}" if n else "–"


def cmd_report(a):
    t = dict(np.load(a.truth))
    cls, st = t["cls"], t["st"]
    truth = cls == CLS_TRUTH
    runs = {}
    for spec in a.runs:
        name, path = spec.split("=", 1)
        runs[name] = dict(np.load(path))
    nt = int(truth.sum())
    print(f"Mates: {cls.size:,}; with a GRCh38 truth {nt:,}; no linear truth {int((cls == CLS_NOLIN).sum()):,}; "
          f"other (excluded or not labelled) {int((cls > CLS_NOLIN).sum()):,}.\n")

    print("**Placement (mates with a truth).** Misplaced: in the SAM file, not at the truth.\n")
    print("| Aligner | in SAM % | at truth % | misplaced | MAPQ ≥ 1 | ≥ 20 | ≥ 40 | ≥ 60 |")
    print("|---|---|---|---|---|---|---|---|")
    for name, r in runs.items():
        pl, ok, mq = r["placed"] & truth, r["correct"], r["mapq"]
        bad = pl & ~ok
        row = [pct(int(pl.sum()), nt), pct(int(ok.sum()), nt), f"{int(bad.sum()):,}"]
        row += [f"{int((bad & (mq >= q)).sum()):,}" for q in (1, 20, 40, 60)]
        print(f"| {name} | " + " | ".join(row) + " |")

    print("\n**Avoidable errors**: misplaced mates that some other aligner placed at the truth "
          "(the rest no aligner placed there), by MAPQ.\n")
    print("| Aligner | misplaced | MAPQ ≥ 1 | ≥ 20 | ≥ 40 | ≥ 60 |")
    print("|---|---|---|---|---|---|")
    for name, r in runs.items():
        others = np.zeros(cls.shape, bool)
        for o, ro in runs.items():
            if o != name:
                others |= ro["correct"]
        bad = r["placed"] & truth & ~r["correct"] & others
        mq = r["mapq"]
        print(f"| {name} | {int(bad.sum()):,} | " + " | ".join(f"{int((bad & (mq >= q)).sum()):,}" for q in (1, 20, 40, 60)) + " |")

    print("\n**Placed against misplaced as the MAPQ threshold falls** (mates with a truth; placed / misplaced).\n")
    qs = (60, 50, 40, 30, 20, 10, 1, 0)
    print("| Aligner | " + " | ".join(f"MAPQ ≥ {q}" for q in qs) + " |")
    print("|---|" + "---|" * len(qs))
    for name, r in runs.items():
        pl, ok, mq = r["placed"] & truth, r["correct"], r["mapq"]
        cells = [f"{int((pl & (mq >= q)).sum()):,} / {int((pl & ~ok & (mq >= q)).sum()):,}" for q in qs]
        print(f"| {name} | " + " | ".join(cells) + " |")

    print("\n**Strata** (mates with a truth in each; at truth % / misplaced with MAPQ ≥ 20).\n")
    print("| Stratum | mates | " + " | ".join(runs) + " |")
    print("|---|---|" + "---|" * len(runs))
    for bit, (sname, _, _) in enumerate(STRATA):
        sm = truth & (((st >> bit) & 1) == 1)
        k = int(sm.sum())
        if not k:
            continue
        cells = []
        for r in runs.values():
            ok = r["correct"] & sm
            bad = r["placed"] & sm & ~r["correct"] & (r["mapq"] >= 20)
            cells.append(f"{pct(int(ok.sum()), k)} / {int(bad.sum()):,}")
        print(f"| {sname} | {k:,} | " + " | ".join(cells) + " |")

    print("\n**Calibration** (mates with a truth): misplaced / expected from MAPQ (Σ 10^(−MAPQ/10)).\n")
    bands = ((0, 9), (10, 19), (20, 29), (30, 39), (40, 49), (50, 60))
    print("| Aligner | " + " | ".join(f"MAPQ {lo}–{hi}" for lo, hi in bands) + " |")
    print("|---|" + "---|" * len(bands))
    for name, r in runs.items():
        pl, ok, mq = r["placed"] & truth, r["correct"], r["mapq"].astype(np.int64)
        cells = []
        for lo, hi in bands:
            b = pl & (mq >= lo) & (mq <= hi)
            exp = float((10.0 ** (-mq[b] / 10.0)).sum())
            cells.append(f"{int(b.sum()):,}: {int((b & ~ok).sum()):,} / {exp:,.1f}")
        print(f"| {name} | " + " | ".join(cells) + " |")

    if a.bio in runs and a.default in runs:
        b, d = runs[a.bio], runs[a.default]
        bo, do = b["correct"], d["correct"]
        bw, dw = b["placed"] & truth & ~bo, d["placed"] & truth & ~do
        print(f"\n**{a.bio} against {a.default}** (mates with a truth).\n")
        print("| | mates |")
        print("|---|---|")
        print(f"| at truth in both | {int((bo & do).sum()):,} |")
        print(f"| moved to the truth ({a.default} wrong or absent) | {int((bo & ~do & truth).sum()):,} |")
        print(f"| moved away ({a.bio} wrong or absent) | {int((do & ~bo).sum()):,} |")
        print(f"| ({a.bio} wrong, {a.default} at truth) | {int((do & bw).sum()):,} |")
        print(f"| wrong in both | {int((bw & dw).sum()):,} |")
        ya = b["ya"]
        for q in (20, 40):
            w = bw & (b["mapq"] >= q)
            print(f"| {a.bio} wrong with MAPQ ≥ {q}: of them YA > 0 (tandem offset) | {int(w.sum()):,}: {int((w & (ya > 0)).sum()):,} |")

    nolin = cls == CLS_NOLIN
    print(f"\n**No linear truth** ({int(nolin.sum()):,} mates: origin not 1:1 in GRCh38): placed, and with MAPQ ≥ 20.\n")
    print("| Aligner | placed | MAPQ ≥ 20 |")
    print("|---|---|---|")
    for name, r in runs.items():
        pl = r["placed"] & nolin
        print(f"| {name} | {pct(int(pl.sum()), int(nolin.sum()))} % | {int((pl & (r['mapq'] >= 20)).sum()):,} |")

    if a.scores and a.bio in runs:
        b = runs[a.bio]
        p1 = b["p1"].astype(np.int64)
        cert = (b["tier"] == TIER["BM"]) & (p1 >= 0)
        n = cls.shape[0]
        print(f"\n**Soundness**: other aligners' end-to-end alignments scored under the {a.bio} model "
              f"against its certified optimum p1 ({int(cert.sum()):,} certified mates). Alignments with an "
              f"insertion next to a deletion are outside the model (THEORY §12) and counted apart.\n")
        print("| Aligner | scored | on certified mates | cheaper than p1 (violations) | equal to p1 | above p1 "
              "| I next to D (outside the model) |")
        print("|---|---|---|---|---|---|---|")
        # The cheapest alignment at the truth any other aligner reported, per mate.
        at_truth = np.full((n, 2), np.iinfo(np.int64).max, np.int64)
        for spec in a.scores:
            name, path = spec.split("=", 1)
            pen = np.full((n, 2), -1, np.int64)
            with open(path) as f:
                for line in f:
                    q, mate, v = line.split("\t")
                    pen[int(q[1:]), max(int(mate) - 1, 0)] = int(v)
            sc = pen >= 0
            on = sc & cert
            print(f"| {name} | {int(sc.sum()):,} | {int(on.sum()):,} | {int((on & (pen < p1)).sum()):,} | "
                  f"{int((on & (pen == p1)).sum()):,} | {int((on & (pen > p1)).sum()):,} | {int((pen == -2).sum()):,} |")
            if name in runs:
                ok = sc & runs[name]["correct"]
                at_truth = np.where(ok, np.minimum(at_truth, pen), at_truth)

        wrong = b["placed"] & truth & ~b["correct"] & cert
        print(f"\n**Why {a.bio} misplaces**: for its misplaced, certified mates, the cheapest alignment at the "
              f"truth that another aligner reported, against {a.bio}'s optimum p1 (model penalty, Phred).\n")
        bins = (("no other aligner at the truth", None, None), ("equal to p1 (a tie)", 0, 0), ("1–9 above p1", 1, 9),
                ("10–29 above", 10, 29), ("30–59 above", 30, 59), ("60 or more above", 60, None))
        print("| MAPQ of the misplaced mate | " + " | ".join(x[0] for x in bins) + " |")
        print("|---|" + "---|" * len(bins))
        found = at_truth < np.iinfo(np.int64).max
        diff = at_truth - p1
        for lo, hi in ((0, 19), (20, 39), (40, 255)):
            w = wrong & (b["mapq"] >= lo) & (b["mapq"] <= hi)
            cells = []
            for _label, x0, x1 in bins:
                if x0 is None:
                    sel = w & ~found
                else:
                    sel = w & found & (diff >= x0) & ((diff <= x1) if x1 is not None else True)
                cells.append(f"{int(sel.sum()):,}")
            print(f"| {lo}–{hi if hi < 255 else 60} | " + " | ".join(cells) + " |")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("score")
    s.add_argument("truth")
    s.add_argument("fai")
    s.add_argument("out")
    r = sub.add_parser("report")
    r.add_argument("truth")
    r.add_argument("runs", nargs="+")
    r.add_argument("--bio", default="certa --bio")
    r.add_argument("--default", default="certa")
    r.add_argument("--scores", nargs="*", default=[])
    a = ap.parse_args()
    {"score": cmd_score, "report": cmd_report}[a.cmd](a)


if __name__ == "__main__":
    main()
