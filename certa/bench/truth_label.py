#!/usr/bin/env python3
"""Test B truth: where real HG002 reads come from, read off the sample's own
genome (bench/VALIDATION_PLAN.md).

  truth_label.py MAT.sam PAT.sam N REF.fai blocks.npz STRAT_DIR BENCH_BED truth.npz
                 [--margin 2] [--check MAT_minibwa.sam PAT_minibwa.sam]

MAT.sam / PAT.sam: CERTA's default mode (edit certificate) on each haplotype
genome (maternal + chrY, paternal + chrX), reads named b<pair index>; N:
pairs. A mate's origin is known when its placement on one haplotype is
certified unique (S0/S1, one locus at the best distance d, none other within
d + margin) and either
  - the other haplotype's certified best distance (S0/S1/SR) is at least
    d + margin, or
  - the other haplotype's placement is certified unique too, at a distance
    within the margin, and both project to the same GRCh38 locus
    (homologous sequence, e.g. chrX/chrY, which both genomes carry).
The origin is projected to GRCh38 as in test A (truth_lift.py). --check: a
second aligner's placement on the origin's haplotype must agree (same contig
and strand, within 20 bp), or the mate is left unlabelled; mates it did not
place keep their label.
"""
import argparse
import re
import sys

import numpy as np

sys.path.insert(0, __import__("os").path.dirname(__file__))
from truth_eval import TIER, WINDOW, read_sam  # noqa: E402
from truth_lift import (CLS_EXCL, CLS_NOLIN, CLS_TRUTH, CLS_UNLAB, STRATA, Lift, Ref,  # noqa: E402
                        Strata, summary, truth_arrays)


def edit_tags(stream, n, contigs):
    """CERTA SAM on a haplotype: placement plus XE (d), the margin to the
    second locus ((XD if >= 0 else XR + 1) - XE) and XB (ties)."""
    tagre = re.compile(r"\t(XT:Z|XE:i|XR:i|XD:i|XB:i):([^\t\n]+)")
    out = dict(placed=np.zeros((n, 2), bool), rc=np.full((n, 2), -1, np.int64), left=np.zeros((n, 2), np.int64),
               span=np.zeros((n, 2), np.int64), strand=np.zeros((n, 2), np.int8), tier=np.full((n, 2), -1, np.int8),
               d=np.full((n, 2), -1, np.int16), margin=np.zeros((n, 2), np.int16), ties=np.zeros((n, 2), np.int16))
    for line in stream:
        if line[0] == "@":
            continue
        t = line.split("\t", 11)
        flag = int(t[1])
        if flag & (4 | 256 | 2048):
            continue
        c = contigs.get(t[2])
        if c is None:
            continue
        i, m = int(t[0][1:]), 1 if flag & 128 else 0
        tg = dict(tagre.findall("\t" + t[11])) if len(t) > 11 else {}
        span = sum(int(k) for k, o in re.findall(r"(\d+)([MD=X])", t[5]))
        out["placed"][i, m] = True
        out["rc"][i, m] = c
        out["left"][i, m] = int(t[3]) - 1
        out["span"][i, m] = span
        out["strand"][i, m] = 1 if flag & 16 else 0
        out["tier"][i, m] = TIER.get(tg.get("XT:Z", ""), -1)
        d = int(tg.get("XE:i", "-1"))
        out["d"][i, m] = d
        if d >= 0:
            xd, xr = int(tg.get("XD:i", "-1")), int(tg.get("XR:i", "-1"))
            out["margin"][i, m] = (xd if xd >= 0 else xr + 1) - d
            out["ties"][i, m] = int(tg.get("XB:i", "1"))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for x in ("mat", "pat", "n", "fai", "blocks", "strat", "bench", "out"):
        ap.add_argument(x)
    ap.add_argument("--margin", type=int, default=2)
    ap.add_argument("--check", nargs=2, default=None)
    a = ap.parse_args()
    n, M = int(a.n), a.margin
    ref, lift = Ref(a.fai), Lift(a.blocks)
    strata = Strata(ref, a.strat, a.bench)
    contigs = dict(lift.hidx)
    for extra in ("chrM", "chrEBV"):
        contigs.setdefault(extra, -2 - len(contigs))  # placed, but outside the projection
    h = []
    for path in (a.mat, a.pat):
        with open(path) as f:
            h.append(edit_tags(f, n, contigs))
    mat, pat = h
    unique = [((x["tier"] == TIER["S0"]) | (x["tier"] == TIER["S1"])) & (x["ties"] == 1) & (x["margin"] >= M) & (x["d"] >= 0)
              for x in h]
    cert = [x["placed"] & (x["d"] >= 0) & np.isin(x["tier"], [TIER["S0"], TIER["S1"], TIER["SR"]]) for x in h]
    only_mat = unique[0] & cert[1] & (pat["d"] >= mat["d"] + M)
    only_pat = unique[1] & cert[0] & (mat["d"] >= pat["d"] + M)
    both = unique[0] & unique[1] & (np.abs(mat["d"].astype(int) - pat["d"]) < M)

    t = truth_arrays(n)
    hap = np.full((n, 2), -1, np.int8)  # origin haplotype: 0 maternal, 1 paternal, 2 both
    counts = {}
    for m in range(2):
        # Projections of every unique placement on each haplotype.
        proj = []
        for x, u in zip(h, unique):
            rows = np.nonzero(u[:, m] & (x["rc"][:, m] >= 0))[0]
            hh = x["rc"][rows, m]
            s = x["left"][rows, m]
            v, rc, left, rs = lift.project(hh, s, s + x["span"][rows, m], x["strand"][rows, m])
            rc, left = rc.astype(np.int64), left.astype(np.int64)
            ref.par_to_x(rc, left)
            full = dict(valid=np.zeros(n, bool), rc=np.full(n, -1, np.int64), left=np.zeros(n, np.int64),
                        strand=np.zeros(n, np.int8))
            full["valid"][rows], full["rc"][rows], full["left"][rows], full["strand"][rows] = v, rc, left, rs
            proj.append(full)
        same = (proj[0]["valid"] & proj[1]["valid"] & (proj[0]["rc"] == proj[1]["rc"]) &
                (proj[0]["strand"] == proj[1]["strand"]) & (np.abs(proj[0]["left"] - proj[1]["left"]) <= WINDOW))
        homo = both[:, m] & same
        # chrM / chrEBV: both genomes carry them, and they are not projected.
        extra = both[:, m] & ((mat["rc"][:, m] < 0) | (pat["rc"][:, m] < 0))
        t["cls"][extra, m] = CLS_EXCL
        counts[f"mate {m + 1}: origin maternal only"] = int(only_mat[:, m].sum())
        counts[f"mate {m + 1}: origin paternal only"] = int(only_pat[:, m].sum())
        counts[f"mate {m + 1}: homologous, same GRCh38 locus"] = int(homo.sum())
        counts[f"mate {m + 1}: on chrM / chrEBV in both (excluded)"] = int(extra.sum())
        counts[f"mate {m + 1}: both haplotypes, different GRCh38 loci (unlabelled)"] = int((both[:, m] & ~same & ~extra).sum())
        for sel, k, code in ((only_mat[:, m], 0, 0), (only_pat[:, m], 1, 1), (homo, 0, 2)):
            hap[sel, m] = code
            x = h[k]
            on_extra = sel & (x["rc"][:, m] < 0)
            t["cls"][on_extra, m] = CLS_EXCL
            p = proj[k]
            ok = sel & p["valid"]
            t["cls"][sel & ~on_extra & ~p["valid"], m] = CLS_NOLIN
            t["cls"][ok, m] = CLS_TRUTH
            t["tc"][ok, m] = p["rc"][ok]
            t["tp"][ok, m] = p["left"][ok] + 1
            t["ts"][ok, m] = p["strand"][ok]
            t["st"][ok, m] = strata.bits(p["rc"][ok], p["left"][ok], x["span"][ok, m])

    if a.check:
        chk = []
        for path in a.check:
            with open(path) as f:
                chk.append(read_sam(f, n, contigs))
        for m in range(2):
            labelled = hap[:, m] >= 0
            agree = np.zeros(n, bool)
            placed = np.zeros(n, bool)
            for k in (0, 1):
                x, c = h[k], chk[k]
                on = labelled & ((hap[:, m] == k) | (hap[:, m] == 2))
                pl = on & c["placed"][:, m]
                ag = pl & (c["rc"][:, m] == x["rc"][:, m]) & (c["strand"][:, m] == x["strand"][:, m]) & \
                    (np.abs(c["left"][:, m] - x["left"][:, m]) <= WINDOW)
                placed |= pl
                agree |= ag
            dis = placed & ~agree
            counts[f"mate {m + 1}: second aligner agrees"] = int(agree.sum())
            counts[f"mate {m + 1}: second aligner disagrees (unlabelled)"] = int(dis.sum())
            counts[f"mate {m + 1}: second aligner did not place"] = int((labelled & ~placed).sum())
            t["cls"][dis, m] = CLS_UNLAB
            hap[dis, m] = -1
    for k, v in counts.items():
        print(f"  {k}: {v:,}", file=sys.stderr)
    np.savez(a.out, names=np.array([x for x, _, _ in STRATA]), ref=np.array(ref.names), hap=hap, **t)
    summary(t, file=sys.stderr)


if __name__ == "__main__":
    main()
