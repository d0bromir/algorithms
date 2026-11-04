#!/usr/bin/env python3
"""A MAPQ informed by T2T-CHM13 (bench/VALIDATION_RESULTS.md, "A MAPQ from
T2T-CHM13"): a GRCh38 placement whose read an aligner places confidently on
T2T-CHM13 at another locus gets MAPQ 0.

  samtools view X.hs1.bam | truth_t2t.py place --chain hs1ToHg38.over.chain.gz truth.npz FAI t2t.npz
  samtools view X.bam ...  | truth_t2t.py place truth.npz FAI X.npz
  truth_t2t.py eval [--second RUN] truth.npz DIR > t2t.md

`place` keeps each mate's primary placement: contig, unclipped leftmost base,
strand and MAPQ, with chrY PAR placements moved to chrX. With --chain the
records are on hs1 and are lifted to GRCh38 by lift_sam.py's rule; a mapped
record that does not lift is kept as such (state 2).

`eval` reads DIR's placements of truth_headroom.py's runs (ivcert,
def.complete with def, bio.complete with bio, minibwa, bwa-mem2,
strobealign, minimap2; the fallback mates of the minibwa pipelines get
MAPQ / 4) and t2t.npz, and scores them at the locus (150 bp) before and
after the rule, at matched sensitivity: each run's N most confident mates
with a truth, N fixed beforehand (the confident shares, MAPQ >= 40, of
interval mode, --bio + minibwa and minibwa as they are).

The rule, fixed before the first evaluation: when the T2T record has
MAPQ >= 40 and lifts more than 150 bp from the GRCh38 placement (or to
another contig or strand), the GRCh38 MAPQ becomes 0. A grid varies the T2T
threshold, counts T2T records that do not lift as disagreeing, and raises
the MAPQ of a placement T2T confirms to T2T's. --second RUN is the control:
the same rule with a GRCh38 run of DIR as the second opinion, so that only
the reference differs from BWA-MEM2 on T2T.
"""
import argparse
import sys

import numpy as np

sys.path.insert(0, __import__("os").path.dirname(__file__))
from lift_sam import lift_record, load_chain  # noqa: E402
from truth_eval import read_sam  # noqa: E402
from truth_lift import CLS_TRUTH, Ref  # noqa: E402

LOCUS = 150


def cmd_place(a):
    n = np.load(a.truth)["cls"].shape[0]
    ref = Ref(a.fai)
    if not a.chain:
        r = read_sam(sys.stdin, n, ref.index)
        ref.par_to_x(r["rc"], r["left"])
        np.savez(a.out, placed=r["placed"], rc=r["rc"].astype(np.int16), left=r["left"].astype(np.int32),
                 strand=r["strand"], mapq=r["mapq"])
        return
    index = load_chain(a.chain)
    state = np.zeros((n, 2), np.int8)  # 0 unmapped, 1 lifted, 2 mapped but not lifted
    rc = np.full((n, 2), -1, np.int64)
    left = np.zeros((n, 2), np.int64)
    strand = np.zeros((n, 2), np.int8)
    mapq = np.zeros((n, 2), np.uint8)
    for line in sys.stdin:
        if line[0] == "@":
            continue
        f = line.split("\t", 6)
        flag = int(f[1])
        if flag & (4 | 256 | 2048):
            continue
        i, m = int(f[0][1:]), 1 if flag & 128 else 0
        mapq[i, m] = min(int(f[4]), 255)
        lr = lift_record(index, f[2], int(f[3]), f[5])
        if lr is None or lr[0] not in ref.index:
            state[i, m] = 2
            continue
        state[i, m] = 1
        rc[i, m], left[i, m], strand[i, m] = ref.index[lr[0]], lr[1], ((flag >> 4) & 1) ^ lr[2]
    ref.par_to_x(rc, left)
    print(f"[t2t place] mapped {int((state > 0).sum()):,}, lifted {int((state == 1).sum()):,}", file=sys.stderr)
    np.savez(a.out, placed=state == 1, state=state, rc=rc.astype(np.int16), left=left.astype(np.int32),
             strand=strand, mapq=mapq)


def cmd_eval(a):
    t = np.load(a.truth)
    truth = t["cls"] == CLS_TRUTH
    nt = int(truth.sum())

    def load(name):
        z = dict(np.load(f"{a.dir}/{name}.npz"))
        z["mapq"] = z["mapq"].astype(np.int64)
        z["left"] = z["left"].astype(np.int64)
        return z

    def at(r, s):  # r's placements within LOCUS of s's (arrays of mates)
        return (r["placed"] & s["placed"] & (r["rc"] == s["rc"]) & (r["strand"] == s["strand"])
                & (np.abs(r["left"] - s["left"]) <= LOCUS))

    tr = {"placed": truth, "rc": t["tc"], "strand": t["ts"], "left": t["tp"].astype(np.int64) - 1}
    CERT, DEF, BIO = "certa interval mode, --iv-mapq cert", "certa + minibwa (MAPQ / 4)", \
        "certa --bio + minibwa (MAPQ / 4)"
    runs = {CERT: load("ivcert"), DEF: load("def.complete"), BIO: load("bio.complete")}
    for name, alone in ((DEF, "def"), (BIO, "bio")):
        r = runs[name]
        fb = r["placed"] & ~load(alone)["placed"]
        r["mapq"] = np.where(fb, r["mapq"] // 4, r["mapq"])
    for x in ("minibwa", "bwa-mem2", "strobealign", "minimap2"):
        runs[x] = load(x)
    ok = {name: at(r, tr) for name, r in runs.items()}
    if a.second:  # the control: a GRCh38 run as the second opinion
        t2 = dict(runs[a.second], state=runs[a.second]["placed"].astype(np.int8))
    else:
        t2 = load("t2t")
    unl = t2["state"] == 2

    def rule(r, qt=40, unlifted=False, confirm=False):
        conf = t2["mapq"] >= qt
        same = at(r, t2)
        away = r["placed"] & conf & ((t2["placed"] & ~same) | (unl if unlifted else False))
        mq = np.where(away, 0, r["mapq"])
        if confirm:
            mq = np.where(same & conf, np.maximum(mq, t2["mapq"]), mq)
        return mq, away

    def at_n(mq, pl, sel, n):
        """Selected mates among the n most confident placed mates (boundary level pro rata)."""
        m, s = mq[pl], sel[pl]
        order = np.argsort(-m, kind="stable")
        m, s = m[order], s[order]
        if n >= len(m):
            return float(s.sum())
        q = m[n - 1]
        above = m > q
        lvl = m == q
        return float(s[above].sum() + (n - int(above.sum())) * (s[lvl].mean() if lvl.any() else 0.0))

    def conf(mq, pl, q=40):
        return int((pl & (mq >= q)).sum())

    refs = {"interval mode": conf(runs[CERT]["mapq"], runs[CERT]["placed"] & truth),
            "--bio + minibwa": conf(runs[BIO]["mapq"], runs[BIO]["placed"] & truth),
            "minibwa": conf(runs["minibwa"]["mapq"], runs["minibwa"]["placed"] & truth)}
    t2ok = at(t2, tr)
    src = a.second or "T2T-CHM13"
    print(f"Mates with a truth: {nt:,}. {src} records: {int((t2['state'] > 0).sum()):,} mapped, "
          f"{int(t2['placed'].sum()):,} on GRCh38; with a truth, {int((t2ok & truth).sum()):,} "
          f"({100 * (t2ok & truth).sum() / nt:.3f} %) at the locus.\n")
    hdr = " | ".join(f"N = {k}'s ({100 * n / nt:.2f} %)" for k, n in refs.items())
    print(f"**The rule** ({src} MAPQ ≥ 40, elsewhere on GRCh38: MAPQ 0). Misplaced per million among each run's N most "
          "confident mates with a truth, as is → with the rule; confident % (MAPQ ≥ 40); and the confident mates "
          "the rule sets to 0: how many, and the share misplaced.\n")
    print(f"| run | confident % | {hdr} | confident mates set to 0 | of them misplaced |")
    print("|---|---|" + "---|" * len(refs) + "---|---|")
    for name, r in runs.items():
        pl = r["placed"] & truth
        mis = pl & ~ok[name]
        mq, away = rule(r)
        cells = [f"{1e6 * at_n(r['mapq'], pl, mis, n) / n:,.0f} → {1e6 * at_n(mq, pl, mis, n) / n:,.0f}"
                 for n in refs.values()]
        hit = away & pl & (r["mapq"] >= 40)
        print(f"| {name} | {100 * conf(r['mapq'], pl) / nt:.2f} → {100 * conf(mq, pl) / nt:.2f} | "
              + " | ".join(cells) + f" | {int(hit.sum()):,} | {100 * (hit & mis).sum() / max(1, hit.sum()):.1f} % |")

    print("\n**Variants**, misplaced per million at N = interval mode's (and confident %): the T2T MAPQ threshold; "
          "T2T records that do not lift counted as disagreeing; and a placement T2T confirms raised to T2T's MAPQ.\n")
    variants = [(qt, u, c) for qt in (20, 40, 60) for u in (False, True) for c in (False, True)]
    print("| run | as is | " + " | ".join(f"≥ {qt}{', unlifted' if u else ''}{', confirm' if c else ''}"
                                         for qt, u, c in variants) + " |")
    print("|---|---|" + "---|" * len(variants))
    n = refs["interval mode"]
    for name in (CERT, "minibwa", "bwa-mem2", "strobealign"):
        r = runs[name]
        pl = r["placed"] & truth
        mis = pl & ~ok[name]
        cells = [f"{1e6 * at_n(r['mapq'], pl, mis, n) / n:,.0f} ({100 * conf(r['mapq'], pl) / nt:.2f} %)"]
        for qt, u, c in variants:
            mq, _ = rule(r, qt, u, c)
            cells.append(f"{1e6 * at_n(mq, pl, mis, n) / n:,.0f} ({100 * conf(mq, pl) / nt:.2f} %)")
        print(f"| {name} | " + " | ".join(cells) + " |")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("place")
    p.add_argument("--chain")
    p.add_argument("truth")
    p.add_argument("fai")
    p.add_argument("out")
    e = sub.add_parser("eval")
    e.add_argument("--second", help="a GRCh38 run (minibwa, bwa-mem2, ...) as the second opinion instead of T2T")
    e.add_argument("truth")
    e.add_argument("dir")
    a = ap.parse_args()
    {"place": cmd_place, "eval": cmd_eval}[a.cmd](a)


if __name__ == "__main__":
    main()
