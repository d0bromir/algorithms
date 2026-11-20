#!/usr/bin/env python3
"""Why a mate is not at its origin: the algorithm, a tie, or the reference
(bench/VALIDATION_RESULTS.md, "How far from the limit").

  truth_oracle.py jobs truth.npz PREFIX NAME=RUN.npz ...
  tools/truth_oracle REF.fa R1.fq R2.fq PREFIX.jobs.bin PREFIX.out.bin -t 64
  truth_oracle.py report truth.npz PREFIX SHARES NAME=RUN.npz ... [--t2t NAME=T2T.npz ...] [--save OUT]

The runs are truth_t2t.py place files. `jobs` takes every pair in which a
run leaves a mate with a truth away from its locus (150 bp) and lists, for
both mates, the loci to score: the truth and each run's placement (loci of a
mate within 150 bp of one another are one job). tools/truth_oracle scores
them: the best local score of the read in the locus's window under CERTA's
objective (BWA-MEM's defaults), by a full DP.

`report` compares, for each such mate of each run, the best pair with the
mate at its truth against the pair the run reported. A pair scores the sum
of its mates' scores, minus BWA-MEM's unpaired penalty (17) when the two are
not a proper pair (same contig, opposite strands, the reverse mate starting
-150..1,000 bases after the forward one); an unplaced mate scores 0. With
the mate at its truth, the other mate is where the run put it or at its own
truth, whichever scores more. The classes:
  missed     the truth scores more: a max-score aligner places the mate there;
  tie        the same score: a coin decides;
  reference  the truth scores less: no aligner that maximizes this score on
             this reference places the mate at its origin.
SHARES: the operating points (shares of the mates with a truth) at which
the classes of each run's N most confident mates are counted.
--t2t: of each class, the mates a T2T-CHM13 run lifts to the truth.
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from truth_lift import CLS_TRUTH  # noqa: E402

LOCUS = 150
PEN_UNPAIRED = 17
JOB = np.dtype([("pair", "<i4"), ("pos", "<i4"), ("contig", "<i2"), ("mate", "i1"), ("strand", "i1"),
                ("unused", "<i4")])
OUT = np.dtype([("score", "<i2"), ("aligned", "<i2"), ("start", "<i4")])


def load_run(spec):
    name, p = spec.split("=", 1)
    z = np.load(p)
    return name, dict(placed=z["placed"], rc=z["rc"].astype(np.int64), left=z["left"].astype(np.int64),
                      strand=z["strand"].astype(np.int64), mapq=z["mapq"].astype(np.int64))


def truth_run(t):
    return dict(placed=t["cls"] == CLS_TRUTH, rc=t["tc"].astype(np.int64), left=t["tp"].astype(np.int64) - 1,
                strand=t["ts"].astype(np.int64))


def at(r, s):
    return (r["placed"] & s["placed"] & (r["rc"] == s["rc"]) & (r["strand"] == s["strand"])
            & (np.abs(r["left"] - s["left"]) <= LOCUS))


def loci_of(tr, runs, sel):
    """Arrays (pairs, 2, K) of the loci of the selected pairs: the truth, then each run's."""
    srcs = [tr] + [r for _, r in runs]
    return {k: np.stack([s[k][sel] for s in srcs], axis=2) for k in ("placed", "rc", "left", "strand")}


def job_map(lc):
    """Job number (per mate, 0-based among that mate's loci) of every locus;
    -1 for none. A locus within LOCUS of an earlier one of the mate shares its job."""
    K = lc["placed"].shape[2]
    jm = np.full(lc["placed"].shape, -1, np.int64)
    njob = np.zeros(lc["placed"].shape[:2], np.int64)
    for k in range(K):
        new = lc["placed"][:, :, k].copy()
        for j in range(k):
            same = (lc["placed"][:, :, k] & lc["placed"][:, :, j] & (lc["rc"][:, :, k] == lc["rc"][:, :, j])
                    & (lc["strand"][:, :, k] == lc["strand"][:, :, j])
                    & (np.abs(lc["left"][:, :, k] - lc["left"][:, :, j]) <= LOCUS) & (jm[:, :, k] < 0))
            jm[:, :, k] = np.where(same, jm[:, :, j], jm[:, :, k])
            new &= ~same
        jm[:, :, k] = np.where(new, njob, jm[:, :, k])
        njob += new
    return jm, njob


def select(t, tr, runs):
    truth = tr["placed"]
    bad = np.zeros(truth.shape, bool)
    for _, r in runs:
        bad |= truth & ~at(r, tr)
    return np.flatnonzero(bad.any(axis=1))


def cmd_jobs(a):
    t = np.load(a.truth)
    tr = truth_run(t)
    runs = [load_run(s) for s in a.runs]
    sel = select(t, tr, runs)
    lc = loci_of(tr, runs, sel)
    jm, njob = job_map(lc)
    # Jobs in (pair, mate, job) order; the center is the first locus of the job.
    base = np.zeros(njob.size + 1, np.int64)
    np.cumsum(njob.ravel(), out=base[1:])
    jobs = np.zeros(int(base[-1]), JOB)
    K = lc["placed"].shape[2]
    first = np.ones(jm.shape, bool)
    for k in range(K):
        for j in range(k):
            first[:, :, k] &= jm[:, :, k] != jm[:, :, j]
        first[:, :, k] &= jm[:, :, k] >= 0
    pi, mi, ki = np.nonzero(first)
    g = base[:-1].reshape(njob.shape)[pi, mi] + jm[pi, mi, ki]
    jobs["pair"][g] = sel[pi]
    jobs["mate"][g] = mi
    jobs["contig"][g] = lc["rc"][pi, mi, ki]
    jobs["pos"][g] = lc["left"][pi, mi, ki]
    jobs["strand"][g] = lc["strand"][pi, mi, ki]
    assert (np.diff(jobs["pair"]) >= 0).all()
    jobs.tofile(a.prefix + ".jobs.bin")
    np.savez(a.prefix + ".jobs.npz", sel=sel, runs=np.array([n for n, _ in runs]))
    print(f"[jobs] {len(sel):,} pairs, {len(jobs):,} jobs", file=sys.stderr)


def cmd_report(a):
    t = np.load(a.truth)
    tr = truth_run(t)
    truth = tr["placed"]
    nt = int(truth.sum())
    runs = [load_run(s) for s in a.runs]
    z = np.load(a.prefix + ".jobs.npz")
    sel = z["sel"]
    assert list(z["runs"]) == [n for n, _ in runs], "the runs differ from those of `jobs`"
    lc = loci_of(tr, runs, sel)
    jm, njob = job_map(lc)
    base = np.zeros(njob.size + 1, np.int64)
    np.cumsum(njob.ravel(), out=base[1:])
    out = np.fromfile(a.prefix + ".out.bin", OUT)
    assert len(out) == base[-1], "the scores do not match the jobs"
    b2 = base[:-1].reshape(njob.shape)
    # score[p, m, k]: the score of mate m at its locus k (0 for none)
    score = np.where(jm >= 0, out["score"][np.clip(b2[:, :, None] + jm, 0, max(len(out) - 1, 0))].astype(np.int64), 0)
    score = np.where(score < 0, 0, score)

    def proper(k0, k1):
        """Mate 0 at its locus k0 and mate 1 at its locus k1 are a proper pair."""
        p0, p1 = lc["placed"][:, 0, k0], lc["placed"][:, 1, k1]
        s0, s1 = lc["strand"][:, 0, k0], lc["strand"][:, 1, k1]
        l0, l1 = lc["left"][:, 0, k0], lc["left"][:, 1, k1]
        d = np.where(s0 == 0, l1 - l0, l0 - l1)  # reverse start minus forward start
        return p0 & p1 & (lc["rc"][:, 0, k0] == lc["rc"][:, 1, k1]) & (s0 != s1) & (d >= -150) & (d <= 1000)

    def pair_score(k0, k1):
        both = lc["placed"][:, 0, k0] & lc["placed"][:, 1, k1]
        return score[:, 0, k0] + score[:, 1, k1] - PEN_UNPAIRED * (both & ~proper(k0, k1))

    t2 = [load_run(s) for s in a.t2t]
    shares = [float(s) for s in a.shares.split(",")]
    print(f"Mates with a truth: {nt:,}. Pairs scored: {len(sel):,} (a mate away from its locus in some run).\n")
    for ri, (name, r) in enumerate(runs):
        k = ri + 1
        ok = at(r, tr)
        mis = truth & ~ok  # misplaced or unplaced
        cls = np.zeros(truth.shape, np.int8)  # 0 at the locus / no truth, 1 missed, 2 tie, 3 reference
        delta = np.zeros(truth.shape, np.int64)
        for m in range(2):
            kk = (lambda km, kp: (km, kp)) if m == 0 else (lambda km, kp: (kp, km))
            c_r = pair_score(*kk(k, k))
            c_t = pair_score(*kk(0, k))
            both_truth = lc["placed"][:, 1 - m, 0]
            c_t = np.where(both_truth, np.maximum(c_t, pair_score(0, 0)), c_t)
            d = c_t - c_r
            mm = mis[sel, m]
            delta[sel, m] = np.where(mm, d, 0)
            cls[sel, m] = np.where(mm, np.where(d > 0, 1, np.where(d == 0, 2, 3)), 0)
        assert (cls[mis] > 0).all()
        if a.save:
            np.savez(f"{a.save}.{ri}.npz", cls=cls, delta=delta.astype(np.int16))
        placed = r["placed"] & truth
        print(f"### {name}\n")
        print(f"At the locus: {int((ok & truth).sum()):,} ({100 * (ok & truth).sum() / nt:.3f} %); "
              f"misplaced {int((mis & r['placed']).sum()):,}, unplaced {int((mis & ~r['placed']).sum()):,}.\n")
        names = ("missed (the truth scores more)", "tie", "reference (the truth scores less)")
        print("| class | mates | per million | MAPQ ≥ 40 | 20–39 | 1–19 | 0 | unplaced |"
              + "".join(f" lifted to the truth by {n} (MAPQ ≥ 40 there) |" for n, _ in t2))
        print("|---|---|---|---|---|---|---|---|" + "---|" * len(t2))
        for c in (1, 2, 3):
            x = cls == c
            q = r["mapq"]
            pl = r["placed"]
            cells = [f"{int(x.sum()):,}", f"{1e6 * x.sum() / nt:,.0f}", f"{int((x & pl & (q >= 40)).sum()):,}",
                     f"{int((x & pl & (q >= 20) & (q < 40)).sum()):,}", f"{int((x & pl & (q >= 1) & (q < 20)).sum()):,}",
                     f"{int((x & pl & (q == 0)).sum()):,}", f"{int((x & ~pl).sum()):,}"]
            for _, s in t2:
                hit = x & at(s, tr)
                cells.append(f"{int(hit.sum()):,} ({int((hit & (s['mapq'] >= 40)).sum()):,})")
            print(f"| {names[c - 1]} | " + " | ".join(cells) + " |")
        ref = cls == 3
        if ref.any():
            d = -delta[ref]
            print(f"\nReference class, how much less the truth scores: 1–4 {int((d <= 4).sum()):,}, "
                  f"5–9 {int(((d >= 5) & (d <= 9)).sum()):,}, 10–16 {int(((d >= 10) & (d <= 16)).sum()):,}, "
                  f"17–29 {int(((d >= 17) & (d <= 29)).sum()):,}, ≥ 30 {int((d >= 30).sum()):,}.")
        # The classes among the N most confident mates (boundary level pro rata).
        print("\n| share of the mates with a truth | misplaced per million | missed | tie | reference |")
        print("|---|---|---|---|---|")
        mq = r["mapq"][placed]
        cl = cls[placed]
        order = np.argsort(-mq, kind="stable")
        mq, cl = mq[order], cl[order]
        for s in shares:
            n = int(round(s * nt))
            if n >= len(mq):
                cnt = [float((cl == c).sum()) for c in (1, 2, 3)]
            else:
                qb = mq[n - 1]
                above, lvl = mq > qb, mq == qb
                w = (n - int(above.sum())) / max(1, int(lvl.sum()))
                cnt = [float((cl[above] == c).sum() + w * (cl[lvl] == c).sum()) for c in (1, 2, 3)]
            print(f"| {100 * s:.2f} % | {1e6 * sum(cnt) / n:,.1f} | " + " | ".join(f"{1e6 * c / n:,.1f}" for c in cnt) + " |")
        print()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    p = sp.add_parser("jobs")
    p.add_argument("truth"), p.add_argument("prefix"), p.add_argument("runs", nargs="+")
    p = sp.add_parser("report")
    p.add_argument("truth"), p.add_argument("prefix"), p.add_argument("shares"), p.add_argument("runs", nargs="+")
    p.add_argument("--t2t", action="append", default=[])
    p.add_argument("--save", help="write each run's classes and score differences to SAVE.<run index>.npz")
    a = ap.parse_args()
    {"jobs": cmd_jobs, "report": cmd_report}[a.cmd](a)


if __name__ == "__main__":
    main()
