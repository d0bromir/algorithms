#!/usr/bin/env python3
"""Pass V's attempts (reads uncertified with --full-cert --adaptive-parts 0)
split into successes and failures (reads still uncertified with --full-cert),
by cheap per-read features from the FASTQ: how many bases are below Q20, in
the last 30 bases, the longest run of low-quality bases, N, poly-G, length.
Prints, per feature threshold, how many attempts it would skip and how many
of those pass V certifies (the cost of skipping them).
  passv_features.py attempts.fq failures.fq
"""
import sys
from collections import Counter


def records(path):
    with open(path) as f:
        while True:
            h = f.readline()
            if not h:
                return
            s = f.readline().rstrip("\n")
            f.readline()
            q = f.readline().rstrip("\n")
            yield h[1:].split()[0], s, q


def features(s, q):
    qs = [ord(c) - 33 for c in q]
    low = [v < 20 for v in qs]
    run = best = 0
    for x in low:
        run = run + 1 if x else 0
        best = max(best, run)
    return {
        "lowq": sum(low),
        "lowq_tail30": sum(low[-30:]),
        "lowq_run": best,
        "n": s.count("N"),
        "polyG": len(s) - len(s.rstrip("G")),
        "short": 150 - len(s),
    }


def main():
    failures = {name for name, _, _ in records(sys.argv[2])}
    rows = [(features(s, q), name in failures) for name, s, q in records(sys.argv[1])]
    n = len(rows)
    nf = sum(f for _, f in rows)
    print(f"pass V attempts {n}, failures {nf} ({100 * nf / n:.1f} %), successes {n - nf}")
    for feat in ["lowq", "lowq_tail30", "lowq_run", "n", "polyG", "short"]:
        print(f"== {feat}: skip if >= t  ->  attempts skipped, of which failures, successes lost (share of all successes)")
        for t in [1, 2, 3, 5, 8, 10, 15, 20, 30]:
            sk = [f for x, f in rows if x[feat] >= t]
            if not sk:
                continue
            fails = sum(sk)
            lost = len(sk) - fails
            print(f"  t={t:<3d} skipped {len(sk):6d} ({100 * len(sk) / n:5.1f} % of attempts), failures {fails:6d} "
                  f"({100 * fails / max(1, len(sk)):5.1f} %), successes lost {lost:5d} ({100 * lost / max(1, n - nf):5.2f} %)")


if __name__ == "__main__":
    main()
