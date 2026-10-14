"""Summarizes an Aim 6 run of tools/lr_ceiling (proposal §6.6).

  aim6_summary.py run_tau16.tsv [run_tau256.tsv]

Prints the pre-registered criterion (share of all reads with sigma > U1),
then descriptive statistics that are not part of it: the margin, why reads
miss, the ceiling's parts and the enumerated occurrences per read, the
trimmed-path variant and a second tau. For a run with --certify, also the
completed certificates: OPT proven, how often it beats the best known
alignment, and the cost (DP cells and time per read).

Miss classes (descriptive, defined here):
  unmapped        no primary alignment;
  clipped         the unaligned read bases plus clip penalties alone reach
                  the ceiling L - U1;
  repeat-limited  the loss L - sigma would fit the ideal ceiling of
                  5 floor(L / 22) (22-base parts tiling the read), so the
                  read misses because its parts are longer or skipped;
  high loss       the rest: even ideal parts would not cover the loss.
"""
import math
import sys


def load(path):
    rows = []
    with open(path) as f:
        head = next(f).rstrip("\n").split("\t")
        for line in f:
            v = line.rstrip("\n").split("\t")
            r = dict(zip(head, v))
            for k in head[1:]:
                r[k] = int(r[k])
            rows.append(r)
    return rows


def wilson(k, n, z=1.96):
    if n == 0:
        return 0.0, 0.0
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def q(xs, p):
    if not xs:
        return float("nan")
    xs = sorted(xs)
    return xs[min(len(xs) - 1, int(p * len(xs)))]


def classify(r):
    if not r["mapped"]:
        return "unmapped"
    L = r["len"]
    ceiling = L - r["U1"]
    unaligned = r["qs"] + (L - r["qe"]) + 5 * (r["qs"] > 0) + 5 * (r["qe"] < L)
    if unaligned >= ceiling:
        return "clipped"
    if L - r["sigma"] < 5 * (L // 22):
        return "repeat-limited"
    return "high loss"


def summary(rows, label):
    n = len(rows)
    fit = sum(r["fits"] for r in rows)
    lo, hi = wilson(fit, n)
    print(f"== {label}")
    print(f"reads {n}; fit {fit} = {100 * fit / n:.2f} % (95 % CI {100 * lo:.2f}-{100 * hi:.2f} %); "
          f"go criterion (>= 90 %): {'met' if fit >= 0.9 * n else 'NOT met'}")
    mapped = [r for r in rows if r["mapped"]]
    print(f"mapped {len(mapped)} ({100 * len(mapped) / n:.2f} %); "
          f"fit with the trimmed path {sum(r['mapped'] and r['sigma_trim'] > r['U1'] for r in rows)}")
    lens = [r["len"] for r in rows]
    print(f"read length: median {q(lens, 0.5)}, 10-90 % {q(lens, 0.1)}-{q(lens, 0.9)}, max {max(lens)}")
    m = [r["margin"] for r in mapped]
    print("margin sigma - U1 (mapped): " + ", ".join(f"{p}% {q(m, p / 100)}" for p in (1, 5, 10, 50, 90)))
    loss = [r["len"] - r["sigma"] for r in mapped]
    ceil = [r["len"] - r["U1"] for r in rows]
    print("loss L - sigma (mapped): " + ", ".join(f"{p}% {q(loss, p / 100)}" for p in (10, 50, 90, 99)))
    print("ceiling L - U1: " + ", ".join(f"{p}% {q(ceil, p / 100)}" for p in (1, 10, 50, 90)))
    parts = [r["parts_rev"] if r["mapped"] and r["strand"] else r["parts_fwd"] for r in rows]
    per_kb = [1000 * p / r["len"] for p, r in zip(parts, rows)]
    print(f"parts on the alignment strand: median {q(parts, 0.5)} ({q(per_kb, 0.5):.1f} per kb; "
          f"22-base tiling gives {1000 / 22:.1f}); skipped bases median {q([r['skipped'] for r in rows], 0.5)}, "
          f"99% {q([r['skipped'] for r in rows], 0.99)}")
    occ = [r["occ_sum"] for r in rows]
    print("enumerated occurrences per read (both strands): " +
          ", ".join(f"{p}% {q(occ, p / 100)}" for p in (50, 90, 99)) + f", max {max(occ)}, "
          f"mean {sum(occ) / n:.0f}")
    misses = [r for r in rows if not r["fits"]]
    cls = {}
    for r in misses:
        c = classify(r)
        cls[c] = cls.get(c, 0) + 1
    print("misses: " + "; ".join(f"{c} {k} ({100 * k / n:.2f} % of reads)" for c, k in
                                 sorted(cls.items(), key=lambda x: -x[1])))
    if "cert" not in rows[0]:
        return
    # --certify: OPT proven by banded DP over the candidate chains.
    cert = [r for r in rows if r["cert"] == 1]
    capped = sum(r["cert"] == 2 for r in rows)
    print(f"certified (OPT proven): {len(cert)} = {100 * len(cert) / n:.2f} % of reads; over the cell cap {capped}")
    if not cert:
        return
    gain = [r["opt"] - max(r["sigma"], r["sigma_trim"]) for r in cert]
    up = [g for g in gain if g > 0]
    print(f"OPT above the best known alignment: {len(up)} reads ({100 * len(up) / len(cert):.1f} %), by "
          + (", ".join(f"{p}% {q(up, p / 100)}" for p in (50, 90, 99)) + f", max {max(up)}" if up else "-"))
    cells = [r["cells"] for r in cert]
    us = [r["us"] for r in cert]
    print("DP cells per read: " + ", ".join(f"{p}% {q(cells, p / 100):.3g}" for p in (50, 90, 99)) +
          f", max {max(cells):.3g}, mean {sum(cells) / len(cells):.3g}")
    print("DP time per read (one thread): " + ", ".join(f"{p}% {q(us, p / 100) / 1000:.1f} ms" for p in (50, 90, 99)) +
          f", max {max(us) / 1000:.0f} ms, mean {sum(us) / len(us) / 1000:.1f} ms")
    s = [r["searched"] for r in cert]
    c = [r["chains"] for r in cert]
    print(f"chains searched by DP: median {q(s, 0.5)}, 99% {q(s, 0.99)}, max {max(s)} (of chains built: median "
          f"{q(c, 0.5)}, max {max(c)})")


def main():
    a = load(sys.argv[1])
    summary(a, sys.argv[1])
    if len(sys.argv) > 2:
        b = load(sys.argv[2])
        summary(b, sys.argv[2])
        fa = {r["read"] for r in a if r["fits"]}
        fb = {r["read"] for r in b if r["fits"]}
        print(f"== {sys.argv[2]} vs {sys.argv[1]}: fit in both {len(fa & fb)}, "
              f"only the first {len(fa - fb)}, only the second {len(fb - fa)}")


if __name__ == "__main__":
    main()
