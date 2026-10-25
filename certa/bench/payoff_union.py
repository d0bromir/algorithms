#!/usr/bin/env python3
"""Combine tools/payoff runs: what techniques certify together, run as a
cascade (each read goes on to the next technique only while uncertified),
and what the cascade costs.

  payoff_union.py N_READS CERTA_CPU_S LABEL=per_read.tsv [LABEL=per_read.tsv ...]

per_read.tsv: `payoff ... --out` on the reads CERTA leaves uncertified
(columns: read, then per mode its certified distance (-1: none) and
microseconds). A technique is LABEL:mode (e.g. q22:rare). N_READS: all reads
of the sample; CERTA_CPU_S: CERTA's CPU time for all of them, to set the
cascade's cost against.
"""
import sys

CASCADES = [
    ['q22:fixed'],
    ['q22:rare'],
    ['q22:err1r'],
    ['q22:err1'],
    ['q22:long'],
    ['q15:rare'],
    ['q22:rare', 'q22:err1r'],
    ['q22:rare', 'q22:err1r', 'q22:err1'],
    ['q22:rare', 'q15:rare', 'q22:err1r'],
    ['q22:rare', 'q15:rare', 'q22:err1r', 'q22:err1'],
    ['q22:rare', 'q15:rare', 'q22:err1r', 'q22:err1', 'q22:long', 'q15:long', 'q15:err1r', 'q15:err1'],
]


def load(path):
    rows = {}
    with open(path) as f:
        head = f.readline().rstrip('\n').split('\t')
        for line in f:
            t = line.rstrip('\n').split('\t')
            rows[t[0]] = {head[i]: int(float(t[i])) for i in range(1, len(t))}
    return rows


def main():
    n_all, certa_cpu = int(sys.argv[1]), float(sys.argv[2])
    runs = {}
    for spec in sys.argv[3:]:
        label, path = spec.split('=', 1)
        runs[label] = load(path)
    reads = sorted(set.intersection(*(set(r) for r in runs.values())))
    print(f'{len(reads)} uncertified reads ({100 * len(reads) / n_all:.2f} % of {n_all}); '
          f"CERTA's CPU for all reads {certa_cpu:.0f} s")
    # Every certificate proves the same minimum, so all must agree.
    conflicts = 0
    for r in reads:
        ds = {row[m] for row in (runs[label][r] for label in runs) for m in row if not m.endswith('_us') and row[m] >= 0}
        conflicts += len(ds) > 1
    print(f'reads whose certified distances disagree between techniques: {conflicts}')
    print(f"{'cascade':58s} {'certified':>10s} {'% of all':>9s} {'<= 6 edits':>11s} {'7+ edits':>9s} "
          f"{'CPU (s)':>9s} {'vs CERTA':>9s}")
    for cas in CASCADES:
        if any(c.split(':')[0] not in runs for c in cas):
            continue
        cert = le6 = gt6 = 0
        us = 0
        for r in reads:
            for c in cas:
                label, mode = c.split(':')
                row = runs[label][r]
                us += row[mode + '_us']
                d = row[mode]
                if d >= 0:
                    cert += 1
                    le6 += d <= 6
                    gt6 += d > 6
                    break
        cpu = us / 1e6
        print(f"{' -> '.join(cas):58s} {cert:>10d} {100 * cert / n_all:>8.2f}% {le6:>11d} {gt6:>9d} "
              f"{cpu:>9.1f} {100 * cpu / certa_cpu:>8.0f}%")


if __name__ == '__main__':
    main()
