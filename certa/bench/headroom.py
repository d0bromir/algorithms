#!/usr/bin/env python3
"""What a full-sensitivity mapper finds for the reads CERTA leaves uncertified:
the headroom of patterns taken from the slow mappers.

  headroom.py yara.sam N_READS LABEL=uncertified.fq [LABEL=uncertified.fq ...]

uncertified.fq: `certa map ... --unc-tags -u` (each name carries parts=,
used= (|S|, the parts that fit the hit budget) and reason=). yara.sam: Yara
on the same reads, every alignment with <= 4 % errors (6 edits at 150 bp),
end to end. Per uncertified read: how many parts were usable (all of them,
or fewer because some were too repetitive) against Yara's best edit count
(none within 6 edits: the read needs more edits or clipping). A read whose
best alignment has d edits is within reach of a certificate of radius d: of
the exact-part certificate once enough parts are usable (R = |S| - 1 <= 5),
or of parts matched with up to one error (R = 2 |S| - 1).
"""
import sys
from collections import Counter

REASON = {2: 'radius_negative', 3: 'not_found', 4: 'cluster_too_wide', 5: 'cross_contig', 1: 'bad_length'}


def yara(path):
    best = {}
    with open(path) as f:
        for line in f:
            if line[0] == '@':
                continue
            t = line.rstrip('\n').split('\t')
            flag = int(t[1])
            if flag & 0x900:
                continue
            best[t[0]] = None if flag & 4 else next(int(x[5:]) for x in t[11:] if x.startswith('NM:i:'))
    return best


def uncertified(path):
    out = []
    with open(path) as f:
        for i, line in enumerate(f):
            if i % 4:
                continue
            name, *tags = line[1:].split()
            kv = dict(x.split('=') for x in tags if '=' in x)
            out.append((name, int(kv['parts']), int(kv['used']), int(kv['reason'])))
    return out


def main():
    y = yara(sys.argv[1])
    n_reads = int(sys.argv[2])  # Yara may leave unmapped reads out of its SAM
    print(f'reads {n_reads}; Yara (-e 4 -y full) aligns {sum(v is not None for v in y.values())} within 6 edits')
    for spec in sys.argv[3:]:
        label, path = spec.split('=', 1)
        unc = uncertified(path)
        n = len(unc)
        print(f'\n== {label}: {n} uncertified reads ({100 * n / n_reads:.2f} % of all)')
        by_reason = Counter(REASON.get(r, str(r)) for _, _, _, r in unc)
        print('reasons: ' + ', '.join(f'{k} {v} ({100 * v / n_reads:.2f} %)' for k, v in by_reason.most_common()))
        # Groups: usable parts (|S|) against Yara's best edit count.
        grid = Counter()
        for name, parts, used, reason in unc:
            d = y.get(name)
            col = 'none <= 6' if d is None else str(d)
            row = f'{used} of {parts}' if used < parts else f'all {parts}'
            grid[(row, col)] += 1
        rows = sorted({r for r, _ in grid}, key=lambda r: (r.startswith('all'), int(r.split()[0]) if not r.startswith('all') else 99))
        cols = [str(d) for d in range(7)] + ['none <= 6']
        print(f"{'usable parts':>14s} " + ' '.join(f'{c:>9s}' for c in cols) + f" {'total':>9s}")
        for r in rows:
            tot = sum(grid[(r, c)] for c in cols)
            print(f'{r:>14s} ' + ' '.join(f'{grid[(r, c)]:>9d}' for c in cols) + f' {tot:>9d}')
        print('  (columns: Yara best edit count; "none <= 6": no end-to-end alignment within 6 edits)')
        # Headroom (as shares of all reads).
        rep_le5 = sum(v for (r, c), v in grid.items() if not r.startswith('all') and c not in ('none <= 6', '6'))
        rep_6 = sum(v for (r, c), v in grid.items() if not r.startswith('all') and c == '6')
        all_le5 = sum(v for (r, c), v in grid.items() if r.startswith('all') and c not in ('none <= 6', '6'))
        all_6 = sum(v for (r, c), v in grid.items() if r.startswith('all') and c == '6')
        none = sum(v for (r, c), v in grid.items() if c == 'none <= 6')
        pct = lambda v: f'{v} ({100 * v / n_reads:.2f} % of all reads)'
        print(f'repeat-limited (some parts unusable), best <= 5 edits: {pct(rep_le5)}')
        print(f'repeat-limited, best 6 edits:                        {pct(rep_6)}')
        print(f'all parts usable, best <= 5 edits (ties, wide clusters): {pct(all_le5)}')
        print(f'all parts usable, best 6 edits:                       {pct(all_6)}')
        print(f'no end-to-end alignment within 6 edits:               {pct(none)}')


if __name__ == '__main__':
    main()
