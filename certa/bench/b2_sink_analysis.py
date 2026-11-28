#!/usr/bin/env python3
"""B2 sink analysis: of the mates placed on sink_* contigs by the augmented
index, how many are stolen real reads (a GRCh38 truth mate attracted to a sink)
vs correct missing-copy abstentions (no-linear-truth mates)? A good sink set
abstains on missing-copy reads without stealing real ones (FixItFelix lesson).

  python3 b2_sink_analysis.py AUG.sam TRUTH.npz
"""
import sys

import numpy as np


def main():
    sam, truth_path = sys.argv[1], sys.argv[2]
    t = np.load(truth_path)
    cls = t["cls"]
    n = cls.shape[0]
    sink = np.zeros((n, 2), bool)
    smq = np.zeros((n, 2), np.uint8)
    nrec = 0
    with open(sam) as f:
        for line in f:
            if line[0] == "@":
                continue
            p = line.split("\t", 11)
            flag = int(p[1])
            if flag & (256 | 2048 | 4):
                continue
            nrec += 1
            if not p[2].startswith("sink_"):
                continue
            i = int(p[0][1:])
            m = 1 if flag & 128 else 0
            sink[i, m] = True
            smq[i, m] = min(int(p[4]), 255)
    truth = cls == 0   # CLS_TRUTH: a GRCh38 truth locus exists
    nolin = cls == 1   # CLS_NOLIN: missing-copy class, no GRCh38 truth
    ns = int(sink.sum())
    print(f"primary records: {nrec:,}; sink placements: {ns:,}")
    if ns:
        print(f"  stolen real reads (GRCh38-truth mate on a sink): {int((sink & truth).sum()):,}")
        print(f"  missing-copy abstentions (no-linear-truth on a sink): {int((sink & nolin).sum()):,}")
        print(f"  other/unlabelled on a sink: {int((sink & ~truth & ~nolin).sum()):,}")
        print(f"  sink placements with MAPQ>=40: {int((sink & (smq >= 40)).sum()):,}")
        print(f"  stolen real reads with MAPQ>=40: {int((sink & truth & (smq >= 40)).sum()):,}")


if __name__ == "__main__":
    main()
