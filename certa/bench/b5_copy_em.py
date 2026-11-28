#!/usr/bin/env python3
"""B5: copy-level joint assignment (prototype).

The missing-copy class is a paralog *collapsed onto* a GRCh38 locus (B2 finding),
so it cannot be recovered by adding gap/decoy sequence. But reads originating from
an unseen copy carry mismatches vs the collapsed reference that *recur at the same
reference positions across reads* -- a signature a per-read aligner cannot see but
a copy-level model can. This prototype, for reads mapped into segmental
duplications, estimates per-locus divergent positions (where a fraction of reads
consistently mismatch) and assigns each read to the reference copy or a phantom
copy. Phantom-assigned reads are the missing-copy candidates (abstain / name the
copy). Validates against the assembly truth: phantom assignment should enrich for
the no-linear-truth (missing-copy) class and for misplacements.

  python3 b5_copy_em.py REF.fa SEGDUPS.bed.gz SAM TRUTH.npz SCORED.npz [--min-reads 10]

REF.fa GRCh38 fasta; SEGDUPS.bed.gz GRCh38 segdups; SAM the GRCh38 mapping;
TRUTH.npz cls (0 truth, 1 no-linear-truth); SCORED.npz correct/mapq (truth_eval).
"""
import argparse
import gzip
from collections import defaultdict


def read_fasta(path):
    seqs = {}
    name = None
    chunks = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def read_bed(path):
    iv = defaultdict(list)
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as f:
        for line in f:
            p = line.split()
            if len(p) < 3 or p[0].startswith(("#", "track", "browser")):
                continue
            c = p[0][3:] if p[0].startswith("chr") else p[0]  # UCSC BED vs Ensembl SAM/ref
            iv[c].append((int(p[1]), int(p[2])))
    for c in iv:
        iv[c].sort()
    return iv


def locus_of(iv, c, pos):
    import bisect
    arr = iv.get(c, [])
    starts = [a for a, b in arr]
    j = bisect.bisect_right(starts, pos) - 1
    if j >= 0 and arr[j][0] <= pos < arr[j][1]:
        return (c, arr[j][0], arr[j][1])
    return None


def mismatches(seq, ref, rpos, cigar):
    """Reference positions where the read mismatches (walk CIGAR M ops)."""
    out = []
    q = 0
    r = rpos
    for op in cigar:
        ln, code = op >> 4, op & 0xF
        if code == 0:  # M
            for k in range(ln):
                if q + k < len(seq) and r + k < len(ref) and seq[q + k] != ref[r + k] and seq[q + k] != "N":
                    out.append(r + k)
            q += ln
            r += ln
        elif code == 1:  # I
            q += ln
        elif code == 2:  # D
            r += ln
        elif code == 4:  # S
            q += ln
    return out


def parse_cigar(s):
    import re
    return [(int(n) << 4) | "MIDSH".index(c) for n, c in re.findall(r"(\d+)([MIDSH])", s)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ref")
    ap.add_argument("segdups")
    ap.add_argument("sam")
    ap.add_argument("truth")
    ap.add_argument("scored")
    ap.add_argument("--min-reads", type=int, default=10)
    ap.add_argument("--div-freq", type=float, default=0.2)
    a = ap.parse_args()
    import numpy as np
    print("loading reference + segdups...", flush=True)
    ref = read_fasta(a.ref)
    sd = read_bed(a.segdups)
    cls = np.load(a.truth)["cls"]
    sc = np.load(a.scored)
    correct, mapq = sc["correct"], sc["mapq"]

    # collect reads per segdup locus: (name_idx, mate, mismatch_positions, nmm)
    loci = defaultdict(list)
    n_sd_reads = 0
    with open(a.sam) as f:
        for line in f:
            if line[0] == "@":
                continue
            p = line.split("\t", 11)
            flag = int(p[1])
            if flag & (256 | 2048 | 4):
                continue
            c = p[2]
            if c not in ref:
                continue
            rpos = int(p[3]) - 1
            loc = locus_of(sd, c, rpos)
            if loc is None:
                continue
            cigar = parse_cigar(p[5])
            mm = mismatches(p[9], ref[c], rpos, cigar)
            i = int(p[0][1:])
            m = 1 if flag & 128 else 0
            loci[loc].append((i, m, mm))
            n_sd_reads += 1
    print(f"segdup reads: {n_sd_reads:,}; loci: {len(loci):,}", flush=True)

    # per-locus copy assignment
    phantom_reads = []   # (i, m)
    ref_reads = []
    n_loci_used = 0
    for loc, reads in loci.items():
        if len(reads) < a.min_reads:
            continue
        n_loci_used += 1
        c, ls, le = loc
        # per-position mismatch count + coverage
        mmcount = defaultdict(int)
        cov = defaultdict(int)
        for i, m, mm in reads:
            for pos in mm:
                mmcount[pos] += 1
        # coverage: reads spanning each position (approx via read start..end)
        # use mismatch-bearing positions only; divergent = freq among reads covering
        # approximate coverage by number of reads whose span includes pos
        spans = []
        for i, m, mm in reads:
            spans.append((mm[0] if mm else ls, mm[-1] if mm else ls))
        # divergent positions: mismatched by >= div-freq of reads (and >=3 reads)
        div = set()
        for pos, cnt in mmcount.items():
            covering = sum(1 for s, e in spans if s <= pos <= e)
            if covering >= 3 and cnt / max(covering, 1) >= a.div_freq and cnt >= 3:
                div.add(pos)
        if not div:
            for i, m, mm in reads:
                ref_reads.append((i, m))
            continue
        for i, m, mm in reads:
            at_div = sum(1 for pos in mm if pos in div)
            # phantom if most mismatches are at divergent positions and >=2 there
            if at_div >= 2 and at_div >= 0.5 * max(len(mm), 1):
                phantom_reads.append((i, m))
            else:
                ref_reads.append((i, m))

    def enrich(rs, label):
        if not rs:
            print(f"{label}: 0 reads")
            return
        n = len(rs)
        nolin = sum(1 for i, m in rs if cls[i, m] == 1)
        truth = sum(1 for i, m in rs if cls[i, m] == 0)
        mis = sum(1 for i, m in rs if cls[i, m] == 0 and not correct[i, m])
        conf_mis = sum(1 for i, m in rs if cls[i, m] == 0 and not correct[i, m] and mapq[i, m] >= 40)
        print(f"{label}: {n:,} reads | no-linear-truth {nolin:,} ({100*nolin/n:.1f}%) | "
              f"misplaced {mis:,} ({100*mis/n:.2f}%) | confident-misplaced {conf_mis:,}")

    print(f"\nloci used (>={a.min_reads} reads): {n_loci_used:,}; divergent-position loci drive assignment")
    enrich(phantom_reads, "phantom-assigned (missing-copy candidates)")
    enrich(ref_reads, "reference-assigned")


if __name__ == "__main__":
    main()
