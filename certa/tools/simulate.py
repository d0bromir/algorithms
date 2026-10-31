#!/usr/bin/env python3
"""Synthetic data for exercising CERTA (standard library only).

  simulate.py genome OUT.fa  [--length 20000000] [--contigs 4] [--seed 1]
  simulate.py reads  REF.fa OUT.fq [--n 1000000] [--length 150]
                     [--error 0.002] [--snv 0.001] [--indel 0.0001] [--seed 2]
                     [--quality flat|novaseqx|novaseq6000] [--titv 2.1]
                     [--mates OUT_2.fq] [--insert 400,100]
  simulate.py eval   OUT.sam [--window 10]

The genome contains repeat families (diverged copies), segmental
duplications, tandem repeats and N gaps, so that repeats and caps are
exercised. Read names carry the truth: r<i>_<contig>_<pos1>_<strand>_<edits>.
Real-data validation (GIAB, NovaSeq X, ONT) is described in
RESEARCH_ROADMAP.md; synthetic data only checks mechanics and speed.

--quality flat (the default) writes Q37 everywhere and draws errors at
--error, as before. novaseqx and novaseq6000 draw each base's quality from
the bins measured on GIAB HG002 NovaSeq X (2/9/24/40, ~95 % in the top bin)
and HG005 NovaSeq 6000 (2/11/25/37, ~92 %; bench/results/hg002_2026-10/
recent_illumina_costs.txt), lower towards the 3' end, with a low-quality
tail on ~1.5 % of reads; each base is then miscalled with the probability
its quality states. Variants are shared by both mates of a fragment: SNVs at
--snv (transitions --titv times as frequent as transversions) and indels at
--indel (ten times as frequent next to a homopolymer of 4 or more, geometric
lengths). With --mates, pairs: mate 1 forward from the fragment's start,
mate 2 the reverse complement of its end, fragment length ~ N(--insert);
names r<i>_<contig>_<pos1>_<strand1>_<pos2>_<strand2>_<edits>, identical
for both mates.
"""
import argparse
import random
import sys

COMP = str.maketrans("ACGTN", "TGCAN")
TRANSITION = {"A": "G", "G": "A", "C": "T", "T": "C"}

# Quality bins (value, share of bases at the read's start), by instrument.
PROFILES = {
    "novaseqx": [(40, 0.965), (24, 0.025), (9, 0.008), (2, 0.002)],
    "novaseq6000": [(37, 0.935), (25, 0.04), (11, 0.02), (2, 0.005)],
}


def mutate(seq, rate, rng):
    out = []
    for c in seq:
        r = rng.random()
        if r < rate * 0.8:
            out.append(rng.choice("ACGT".replace(c, "") or "ACGT"))
        elif r < rate * 0.9:
            out.append(c)
            out.append(rng.choice("ACGT"))
        elif r < rate:
            continue
        else:
            out.append(c)
    return "".join(out)


def cmd_genome(a):
    rng = random.Random(a.seed)
    per = a.length // a.contigs
    families = [("".join(rng.choice("ACGT") for _ in range(rng.randint(300, 6000))), rng.uniform(0.0, 0.15))
                for _ in range(40)]
    with open(a.out, "w") as f:
        for c in range(a.contigs):
            parts, size = [], 0
            while size < per:
                r = rng.random()
                if r < 0.60:
                    s = "".join(rng.choices("ACGT", k=rng.randint(2000, 20000)))
                elif r < 0.85:  # interspersed repeat copy, 0-15% diverged
                    seq, div = rng.choice(families)
                    s = mutate(seq, rng.uniform(0, div), rng)
                elif r < 0.93:  # tandem repeat
                    unit = "".join(rng.choices("ACGT", k=rng.randint(2, 60)))
                    s = unit * rng.randint(5, 200)
                elif r < 0.97 and size > 0:  # segmental duplication of earlier sequence
                    src = "".join(parts)
                    st = rng.randint(0, max(0, len(src) - 1))
                    s = mutate(src[st:st + rng.randint(1000, 30000)], rng.uniform(0, 0.02), rng)
                else:
                    s = "N" * rng.randint(100, 5000)
                parts.append(s)
                size += len(s)
            seq = "".join(parts)[:per]
            f.write(f">chr{c + 1}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + "\n")
    print(f"wrote {a.out}: {a.contigs} contigs x {per} bp", file=sys.stderr)


def read_fasta(path):
    names, seqs, cur = [], [], []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                names.append(line[1:].split()[0])
                cur = []
            else:
                cur.append(line.strip().upper())
    if cur:
        seqs.append("".join(cur))
    return names, seqs


def cmd_reads(a):
    if a.quality != "flat" or a.mates:
        return cmd_reads_profile(a)
    rng = random.Random(a.seed)
    names, seqs = read_fasta(a.ref)
    weights = [len(s) for s in seqs]
    L = a.length
    qual = "F" * L  # NovaSeq-style binned Q37
    n = 0
    with open(a.out, "w") as f:
        while n < a.n:
            c = rng.choices(range(len(seqs)), weights)[0]
            s = seqs[c]
            if len(s) < L + 20:
                continue
            p = rng.randint(0, len(s) - L - 20)
            frag = s[p:p + L + 20]
            if frag.count("N") > L // 2:
                continue
            out, edits, i = [], 0, 0
            while len(out) < L and i < len(frag):
                c0, r = frag[i], rng.random()
                if r < a.indel / 2:
                    edits += 1
                    i += 1  # deletion
                    continue
                if r < a.indel:
                    edits += 1
                    out.append(rng.choice("ACGT"))  # insertion
                    continue
                if rng.random() < a.snv + a.error:
                    edits += 1
                    out.append(rng.choice([b for b in "ACGT" if b != c0]))
                else:
                    out.append(c0)
                i += 1
            read = "".join(out)[:L]
            if len(read) < L:
                continue
            strand = "+"
            if rng.random() < 0.5:
                read = read.translate(COMP)[::-1]
                strand = "-"
            f.write(f"@r{n}_{names[c]}_{p + 1}_{strand}_{edits}\n{read}\n+\n{qual}\n")
            n += 1
    print(f"wrote {a.out}: {n} reads of {L} bp", file=sys.stderr)


def sample_haplotype(frag, a, rng):
    """The sample's sequence over a reference stretch: SNVs (ti/tv) and indels
    (10x next to a homopolymer of >= 4). Returns (sequence, variant count)."""
    out, n, i = [], 0, 0
    pti = a.titv / (1 + a.titv)
    while i < len(frag):
        c = frag[i]
        run = i >= 3 and frag[i - 3:i + 1] == c * 4
        if rng.random() < a.indel * (10 if run else 1):
            n += 1
            k = 1
            while k < 10 and rng.random() < 0.5:
                k += 1
            if rng.random() < 0.5:
                i += k  # deletion
                continue
            out.extend(c if run else rng.choice("ACGT") for _ in range(k))  # insertion
        if c in TRANSITION and rng.random() < a.snv:
            n += 1
            c = TRANSITION[c] if rng.random() < pti else rng.choice("AC" if c in "GT" else "GT")
        out.append(c)
        i += 1
    return "".join(out), n


def read_quals(L, prof, rng):
    """Binned qualities: the low bins grow towards the 3' end; ~1.5 % of reads
    end in a low-quality run."""
    q = []
    for i in range(L):
        t = i / L
        u = rng.random()
        boost = 1 + 4 * t ** 3
        acc, v = 0.0, prof[0][0]
        for val, share in reversed(prof[1:]):
            acc += share * boost
            if u < acc:
                v = val
                break
        q.append(v)
    if rng.random() < 0.015:
        tail = rng.randint(5, 30)
        for i in range(L - tail, L):
            q[i] = prof[-1][0] if rng.random() < 0.6 else prof[-2][0]
    return q


def sequence(read, quals, rng):
    """Miscalls each base with the probability its quality states."""
    out, n = [], 0
    for c, q in zip(read, quals):
        if c != "N" and rng.random() < 10 ** (-q / 10):
            c = rng.choice("ACGT".replace(c, ""))
            n += 1
        out.append(c)
    return "".join(out), n


def cmd_reads_profile(a):
    rng = random.Random(a.seed)
    names, seqs = read_fasta(a.ref)
    weights = [len(s) for s in seqs]
    L = a.length
    prof = PROFILES.get(a.quality, PROFILES["novaseqx"])
    mean, sd = (float(x) for x in a.insert.split(","))
    f2 = open(a.mates, "w") if a.mates else None
    n = 0
    with open(a.out, "w") as f1:
        while n < a.n:
            c = rng.choices(range(len(seqs)), weights)[0]
            s = seqs[c]
            F = max(L + 20, int(rng.gauss(mean, sd))) if f2 else L
            if len(s) < F + 60:
                continue
            p = rng.randint(0, len(s) - F - 60)
            frag = s[p:p + F + 60]
            if frag[:F].count("N") > L // 2:
                continue
            hap, nv = sample_haplotype(frag, a, rng)
            if len(hap) < F:
                continue
            q1 = read_quals(L, prof, rng)
            if not f2:
                r1, ne = sequence(hap[:L], q1, rng)
                strand = "+"
                if rng.random() < 0.5:
                    r1, q1, strand = r1.translate(COMP)[::-1], q1[::-1], "-"
                qs = "".join(chr(33 + q) for q in q1)
                f1.write(f"@r{n}_{names[c]}_{p + 1}_{strand}_{nv + ne}\n{r1}\n+\n{qs}\n")
            else:
                q2 = read_quals(L, prof, rng)
                left, e1 = sequence(hap[:L], q1, rng)
                right, e2 = sequence(hap[F - L:F].translate(COMP)[::-1], q2, rng)
                p1, p2 = p + 1, p + F - L + 1
                if rng.random() < 0.5:
                    m1, m2, t = left, right, f"{p1}_+_{p2}_-"
                    s1, s2 = q1, q2
                else:
                    m1, m2, t = right, left, f"{p2}_-_{p1}_+"
                    s1, s2 = q2, q1
                name = f"r{n}_{names[c]}_{t}_{nv + e1 + e2}"
                f1.write(f"@{name}\n{m1}\n+\n{''.join(chr(33 + q) for q in s1)}\n")
                f2.write(f"@{name}\n{m2}\n+\n{''.join(chr(33 + q) for q in s2)}\n")
            n += 1
    if f2:
        f2.close()
    print(f"wrote {a.out}{' and ' + a.mates if f2 else ''}: {n} {'pairs' if f2 else 'reads'} of {L} bp "
          f"({a.quality} qualities)", file=sys.stderr)


def cmd_eval(a):
    total, correct, by_mapq, by_tier = 0, 0, {}, {}
    with open(a.sam) as f:
        for line in f:
            if line.startswith("@"):
                continue
            t = line.split("\t")
            fields = t[0].split("_")
            if len(fields) >= 7:  # pair: r<i>_<contig>_<pos1>_<s1>_<pos2>_<s2>_<edits>
                chrom = "_".join(fields[1:-5])
                mate2 = int(t[1]) & 128
                pos = fields[-3] if mate2 else fields[-5]
            else:
                _, chrom, pos, _strand, _edits = t[0].rsplit("_", 4)
            ok = t[2] == chrom and abs(int(t[3]) - int(pos)) <= a.window
            mq = int(t[4])
            b = by_mapq.setdefault(mq, [0, 0])
            b[0] += 1
            b[1] += ok
            tier = next((x[5:] for x in t[11:] if x.startswith("XT:Z:")), "?").strip()
            bt = by_tier.setdefault(tier, [0, 0])
            bt[0] += 1
            bt[1] += ok
            total += 1
            correct += ok
    print(f"certified reads: {total}, placed at truth (+-{a.window} bp): {correct} "
          f"({100.0 * correct / max(total, 1):.4f}%)")
    for tier in sorted(by_tier):
        n, ok = by_tier[tier]
        print(f"  tier {tier:2s}: {n:9d} reads, {100.0 * ok / n:8.4f}% at truth")
    for mq in sorted(by_mapq, reverse=True):
        n, ok = by_mapq[mq]
        print(f"  MAPQ {mq:2d}: {n:9d} reads, {100.0 * ok / n:8.4f}% at truth")
    # Calibration: misplaced reads expected from MAPQ (sum of 10^(-MAPQ/10))
    # against those observed, in MAPQ bands.
    for lo, hi in ((0, 9), (10, 19), (20, 29), (30, 39), (40, 49), (50, 60)):
        n = sum(by_mapq[q][0] for q in by_mapq if lo <= q <= hi)
        if not n:
            continue
        bad = sum(by_mapq[q][0] - by_mapq[q][1] for q in by_mapq if lo <= q <= hi)
        exp = sum(by_mapq[q][0] * 10 ** (-q / 10) for q in by_mapq if lo <= q <= hi)
        print(f"  MAPQ {lo:2d}-{hi:2d}: {n:9d} reads, misplaced {bad:7d}, expected from MAPQ {exp:9.1f}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("genome")
    g.add_argument("out")
    g.add_argument("--length", type=int, default=20_000_000)
    g.add_argument("--contigs", type=int, default=4)
    g.add_argument("--seed", type=int, default=1)
    r = sub.add_parser("reads")
    r.add_argument("ref")
    r.add_argument("out")
    r.add_argument("--n", type=int, default=1_000_000)
    r.add_argument("--length", type=int, default=150)
    r.add_argument("--error", type=float, default=0.002)
    r.add_argument("--snv", type=float, default=0.001)
    r.add_argument("--indel", type=float, default=0.0001)
    r.add_argument("--seed", type=int, default=2)
    r.add_argument("--quality", choices=["flat", "novaseqx", "novaseq6000"], default="flat")
    r.add_argument("--titv", type=float, default=2.1)
    r.add_argument("--mates", default="")
    r.add_argument("--insert", default="400,100")
    e = sub.add_parser("eval")
    e.add_argument("sam")
    e.add_argument("--window", type=int, default=10)
    a = ap.parse_args()
    {"genome": cmd_genome, "reads": cmd_reads, "eval": cmd_eval}[a.cmd](a)


if __name__ == "__main__":
    main()
