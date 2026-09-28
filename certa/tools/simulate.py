#!/usr/bin/env python3
"""Synthetic data for exercising CERTA (standard library only).

  simulate.py genome OUT.fa  [--length 20000000] [--contigs 4] [--seed 1]
  simulate.py reads  REF.fa OUT.fq [--n 1000000] [--length 150]
                     [--error 0.002] [--snv 0.001] [--indel 0.0001] [--seed 2]
  simulate.py eval   OUT.sam [--window 10]

The genome contains repeat families (diverged copies), segmental
duplications, tandem repeats and N gaps, so that repeats and caps are
exercised. Read names carry the truth: r<i>_<contig>_<pos1>_<strand>_<edits>.
Real-data validation (GIAB, NovaSeq X, ONT) is described in
RESEARCH_ROADMAP.md; synthetic data only checks mechanics and speed.
"""
import argparse
import random
import sys

COMP = str.maketrans("ACGTN", "TGCAN")


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


def cmd_eval(a):
    total, correct, by_mapq = 0, 0, {}
    with open(a.sam) as f:
        for line in f:
            if line.startswith("@"):
                continue
            t = line.split("\t")
            _, chrom, pos, _strand, _edits = t[0].rsplit("_", 4)
            ok = t[2] == chrom and abs(int(t[3]) - int(pos)) <= a.window
            mq = int(t[4])
            b = by_mapq.setdefault(mq, [0, 0])
            b[0] += 1
            b[1] += ok
            total += 1
            correct += ok
    print(f"certified reads: {total}, placed at truth (+-{a.window} bp): {correct} "
          f"({100.0 * correct / max(total, 1):.4f}%)")
    for mq in sorted(by_mapq, reverse=True):
        n, ok = by_mapq[mq]
        print(f"  MAPQ {mq:2d}: {n:9d} reads, {100.0 * ok / n:8.4f}% at truth")


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
    e = sub.add_parser("eval")
    e.add_argument("sam")
    e.add_argument("--window", type=int, default=10)
    a = ap.parse_args()
    {"genome": cmd_genome, "reads": cmd_reads, "eval": cmd_eval}[a.cmd](a)


if __name__ == "__main__":
    main()
