"""Simulated paired-end reads with known origin, for MAPQ calibration (Aim 3).

  simulate_pairs.py ref.fa N out_1.fq out_2.fq [seed]

Fragments are drawn uniformly from the primary chromosomes (1-22, X, Y, or
chr1-chr22, chrX, chrY), length ~ N(400, 80) within [200, 800]; fragments
containing N are drawn again. Each 150-base read differs from the reference
independently (no shared haplotype):
  - variants: SNVs at 0.1 % of bases, indels at 0.01 % (length 1-10,
    geometric);
  - sequencing errors: substitutions at 0.2 % of bases; 3 % of reads at 2 %.
Mate 1 is the fragment's left end on the forward strand, or (half of the
fragments) its right end reverse-complemented; mate 2 is the other end.
Read names: sim<i>:<contig>:<start1>:<strand1>:<start2>:<strand2>, with the
1-based leftmost reference position of each mate's span (indels near a read's
start can shift its alignment by a few bases).
"""
import random
import sys

L = 150
COMP = bytes.maketrans(b"ACGTN", b"TGCAN")
PRIMARY = {str(i) for i in range(1, 23)} | {"X", "Y"}


def load(path):
    seqs, name, chunks = {}, None, []
    with open(path, "rb") as f:
        for line in f:
            if line.startswith(b">"):
                if name in PRIMARY:
                    seqs[name] = b"".join(chunks).upper()
                name = line[1:].split()[0].decode()
                name = name[3:] if name.startswith("chr") else name
                chunks = []
            elif name in PRIMARY:
                chunks.append(line.rstrip())
    if name in PRIMARY:
        seqs[name] = b"".join(chunks).upper()
    return seqs


def mutate(rng, s, err):
    """Indels (variants), then substitutions (variants + errors), on a list of bases."""
    i = 0
    out = bytearray()
    while i < len(s):
        if rng.random() < 1e-4:
            k = 1
            while k < 10 and rng.random() < 0.5:
                k += 1
            if rng.random() < 0.5:
                i += k  # deletion
                continue
            out.extend(rng.choice(b"ACGT") for _ in range(k))  # insertion
        out.append(s[i])
        i += 1
    rate = 0.001 + err
    for j in range(len(out)):
        if rng.random() < rate:
            out[j] = rng.choice([b for b in b"ACGT" if b != out[j]])
    return bytes(out)


def main():
    ref, n, o1, o2 = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4]
    rng = random.Random(int(sys.argv[5]) if len(sys.argv) > 5 else 7)
    seqs = load(ref)
    names = sorted(seqs)
    weights = [len(seqs[c]) for c in names]
    qual = b"I" * L
    with open(o1, "wb") as f1, open(o2, "wb") as f2:
        i = 0
        while i < n:
            c = rng.choices(names, weights)[0]
            s = seqs[c]
            F = max(200, min(800, int(rng.gauss(400, 80))))
            a = rng.randrange(0, len(s) - F - 40)
            frag = s[a:a + F + 40]
            if b"N" in frag:
                continue
            err = 0.02 if rng.random() < 0.03 else 0.002
            left = mutate(rng, frag[:L + 20], err)[:L]
            err = 0.02 if rng.random() < 0.03 else 0.002
            # The right end: mutate the reversed template so indels shift its far end.
            right = mutate(rng, frag[F - L - 20:F][::-1], err)[:L][::-1]
            right_rc = right.translate(COMP)[::-1]
            p1, p2 = a + 1, a + F - L + 1
            if rng.random() < 0.5:
                r1, r2, t = left, right_rc, f"{c}:{p1}:+:{p2}:-"
            else:
                r1, r2, t = right_rc, left, f"{c}:{p2}:-:{p1}:+"
            name = f"sim{i}:{t}".encode()
            f1.write(b"@" + name + b"\n" + r1 + b"\n+\n" + qual + b"\n")
            f2.write(b"@" + name + b"\n" + r2 + b"\n+\n" + qual + b"\n")
            i += 1


if __name__ == "__main__":
    main()
