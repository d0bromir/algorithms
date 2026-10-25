#!/usr/bin/env python3
"""Splits single-end reads into classes and writes one FASTQ of `size` reads
per class (the class's reads repeated under unique names when it has fewer),
for timing CERTA per class.
  read_classes.py READS.fq OUT_DIR SIZE [uncertified.fq ...]
Classes: random (every k-th read), polyG (>= 10 3' G), adapter (AGATCGGAAGAGC),
lowq (>= 5 of the last 10 bases below Q20), n (an N), short (< 146 bp), and
one class per uncertified FASTQ given (named after the file).
"""
import os
import sys

ADAPTER = "AGATCGGAAGAGC"


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


def classes(s, q):
    out = []
    if len(s) - len(s.rstrip("G")) >= 10:
        out.append("polyG")
    if ADAPTER in s:
        out.append("adapter")
    if len(q) >= 10 and sum(ord(c) - 33 < 20 for c in q[-10:]) >= 5:
        out.append("lowq")
    if "N" in s:
        out.append("n")
    if len(s) < 146:
        out.append("short")
    return out


def write(path, recs, size):
    with open(path, "w") as f:
        for i in range(size):
            name, s, q = recs[i % len(recs)]
            f.write(f"@{name}_{i // len(recs)}\n{s}\n+\n{q}\n")


def main():
    reads, out, size = sys.argv[1], sys.argv[2], int(sys.argv[3])
    os.makedirs(out, exist_ok=True)
    total = 0
    by = {"random": []}
    for name, s, q in records(reads):
        total += 1
        for c in classes(s, q):
            by.setdefault(c, []).append((name, s, q))
    step = max(1, total // size)
    by["random"] = [r for i, r in enumerate(records(reads)) if i % step == 0][:size]
    for u in sys.argv[4:]:
        by[os.path.basename(u).rsplit(".", 1)[0]] = list(records(u))
    for c, recs in sorted(by.items()):
        if recs:
            write(os.path.join(out, c + ".fq"), recs, size)
        print(f"{c}\t{len(recs)}\t{100 * len(recs) / total:.3f} %")
    print(f"total\t{total}")


if __name__ == "__main__":
    main()
