#!/usr/bin/env python3
"""Measure, against CERTA certificates, where other aligners provably fail.

  weaknesses.py CERTA.sam TOOL=tool.sam [TOOL=tool.sam ...]

For every read CERTA certifies (tiers S0/S1/SR) the best end-to-end locus
is proven, together with whether it is unique or tied. Per tool it counts:

  missed better locus  the tool reports a different locus (>10 bp away, or
                       another contig) whose alignment scores at least one
                       edit (5 points) worse, under bwa-mem's own default
                       scheme (+1/-4, gap -6-1k, clip -5 per end), than the
                       alignment CERTA reports at the certified locus;
  unmapped             the tool reports nothing although an alignment with
                       at most R <= 5 edits provably exists;
  overconfident ties   the read has >= 2 identical exact copies (a proven
                       tie under any scoring), yet the tool gives MAPQ >= 20.

Scores for every tool are recomputed from CIGAR and NM with the same scheme,
so tools with different native scores are compared on one scale.
"""
import re
import sys

CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def bwa_score(cigar, nm):
    ops = [(int(n), o) for n, o in CIG.findall(cigar)]
    m = sum(n for n, o in ops if o in "M=X")
    gaps = [n for n, o in ops if o in "ID"]
    mism = max(0, nm - sum(gaps))
    clips = sum(1 for n, o in ops if o in "SH")
    return (m - mism) - 4 * mism - sum(6 + g for g in gaps) - 5 * clips


def load(path):
    d = {}
    with open(path) as f:
        for line in f:
            if line[0] == "@":
                continue
            t = line.rstrip("\n").split("\t")
            flag = int(t[1])
            if flag & 0x900:
                continue
            if flag & 4:
                d[t[0]] = None
                continue
            tags = {x[:2]: x[5:] for x in t[11:]}
            nm = int(tags.get("NM", 0))
            d[t[0]] = {
                "chrom": t[2], "pos": int(t[3]), "mapq": int(t[4]),
                "score": bwa_score(t[5], nm), "xt": tags.get("XT"),
                "xb": int(tags.get("XB", 1)), "xe": int(tags.get("XE", nm)),
            }
    return d


def main():
    certa = load(sys.argv[1])
    cert = {k: v for k, v in certa.items() if v and v["xt"] in ("S0", "S1", "SR")}
    ties = {k for k, v in cert.items() if v["xe"] == 0 and (v["xt"] == "SR" or v["xb"] > 1)}
    n_reads = len(certa)
    print(f"reads {n_reads}; CERTA-certified {len(cert)} ({100 * len(cert) / n_reads:.2f}%); "
          f"proven exact multi-copy ties {len(ties)}")
    print(f"{'tool':12s} {'missed better locus':>24s} {'of which MAPQ>=20':>18s} "
          f"{'unmapped (<=5 edits exist)':>28s} {'overconfident ties (MAPQ>=20)':>32s}")
    for spec in sys.argv[2:]:
        name, path = spec.split("=", 1)
        tool = load(path)
        worse = worse_conf = unmapped = over = 0
        for k, c in cert.items():
            t = tool.get(k)
            if t is None:
                unmapped += 1
                continue
            same = t["chrom"] == c["chrom"] and abs(t["pos"] - c["pos"]) <= 10
            if not same and t["score"] <= c["score"] - 5:
                worse += 1
                worse_conf += t["mapq"] >= 20
        for k in ties:
            t = tool.get(k)
            if t and t["mapq"] >= 20:
                over += 1
        n = len(cert)
        print(f"{name:12s} {worse:>12d} ({100 * worse / n:.3f}%) {worse_conf:>18d} "
              f"{unmapped:>14d} ({100 * unmapped / n:.3f}%) "
              f"{over:>16d} ({100 * over / max(len(ties), 1):.2f}% of ties)")


if __name__ == "__main__":
    main()
