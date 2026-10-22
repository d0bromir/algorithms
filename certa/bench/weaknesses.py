#!/usr/bin/env python3
"""Measure, against CERTA certificates, where other aligners provably fail.

  weaknesses.py CERTA.sam TOOL=tool.sam [TOOL=tool.sam ...]

A tool's primary record counts. TOOL:all=tool.sam takes instead the
best-scoring of each read's records, for all-mappers whose primary flag does
not mark their best (RazerS 3 flags its first match).

For every read CERTA certifies (tiers S0/S1/SR) the best end-to-end locus
is proven, together with whether it is unique or tied; for tier SL the
maximum local (clipped) score is proven. Per tool it counts:

  missed better locus  the tool reports a different locus (>10 bp away, or
                       another contig) whose alignment scores at least one
                       edit (5 points) worse, under bwa-mem's own default
                       scheme (+1/-4, gap -6-1k, clip -5 per end), than the
                       alignment CERTA reports at the certified locus;
  unmapped             the tool reports nothing although an alignment with
                       at most R edits (R <= 5; pass V: R <= 12; SL: with the
                       certified score) provably exists, split below by the
                       proven edit count (a tool searching up to e edits
                       cannot report a read whose optimum has more);
  overconfident ties   the read has >= 2 identical exact copies (a proven
                       tie under any scoring), yet the tool gives MAPQ >= 20
                       (MAPQ 255, "unavailable", never counts: RazerS 3).

Scores for every tool are recomputed from CIGAR and NM with the same scheme,
so tools with different native scores are compared on one scale. A second
table splits the first two counts by certificate: end-to-end (S0, S1, SR: the
minimum edit distance XE <= 5 proven) or local (SL, PR: the best clipped
score proven, which an end-to-end mapper may not reach). For end-to-end
certificates it also counts, independently of any scoring, the reads a tool
places elsewhere with an unclipped alignment of more edits than XE.
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


def load(path, best_of_all=False):
    """The primary record of each read; with best_of_all, the best-scoring of
    its primary and secondary records."""
    d = {}
    with open(path) as f:
        for line in f:
            if line[0] == "@":
                continue
            t = line.rstrip("\n").split("\t")
            flag = int(t[1])
            if flag & (0x800 if best_of_all else 0x900):
                continue
            key = t[0] + ("/1" if flag & 64 else "/2" if flag & 128 else "")  # mates share a name
            if flag & 4:
                d.setdefault(key, None)
                continue
            tags = {x[:2]: x[5:] for x in t[11:]}
            nm = int(tags.get("NM", 0))
            mapq = int(t[4])
            r = {
                "chrom": t[2], "pos": int(t[3]), "mapq": -1 if mapq == 255 else mapq,  # 255: unavailable
                "score": bwa_score(t[5], nm), "xt": tags.get("XT"),
                "xb": int(tags.get("XB", 1)), "xe": int(tags.get("XE", nm)),
                "nm": nm, "clipped": "S" in t[5] or "H" in t[5],
            }
            old = d.get(key)
            if old is None or (best_of_all and r["score"] > old["score"]):
                d[key] = r
    return d


def main():
    certa = load(sys.argv[1])
    cert = {k: v for k, v in certa.items() if v and v["xt"] in ("S0", "S1", "SR", "SL", "PR")}
    ties = {k for k, v in cert.items() if v["xe"] == 0 and (v["xt"] == "SR" or v["xb"] > 1)}
    n_reads = len(certa)
    print(f"reads {n_reads}; CERTA-certified {len(cert)} ({100 * len(cert) / n_reads:.2f}%); "
          f"proven exact multi-copy ties {len(ties)}")
    print(f"{'tool':12s} {'missed better locus':>24s} {'of which MAPQ>=20':>18s} "
          f"{'unmapped (<=R edits exist)':>28s} {'overconfident ties (MAPQ>=20)':>32s}")
    # The same counts split by certificate: end-to-end (S0, S1, SR) or local
    # (SL, PR); see the docstring.
    local_tiers = ("SL", "PR")
    n_local = sum(1 for v in cert.values() if v["xt"] in local_tiers)
    by_tier = []
    for spec in sys.argv[2:]:
        name, path = spec.split("=", 1)
        tool = load(path, name.endswith(":all"))
        worse = worse_conf = unmapped = over = 0
        local_worse = local_unmapped = more_edits = 0
        unm_by_xe = [0, 0, 0]  # end-to-end unmapped with a proven alignment of <= 5, 6, >= 7 edits
        for k, c in cert.items():
            t = tool.get(k)
            local = c["xt"] in local_tiers
            if t is None:
                unmapped += 1
                local_unmapped += local
                if not local:
                    unm_by_xe[0 if c["xe"] <= 5 else 1 if c["xe"] == 6 else 2] += 1
                continue
            same = t["chrom"] == c["chrom"] and abs(t["pos"] - c["pos"]) <= 10
            if not same and t["score"] <= c["score"] - 5:
                worse += 1
                worse_conf += t["mapq"] >= 20
                local_worse += local
            more_edits += not local and not same and not t["clipped"] and t["nm"] > c["xe"]
        for k in ties:
            t = tool.get(k)
            if t and t["mapq"] >= 20:
                over += 1
        n = len(cert)
        print(f"{name:12s} {worse:>12d} ({100 * worse / n:.3f}%) {worse_conf:>18d} "
              f"{unmapped:>14d} ({100 * unmapped / n:.3f}%) "
              f"{over:>16d} ({100 * over / max(len(ties), 1):.2f}% of ties)")
        by_tier.append((name, worse - local_worse, more_edits, unmapped - local_unmapped, local_worse,
                        local_unmapped, unm_by_xe))
    n_e2e = len(cert) - n_local
    xe_n = [sum(1 for v in cert.values() if v["xt"] not in local_tiers and f(v["xe"]))
            for f in (lambda x: x <= 5, lambda x: x == 6, lambda x: x >= 7)]
    print(f"by certificate: end-to-end (S0 S1 SR) {n_e2e} (proven edits <= 5: {xe_n[0]}, 6: {xe_n[1]}, "
          f">= 7: {xe_n[2]}), local (SL PR) {n_local}")
    print(f"{'tool':12s} {'e2e: missed better':>22s} {'e2e: elsewhere, NM > XE':>24s} {'e2e: unmapped':>20s} "
          f"{'local: missed better':>22s} {'local: unmapped':>22s}")
    pc = lambda x, n, d: f"{x:>12d} ({100 * x / max(n, 1):.{d}f}%)"
    for name, ew, em, eu, lw, lu, ux in by_tier:
        print(f"{name:12s} {pc(ew, n_e2e, 3):>22s} {pc(em, n_e2e, 3):>24s} {pc(eu, n_e2e, 3):>20s} "
              f"{pc(lw, n_local, 2):>22s} {pc(lu, n_local, 2):>22s}")
    print("end-to-end unmapped by the proven edit count: tool, <= 5, 6, >= 7")
    for name, ew, em, eu, lw, lu, ux in by_tier:
        print(f"{name:12s} {pc(ux[0], xe_n[0], 4):>22s} {pc(ux[1], xe_n[1], 2):>22s} {pc(ux[2], xe_n[2], 2):>22s}")


if __name__ == "__main__":
    main()
