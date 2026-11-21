#!/usr/bin/env python3
"""Per-base gnomAD context for the clinical report, reduced to runs inside the exome regions.

    build_gnomad_sites.py coverage REGIONS.bed OUT.tsv   < gnomad.exomes.v4.0.coverage.summary.tsv (decompressed)
    build_gnomad_sites.py pext     REGIONS.bed OUT.tsv   < gnomad.pext.gtex_v10.base_level.tsv (decompressed)

coverage  gnomAD v4.0 exome coverage summary (730,947 exomes; gs://gcp-public-data--gnomad/release/4.0/coverage/
          exomes/gnomad.exomes.v4.0.coverage.summary.tsv.bgz, 4 GB, CC0): mean depth (rounded) and the share of
          exomes with >= 20 reads (rounded to 0.05), as runs: chrom start end mean over_20
pext      gnomAD v4.1 base-level pext from GTEx v10 (release/4.1/pext/gnomad.pext.gtex_v10.base_level.tsv.gz, 198 MB,
          CC0): the share of the gene's expression that includes the base, averaged over tissues (rounded to 0.05),
          as runs per gene: chrom start end gene pext pext_max_tissue (the highest share in any single tissue);
          plus '# gene_max' lines with each gene's highest mean value
REGIONS.bed is the workflow's regions BED (targets plus coding exons; Ensembl or chr names), merged or not.
Positions are 1-based inclusive. Standard library only; streams the input once.
"""
import bisect
import collections
import math
import sys


def norm(c):
    return c[3:] if c.startswith("chr") else c


def load_regions(path):
    raw = collections.defaultdict(list)
    for line in open(path):
        t = line.split("\t")
        if len(t) >= 3 and not line.startswith(("#", "track", "browser")):
            raw[norm(t[0])].append((int(t[1]) + 1, int(t[2])))
    out = {}
    for c, iv in raw.items():
        iv.sort()
        merged = []
        for s, e in iv:
            if merged and s <= merged[-1][1] + 1:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        out[c] = ([s for s, _ in merged], [e for _, e in merged])
    return out


def inside(regions, c, p):
    r = regions.get(c)
    if not r:
        return False
    i = bisect.bisect_right(r[0], p) - 1
    return i >= 0 and p <= r[1][i]


def r05(x):
    return round(round(x * 20) / 20, 2)


def main(kind, regions_path, out_path):
    regions = load_regions(regions_path)
    head = sys.stdin.readline().rstrip("\n").split("\t")
    ix = {h: i for i, h in enumerate(head)}
    run = None  # [chrom, start, end, values...]
    gene_max = collections.defaultdict(float)
    n_in = n_out = 0
    with open(out_path, "w", newline="\n") as o:
        o.write(f"# gnomAD {kind} runs inside {regions_path.rsplit('/', 1)[-1]}\n")

        def flush():
            nonlocal n_out
            if run:
                o.write("\t".join(map(str, run)) + "\n")
                n_out += 1

        for line in sys.stdin:
            t = line.rstrip("\n").split("\t")
            c, _, p = t[0].partition(":")
            c, p = norm(c), int(p)
            if not inside(regions, c, p):
                continue
            n_in += 1
            if kind == "coverage":
                vals = [int(round(float(t[ix["mean"]]))), r05(float(t[ix["over_20"]]))]
            else:
                x = t[ix["exp_prop_mean"]]
                x = float(x) if x not in ("NaN", "nan", "") else float("nan")
                if math.isnan(x):
                    continue
                g = t[ix["gene_symbol"]]
                tissues = [float(y) for y in t[ix["exp_prop_mean"] + 1:] if y not in ("NaN", "nan", "")]
                mx = max(tissues) if tissues else x
                gene_max[g] = max(gene_max[g], x)
                vals = [g, r05(x), r05(mx)]
            if run and run[0] == c and run[2] == p - 1 and run[3:] == vals:
                run[2] = p
            else:
                flush()
                run = [c, p, p] + vals
        flush()
        for g, x in sorted(gene_max.items()):
            o.write(f"# gene_max\t{g}\t{r05(x)}\n")
    print(f"{n_in:,} bases in the regions -> {n_out:,} runs", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) != 4 or sys.argv[1] not in ("coverage", "pext"):
        raise SystemExit(__doc__)
    main(*sys.argv[1:])
