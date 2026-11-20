#!/usr/bin/env python3
"""Mitochondrial DNA variants and heteroplasmy from the reads of an exome (or genome).

Exome captures pull in many mtDNA reads off target (about 50x in the GIAB IDT exomes), enough to call homoplasmic
variants and heteroplasmy down to ~10 % with read support. This tool piles up the reads on the mitochondrial
chromosome (MAPQ >= 20 removes most NUMT reads, base quality >= 20), calls every allele above the heteroplasmy cut-off
with enough supporting reads on both strands, and annotates it with:
  ClinVar (CLNSIG, review status, condition), gnomAD v3.1 chrM (homoplasmic and heteroplasmic frequency among 56,434
  genomes, maximum heteroplasmy seen, haplogroup-defining flag, common low-heteroplasmy flag, MitoTIP and PON-mt-tRNA
  predictions for tRNA variants) and the rCRS gene.

Outputs: a TSV of calls (one row per allele) and a JSON summary (coverage, detection limit, counts).
Standard library plus samtools.

    wes_mito.py --bam reads.bam --reference ref.fa [--gnomad-mt gnomad.chrM.vcf.bgz] [--clinvar clinvar.vcf.gz]
                --tsv calls.tsv --json summary.json [--min-hl 0.10] [--min-alt 3] [--min-depth 20]
"""
import argparse
import gzip
import json
import re
import statistics
import subprocess
import sys

VERSION = "0.1.0"
# rCRS (NC_012920.1) genes: start, end (1-based, inclusive), symbol
MT_GENES = [(577, 647, "MT-TF"), (648, 1601, "MT-RNR1"), (1602, 1670, "MT-TV"), (1671, 3229, "MT-RNR2"),
            (3230, 3304, "MT-TL1"), (3307, 4262, "MT-ND1"), (4263, 4331, "MT-TI"), (4329, 4400, "MT-TQ"),
            (4402, 4469, "MT-TM"), (4470, 5511, "MT-ND2"), (5512, 5579, "MT-TW"), (5587, 5655, "MT-TA"),
            (5657, 5729, "MT-TN"), (5761, 5826, "MT-TC"), (5826, 5891, "MT-TY"), (5904, 7445, "MT-CO1"),
            (7446, 7514, "MT-TS1"), (7518, 7585, "MT-TD"), (7586, 8269, "MT-CO2"), (8295, 8364, "MT-TK"),
            (8366, 8572, "MT-ATP8"), (8527, 9207, "MT-ATP6"), (9207, 9990, "MT-CO3"), (9991, 10058, "MT-TG"),
            (10059, 10404, "MT-ND3"), (10405, 10469, "MT-TR"), (10470, 10766, "MT-ND4L"), (10760, 12137, "MT-ND4"),
            (12138, 12206, "MT-TH"), (12207, 12265, "MT-TS2"), (12266, 12336, "MT-TL2"), (12337, 14148, "MT-ND5"),
            (14149, 14673, "MT-ND6"), (14674, 14742, "MT-TE"), (14747, 15887, "MT-CYB"), (15888, 15953, "MT-TT"),
            (15956, 16023, "MT-TP")]
COLS = ["pos", "ref", "alt", "gene", "depth", "alt_reads", "alt_fwd", "alt_rev", "heteroplasmy", "call",
        "gnomad_af_hom", "gnomad_af_het", "gnomad_max_hl", "haplogroup_defining", "common_low_heteroplasmy",
        "trna_prediction", "clinvar", "clinvar_review", "clinvar_condition", "flags"]
INDEL = re.compile(r"[+-](\d+)")
# MT-RNR1 variants of the CPIC aminoglycoside guideline (McDermott et al., Clin Pharmacol Ther 2022)
PGX_SITES = {1095: "T>C", 1494: "C>T", 1555: "A>G"}
MIXTURE_SITES = 3  # heteroplasmic calls at common haplogroup variants that flag a second lineage


def gene_of(pos):
    hits = [g for s, e, g in MT_GENES if s <= pos <= e]
    return ",".join(hits) if hits else ("control region" if pos >= 16024 or pos <= 576 else "intergenic")


def opener(path):
    with open(path, "rb") as f:
        gz = f.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if gz else open(path)


def mt_contig(bam):
    out = subprocess.run(["samtools", "idxstats", bam], capture_output=True, text=True, check=True).stdout
    for line in out.splitlines():
        c = line.split("\t")[0]
        if c in ("MT", "chrM", "M", "chrMT", "NC_012920.1"):
            return c
    raise SystemExit("no mitochondrial contig (MT/chrM) in the BAM")


def pileup(bam, ref, contig, min_mapq, min_bq):
    """Per position: depth, base counts by strand, indel counts by strand (allele as VCF ref/alt)."""
    cmd = ["samtools", "mpileup", "-r", contig, "-f", ref, "-q", str(min_mapq), "-Q", str(min_bq), "-B", "-d", "0",
           "-a", bam]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    for line in p.stdout:
        t = line.rstrip("\n").split("\t")
        pos, rb, depth, s = int(t[1]), t[2].upper(), int(t[3]), t[4] if len(t) > 4 else ""
        counts, indels, i = {}, {}, 0
        while i < len(s):
            c = s[i]
            if c == "^":
                i += 2
                continue
            if c == "$":
                i += 1
                continue
            if c in "+-":
                m = INDEL.match(s, i)
                n = int(m.group(1))
                seq = s[m.end():m.end() + n]
                i = m.end() + n
                key = (rb + seq.upper(), rb) if c == "-" else (rb, rb + seq.upper())
                fwd = seq.isupper()
                d = indels.setdefault(key, [0, 0])
                d[0 if fwd else 1] += 1
                continue
            if c in ".,":
                b, fwd = rb, c == "."
            elif c in "ACGTNacgtn":
                b, fwd = c.upper(), c.isupper()
            else:  # '*' (deleted base) and others
                i += 1
                continue
            d = counts.setdefault(b, [0, 0])
            d[0 if fwd else 1] += 1
            i += 1
        yield pos, rb, depth, counts, indels
    p.wait()
    if p.returncode:
        raise SystemExit("samtools mpileup failed: " + p.stderr.read()[:500])


def load_vcf_mt(path, fields):
    """(pos, ref, alt) -> {field: value} for the mitochondrial records of a VCF."""
    out = {}
    if not path:
        return out
    with opener(path) as f:
        for line in f:
            if line[0] == "#":
                continue
            t = line.rstrip("\r\n").split("\t", 8)
            if t[0] not in ("MT", "chrM", "M", "NC_012920.1"):
                continue
            info = {}
            for kv in t[7].split(";"):
                k, _, v = kv.partition("=")
                if k in fields:
                    info[k] = v if v else True
            for alt in t[4].split(","):
                out[(int(t[1]), t[3], alt)] = info
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--bam", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--gnomad-mt")
    ap.add_argument("--clinvar")
    ap.add_argument("--min-hl", type=float, default=0.10, help="lowest heteroplasmy reported (default 0.10)")
    ap.add_argument("--min-alt", type=int, default=3, help="supporting reads needed (default 3)")
    ap.add_argument("--min-depth", type=int, default=20, help="depth needed to call heteroplasmy (default 20)")
    ap.add_argument("--min-depth-hom", type=int, default=5,
                    help="depth needed to call a (near-)homoplasmic variant, >= 90 %% of reads (default 5)")
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--min-bq", type=int, default=20)
    ap.add_argument("--tsv", required=True)
    ap.add_argument("--json", required=True)
    a = ap.parse_args(argv)
    contig = mt_contig(a.bam)
    gn = load_vcf_mt(a.gnomad_mt, {"AF_hom", "AF_het", "max_hl", "hap_defining_variant", "common_low_heteroplasmy",
                                  "mitotip_trna_prediction", "pon_mt_trna_prediction", "filters"})
    cv = load_vcf_mt(a.clinvar, {"CLNSIG", "CLNREVSTAT", "CLNDN"})
    depths, rows, pgx_depth = [], [], {}
    for pos, rb, depth, counts, indels in pileup(a.bam, a.reference, contig, a.min_mapq, a.min_bq):
        depths.append(depth)
        if pos in PGX_SITES:
            pgx_depth[str(pos)] = depth
        alleles = [((rb, b), v) for b, v in counts.items() if b != rb and b != "N"] + list(indels.items())
        for (ref, alt), (f, r) in alleles:
            n = f + r
            if "N" in ref or "N" in alt or not depth:
                continue  # rCRS 3107N placeholder
            hl = n / depth
            homo = hl >= 0.9 and depth >= a.min_depth_hom
            if not homo and (depth < a.min_depth or n < a.min_alt or hl < a.min_hl):
                continue
            flags = []
            if min(f, r) == 0 and n >= 6:
                flags.append("one strand only")
            if pos in range(302, 317) or pos in range(513, 526) or pos in range(16182, 16195) or pos == 3107:
                flags.append("homopolymer/artefact-prone (302-316, 513-525, 16182-16194, 3107N)")
            if hl < 0.5 and len(ref) == len(alt):
                flags.append("heteroplasmic: confirm (NUMT reads or contamination can mimic low levels)")
            g = gn.get((pos, ref, alt), {})
            c = cv.get((pos, ref, alt), {})
            trna = "; ".join(x for x in (g.get("mitotip_trna_prediction") and f"MitoTIP {g['mitotip_trna_prediction']}",
                                         g.get("pon_mt_trna_prediction") and f"PON-mt-tRNA {g['pon_mt_trna_prediction']}")
                             if x and "None" not in x)
            rows.append([pos, ref, alt, gene_of(pos), depth, n, f, r, f"{hl:.3f}",
                         "homoplasmic" if hl >= 0.9 else "heteroplasmic", g.get("AF_hom", ""), g.get("AF_het", ""),
                         g.get("max_hl", ""), "yes" if g.get("hap_defining_variant") else "",
                         "yes" if g.get("common_low_heteroplasmy") else "", trna.replace("_", " "),
                         str(c.get("CLNSIG", "")), str(c.get("CLNREVSTAT", "")).replace("_", " "),
                         str(c.get("CLNDN", "")).replace("_", " "), "; ".join(flags)])
    with open(a.tsv, "w") as f:
        f.write("\t".join(COLS) + "\n")
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")
    cov = depths or [0]
    # A second mtDNA lineage (sample contamination, or NUMT reads in a low-coverage exome) shows as heteroplasmy at
    # positions that are common homoplasmic variants of other haplogroups (gnomAD AF_hom >= 1 %), like haplocheck.
    mix = [r for r in rows if r[9] == "heteroplasmic" and r[10] not in ("", None) and float(r[10]) >= 0.01]
    for r in mix:
        r[19] = "; ".join(x for x in (r[19], "common variant of another haplogroup at heteroplasmic level") if x)
    summary = {"tool": f"wes_mito {VERSION}", "contig": contig, "length": len(depths),
               "mean_depth": round(statistics.mean(cov), 1), "median_depth": statistics.median(cov),
               "fraction_ge_5x": round(sum(d >= 5 for d in cov) / max(len(cov), 1), 4),
               "fraction_ge_20x": round(sum(d >= 20 for d in cov) / max(len(cov), 1), 4),
               "fraction_ge_100x": round(sum(d >= 100 for d in cov) / max(len(cov), 1), 4),
               "min_hl": a.min_hl, "min_alt": a.min_alt, "min_depth": a.min_depth,
               "detection_limit_at_median": round(max(a.min_hl, a.min_alt / max(statistics.median(cov), 1)), 3),
               "homoplasmic": sum(r[9] == "homoplasmic" for r in rows),
               "heteroplasmic": sum(r[9] == "heteroplasmic" for r in rows),
               "mixture_sites": len(mix),
               "mixture": "possible second mtDNA lineage (contamination or NUMT reads): heteroplasmy at "
                          f"{len(mix)} common haplogroup variant sites" if len(mix) >= MIXTURE_SITES else "",
               "pgx_depth": pgx_depth, "gnomad_mt": bool(gn), "clinvar": bool(cv)}
    with open(a.json, "w") as f:
        json.dump(summary, f, indent=1)
    print(json.dumps(summary), file=sys.stderr)


if __name__ == "__main__":
    main()
