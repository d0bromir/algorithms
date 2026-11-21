#!/usr/bin/env python3
"""Evaluate KIR gene content from exomes (wes_kir) against genomes and assemblies.

1000 Genomes panel (--panel, --res): the same 104 males as the haplogroup panel (tests/haplogroup_panel_1kg.tsv).
For each, wes_kir ran on the KIR reads of the GRCh38 exome CRAM (--preset kir-wes) and of the 30x genome CRAM of the
same person (--preset kir-wgs), both from chr19:54.6-55.0 Mb and the chr19 alternative haplotypes (KIR haplotypes
are in the alts of GRCh38DH), with --keep (T1K's raw calls in {sample}.{ex,wgs}.raw/t1k_genotype.tsv).
Truth = gene present in the genome: best allele with T1K quality >= 1 and abundance >= --truth-abundance of the
framework genes' mean. The exome presence rule (wes_kir MIN_QUALITY, MIN_ABUNDANCE) is chosen on the tuning half
(even panel lines) by accuracy over all gene calls, and reported on the held-out half (odd lines).

GIAB (--giab DIR with {sample}.raw/t1k_genotype.tsv, --skirt DIR with the SKIRT *.hap.csv of HG002 and HG005, Zenodo
8094803, CC BY 4.0): gene content against the HPRC assemblies, and trio consistency (a child's gene must be in a
parent; a gene in neither parent is an error; a gene absent in the child although a parent is homozygous for it
cannot be checked without copy number, so only the first is counted).

    kir_eval.py --wes-kir tools/wes_kir/wes_kir.py --panel panel.tsv --res DIR [--giab DIR --skirt DIR]
"""
import argparse
import collections
import csv
import glob
import importlib.util
import os
import statistics

TRIOS = [("HG002", "HG003", "HG004"), ("HG005", "HG006", "HG007"), ("HG001", "NA12891", "NA12892")]
COLS = {"3DL3", "2DS2", "2DL2", "2DL3", "2DL5B", "2DS3", "2DS5", "2DP1", "2DL1", "3DP1", "2DL4", "3DL1", "3DS1",
        "2DL5A", "2DS1", "2DS4", "3DL2"}


def load_tool(path):
    spec = importlib.util.spec_from_file_location("wes_kir", path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def present(k, raw, min_q, min_ab):
    if not os.path.exists(raw):
        return None
    rows, _ = k.genotypes(k.parse_genotype(raw), min_q, min_ab)
    return {r["gene"] for r in rows if r["present"]}


def skirt(d, sample):
    """Genes on either haplotype of the SKIRT annotation (any allele name, novel '#' alleles and fusions included)."""
    genes = set()
    for f in glob.glob(os.path.join(d, f"{sample.lower()}.*.hap.csv")):
        rows = list(csv.reader(open(f)))
        head = rows[0]
        for r in rows[1:]:
            for col, v in zip(head, r):
                if col in COLS and v.strip():
                    genes.add("KIR" + col)
                elif col == "central region" and v.strip():  # fusion such as KIR3DS1~KIR3DP1: both parts
                    genes.update("KIR" + x.split("*")[0].split("#")[0].replace("KIR", "")
                                 for x in v.split("~") if x.strip())
    return genes


def score(pairs, genes):
    """pairs: [(called set, truth set)] -> counts per gene and overall."""
    c = collections.Counter()
    per = {g: collections.Counter() for g in genes}
    for call, truth in pairs:
        for g in genes:
            k = ("TP" if g in truth else "FP") if g in call else ("FN" if g in truth else "TN")
            c[k] += 1
            per[g][k] += 1
    return c, per


def fmt(c):
    n = sum(c.values())
    sens = c["TP"] / (c["TP"] + c["FN"]) if c["TP"] + c["FN"] else float("nan")
    spec = c["TN"] / (c["TN"] + c["FP"]) if c["TN"] + c["FP"] else float("nan")
    return (f"accuracy {c['TP'] + c['TN']}/{n} ({100 * (c['TP'] + c['TN']) / n:.1f} %), sensitivity {100 * sens:.1f} %, "
            f"specificity {100 * spec:.1f} % (FP {c['FP']}, FN {c['FN']})")


def times(res, kind):
    t = []
    for f in glob.glob(os.path.join(res, f"*.{kind}.time")):
        sec, kb = open(f).read().split()[-2:]
        t.append((float(sec), int(kb)))
    if t:
        print(f"  {kind}: median {statistics.median(x[0] for x in t):.0f} s (max {max(x[0] for x in t):.0f} s), peak "
              f"memory median {statistics.median(x[1] for x in t) / 1024 ** 2:.2f} GB (max "
              f"{max(x[1] for x in t) / 1024 ** 2:.2f} GB), n={len(t)}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--wes-kir", required=True)
    ap.add_argument("--panel")
    ap.add_argument("--res")
    ap.add_argument("--truth-abundance", type=float, default=0.25)
    ap.add_argument("--giab")
    ap.add_argument("--skirt")
    a = ap.parse_args(argv)
    k = load_tool(a.wes_kir)
    genes = k.GENES

    if a.panel:
        tune, hold = [], []
        for i, r in enumerate(csv.reader(open(a.panel), delimiter="\t")):
            s = r[0]
            truth = present(k, os.path.join(a.res, f"{s}.wgs.raw", "t1k_genotype.tsv"), 1, a.truth_abundance)
            raw = os.path.join(a.res, f"{s}.ex.raw", "t1k_genotype.tsv")
            if truth is None or not os.path.exists(raw):
                continue
            (tune if i % 2 == 0 else hold).append((s, raw, truth))
        print(f"== 1000 Genomes exomes vs genomes: {len(tune)} tuning, {len(hold)} held-out samples")
        grid = [(q, ab) for q in (0, 1, 5, 10, 15, 20, 25, 30) for ab in (0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3)]
        best = None
        for q, ab in grid:
            c, _ = score([(present(k, raw, q, ab), t) for _, raw, t in tune], genes)
            acc = (c["TP"] + c["TN"]) / sum(c.values())
            if best is None or acc > best[0] + 1e-12:
                best = (acc, q, ab)
        _, q, ab = best
        print(f"chosen on the tuning half: min quality {q}, min abundance {ab} (tool defaults {k.MIN_QUALITY}, "
              f"{k.MIN_ABUNDANCE})")
        for name, part in (("tuning", tune), ("held-out", hold)):
            pairs = [(present(k, raw, q, ab), t) for _, raw, t in part]
            c, per = score(pairs, genes)
            print(f"{name}: {fmt(c)}")
            gt = sum(k.haplotypes([{"gene": g, "present": g in call} for g in genes])["genotype"] ==
                     k.haplotypes([{"gene": g, "present": g in t} for g in genes])["genotype"] for call, t in pairs)
            bc = sum(k.haplotypes([{"gene": g, "present": g in call} for g in genes])["b_content_group"] ==
                     k.haplotypes([{"gene": g, "present": g in t} for g in genes])["b_content_group"] for call, t in pairs)
            exact = sum(call == t for call, t in pairs)
            print(f"  gene content identical {exact}/{len(pairs)}; AA/Bx genotype {gt}/{len(pairs)}; "
                  f"B-content group {bc}/{len(pairs)}")
            if name == "held-out":
                print("  per gene (TP/FN/FP/TN): " + "; ".join(
                    f"{g[3:]} {per[g]['TP']}/{per[g]['FN']}/{per[g]['FP']}/{per[g]['TN']}" for g in genes))
        print("time and memory per sample:")
        times(a.res, "ex")
        times(a.res, "wgs")

    if a.giab:
        q, ab = k.MIN_QUALITY, k.MIN_ABUNDANCE
        calls = {}
        for f in glob.glob(os.path.join(a.giab, "*.raw", "t1k_genotype.tsv")):
            s = os.path.basename(os.path.dirname(f)).split(".")[0]
            calls[s] = present(k, f, q, ab)
        print(f"== GIAB IDT exomes (min quality {q}, min abundance {ab})")
        if a.skirt:
            for s in ("HG002", "HG005"):
                if s in calls:
                    t = skirt(a.skirt, s)
                    c, _ = score([(calls[s], t)], genes)
                    print(f"{s} vs SKIRT assembly: {fmt(c)}; missed {sorted(t - calls[s])}, extra {sorted(calls[s] - t)}")
        for child, fa, mo in TRIOS:
            if all(x in calls for x in (child, fa, mo)):
                bad = sorted(calls[child] - calls[fa] - calls[mo])
                print(f"trio {child}: {len(calls[child])} genes, {len(calls[child]) - len(bad)} found in a parent"
                      + (f"; not in either parent: {bad}" if bad else ""))
        for s in sorted(calls):
            hp = k.haplotypes([{"gene": g, "present": g in calls[s]} for g in genes])
            print(f"{s}\t{len(calls[s])} genes\t{hp['genotype']}\tCen {hp['centromeric']}\tTel {hp['telomeric']}\t"
                  f"{hp['b_content_group']}")
        times(a.giab, "giab")


if __name__ == "__main__":
    main()
