#!/usr/bin/env python3
"""Sensitivity of couple carrier screening, by spiking annotated pathogenic alleles into two real exomes.

    couple_spike_eval.py A.annotated.vcf B.annotated.vcf BENCH.annotated.vcf BENCH_truth.tsv carrier_genes.tsv OUTDIR \\
        --report wes_clinical_report.py --couple wes_couple.py [--trials 10] [--genes 20] [--mode clinvar|novel] \\
        [--jobs 16] -- <report knowledge arguments>

BENCH.annotated.vcf is the VCF workflow's annotation of tests/acmg_clinvar_set.py's pseudo-sample (every record a
heterozygous call with SnpEff, gnomAD, ClinVar, AlphaMissense/REVEL and SpliceAI annotation). Each trial picks
--genes recessive genes (carrier_genes.tsv, AR) with >= 2 ClinVar P/LP (>= 2 stars) benchmark alleles on even
chromosomes (the held-out half of the ACMG benchmark), adds one allele to partner A's VCF and another (or the same)
to partner B's, runs the clinical report on both and the couple tool on the two findings tables, and counts the
spiked genes reported as "both carriers". --mode novel removes the ClinVar annotation (CLN* fields) of the spiked
records: detection then rests on the report's automated ACMG classification alone, as for a variant not yet in ClinVar.
Standard library plus the two tools; runs the report jobs in parallel.
"""
import argparse
import collections
import csv
import os
import random
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

TAB = "\t"


def read_vcf(path):
    head, recs = [], []
    with open(path) as f:
        for line in f:
            (head if line.startswith("#") else recs).append(line)
    return head, recs


def chrom_key(c):
    return (0, int(c)) if c.isdigit() else (1, c)


def strip_clinvar(line):
    t = line.split("\t")
    t[7] = ";".join(kv for kv in t[7].split(";") if not kv.startswith("CLN") and not kv.startswith("ALLELEID")) or "."
    return "\t".join(t)


def spiked_vcf(head, recs, extra, sample_col_from, out):
    """Write recs + extra (lines from the benchmark VCF, sample column replaced by a het call) sorted."""
    lines = list(recs)
    for e in extra:
        t = e.rstrip("\n").split("\t")
        t[8:] = ["GT:DP:AD:GQ", "0/1:40:20,20:99"]
        lines.append("\t".join(t) + "\n")
    lines.sort(key=lambda l: (chrom_key(l.split("\t", 2)[0]), int(l.split("\t", 2)[1])))
    with open(out, "w") as f:
        f.writelines(head)
        f.writelines(lines)


def main():
    argv = sys.argv[1:]
    extra_args = argv[argv.index("--") + 1:] if "--" in argv else []
    argv = argv[:argv.index("--")] if "--" in argv else argv
    ap = argparse.ArgumentParser()
    for x in ("a", "b", "bench", "truth", "carrier_genes", "outdir"):
        ap.add_argument(x)
    ap.add_argument("--report", required=True)
    ap.add_argument("--couple", required=True)
    ap.add_argument("--trials", type=int, default=10)
    ap.add_argument("--genes", type=int, default=20)
    ap.add_argument("--mode", choices=("clinvar", "novel"), default="clinvar")
    ap.add_argument("--jobs", type=int, default=16)
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)
    ar = {r["gene"] for r in csv.DictReader(open(a.carrier_genes), delimiter="\t") if "AR" in r["moi"].split(",")}
    truth = {}
    for r in csv.DictReader(open(a.truth), delimiter="\t"):
        truth[(r["chrom"], r["pos"], r["ref"], r["alt"])] = r
    bhead, brecs = read_vcf(a.bench)
    pool = collections.defaultdict(list)
    for line in brecs:
        t = line.split("\t", 5)
        r = truth.get((t[0], t[1], t[3], t[4]))
        if r and r["class"] in ("P", "LP") and int(r["stars"]) >= 2 and r["gene"] in ar and t[0].isdigit() \
                and int(t[0]) % 2 == 0:
            pool[r["gene"]].append(line)
    genes = sorted(g for g, v in pool.items() if len(v) >= 2)
    print(f"{len(genes)} recessive genes with >= 2 held-out ClinVar P/LP alleles", file=sys.stderr)
    ahead, arecs = read_vcf(a.a)
    bhead2, brecs2 = read_vcf(a.b)
    rnd = random.Random(a.seed)
    trials = []
    for i in range(a.trials):
        gs = rnd.sample(genes, a.genes)
        pa, pb = [], []
        for g in gs:
            x, y = rnd.sample(pool[g], 2)
            if a.mode == "novel":
                x, y = strip_clinvar(x), strip_clinvar(y)
            pa.append(x)
            pb.append(y)
        trials.append((i, gs, pa, pb))

    def run(i, gs, pa, pb):
        d = os.path.join(a.outdir, f"t{i}")
        os.makedirs(d, exist_ok=True)
        spiked_vcf(ahead, arecs, pa, None, os.path.join(d, "A.vcf"))
        spiked_vcf(bhead2, brecs2, pb, None, os.path.join(d, "B.vcf"))
        for p in "AB":
            subprocess.run([sys.executable, a.report, "--vcf", os.path.join(d, f"{p}.vcf"), "--sample-name", f"{p}{i}",
                            "--html", os.path.join(d, f"{p}.html"), "--tsv", os.path.join(d, f"{p}.tsv"), *extra_args],
                           check=True, capture_output=True)
        subprocess.run([sys.executable, a.couple, "--a", os.path.join(d, "A.tsv"), "--b", os.path.join(d, "B.tsv"),
                        "--carrier-genes", a.carrier_genes, "--html", os.path.join(d, "couple.html"),
                        "--tsv", os.path.join(d, "couple.tsv"), "--a-sex", "female", "--b-sex", "male"],
                       check=True, capture_output=True)
        cats = collections.defaultdict(set)
        for r in csv.DictReader(open(os.path.join(d, "couple.tsv")), delimiter=TAB):
            cats[r["category"]].add(r["gene"])
        rep = {}
        for p in "AB":
            rows = list(csv.DictReader(open(os.path.join(d, f"{p}.tsv")), delimiter=TAB))
            rep[p] = {r["gene"] for r in rows if r.get("pos")}
        return i, gs, cats, rep

    with ThreadPoolExecutor(a.jobs) as ex:
        res = list(ex.map(lambda x: run(*x), trials))
    n = hit = conf = 0
    miss = collections.Counter()
    extra = extra_conf = 0
    for i, gs, cats, rep in res:
        for g in gs:
            n += 1
            if g in cats["both_carriers"]:
                hit += 1
            elif g in cats["both_carriers_confirm"]:
                conf += 1
            else:
                miss["A not listed" if g not in rep["A"] else ("B not listed" if g not in rep["B"] else "listed, not P/LP")] += 1
        extra += len(cats["both_carriers"] - set(gs))
        extra_conf += len(cats["both_carriers_confirm"] - set(gs))
    print(f"mode {a.mode}: {hit}/{n} spiked couple genes found ({100 * hit / n:.1f} %), {conf} more to confirm "
          f"({100 * (hit + conf) / n:.1f} % in all); missed: {dict(miss)}; other couple genes from the real exomes: "
          f"{extra} found, {extra_conf} to confirm, over {len(res)} trials")


if __name__ == "__main__":
    main()
