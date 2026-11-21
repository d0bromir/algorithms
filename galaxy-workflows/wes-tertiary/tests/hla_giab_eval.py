#!/usr/bin/env python3
"""Accuracy of wes_hla on the GIAB exomes against published HLA types, and Mendelian consistency in the trios.

    hla_giab_eval.py DIR [--ref hla_dna_seq.fa]

DIR holds SAMPLE.json (wes_hla output) and optionally SAMPLE.time (GNU time -v of the run) for the nine GIAB samples.
Truth, two-field resolution:
  * HG001 (NA12878), NA12891, NA12892: 1000 Genomes HLA types (Abi-Rached et al., PLoS One 2018; PolyPheMe on
    exomes, IMGT 3.28): ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HLA_types/
    20181129_HLA_types_full_1000_Genomes_Project_panel.txt (HLA-A, -B, -C, -DQB1, -DRB1).
  * HG002: clinical sequence-based typing (Stanford Blood Center), trio-phased; Chin et al., Nat Commun 2020
    (11:4794), Supplementary Table 4 (HLA-A, -B, -C, -DQA1, -DQB1, -DRB1; G groups).
  * HG003 / HG004: the allele each passed to HG002 (paternal / maternal haplotype of the same table): one known
    allele per gene.
  * HG005 trio: no published types; checked only for Mendelian consistency.
With --ref (T1K's allele reference, built by t1k-build.pl from hla.dat), a call also counts as matching when its
antigen-recognition exons (class I exons 2+3, class II exon 2) equal those of the truth allele: the same G group,
the resolution that exome reads (and the clinical typing of HG002) define.
"""
import argparse
import collections
import importlib.util
import json
import os
import re

TRUTH = {
    "HG001": {"A": ["01:01", "11:01"], "B": ["08:01", "56:01"], "C": ["01:02", "07:01"], "DQB1": ["02:01", "05:01"],
              "DRB1": ["01:01", "03:01"]},
    "NA12891": {"A": ["01:01", "24:02"], "B": ["07:02", "08:01"], "C": ["07:01", "07:02"],
                "DQB1": ["02:01", "06:02"], "DRB1": ["03:01", "15:01"]},
    "NA12892": {"A": ["02:01", "11:01"], "B": ["15:01", "56:01"], "C": ["01:02", "04:01"],
                "DQB1": ["05:01", "05:01"], "DRB1": ["01:01", "01:01"]},
    "HG002": {"A": ["01:01", "26:01"], "B": ["35:08", "38:01"], "C": ["04:01", "12:03"], "DQA1": ["01:01", "03:01"],
              "DQB1": ["05:01", "03:02"], "DRB1": ["10:01", "04:02"]},
    "HG004": {"A": ["01:01"], "B": ["35:08"], "C": ["04:01"], "DQA1": ["01:01"], "DQB1": ["05:01"], "DRB1": ["10:01"]},
    "HG003": {"A": ["26:01"], "B": ["38:01"], "C": ["12:03"], "DQA1": ["03:01"], "DQB1": ["03:02"], "DRB1": ["04:02"]},
}
SAMPLES = ["HG001", "NA12891", "NA12892", "HG002", "HG003", "HG004", "HG005", "HG006", "HG007"]
TRIOS = {"HG001": ("NA12891", "NA12892"), "HG002": ("HG003", "HG004"), "HG005": ("HG006", "HG007")}
GENES = ["A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]


def load(d, s):
    """SAMPLE.json -> {gene: ([allele1, allele2] without the gene, confidence)} of the typed genes."""
    p = os.path.join(d, f"{s}.json")
    if not os.path.exists(p):
        return None
    out = {}
    for r in json.load(open(p))["genotypes"]:
        if r["confidence"] != "none":
            out[r["gene"].replace("HLA-", "")] = ([r["allele1"].split("*")[1], r["allele2"].split("*")[1]],
                                                  r["confidence"])
    return out


def matched(truth, pred, same=lambda a, b: a == b):
    """Number of truth alleles found in the predicted genotype (multiset)."""
    left = list(pred)
    n = 0
    for a in truth:
        hit = next((x for x in left if same(a, x)), None)
        if hit is not None:
            left.remove(hit)
            n += 1
    return n


def ars_equivalence(ref, gene_alleles):
    """same(a, b) on 'GENE*xx:yy' names: equal, or equal antigen-recognition exons (wes_hla.ars_sequences)."""
    here = os.path.dirname(os.path.abspath(__file__))
    spec = importlib.util.spec_from_file_location("wes_hla", os.path.join(here, "..", "tools", "wes_hla", "wes_hla.py"))
    w = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(w)
    ars = w.ars_sequences(ref, ["HLA-" + a for a in gene_alleles])

    def same(a, b):
        x = ars.get("HLA-" + a)
        return a == b or (x is not None and x == ars.get("HLA-" + b))
    return same


def usage(d, s):
    p = os.path.join(d, f"{s}.time")
    if not os.path.exists(p):
        return "", ""
    t = open(p).read()
    w = re.search(r"Elapsed \(wall clock\) time.*: (\S+)", t)
    m = re.search(r"Maximum resident set size \(kbytes\): (\d+)", t)
    return (w.group(1) if w else ""), (f"{int(m.group(1)) / 2 ** 20:.1f} GiB" if m else "")


def main(d, ref=None):
    calls = {s: load(d, s) for s in SAMPLES}
    same = None
    if ref:
        names = {f"{g}*{a}" for t in TRUTH.values() for g, al in t.items() for a in al}
        for p in calls.values():
            names |= {f"{g}*{a}" for g, (al, _c) in (p or {}).items() for a in al}
        same = ars_equivalence(ref, sorted(names))
    # per gene: two-field found, truth alleles, G group found, found among high-confidence calls, their truth alleles
    per_gene = collections.defaultdict(lambda: [0, 0, 0, 0, 0])
    print("sample\ttwo-field match\tsame G group\thigh-confidence calls: match\tlow-confidence genes\ttime\tmemory"
          "\tdifferences")
    for s in SAMPLES:
        p = calls[s]
        if p is None:
            print(f"{s}\t(missing)")
            continue
        f = t = fg = fh = th = 0
        diffs = []
        for g, tr in TRUTH.get(s, {}).items():
            pred, conf = p.get(g, (["", ""], "none"))
            k = matched(tr, pred)
            kg = matched([f"{g}*{a}" for a in tr], [f"{g}*{a}" for a in pred], same) if same else k
            x = per_gene[g]
            x[0] += k
            x[1] += len(tr)
            x[2] += kg
            f, t, fg = f + k, t + len(tr), fg + kg
            if conf == "high":
                x[3] += k
                x[4] += len(tr)
                fh, th = fh + k, th + len(tr)
            if k < len(tr):
                diffs.append(f"{g}: truth {'/'.join(tr)}, called {'/'.join(pred)} ({conf}"
                             + (", same G group" if kg == len(tr) else "") + ")")
        low = sum(1 for g, (_a, c) in p.items() if c == "low" and g in GENES)
        w, m = usage(d, s)
        print(f"{s}\t{f}/{t}\t{fg}/{t}\t{fh}/{th}\t{low}\t{w}\t{m}\t" + "; ".join(diffs))
    print("\ngene\ttwo-field match\tsame G group\thigh-confidence calls: match")
    for g in GENES:
        if g in per_gene:
            x = per_gene[g]
            print(f"{g}\t{x[0]}/{x[1]}\t{x[2]}/{x[1]}\t{x[3]}/{x[4]}")
    tot = [sum(x[i] for x in per_gene.values()) for i in range(5)]

    def pc(a, b):
        return f"{a}/{b} ({100 * a / max(b, 1):.1f} %)"
    print(f"all\t{pc(tot[0], tot[1])}\t{pc(tot[2], tot[1])}\t{pc(tot[3], tot[4])}")
    print("\nMendelian consistency (child shares one allele with each parent; two fields):")
    for child, (fa, mo) in TRIOS.items():
        c, x, y = calls[child], calls[fa], calls[mo]
        if not (c and x and y):
            print(f"  {child}: missing")
            continue
        ok, bad = 0, []
        for g in GENES:
            if g not in c or g not in x or g not in y:
                continue
            (c1, c2), (p1, p2), (m1, m2) = c[g][0], x[g][0], y[g][0]
            if (c1 in (p1, p2) and c2 in (m1, m2)) or (c2 in (p1, p2) and c1 in (m1, m2)):
                ok += 1
            else:
                bad.append(f"{g} child {c1}/{c2}, father {p1}/{p2}, mother {m1}/{m2}")
        print(f"  {child} ({fa} x {mo}): {ok} of {ok + len(bad)} genes consistent" + ("; " + "; ".join(bad) if bad else ""))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dir")
    ap.add_argument("--ref", help="T1K allele reference (hla_dna_seq.fa) for G-group equivalence")
    a = ap.parse_args()
    main(a.dir, a.ref)
