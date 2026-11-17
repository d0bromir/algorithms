#!/usr/bin/env python3
"""Compare the clinical findings of a sequenced sample with those obtained from its GIAB truth VCF.

usage: compare_findings.py calls_findings.tsv truth_findings.tsv eval.bed [label]

Both TSVs are "clinical_findings" outputs of the WES workflows (the truth one from running the VCF workflow on
the GIAB v4.2.1 benchmark VCF). Only variants inside eval.bed (capture targets x GIAB benchmark regions) are
compared. A finding is matched on CHROM/POS/REF/ALT; category agreement (e.g. carrier vs diagnostic) is
reported separately because it also depends on the genotype.
"""
import bisect
import collections
import csv
import sys

GROUPS = collections.OrderedDict([
    ("ClinVar P/LP (diagnostic, review, carrier, low penetrance)",
     ("diagnostic", "pathogenic_review", "carrier", "low_penetrance")),
    ("ACMG SF (reportable or review)", ("acmg_sf_reportable", "acmg_sf_review", "acmg_sf_not")),
    ("Candidates not in ClinVar (LoF, missense, splicing)", ("candidate", "candidate_inheritance_not_met")),
    ("ClinVar drug response", ("clinvar_drug_response",)),
    ("ClinVar risk factor / association", ("clinvar_association",)),
    ("Fixed risk alleles (F5, F2, HFE)", ("risk_allele",)),
])


def load_bed(path):
    iv = collections.defaultdict(list)
    for line in open(path):
        t = line.split("\t")
        if len(t) >= 3:
            iv[t[0].removeprefix("chr")].append((int(t[1]), int(t[2])))
    out = {}
    for c, v in iv.items():
        v.sort()
        out[c] = ([a for a, _ in v], [b for _, b in v])
    return out


def inside(bed, chrom, pos):
    if chrom not in bed:
        return False
    starts, ends = bed[chrom]
    i = bisect.bisect_right(starts, pos - 1) - 1
    return i >= 0 and starts[i] < pos <= ends[i]


def load(path, bed):
    rows = {}
    pgx = {}
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["category"] in ("pgx", "apoe"):
                pgx[r["gene"]] = (r["genotype"], r["note"])
                continue
            if not r["pos"] or not inside(bed, r["chrom"], int(r["pos"])):
                continue
            key = (r["chrom"], int(r["pos"]), r["ref"], r["alt"])
            rows.setdefault(key, {})[r["category"]] = r
    return rows, pgx


def group_of(cats, members):
    return [c for c in cats if c in members or any(c.startswith(m) for m in members)]


def main():
    calls_path, truth_path, bed_path = sys.argv[1:4]
    label = sys.argv[4] if len(sys.argv) > 4 else ""
    bed = load_bed(bed_path)
    calls, pgx_c = load(calls_path, bed)
    truth, pgx_t = load(truth_path, bed)
    print(f"## Clinical concordance {label} (inside {bed_path})\n")
    print("| Finding group | In truth report | In sample report | Matched | Missed (FN) | Extra (FP) | Recall | Precision "
          "| Precision, high-confidence calls | Same category |")
    print("|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|")
    details = []
    for name, members in GROUPS.items():
        t = {k for k, v in truth.items() if group_of(v, members)}
        c = {k for k, v in calls.items() if group_of(v, members)}
        tp, fn, fp = t & c, t - c, c - t
        same = sum(1 for k in tp if set(group_of(truth[k], members)) == set(group_of(calls[k], members)))
        rec = f"{len(tp) / len(t):.3f}" if t else "-"
        prec = f"{len(tp) / len(c):.3f}" if c else "-"
        hc = {k for k in c if not next(iter(calls[k].values())).get("confidence", "high").startswith("low")}
        prec_hc = f"{len(tp & hc) / len(hc):.3f} ({len(hc - t)} FP)" if hc else "-"
        print(f"| {name} | {len(t)} | {len(c)} | {len(tp)} | {len(fn)} | {len(fp)} | {rec} | {prec} | {prec_hc} "
              f"| {same}/{len(tp)} |")
        if not name.startswith("ClinVar drug") and not name.startswith("ClinVar risk"):
            for k in sorted(fn):
                r = next(iter(truth[k].values()))
                details.append(f"- missed {name.split(' (')[0]}: {r['gene']} {r['hgvs_c']} {r['hgvs_p']} "
                               f"({','.join(truth[k])}, truth GT {r['genotype']})")
            for k in sorted(fp):
                r = next(iter(calls[k].values()))
                details.append(f"- extra {name.split(' (')[0]}: {r['gene']} {r['hgvs_c']} {r['hgvs_p']} "
                               f"({','.join(calls[k])}, GT {r['genotype']}, DP {r['depth']}, confidence: {r.get('confidence', '?')}, "
                               f"flags: {r['qc_flags'] or 'none'})")
            for k in sorted(tp):
                a, b = set(group_of(truth[k], members)), set(group_of(calls[k], members))
                if a != b:
                    r = next(iter(calls[k].values()))
                    details.append(f"- category differs {r['gene']} {r['hgvs_c']}: truth {','.join(sorted(a))}, "
                                   f"sample {','.join(sorted(b))} (GT {r['genotype']})")
    if details:
        print("\nDiscordant findings:\n" + "\n".join(details))
    if pgx_c:
        print("\nPharmacogenomics (sample report; truth report assumes reference where the truth VCF has no call):\n")
        print("| Gene | Sample genotype | Sample phenotype | Truth-VCF genotype | Truth-VCF phenotype |")
        print("|---|---|---|---|---|")
        for g in pgx_c:
            a, b = pgx_c[g], pgx_t.get(g, ("", ""))
            print(f"| {g} | {a[0]} | {a[1]} | {b[0]} | {b[1]} |")


if __name__ == "__main__":
    main()
