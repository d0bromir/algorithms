#!/usr/bin/env python3
"""Recessive and X-linked genes for couple carrier screening, with carrier frequencies, treatability and newborn screening.

    bcftools view -i 'INFO/CLNSIG~"athogenic" && INFO/CLNSIG!~"onflicting"' clinvar.vcf.gz \\
      | bcftools annotate -a gnomad.v4.1.joint.bcf -c INFO/gnomAD_AF,INFO/gnomAD_grpmax_AF,INFO/gnomAD_grpmax,INFO/gnomAD_filter \\
      | build_carrier_genes.py gene_disease.tsv actionability_pediatric.json carrier_genes.tsv

Inputs: ClinVar's pathogenic/likely pathogenic records (>= 1 star, not conflicting) annotated with the gnomAD v4.1
joint frequencies (the workflow's extract); the curated gene-disease table (data/build_gene_disease.py); ClinGen's
pediatric actionability curations (https://actionability.clinicalgenome.org/ac/Pediatric/api/summ?flavor=flat).

Carrier frequency of a gene = 2 x the summed allele frequency of its P/LP alleles (Hardy-Weinberg, carriers of one
allele), overall and as an upper bound for the most affected ancestry group (summing each allele's highest group
frequency). Alleles at sites gnomAD filters (not PASS) and ClinVar alleles labelled low penetrance or risk factor
(HFE p.Cys282Tyr, SERPINA1 S and Z) are not counted, nor alleles above 2 % overall unless an expert panel or
practice guideline classified them (CFTR poly-T 5T is out; F508del, HBB sickle and GJB2 p.Val37Ile are in). ACMG 2021 (Gregg et al., Genet Med 23:1793) recommends screening genes with a carrier frequency >= 1/200
in any group (tier 3); this table computes that criterion from the current data instead of a fixed list.

Output columns: gene moi carrier_freq carrier_freq_max_group max_group plp_alleles in_gnomad tier3 newborn_screening
actionability conditions
Standard library only.
"""
import collections
import json
import re
import sys

# Genes of the core conditions of the US Recommended Uniform Screening Panel (HRSA, 2024): screened at birth because
# early treatment changes the outcome. SMA (SMN1) and CAH (CYP21A2) are on the panel but not callable from exome SNVs.
RUSP = {
    "PAH": "phenylketonuria", "BCKDHA": "maple syrup urine disease", "BCKDHB": "maple syrup urine disease",
    "DBT": "maple syrup urine disease", "MMUT": "methylmalonic acidemia", "MMAA": "methylmalonic acidemia (cblA)",
    "MMAB": "methylmalonic acidemia (cblB)", "MMACHC": "methylmalonic acidemia with homocystinuria (cblC)",
    "PCCA": "propionic acidemia", "PCCB": "propionic acidemia", "IVD": "isovaleric acidemia",
    "MCCC1": "3-methylcrotonyl-CoA carboxylase deficiency", "MCCC2": "3-methylcrotonyl-CoA carboxylase deficiency",
    "HLCS": "holocarboxylase synthetase deficiency", "HMGCL": "3-hydroxy-3-methylglutaric aciduria",
    "ACAT1": "beta-ketothiolase deficiency", "GCDH": "glutaric acidemia type I", "ACADM": "MCAD deficiency",
    "ACADVL": "VLCAD deficiency", "HADHA": "LCHAD / trifunctional protein deficiency",
    "HADHB": "trifunctional protein deficiency", "SLC22A5": "carnitine uptake defect", "ASS1": "citrullinemia type I",
    "ASL": "argininosuccinic aciduria", "CBS": "homocystinuria", "FAH": "tyrosinemia type I",
    "BTD": "biotinidase deficiency", "GALT": "classic galactosemia", "CFTR": "cystic fibrosis",
    "HBB": "sickle cell disease, beta thalassemia", "CYP21A2": "congenital adrenal hyperplasia",
    "GAA": "Pompe disease", "IDUA": "mucopolysaccharidosis type I", "IDS": "mucopolysaccharidosis type II",
    "ABCD1": "X-linked adrenoleukodystrophy", "SMN1": "spinal muscular atrophy",
    "IL2RG": "severe combined immunodeficiency", "RAG1": "severe combined immunodeficiency",
    "RAG2": "severe combined immunodeficiency", "ADA": "severe combined immunodeficiency",
    "JAK3": "severe combined immunodeficiency", "IL7R": "severe combined immunodeficiency",
    "DCLRE1C": "severe combined immunodeficiency", "PTPRC": "severe combined immunodeficiency",
}
STARS0 = ("no_assertion", "no_classification", "no_interpretation")
COMMON = 0.02


def main(gd_path, act_path, out_path):
    moi = collections.defaultdict(set)
    cond = collections.defaultdict(list)
    for line in open(gd_path, encoding="utf-8"):
        if line.startswith("#"):
            continue
        t = line.rstrip("\n").split("\t")
        if t[0] == "gene" or len(t) < 6:
            continue
        valid = (t[1] == "Orphanet" and t[5] == "assessed" and "germline" in t[7]) or \
            t[5] in ("definitive", "strong", "moderate")
        mois = t[4].split(",")
        if not valid or (t[1] == "Orphanet" and len(mois) > 1):
            continue
        for m in mois:
            if m in ("AR",):
                moi[t[0]].add("AR")
            elif m in ("XLR", "XL"):
                moi[t[0]].add("XL")
            elif m in ("AD", "XLD", "SD"):
                moi[t[0]].add("dominant")
        if t[2] not in cond[t[0]]:
            cond[t[0]].append(t[2])
    act = collections.defaultdict(list)
    d = json.load(open(act_path))
    for r in (dict(zip(d["columns"], x)) for x in d["rows"]):
        if r.get("status-overall") == "Released" and r.get("effectiveness"):
            act[r["geneOrVariant"]].append((r["effectiveness"], r["overall"], r["disease"]))
    af = collections.defaultdict(float)
    afmax = collections.defaultdict(float)
    grp = collections.defaultdict(collections.Counter)
    n_plp = collections.Counter()
    n_gn = collections.Counter()
    for line in sys.stdin:
        if line[0] == "#":
            continue
        t = line.split("\t", 8)
        info = dict(kv.split("=", 1) for kv in t[7].split(";") if "=" in kv)
        rev = info.get("CLNREVSTAT", "")
        sig = info.get("CLNSIG", "")
        if not rev or any(x in rev for x in STARS0) or "onflicting" in sig or "athogenic" not in sig:
            continue
        genes = [x.split(":")[0] for x in info.get("GENEINFO", "").split("|") if x]
        rec = [x for x in genes if {"AR", "XL"} & moi.get(x, set())]
        if len(rec) != 1:
            continue  # no recessive gene, or a record spanning several of them
        gene = rec[0]
        n_plp[gene] += 1
        x = info.get("gnomAD_AF")
        if info.get("gnomAD_filter", "PASS").strip() != "PASS" or re.search(r"low_penetrance|risk_factor", sig):
            x = None  # gnomAD-filtered sites are unreliable; low-penetrance/risk alleles are not carrier alleles here
        elif x not in (None, ".") and float(x.split(",")[0]) > COMMON and not re.search(r"expert_panel|practice_guideline", rev):
            x = None  # common alleles (e.g. CFTR poly-T 5T) count only with an expert-panel or guideline classification
        if x not in (None, "."):
            n_gn[gene] += 1
            af[gene] += float(x.split(",")[0])
            m = info.get("gnomAD_grpmax_AF")
            if m not in (None, "."):
                afmax[gene] += float(m.split(",")[0])
                grp[gene][info.get("gnomAD_grpmax", "?").split(",")[0]] += float(m.split(",")[0])
    rows = []
    for g, ms in moi.items():
        if not ({"AR", "XL"} & ms):
            continue
        cf, cfm = 2 * af.get(g, 0.0), 2 * afmax.get(g, 0.0)
        rows.append([g, ",".join(sorted(ms & {"AR", "XL"})) + (",+dominant" if "dominant" in ms else ""),
                     f"{cf:.3g}", f"{cfm:.3g}", grp[g].most_common(1)[0][0] if grp[g] else "", n_plp[g], n_gn[g],
                     "yes" if cfm >= 1 / 200 else "", RUSP.get(g, ""),
                     "; ".join(f"{e} ({o}) {dz}" for e, o, dz in act.get(g, [])[:2]),
                     "; ".join(cond[g][:3])])
    rows.sort(key=lambda r: -float(r[3]))
    with open(out_path, "w", encoding="utf-8", newline="\n") as o:
        o.write("gene\tmoi\tcarrier_freq\tcarrier_freq_max_group\tmax_group\tplp_alleles\tin_gnomad\ttier3\t"
                "newborn_screening\tactionability\tconditions\n")
        for r in rows:
            o.write("\t".join(re.sub(r"[\t\n]", " ", str(x)) for x in r) + "\n")
    print(f"{len(rows)} recessive/X-linked genes; {sum(r[7] == 'yes' for r in rows)} with carrier frequency >= 1/200 "
          f"in some group; {sum(bool(r[8]) for r in rows)} newborn-screened; {sum(bool(r[9]) for r in rows)} with "
          "pediatric actionability", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        raise SystemExit(__doc__)
    main(*sys.argv[1:])
