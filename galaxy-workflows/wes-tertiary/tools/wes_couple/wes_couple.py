#!/usr/bin/env python3
"""Couple (reproductive) carrier screening from the findings tables of two exome reports.

Each partner's "Clinical findings (table)" of the clinical report (wes_clinical_report >= 0.13) lists their pathogenic
and likely pathogenic variants and their automated ACMG classes. This tool finds:
  * autosomal recessive genes where both partners carry a pathogenic or likely pathogenic allele: 1 in 4 risk per
    pregnancy of an affected child (1 in 2 if a partner already has two alleles);
  * X-linked genes where the female partner carries one: 1 in 2 sons affected, 1 in 2 daughters carriers;
  * genes where one partner is a carrier, with the other partner's chance of being an undetected carrier: the
    population carrier frequency of the gene (data/build_carrier_genes.py) times the share of carriers an exome
    misses (--miss), so that targeted follow-up can be offered where it matters;
  * for each, whether the condition is treatable: on the US newborn screening panel (RUSP) or scored by ClinGen's
    pediatric actionability curation (effectiveness of intervention 2 or 3).
High-confidence CNV deletions in the reports (exons of a recessive gene) count as carrier alleles too.

    wes_couple.py --a A.tsv --b B.tsv --carrier-genes carrier_genes.tsv --html couple.html --tsv couple.tsv
        [--a-name NAME] [--b-name NAME] [--a-sex male|female] [--b-sex male|female] [--miss 0.1]

Standard library only.
"""
import argparse
import csv
import datetime
import html
import re
import sys

VERSION = "0.1.0"
PLP_ACMG = ("pathogenic", "likely pathogenic")
SKIP_CATS = ("report_meta", "pgx", "pgx_pharmcat", "apoe", "low_penetrance", "clinvar_association",
             "clinvar_drug_response", "mito_candidate", "phenotype_match")


def esc(x):
    return html.escape(str(x), quote=True)


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def load_genes(path):
    out = {}
    with open(path, encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            out[r["gene"]] = r
    return out


def is_plp(r):
    """A pathogenic or likely pathogenic allele: ClinVar P/LP (not conflicting, not low penetrance/risk factor) or the
    report's automated class P/LP. Low-confidence calls count too (the result then says to confirm them): carrier
    screening must not drop an allele in a paralogous gene (HBA2, VWF, CHEK2, NEB...) only because short reads are
    unreliable there."""
    cv = r.get("clinvar", "")
    clin = "athogenic" in cv and "onflicting" not in cv and not re.search(r"low_penetrance|risk_factor", cv)
    auto = r.get("acmg_class", "") in PLP_ACMG
    return clin or auto


def partner_alleles(path, genes):
    """gene -> [allele dict] for one partner; and the report's meta rows."""
    meta, alleles = {}, {}
    with open(path, encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    seen = set()
    for r in rows:
        cat = r.get("category", "")
        if cat == "report_meta":
            meta[r["gene"]] = r.get("note", "")
            continue
        if cat in SKIP_CATS or cat.startswith("change_"):
            continue
        if cat.startswith("cnv_"):
            if cat not in ("cnv_pathogenic", "cnv_likely_pathogenic", "cnv_high_confidence") or "DEL" not in r["genotype"]:
                continue
            for g in r["gene"].split(","):
                if g in genes and {"AR", "XL"} & set(genes[g]["moi"].split(",")):
                    k = (g, "cnv", r["chrom"], r["pos"])
                    if k not in seen:
                        seen.add(k)
                        alleles.setdefault(g, []).append({"name": f"deletion {r['chrom']}:{r['pos']} ({r['genotype']})",
                                                          "zyg": "het", "evidence": "CNV " + cat[4:], "key": k})
            continue
        g = r.get("gene", "")
        if not g or g not in genes or not r.get("pos") or not is_plp(r):
            continue
        k = (g, r["chrom"], r["pos"], r["ref"], r["alt"])
        if k in seen:
            continue
        seen.add(k)
        ev = []
        if "athogenic" in r.get("clinvar", ""):
            ev.append("ClinVar " + r["clinvar"].replace("_", " ") + (f" ({r['clinvar_review']})" if r.get("clinvar_review") else ""))
        if r.get("acmg_class"):
            ev.append(f"automated {r['acmg_class']} ({r.get('acmg_codes', '')})")
        zyg = "hom" if r.get("genotype", "") in ("1/1", "1|1") else "het"
        if r.get("confidence", "high").startswith("low"):
            ev.append("CONFIRM: " + r["confidence"][5:] + " (low-confidence call)")
        alleles.setdefault(g, []).append({"name": (r.get("hgvs_p") or r.get("hgvs_c") or f"{r['chrom']}:{r['pos']}"),
                                          "zyg": zyg, "evidence": "; ".join(ev), "key": k,
                                          "refseq": r.get("hgvs_refseq", ""), "gnomad": r.get("gnomad_af", "")})
    return alleles, meta


def treatable(gene_row):
    nbs = gene_row.get("newborn_screening", "")
    act = gene_row.get("actionability", "")
    eff = [int(m) for m in re.findall(r"(?:^|; )([0-3])[A-E]", act)]
    if nbs:
        return f"yes: newborn screening ({nbs})"
    if eff and max(eff) >= 2:
        return f"yes: ClinGen pediatric actionability, effectiveness {max(eff)} of 3"
    if eff:
        return "limited: ClinGen pediatric actionability, effectiveness " + str(max(eff))
    return ""


def fmt_p(p):
    return "—" if p is None else (f"1 in {round(1 / p):,}" if p > 0 else "0")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--a", required=True, help="findings TSV of partner A")
    ap.add_argument("--b", required=True, help="findings TSV of partner B")
    ap.add_argument("--a-name")
    ap.add_argument("--b-name")
    ap.add_argument("--a-sex", choices=("male", "female"))
    ap.add_argument("--b-sex", choices=("male", "female"))
    ap.add_argument("--carrier-genes", required=True, help="data/build_carrier_genes.py table")
    ap.add_argument("--miss", type=float, default=0.10,
                    help="share of carriers an exome misses (CNVs, deep-intronic, repeats; default 0.10)")
    ap.add_argument("--html", required=True)
    ap.add_argument("--tsv", required=True)
    a = ap.parse_args(argv)
    genes = load_genes(a.carrier_genes)
    A, ma = partner_alleles(a.a, genes)
    B, mb = partner_alleles(a.b, genes)
    name = {"A": a.a_name or ma.get("sample", "partner A"), "B": a.b_name or mb.get("sample", "partner B")}
    sex = {"A": a.a_sex or ma.get("sex", ""), "B": a.b_sex or mb.get("sex", "")}
    for p in "AB":
        if sex[p] not in ("male", "female"):
            sex[p] = ""
    P = {"A": A, "B": B}
    couple, xl, single = [], [], []
    for g in sorted(set(A) | set(B)):
        gr = genes[g]
        moi = set(gr["moi"].split(","))
        if "XL" in moi:
            for p in "AB":
                if g in P[p] and sex[p] == "female":
                    xl.append((g, p))
                elif g in P[p] and sex[p] == "male":
                    xl.append((g, p))
        if "AR" in moi:
            if g in A and g in B:
                two = any(len(P[p][g]) > 1 or P[p][g][0]["zyg"] == "hom" for p in "AB")
                couple.append((g, "1 in 2" if two else "1 in 4"))
            else:
                p = "A" if g in A else "B"
                cf = fnum(gr["carrier_freq_max_group"])
                cfo = fnum(gr["carrier_freq"])
                resid = (cf or 0) * a.miss
                single.append((g, p, cfo, cf, gr["max_group"], resid))
    single.sort(key=lambda x: -x[5])
    today = datetime.date.today().isoformat()
    out = [f"<!doctype html><html><head><meta charset='utf-8'><title>Couple carrier screening: {esc(name['A'])} and "
           f"{esc(name['B'])}</title><style>body{{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;"
           "margin:24px auto;max-width:1100px;color:#1d2430;line-height:1.45;padding:0 16px}table{border-collapse:collapse;"
           "width:100%;font-size:.88em;margin:8px 0 14px}th,td{border:1px solid #d8dee8;padding:4px 6px;text-align:left;"
           "vertical-align:top}th{background:#eef2f7}.box{border:1px solid #d8dee8;border-radius:8px;padding:12px 16px;"
           "background:#f7f9fc}.alert{border-color:#d9534f;background:#fdf1f0}.small{font-size:.85em;color:#5b6576}"
           "@media print{body{max-width:none;font-size:9pt}tr{break-inside:avoid}}</style></head><body>",
           f"<h1>Couple carrier screening</h1><p class='small'>{esc(name['A'])} ({esc(sex['A'] or 'sex unknown')}) and "
           f"{esc(name['B'])} ({esc(sex['B'] or 'sex unknown')}) &middot; {today} &middot; wes_couple {VERSION}</p>"]
    risky = couple + [x for x in xl if sex[x[1]] == "female"]
    out.append(f"<div class='box{' alert' if risky else ''}'><b>{len(couple)}</b> autosomal recessive gene(s) where both "
               f"partners carry a pathogenic allele; <b>{len([x for x in xl if sex[x[1]] == 'female'])}</b> X-linked "
               f"gene(s) carried by the female partner; {len(single)} gene(s) carried by one partner only.</div>")

    def alle(p, g):
        return "<br>".join(esc(f"{x['name']} ({x['zyg']})") + f"<br><span class='small'>{esc(x['evidence'])}</span>"
                           for x in P[p].get(g, [])) or "<span class='small'>none found</span>"
    out.append("<h2>1. Autosomal recessive: both partners carriers</h2>")
    rows = []
    for g, risk in couple:
        gr = genes[g]
        rows.append(f"<tr><td><b>{esc(g)}</b><br><span class='small'>{esc(gr['conditions'])}</span></td><td>{alle('A', g)}</td>"
                    f"<td>{alle('B', g)}</td><td><b>{risk}</b> per pregnancy</td><td>{esc(treatable(gr)) or '—'}</td></tr>")
    out.append("<table><tr><th>Gene, condition</th><th>" + esc(name["A"]) + "</th><th>" + esc(name["B"])
               + "</th><th>Risk of an affected child</th><th>Treatable</th></tr>" + "".join(rows) + "</table>"
               if rows else "<p class='small'>None.</p>")
    out.append("<h2>2. X-linked</h2>")
    rows = []
    for g, p in xl:
        gr = genes[g]
        risk = ("1 in 2 sons affected, 1 in 2 daughters carriers" if sex[p] == "female" else
                "the partner is hemizygous (possibly affected); all daughters will be carriers")
        rows.append(f"<tr><td><b>{esc(g)}</b><br><span class='small'>{esc(gr['conditions'])}</span></td><td>"
                    f"{esc(name[p])} ({esc(sex[p])})</td><td>{alle(p, g)}</td><td>{risk}</td><td>{esc(treatable(gr)) or '—'}</td></tr>")
    out.append("<table><tr><th>Gene, condition</th><th>Carrier</th><th>Allele</th><th>Risk</th><th>Treatable</th></tr>"
               + "".join(rows) + "</table>" if rows else "<p class='small'>None.</p>")
    out.append(f"<h2>3. One partner a carrier</h2><p class='small'>The other partner's chance of being an undetected carrier "
               f"= the gene's carrier frequency in the most affected ancestry group (gnomAD v4.1 and ClinVar; an upper "
               f"bound) &times; the share of carriers an exome misses ({a.miss:.0%}). The couple's risk is then that chance "
               f"&times; 1/4. Genes with a carrier frequency &ge; 1/200 (ACMG 2021 tier 3) first.</p>")
    rows = []
    for g, p, cfo, cf, grp, resid in single:
        gr = genes[g]
        rows.append(f"<tr><td><b>{esc(g)}</b>{' (tier 3)' if gr.get('tier3') else ''}<br><span class='small'>"
                    f"{esc(gr['conditions'])}</span></td><td>{esc(name[p])}</td><td>{alle(p, g)}</td>"
                    f"<td>{fmt_p(cfo)} overall; up to {fmt_p(cf)} ({esc(grp.upper())})</td><td>{fmt_p(resid)}</td>"
                    f"<td>{fmt_p(resid / 4) if resid else '—'}</td><td>{esc(treatable(gr)) or '—'}</td></tr>")
    out.append("<table><tr><th>Gene, condition</th><th>Carrier</th><th>Allele</th><th>Carrier frequency</th><th>Other "
               "partner: residual carrier chance</th><th>Couple risk</th><th>Treatable</th></tr>" + "".join(rows)
               + "</table>" if rows else "<p class='small'>None.</p>")
    out.append("<h2>4. Methods and limits</h2><ul class='small'>"
               "<li>Carrier alleles: ClinVar pathogenic/likely pathogenic (not conflicting, not low-penetrance or risk "
               "factor) or the report's automated P/LP class, high-confidence calls; high-confidence CNV deletions of a "
               "recessive gene's exons.</li>"
               "<li>Recessive and X-linked genes: ClinGen, Gene2Phenotype and Orphanet curations (moderate or stronger). "
               "Carrier frequencies: 2 &times; summed gnomAD v4.1 frequency of the gene's ClinVar P/LP alleles "
               "(data/build_carrier_genes.py).</li>"
               "<li>Not detected by an exome: SMN1 copy number (SMA), HBA1/HBA2 deletions (alpha thalassemia), FMR1 "
               "repeats (fragile X), CYP21A2 and other pseudogene-rich genes, deep-intronic variants. Couples may need "
               "dedicated tests for these.</li>"
               "<li>Risk figures assume the variants are pathogenic and fully penetrant; variants of uncertain "
               "significance are not listed. Genetic counselling is needed before acting on the results.</li></ul>")
    out.append("</body></html>")
    with open(a.html, "w", encoding="utf-8") as f:
        f.write("\n".join(out))
    with open(a.tsv, "w", encoding="utf-8", newline="\n") as f:
        f.write("category\tgene\tpartner\talleles\trisk\tcarrier_freq\tresidual\ttreatable\tconditions\n")
        for g, risk in couple:
            f.write("\t".join(["both_carriers", g, "A+B", " | ".join(x["name"] for p in "AB" for x in P[p][g]), risk, "",
                               "", treatable(genes[g]), genes[g]["conditions"]]) + "\n")
        for g, p in xl:
            f.write("\t".join(["x_linked", g, p, " | ".join(x["name"] for x in P[p][g]),
                               "1 in 2 sons" if sex[p] == "female" else "daughters carriers", "", "",
                               treatable(genes[g]), genes[g]["conditions"]]) + "\n")
        for g, p, cfo, cf, grp, resid in single:
            f.write("\t".join(["one_carrier", g, p, " | ".join(x["name"] for x in P[p][g]), "", f"{cf or 0:.3g}",
                               f"{resid:.3g}", treatable(genes[g]), genes[g]["conditions"]]) + "\n")
    print(f"{len(couple)} shared recessive genes, {len(xl)} X-linked, {len(single)} single-carrier genes", file=sys.stderr)


if __name__ == "__main__":
    main()
