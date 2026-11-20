#!/usr/bin/env python3
"""Curated gene-disease table for the clinical report: one row per gene, disease and source.

Merges three expert curations that the report shows side by side (they often disagree, and the report says so):

  ClinGen gene-disease validity   https://search.clinicalgenome.org/kb/gene-validity/download (CSV)
  Gene2Phenotype (G2P, all panels) https://ftp.ebi.ac.uk/pub/databases/gene2phenotype/G2P_data_downloads/ (allG2P_*.csv.gz)
  Orphanet gene associations       https://www.orphadata.com/data/xml/en_product6.xml
  Orphanet inheritance and onset   https://www.orphadata.com/data/xml/en_product9_ages.xml

Output columns (tab-separated, header line first, '# source' lines record the input files):

    gene  source  disease  disease_id  moi  validity  mechanism  role  onset  hpo  ref

moi        AD, AR, XLD, XLR, XL, SD (semidominant), MT, YL, MF (multigenic/multifactorial), '' (not given); a
           comma-separated list when Orphanet gives several for the disease
validity   definitive, strong, moderate, limited, disputed, refuted, no known (ClinGen and G2P); assessed or
           not yet assessed (Orphanet's own review status)
mechanism  LoF, GoF, DN (dominant negative), undetermined, '' (G2P's molecular mechanism; Orphanet's
           'loss/gain of function' association types)
role       Orphanet's association type (disease-causing germline, major susceptibility factor, modifier, ...)
hpo        G2P's phenotype terms (HP:... separated by ';')
ref        G2P id, ClinGen report URL or ORPHA code
xref       OMIM id of a G2P disease (links it to the HPO disease annotations)

Somatic, fusion-gene and biomarker associations are left out (the report is germline).

Usage: build_gene_disease.py CLINGEN.csv ALLG2P.csv.gz PRODUCT6.xml PRODUCT9_AGES.xml OUT.tsv
Standard library only; takes a few seconds; the output is about 3 MB.
"""
import csv
import gzip
import hashlib
import os
import sys
import xml.etree.ElementTree as ET

COLS = ["gene", "source", "disease", "disease_id", "moi", "validity", "mechanism", "role", "onset", "hpo", "ref", "xref"]

MOI_TEXT = {
    "autosomal dominant": "AD", "autosomal recessive": "AR", "x-linked dominant": "XLD", "x-linked recessive": "XLR",
    "x-linked": "XL", "semidominant": "SD", "semi-dominant": "SD", "mitochondrial": "MT",
    "mitochondrial inheritance": "MT", "y-linked": "YL", "multigenic/multifactorial": "MF",
    "monoallelic_autosomal": "AD", "biallelic_autosomal": "AR", "monoallelic_x_hemizygous": "XLR",
    "monoallelic_x_heterozygous": "XLD", "monoallelic_x": "XL", "biallelic_par": "AR", "monoallelic_par": "AD",
    "monoallelic_y_hemizygous": "YL", "mitochondrial inheritance (mt)": "MT",
}
ORPHA_SKIP = ("somatic", "fusion gene", "biomarker")
ORPHA_MECH = {"loss of function": "LoF", "gain of function": "GoF"}
G2P_MECH = {"loss of function": "LoF", "gain of function": "GoF", "dominant negative": "DN",
            "undetermined": "undetermined", "undetermined non-loss-of-function": "non-LoF"}


def moi_code(text):
    t = text.strip().lower()
    if t in MOI_TEXT:
        return MOI_TEXT[t]
    for k, v in MOI_TEXT.items():  # ClinGen writes e.g. "Autosomal dominant inheritance" / "AD"
        if t.startswith(k):
            return v
    return {"ad": "AD", "ar": "AR", "xl": "XL", "xld": "XLD", "xlr": "XLR", "sd": "SD", "mt": "MT"}.get(t, "")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def clingen(path):
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        lines = [ln for ln in f if not ln.startswith('"+++')]
    r = csv.reader(lines)
    head = None
    for t in r:
        if head is None:
            if t and t[0] == "GENE SYMBOL":
                head = {h: i for i, h in enumerate(t)}
            continue
        if len(t) < len(head):
            continue
        rows.append({"gene": t[head["GENE SYMBOL"]], "source": "ClinGen", "disease": t[head["DISEASE LABEL"]],
                     "disease_id": t[head["DISEASE ID (MONDO)"]], "moi": moi_code(t[head["MOI"]]),
                     "validity": t[head["CLASSIFICATION"]].strip().lower().replace(" evidence", ""),
                     "mechanism": "", "role": "", "onset": "", "hpo": "", "ref": t[head["ONLINE REPORT"]], "xref": ""})
    return rows


def g2p(path):
    rows = []
    with gzip.open(path, "rt", encoding="utf-8", newline="") as f:
        for t in csv.DictReader(f):
            mondo = t.get("disease MONDO", "")
            omim = t.get("disease mim", "")
            rows.append({"gene": t["gene symbol"], "source": "G2P", "disease": t["disease name"],
                         "disease_id": mondo or (f"OMIM:{omim}" if omim else ""),
                         "moi": moi_code(t["allelic requirement"]), "validity": t["confidence"].strip().lower(),
                         "mechanism": G2P_MECH.get(t["molecular mechanism"].strip().lower(), ""),
                         "role": t.get("variant consequence", ""), "onset": "",
                         "hpo": ";".join(x.strip() for x in t.get("phenotypes", "").split(";") if x.strip()),
                         "ref": t["g2p id"], "xref": f"OMIM:{omim}" if omim else ""})
    return rows


def orphanet(p6, p9):
    inh, onset = {}, {}
    for d in ET.parse(p9).getroot().iter("Disorder"):
        code = d.findtext("OrphaCode")
        inh[code] = [moi_code(x.findtext("Name")) for x in d.iter("TypeOfInheritance")]
        onset[code] = [x.findtext("Name") for x in d.iter("AverageAgeOfOnset")]
    rows = []
    for d in ET.parse(p6).getroot().iter("Disorder"):
        code, name = d.findtext("OrphaCode"), d.findtext("Name")
        for a in d.iter("DisorderGeneAssociation"):
            sym = a.findtext("Gene/Symbol")
            typ = (a.findtext("DisorderGeneAssociationType/Name") or "").strip()
            status = (a.findtext("DisorderGeneAssociationStatus/Name") or "").strip().lower()
            if not sym or any(s in typ.lower() for s in ORPHA_SKIP):
                continue
            mech = next((v for k, v in ORPHA_MECH.items() if k in typ.lower()), "")
            mois = [m for m in inh.get(code, []) if m]
            rows.append({"gene": sym, "source": "Orphanet", "disease": name, "disease_id": f"ORPHA:{code}",
                         "moi": ",".join(dict.fromkeys(mois)), "validity": status, "mechanism": mech, "role": typ,
                         "onset": ",".join(onset.get(code, [])), "hpo": "", "ref": f"ORPHA:{code}", "xref": ""})
    return rows


def main(cg, g2, p6, p9, out):
    rows = clingen(cg) + g2p(g2) + orphanet(p6, p9)
    rows.sort(key=lambda r: (r["gene"], r["source"], r["disease"]))
    with open(out, "w", encoding="utf-8", newline="\n") as f:
        for p in (cg, g2, p6, p9):
            f.write(f"# source\t{os.path.basename(p)}\tsha256:{sha256(p)}\n")
        f.write("\t".join(COLS) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]).replace("\t", " ").replace("\n", " ") for c in COLS) + "\n")
    n = {s: sum(r["source"] == s for r in rows) for s in ("ClinGen", "G2P", "Orphanet")}
    print(f"{len(rows)} rows, {len({r['gene'] for r in rows})} genes: {n}", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) != 6:
        raise SystemExit(__doc__)
    main(*sys.argv[1:])
