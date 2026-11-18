#!/usr/bin/env python3
"""Protein-level ClinVar table for the clinical report (ACMG PS1, PM5 and missense hotspots).

Reads ClinVar's variant_summary.txt.gz (https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/) and writes one row
per GRCh38 single-nucleotide missense variant classified pathogenic/likely pathogenic (P) or benign/likely benign (B):

    gene  residue  ref_aa  alt_aa  class  stars  variation_id  chrom  pos  ref  alt  transcript

Variants of uncertain significance and with conflicting classifications are left out. Residues are numbered on
ClinVar's preferred RefSeq transcript (usually MANE Select), so the report checks the reference amino acid before it
compares a variant with them.

Usage: build_clinvar_protein.py variant_summary.txt.gz clinvar_protein.tsv.gz
Standard library only. Takes about 1 minute; the output is about 4 MB.
"""
import gzip
import re
import sys

NAME = re.compile(r"^(N[MR]_\d+(?:\.\d+)?)\(([^)]+)\):c\.[^ ]+ \(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)$")
STARS = {"practice guideline": 4, "reviewed by expert panel": 3,
         "criteria provided, multiple submitters, no conflicts": 2, "criteria provided, single submitter": 1,
         "criteria provided, conflicting classifications": 1, "criteria provided, conflicting interpretations": 1}


def klass(sig):
    s = sig.lower()
    if "conflicting" in s:
        return "CONF"
    if s in ("pathogenic", "likely pathogenic", "pathogenic/likely pathogenic") or s.startswith(
            ("pathogenic;", "likely pathogenic;", "pathogenic/likely pathogenic;")):
        return "P"
    if s in ("benign", "likely benign", "benign/likely benign") or s.startswith(
            ("benign;", "likely benign;", "benign/likely benign;")):
        return "B"
    if s.startswith("uncertain significance"):
        return "VUS"
    return "OTHER"


def main(src, dst):
    n = 0
    seen = set()
    with gzip.open(src, "rt", encoding="utf-8", errors="replace") as f, gzip.open(dst, "wt") as out:
        ix = {h: i for i, h in enumerate(f.readline().lstrip("#").rstrip("\n").split("\t"))}
        out.write("#gene\tresidue\tref_aa\talt_aa\tclass\tstars\tvariation_id\tchrom\tpos\tref\talt\ttranscript\n")
        for line in f:
            t = line.rstrip("\n").split("\t")
            if t[ix["Assembly"]] != "GRCh38" or t[ix["Type"]] != "single nucleotide variant":
                continue
            m = NAME.match(t[ix["Name"]])
            if not m or m.group(5) in ("Ter", "Xaa") or m.group(3) == m.group(5):
                continue
            cls = klass(t[ix["ClinicalSignificance"]])
            if cls not in ("P", "B"):
                continue
            vid = t[ix["VariationID"]]
            if vid in seen:
                continue
            seen.add(vid)
            tx, gene, ref_aa, res, alt_aa = m.groups()
            ref, alt = t[ix["ReferenceAlleleVCF"]], t[ix["AlternateAlleleVCF"]]
            if len(ref) != 1 or len(alt) != 1 or ref == "na":
                continue
            out.write("\t".join((gene.split(":")[0], res, ref_aa, alt_aa, cls,
                                 str(STARS.get(t[ix["ReviewStatus"]], 0)), vid, t[ix["Chromosome"]],
                                 t[ix["PositionVCF"]], ref, alt, tx)) + "\n")
            n += 1
    print(f"{n} pathogenic or benign missense variants written to {dst}", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
