#!/usr/bin/env python3
"""Tests of the clinical report tool (tools/wes_clinical_report) on small synthetic inputs.

    python -m pytest tests/test_clinical_report.py      (or run this file directly)

Only the standard library is needed; the BAM/reference features (samtools) are not exercised here.
"""
import csv
import importlib.util
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
TOOL = os.path.join(HERE, "..", "tools", "wes_clinical_report")
DATA = os.path.join(TOOL, "test-data")

spec = importlib.util.spec_from_file_location("wes_clinical_report", os.path.join(TOOL, "wes_clinical_report.py"))
rep = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rep)


def run_report(extra=()):
    out = tempfile.mkdtemp()
    html, tsv = os.path.join(out, "r.html"), os.path.join(out, "r.tsv")
    rep.main(["--vcf", os.path.join(DATA, "annotated.vcf"), "--gencc", os.path.join(DATA, "gencc.tsv"),
              "--constraint", os.path.join(DATA, "constraint.tsv"), "--html", html, "--tsv", tsv, *extra])
    with open(tsv) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    with open(html) as f:
        return rows, f.read()


def categories(rows, gene):
    return {r["category"] for r in rows if r["gene"] == gene}


def test_acmg_sf_list():
    assert len(rep.ACMG_SF) == 84, "ACMG SF v3.3"
    assert rep.ACMG_SF["MUTYH"] == "biallelic" and rep.ACMG_SF["HFE"] == "hfe" and rep.ACMG_SF["TTN"] == "ttn"
    assert rep.ACMG_SF["PLN"] == "het" and rep.ACMG_SF["CYP27A1"] == "biallelic" and rep.ACMG_SF["ABCD1"] == "xl"
    assert set(rep.ACMG_SF) == set(rep.ACMG_SF_CATEGORY)


def test_sf_v33_genes():
    rows, html = run_report()
    assert categories(rows, "PLN") == {"diagnostic", "acmg_sf_reportable"}, "PLN P/LP het: reportable (dominant)"
    assert "acmg_sf_not" in categories(rows, "ABCD1"), "a single heterozygous ABCD1 variant is not an SF"
    assert "acmg_sf_reportable" not in categories(rows, "ABCD1")
    assert "ACMG SF v3.3 (84 genes)" in html


def test_gnomad_v41_frequency_rules():
    rows, html = run_report()
    assert categories(rows, "GENEF") == {"pathogenic_review"}, \
        "overall AF 0.0004 but grpmax FAF95 0.004: too common for a dominant disorder"
    f = next(r for r in rows if r["gene"] == "GENEF")
    assert f["gnomad_faf95"] == "0.004" and f["gnomad_grpmax_af"] == "0.006"
    assert "gnomAD site filtered (EXOMES_FILTERED)" in f["qc_flags"] and f["confidence"] == "high"
    assert categories(rows, "GENER") == {"candidate_inheritance_not_met"}, "6 homozygotes in gnomAD"
    assert "6 homozygotes in gnomAD" in next(r["note"] for r in rows if r["gene"] == "GENER")
    assert "FAF95 0.004 (afr)" in html


def test_freq():
    v = rep.Variant()
    v.af, v.faf = None, None
    assert rep.freq(v) is None
    v.af, v.faf = 0.02, 0.001
    assert rep.freq(v) == 0.02, "founder variants: FAF excludes bottlenecked groups"
    v.af, v.faf = 0.0004, 0.004
    assert rep.freq(v) == 0.004


def test_findings_by_inheritance():
    rows, html = run_report()
    assert categories(rows, "BRCA2") == {"diagnostic", "acmg_sf_reportable"}, "dominant P/LP: diagnostic + SF"
    assert categories(rows, "PAH") == {"carrier"}, "single P/LP allele in a recessive gene: carrier"
    assert "carrier" in categories(rows, "MUTYH") and "acmg_sf_reportable" not in categories(rows, "MUTYH"), \
        "MUTYH is reportable as a secondary finding only when biallelic"
    assert categories(rows, "GENEX") == {"candidate"}, "rare stop-gain in a LoF-intolerant dominant gene"
    assert not categories(rows, "TP53"), "benign common variant must not be reported"
    assert "Research use / decision support only" in html


def test_common_pathogenic_is_reviewed():
    rows, _ = run_report()
    assert categories(rows, "GENED") == {"pathogenic_review"}, "ClinVar P but gnomAD AF 0.2 must not be diagnostic"


def test_mixed_inheritance_homozygous_lof():
    rows, _ = run_report()
    assert categories(rows, "GENEM") == {"candidate"}, \
        "homozygous LoF in a gene with a recessive disorder fits, even if the gene is LoF-tolerant for dominant disease"


def test_artefact_list():
    rows, html = run_report()
    g = next(r for r in rows if r["gene"] == "GENEX")
    assert g["confidence"] == "high" and "recurrent artefact" not in g["qc_flags"]
    assert "No recurrent-artefact list was provided" in html
    rows, html = run_report(["--artefacts", os.path.join(DATA, "artefacts.vcf")])
    g = next(r for r in rows if r["gene"] == "GENEX")
    assert "recurrent artefact: called in 12/40" in g["qc_flags"], "matched on CHROM (chr prefix ignored)/POS/REF/ALT"
    assert g["confidence"].startswith("low") and "recurrent assay artefact" in g["confidence"]
    assert categories(rows, "GENEX") == {"candidate"}, "flagged, not removed"
    assert "2 alleles called in many" in html and "1 calls of this sample match it" in html


def test_candidate_evidence_codes():
    rows, _ = run_report()
    note = next(r["note"] for r in rows if r["gene"] == "GENEX")
    assert "PVS1" in note and "PM2_supporting" in note


def test_pharmacogenomics_vcf_only():
    rows, html = run_report()
    pgx = {r["gene"]: (r["genotype"], r["note"]) for r in rows if r["category"] == "pgx"}
    assert pgx["CYP2C19"] == ("*2/*2", "Poor metabolizer")
    assert pgx["CYP2C9"][1].startswith("Normal metabolizer")
    assert all(len(r) == 23 and None not in r.values() for r in rows)
    assert "absence of a variant is assumed to be the reference allele" in html


def test_hpo_ranking(tmp_path=None):
    d = tmp_path or tempfile.mkdtemp()
    p = os.path.join(str(d), "p2g.txt")
    with open(p, "w") as f:
        f.write("hpo_id\thpo_name\tncbi_gene_id\tgene_symbol\tdisease_id\nHP:0000001\tExample\t1\tGENEX\tOMIM:1\n")
    rows, html = run_report(["--hpo-genes", p, "--hpo-terms", "HP:0000001"])
    assert "phenotype match: Example" in html
    note = next(r["note"] for r in rows if r["gene"] == "GENEX")
    assert "PP4" in note


def test_low_confidence():
    def var(ref, alt, zyg="het", ab=0.5, qual=200.0, flags=()):
        v = rep.Variant()
        v.ref, v.alt, v.zyg, v.ab, v.qual, v.flags = ref, alt, zyg, ab, qual, list(flags)
        return v
    assert rep.low_confidence(var("C", "T")) == []
    assert rep.low_confidence(var("C", "T", ab=0.30, qual=36, flags=["low depth (10x)"])) == [], "TPP1-like true call"
    assert "alt fraction < 0.25" in rep.low_confidence(var("C", "T", ab=0.17))
    assert "QUAL < 30" in rep.low_confidence(var("C", "T", qual=20))
    assert "indel in a repeat" in rep.low_confidence(var("CT", "C", ab=0.28, flags=["homopolymer (T x10)"]))
    assert rep.low_confidence(var("CT", "C", flags=["homopolymer (T x6)"])) == []
    assert rep.low_confidence(var("C", "T", flags=["segmental duplication / pseudogene: short-read calls unreliable"]))


def test_repeat_flag():
    assert rep.repeat_flag("ACGTTTTTTTTGCA", 3, 2) == "homopolymer (T x8)"
    assert rep.repeat_flag("GGACACACACACTT", 3, 1).startswith("tandem repeat")
    assert rep.repeat_flag("ACGTACGGTCAT", 5, 1) is None


def test_pgx_rules():
    # DPYD: one no-function and one decreased-function allele (assumed in trans) -> activity score 0.5
    dpyd = [("*2A (c.1905+1G>A)", "", "", 1, None, ""), ("*13 (c.1679T>G)", "", "", 0, None, ""),
            ("c.2846A>T", "", "", 1, None, ""), ("HapB3 (c.1236G>A)", "", "", 0, None, "")]
    assert rep.pgx_pheno("DPYD", dpyd) == ("activity score 0.5", "Poor metabolizer (two variants assumed in trans)")
    cyp2c19 = [("*2", "", "", 1, None, ""), ("*3", "", "", 0, None, ""), ("*17", "", "", 1, None, "")]
    assert rep.pgx_pheno("CYP2C19", cyp2c19) == ("*2/*17", "Intermediate metabolizer")
    cyp2c9 = [("*2", "", "", 1, None, ""), ("*3", "", "", 1, None, "")]
    assert rep.pgx_pheno("CYP2C9", cyp2c9) == ("*2/*3", "Poor metabolizer (activity score 0.5)")
    tpmt = [("*2 (c.238G>C)", "", "", 0, None, ""), ("*3B (c.460G>A)", "", "", 1, None, ""), ("*3C (c.719A>G)", "", "", 1, None, "")]
    assert rep.pgx_pheno("TPMT", tpmt) == ("*3A/*1", "Intermediate metabolizer")
    unknown = [("*2", "", "", None, None, ""), ("*3", "", "", 0, None, ""), ("*17", "", "", 0, None, "")]
    assert rep.pgx_pheno("CYP2C19", unknown)[0] == "indeterminate"


if __name__ == "__main__":
    failed = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"ok   {name}")
            except AssertionError as e:
                failed += 1
                print(f"FAIL {name}: {e}")
    sys.exit(1 if failed else 0)
