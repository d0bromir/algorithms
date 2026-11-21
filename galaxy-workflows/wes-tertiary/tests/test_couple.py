#!/usr/bin/env python3
"""Tests of the couple carrier screening tool (tools/wes_couple) on synthetic findings tables."""
import csv
import importlib.util
import os
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
TOOL = os.path.join(HERE, "..", "tools", "wes_couple")
DATA = os.path.join(TOOL, "test-data")
spec = importlib.util.spec_from_file_location("wes_couple", os.path.join(TOOL, "wes_couple.py"))
couple = importlib.util.module_from_spec(spec)
spec.loader.exec_module(couple)


def run(extra=()):
    d = tempfile.mkdtemp()
    h, t = os.path.join(d, "c.html"), os.path.join(d, "c.tsv")
    couple.main(["--a", os.path.join(DATA, "partnerA.tsv"), "--b", os.path.join(DATA, "partnerB.tsv"),
                 "--carrier-genes", os.path.join(DATA, "carrier_genes.tsv"), "--html", h, "--tsv", t, *extra])
    with open(t) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    with open(h) as f:
        return rows, f.read()


def test_couple_risks():
    rows, html = run()
    by = {(r["category"], r["gene"]): r for r in rows}
    assert by[("both_carriers", "CFTR")]["risk"] == "1 in 4", "both partners carry a CFTR P allele"
    assert "newborn screening" in by[("both_carriers", "CFTR")]["treatable"]
    assert ("x_linked", "DMD") in by and by[("x_linked", "DMD")]["risk"] == "1 in 2 sons", "female carrier (automated P)"
    assert ("both_carriers", "PAH") not in by, "partner A's PAH variant is a VUS"
    assert by[("one_carrier", "PAH")]["partner"] == "B", "partner B's PAH deletion (high-confidence CNV) counts"
    assert ("both_carriers", "GJB2") not in by and by[("one_carrier", "GJB2")]["partner"] == "A", \
        "partner B's GJB2 allele is labelled risk factor in ClinVar: not a carrier allele"
    assert by[("one_carrier", "BRCA2")]["partner"] == "A", "low-confidence calls are kept for confirmation"
    assert "CONFIRM: alt fraction &lt; 0.25 (low-confidence call)" in html
    assert "ANNA (female)" in html and "BEN (male)" in html


def test_sex_override_and_residual():
    rows, _ = run(["--a-sex", "male", "--miss", "0.2"])
    by = {(r["category"], r["gene"]): r for r in rows}
    assert by[("x_linked", "DMD")]["risk"] == "daughters carriers"
    assert abs(float(by[("one_carrier", "GJB2")]["residual"]) - 0.04) < 1e-9, "0.2 carrier frequency x 0.2 missed"
