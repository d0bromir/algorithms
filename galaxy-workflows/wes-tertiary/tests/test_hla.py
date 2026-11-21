#!/usr/bin/env python3
"""Tests of the HLA typing tool (tools/wes_hla) without T1K: allele names, G-group grouping, associations, outputs."""
import importlib.util
import json
import os
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("wes_hla", os.path.join(HERE, "..", "tools", "wes_hla", "wes_hla.py"))
hla = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hla)


def call(n, a1, q1, a2=".", q2=-1):
    return {"n": n, "a1": a1, "ab1": 50.0, "q1": q1, "a2": a2, "ab2": 40.0, "q2": q2}


def test_names_and_order():
    assert hla.two_field("HLA-A*01:01:01:01") == "A*01:01"
    assert hla.two_field("HLA-DQB1*02:163N") == "DQB1*02:163N"
    assert hla.two_field("HLA-C*04:09:01N") == "C*04:09N", "expression suffix kept on the two-field name"
    assert sorted(["A*11:01", "A*02:01", "A*101:01", "A*11:01N"], key=hla.allele_key) == \
        ["A*02:01", "A*11:01", "A*11:01N", "A*101:01"]


def test_ars_sequences_and_grouping():
    d = tempfile.mkdtemp()
    fa = os.path.join(d, "ref.fa")
    # exon 2 = positions 4-7 (0-based, inclusive end as in t1k-build.pl headers): 02:01 and 02:198 share it
    with open(fa, "w") as f:
        f.write(">HLA-DQB1*02:01:01:01 3 0 3 4 7 8 11\nAAAACCCCGGGG\n"
                ">HLA-DQB1*02:198 3 0 3 4 7 8 11\nTTTTCCCCGGGG\n"
                ">HLA-DQB1*02:163N 3 0 3 4 7 8 11\nAAAACCCCTTTT\n"
                ">HLA-DQB1*05:01:01:01 3 0 3 4 7 8 11\nAAAAGGGGGGGG\n")
    ars = hla.ars_sequences(fa, ["HLA-DQB1*02:01:01", "HLA-DQB1*02:198", "HLA-DQB1*02:163N", "HLA-DQB1*05:01:01"])
    assert ars["HLA-DQB1*02:01:01"] == ars["HLA-DQB1*02:198"] == "CCCC" != ars["HLA-DQB1*05:01:01"]
    call1, cands, groups, eq = hla.resolve("HLA-DQB1*02:163N,HLA-DQB1*02:198,HLA-DQB1*02:01:01", ars)
    assert (call1, groups) == ("DQB1*02:01", 1), "one G group: the lowest expressed allele is called"
    assert eq == ["DQB1*02:163N", "DQB1*02:198"]
    calls = {"HLA-DQB1": call(2, "HLA-DQB1*02:163N,HLA-DQB1*02:198,HLA-DQB1*02:01:01", 40, "HLA-DQB1*05:01:01", 45)}
    r = {x["gene"]: x for x in hla.genotypes(calls, 20, ars)}["HLA-DQB1"]
    assert (r["allele1"], r["allele2"], r["confidence"]) == ("DQB1*02:01", "DQB1*05:01", "high")
    assert "G group" in r["note"]
    calls["HLA-DQB1"]["a2"] = "HLA-DQB1*05:01:01,HLA-DQB1*02:01:01"  # candidates differ in exon 2
    r = {x["gene"]: x for x in hla.genotypes(calls, 20, ars)}["HLA-DQB1"]
    assert r["confidence"] == "low" and "differ in the antigen-recognition exons" in r["note"]


def test_quality_homozygous_and_untyped():
    rows = {x["gene"]: x for x in hla.genotypes({"HLA-A": call(1, "HLA-A*02:01:01", 60),
                                                 "HLA-B": call(2, "HLA-B*08:01:01", 12, "HLA-B*07:02:01", 50)}, 20)}
    assert (rows["HLA-A"]["allele1"], rows["HLA-A"]["allele2"], rows["HLA-A"]["confidence"]) == \
        ("A*02:01", "A*02:01", "high")
    assert "homozygous" in rows["HLA-A"]["note"]
    assert rows["HLA-B"]["allele1"] == "B*07:02" and rows["HLA-B"]["quality1"] == 50, "alleles sorted with quality"
    assert rows["HLA-B"]["confidence"] == "low"
    assert rows["HLA-C"]["confidence"] == "none" and rows["HLA-DRB4"]["note"].startswith("not typed (no reads, or")


def test_associations_and_pharmcat():
    calls = {"HLA-A": call(2, "HLA-A*31:01:02", 60, "HLA-A*02:01:01", 60),
             "HLA-B": call(2, "HLA-B*57:01:01", 60, "HLA-B*27:05:02", 60),
             "HLA-DQA1": call(2, "HLA-DQA1*03:01:01", 50, "HLA-DQA1*01:01:01", 50),
             "HLA-DQB1": call(2, "HLA-DQB1*03:02:01", 50, "HLA-DQB1*05:01:01", 50)}
    rows = hla.genotypes(calls, 20)
    a = {x["allele"]: x for x in hla.associations(rows)}
    assert a["B*57:01"]["status"] == "present" and "abacavir" in a["B*57:01"]["condition"]
    assert a["A*31:01"]["status"] == "present" and a["B*15:02"]["status"] == "absent"
    assert a["B*27"]["status"] == "present" and a["B*27"]["copies"] == 1
    assert "B*51" not in a, "an absent disease allele of a confidently typed gene is not listed"
    assert a["DQ8"]["status"] == "present" and "celiac" in a["DQ8"]["condition"]
    assert hla.pharmcat_lines(rows) == ["HLA-A\t*02:01/*31:01\n", "HLA-B\t*27:05/*57:01\n"]
    rows2 = hla.genotypes({"HLA-DQA1": call(2, "HLA-DQA1*01:01:01", 50, "HLA-DQA1*01:02:01", 50),
                           "HLA-DQB1": call(2, "HLA-DQB1*05:01:01", 50, "HLA-DQB1*06:02:01", 50)}, 20)
    a2 = {x["condition"]: x for x in hla.associations(rows2)}
    assert a2["celiac disease"]["allele"] == "no DQ2/DQ8" and a2["celiac disease"]["status"] == "absent"
    assert a2["narcolepsy type 1"]["status"] == "present"
    assert a2["abacavir: hypersensitivity reaction"]["status"] == "not determined", "HLA-B not typed"
    assert hla.pharmcat_lines(rows2) == []


def test_without_database():
    d = tempfile.mkdtemp()
    out = {k: os.path.join(d, k) for k in ("t.tsv", "s.json", "pc.tsv")}
    hla.main(["--bam", os.path.join(d, "absent.bam"), "--sample", "S", "--tsv", out["t.tsv"], "--json", out["s.json"],
              "--pharmcat", out["pc.tsv"]])
    j = json.load(open(out["s.json"]))
    assert j["database"] == "not given" and all(r["confidence"] == "none" for r in j["genotypes"])
    assert "no HLA allele database" in j["genotypes"][0]["note"]
    assert open(out["pc.tsv"]).read() == ""


def test_dr_haplotypes_and_celiac_from_dqb1():
    calls = {"HLA-DRB1": call(2, "HLA-DRB1*04:05:01", 50, "HLA-DRB1*09:01:02", 50),
             "HLA-DRB3": call(1, "HLA-DRB3*02:02:01", 30),
             "HLA-DRB4": call(1, "HLA-DRB4*01:03:01", 40),
             "HLA-DQA1": call(2, "HLA-DQA1*03:02", 12, "HLA-DQA1*03:03:01", 30),
             "HLA-DQB1": call(2, "HLA-DQB1*03:03:02", 50, "HLA-DQB1*04:01:01", 45)}
    rows = {r["gene"]: r for r in hla.genotypes(calls, 20)}
    assert rows["HLA-DRB4"]["confidence"] == "high", "DR4 and DR9 carry DRB4"
    assert rows["HLA-DRB3"]["confidence"] == "low" and "not expected with DRB1*04:05/DRB1*09:01" in rows["HLA-DRB3"]["note"]
    assert rows["HLA-DQA1"]["confidence"] == "low"
    a = {x["condition"]: x for x in hla.associations(list(rows.values()))}
    assert a["celiac disease"]["status"] == "absent", "no DQB1*02/03:02: DQ2/DQ8 excluded although DQA1 is uncertain"
    calls["HLA-DQB1"]["a2"] = "HLA-DQB1*03:02:01"
    a = {x["condition"]: x for x in hla.associations(hla.genotypes(calls, 20))}
    assert a["celiac disease"]["status"] == "not determined", "DQB1*03:02 with an uncertain DQA1"
