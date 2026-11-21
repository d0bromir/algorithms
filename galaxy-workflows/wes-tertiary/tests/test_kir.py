#!/usr/bin/env python3
"""Tests of the KIR tool (tools/wes_kir) without T1K: presence rule, A/B genotype, B content, ligands, outputs."""
import importlib.util
import json
import os
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("wes_kir", os.path.join(HERE, "..", "tools", "wes_kir", "wes_kir.py"))
kir = importlib.util.module_from_spec(spec)
spec.loader.exec_module(kir)


def t1k(**genes):
    """gene=(allele, abundance, quality) -> T1K calls; framework genes at abundance 100, quality 60."""
    calls = {g: {"n": 1, "a1": g + "*001", "ab1": 100.0, "q1": 60, "a2": ".", "ab2": 0.0, "q2": -1}
             for g in kir.FRAMEWORK}
    for g, (al, ab, q) in genes.items():
        calls[g] = {"n": 1, "a1": al, "ab1": ab, "q1": q, "a2": ".", "ab2": 0.0, "q2": -1}
    return calls


def present(rows):
    return {r["gene"] for r in rows if r["present"]}


def test_presence_rule():
    rows, fw = kir.genotypes(t1k(KIR2DL1=("KIR2DL1*003", 90.0, 40), KIR2DS3=("KIR2DS3*001", 2.0, 0),
                                 KIR2DL5A=("KIR2DL5A*001", 30.0, 15), KIR2DS2=("KIR2DS2*001", 10.0, 35)), 20, 0.15)
    assert fw == 100.0
    p = present(rows)
    assert "KIR2DL1" in p and "KIR2DS3" not in p
    r = {x["gene"]: x for x in rows}
    assert "quality 15 < 20" in r["KIR2DL5A"]["note"], "low quality: absent"
    assert "abundance 0.10" in r["KIR2DS2"]["note"], "quality 35 but 10 % of the framework: absent"
    assert r["KIR2DL1"]["alleles"] == "KIR2DL1*003" and r["KIR3DL1"]["note"] == "no reads assigned"
    rows, _ = kir.genotypes(t1k(KIR2DL5A=("KIR2DL5A*001", 30.0, 15)), 10, 0.15)
    assert "KIR2DL5A" in present(rows), "the cut-offs are parameters"


def hap(*genes):
    return kir.haplotypes([{"gene": g, "present": g in genes} for g in kir.GENES])


def test_haplotypes():
    a = hap("KIR3DL3", "KIR2DL3", "KIR2DP1", "KIR2DL1", "KIR3DP1", "KIR2DL4", "KIR3DL1", "KIR2DS4", "KIR3DL2")
    assert (a["genotype"], a["centromeric"], a["telomeric"], a["b_content_group"], a["b_content_score"]) == \
        ("AA", "A/A", "A/A", "Neutral", 0)
    assert a["framework_complete"]
    b = hap("KIR3DL3", "KIR2DS2", "KIR2DL2", "KIR2DP1", "KIR2DL1", "KIR3DP1", "KIR2DL4", "KIR3DL1", "KIR2DS4",
            "KIR3DL2")
    assert (b["genotype"], b["centromeric"], b["b_content_group"]) == ("Bx", "B/B", "Best"), "Cen-B/B: Best"
    c = hap("KIR3DL3", "KIR2DL3", "KIR2DS2", "KIR3DP1", "KIR2DL4", "KIR3DL1", "KIR3DS1", "KIR2DS1", "KIR3DL2")
    assert (c["centromeric"], c["telomeric"], c["b_content_score"], c["b_content_group"]) == \
        ("A/B", "A/B", 2, "Better")
    d = hap("KIR3DL3", "KIR3DP1", "KIR2DL4", "KIR3DL2")
    assert d["centromeric"] == "undetermined" and d["b_content_group"] == "undetermined"


def test_ligands():
    rows = [{"gene": g, "present": g in ("KIR2DL1", "KIR2DL3", "KIR3DL1")} for g in kir.GENES]
    hla = {"genotypes": [
        {"gene": "HLA-C", "allele1": "C*06:02", "allele2": "C*07:01", "kir_ligand1": "C2", "kir_ligand2": "C1",
         "confidence": "high"},
        {"gene": "HLA-B", "allele1": "B*08:01", "allele2": "B*57:01", "kir_ligand1": "Bw6", "kir_ligand2": "Bw4-80I",
         "confidence": "high"},
        {"gene": "HLA-A", "allele1": "A*24:02", "allele2": "A*24:02", "kir_ligand1": "Bw4-80I",
         "kir_ligand2": "Bw4-80I", "confidence": "low"}]}
    p = {(x["kir"], x["ligand"]): x for x in kir.ligands(rows, hla)}
    assert p[("KIR2DL1", "C2")]["status"] == "pair present" and p[("KIR2DL1", "C2")]["ligand_alleles"] == ["C*06:02"]
    assert p[("KIR3DL1", "Bw4")]["ligand_alleles"] == ["B*57:01"], "low-confidence HLA-A not used"
    assert p[("KIR2DL2", "C1")]["status"] == "KIR absent"
    assert p[("KIR3DS1", "Bw4-80I")]["status"] == "KIR absent"
    assert all(x["status"] == "ligand not typed" for x in kir.ligands(rows, None))


def test_without_database():
    d = tempfile.mkdtemp()
    out = {k: os.path.join(d, k) for k in ("k.tsv", "k.json")}
    kir.main(["--bam", os.path.join(d, "absent.bam"), "--sample", "S", "--tsv", out["k.tsv"], "--json", out["k.json"]])
    j = json.load(open(out["k.json"]))
    assert j["database"] == "not given" and j["haplotypes"]["genotype"] == "not typed"
    assert all(not r["present"] and "kir.dat" in r["note"] for r in j["genes"])


def test_db_release():
    d = tempfile.mkdtemp()
    f = os.path.join(d, "kir.dat")
    with open(f, "w") as o:
        o.write("ID   KIR00001; SV 1; standard; DNA; HUM; 14800 BP.\nCC   IPD-KIR Release Version 2.15.0\nSQ\n")
    assert kir.db_release(f) == "2.15.0"
