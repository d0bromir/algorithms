#!/usr/bin/env python3
"""Tests of the haplogroup code without HaploGrep3/Yleaf: wes_mito (HSD ranges and notation, parsing of the HaploGrep3
output, confidence) and wes_yhaplo (Y check, parsing of Yleaf's report). The programs are replaced by stubs."""
import importlib.util
import json
import os
import subprocess
import tempfile
import types

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, "..", "tools", *path))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


mito = load("wes_mito", ("wes_mito", "wes_mito.py"))
yh = load("wes_yhaplo", ("wes_yhaplo", "wes_yhaplo.py"))

HG_OUT = ('"SampleID"\t"Haplogroup"\t"Rank"\t"Quality"\t"Range"\t"Not_Found_Polys"\t"Found_Polys"\t"Remaining_Polys"\t'
          '"AAC_In_Remainings"\t"Input_Sample"\n'
          '"s"\t"K1a1b1a"\t"1"\t"0.9830"\t"1-16569"\t"16093C"\t"73G"\t"152C (localPrivateMutation) 16519C (hotspot)"'
          '\t""\t""\n'
          '"s"\t"K1a1b1"\t"2"\t"0.9477"\t"1-16569"\t""\t""\t""\t""\t""\n')


def test_ranges_and_notation():
    assert mito.ranges([1, 2, 3, 7, 8, 10]) == "1-3;7-8;10-10", "single positions as ranges"
    assert mito.ranges([]) == ""
    assert mito.hsd_variants(73, "A", "G") == ["73G"]
    assert mito.hsd_variants(522, "ACA", "A") == ["523d", "524d"], "deletion after the anchor base"
    assert mito.hsd_variants(315, "C", "CCT") == ["315.1C", "315.2T"]


def fake_haplogrep(monkeypatch_out):
    """Replace shutil.which / subprocess.run in wes_mito by a HaploGrep3 stub that writes HG_OUT and keeps the input."""
    seen = {}

    def run(cmd, **kw):
        seen["hsd"] = open(cmd[cmd.index("--in") + 1]).read()
        with open(cmd[cmd.index("--out") + 1], "w") as f:
            f.write(monkeypatch_out)
        return subprocess.CompletedProcess(cmd, 0, "", "")
    mito.shutil = types.SimpleNamespace(which=lambda x: "/bin/haplogrep3", rmtree=lambda *a, **k: None)
    mito.subprocess = types.SimpleNamespace(run=run, Popen=subprocess.Popen, PIPE=subprocess.PIPE)
    return seen


def restore():
    import shutil
    mito.shutil, mito.subprocess = shutil, subprocess


def test_haplogroup_parsing_and_confidence():
    try:
        seen = fake_haplogrep(HG_OUT)
        covered = list(range(1, 16570))
        r = mito.haplogroup("s", covered, [(73, "A", "G"), (3106, "CN", "C"), (310, "T", "TC")], mito.TREE)
        assert seen["hsd"].startswith("s\t1-16569\t?\t73G\t") and "3107d" not in seen["hsd"], "rCRS 3107N skipped"
        assert "310.1C" in seen["hsd"]
        assert (r["haplogroup"], r["quality"], r["confidence"]) == ("K1a1b1a", 0.983, "high")
        assert r["expected_not_found"] == ["16093C"]
        assert r["extra_variants"] == ["152C local private mutation", "16519C hotspot"]
        assert r["alternatives"] == [{"haplogroup": "K1a1b1", "quality": 0.9477}]
        low = mito.haplogroup("s", covered[:3000], [], mito.TREE)
        assert low["confidence"] == "low" and "18% of the mtDNA covered" in low["note"], \
            "quality 0.98 but < 40 % covered: low"
        assert mito.haplogroup("s", [], [], mito.TREE)["status"].startswith("not determined")
    finally:
        restore()
    mito.shutil = types.SimpleNamespace(which=lambda x: None)
    try:
        assert mito.haplogroup("s", [1], [], mito.TREE)["status"] == "not run: haplogrep3 not found"
    finally:
        restore()


def test_yhaplo(monkeypatch=None):
    d = tempfile.mkdtemp()
    out = {k: os.path.join(d, k) for k in ("y.json", "y.tsv")}
    args = ["--bam", "x.bam", "--reference", "r.fa", "--json", out["y.json"], "--tsv", out["y.tsv"], "--sample", "S"]
    orig_idx, orig_yleaf = yh.idxstats, yh.yleaf
    try:
        yh.idxstats = lambda bam: {"X": 1_300_000, "Y": 16_000, "1": 9}
        yh.yleaf = lambda *a: (_ for _ in ()).throw(AssertionError("Yleaf must not run for a female"))
        yh.main(args)
        j = json.load(open(out["y.json"]))
        assert j["status"] == "not typed" and j["y_x_ratio"] == 0.0123 and "female" in j["note"]
        yh.idxstats = lambda bam: {"X": 700_000, "Y": 185_000}
        report = {"yfull": {"haplogroup": "I1*(xI-Y33791,I-Y233039)", "haplogroup_base": "I1",
                            "excluded_subclades": ["I-Y33791", "I-Y233039"], "markers": ["CTS2970"],
                            "valid_markers": 1235, "qc_score": 1.0, "qc1": 1.0, "qc2": 1.0, "qc3": 1.0},
                  "isogg": {"haplogroup": "I1*(xI1a1a1a1a)", "markers": ["Z2760"]}}
        path = [{"haplogroup": "ROOT", "marker": None, "state": "NC"}, {"haplogroup": "F", "marker": "M235", "state": "D"},
                {"haplogroup": "I1", "marker": "CTS2970", "state": "D"}]
        yh.yleaf = lambda *a: (report, path)
        yh.main(args)
        j = json.load(open(out["y.json"]))
        assert (j["haplogroup"], j["isogg"], j["confidence"], j["excluded_subclades"]) == ("I1*", "I1*", "high", 2)
        assert j["lineage"] == "F > I1" and "paragroup" in j["note"]
        assert open(out["y.tsv"]).read() == "haplogroup\tmarker\nF\tM235\nI1\tCTS2970\n"
        report["yfull"] = {"haplogroup": "NA", "valid_markers": 13}
        yh.main(args)
        assert json.load(open(out["y.json"]))["status"] == "not determined"
        yh.idxstats = lambda bam: {"Y": 500}  # no X in the BAM (a Y-only extract): the check is off with 0
        report["yfull"] = {"haplogroup": "J-L816", "valid_markers": 900, "qc_score": 0.8}
        yh.main(args + ["--min-y-ratio", "0"])
        j = json.load(open(out["y.json"]))
        assert j["haplogroup"] == "J-L816" and j["confidence"] == "low", "QC 0.8 < 0.95"
    finally:
        yh.idxstats, yh.yleaf = orig_idx, orig_yleaf
