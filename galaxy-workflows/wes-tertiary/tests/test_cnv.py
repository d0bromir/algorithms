#!/usr/bin/env python3
"""Tests of the exome CNV caller (tools/wes_cnv) on simulated per-target depths.

    python -m pytest tests/test_cnv.py

Needs numpy. Depth per target = Poisson(reads) with a capture efficiency per target shared by all samples (the
dominant source of depth variation in exomes) and 150 bp targets.
"""
import importlib.util
import json
import os
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("wes_cnv", os.path.join(HERE, "..", "tools", "wes_cnv", "wes_cnv.py"))
cnv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cnv)

CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y"]
PER_CHROM = 400


def targets():
    out = []
    for c in CHROMS:
        n = PER_CHROM if c != "Y" else 60
        for i in range(n):
            s = 3_000_000 + i * 5000  # outside the PARs
            out.append((c, s, s + 150))
    return out


T = targets()
EFF = np.random.default_rng(1).lognormal(0, 0.35, len(T))


def write_sample(path, sex, events=(), seed=0, depth=80.0):
    """events: [(chrom, first_target_index_in_chrom, n_targets, copy_factor)] relative to two copies."""
    rng = np.random.default_rng(seed)
    with open(path, "w") as f:
        idx_in_chrom = {}
        for t, (c, s, e) in enumerate(T):
            i = idx_in_chrom[c] = idx_in_chrom.get(c, -1) + 1
            copies = 2.0
            if c == "X" and sex == "male":
                copies = 1.0
            if c == "Y":
                copies = 1.0 if sex == "male" else 0.0
            for ec, first, n, factor in events:
                if ec == c and first <= i < first + n:
                    copies *= factor
            lam = depth * EFF[t] * copies / 2 * 150 / 150
            f.write(f"{c}\t{s}\t{e}\t{rng.poisson(lam):.2f}\n")


def run(sample_events, sex="female", n_panel=6, extra=()):
    d = tempfile.mkdtemp()
    panel = []
    for k in range(n_panel):
        p = os.path.join(d, f"P{k}.bed")
        write_sample(p, "male" if k % 2 else "female", seed=100 + k)
        panel.append(p)
    s = os.path.join(d, "S.bed")
    write_sample(s, sex, sample_events, seed=7)
    calls, summ = os.path.join(d, "calls.tsv"), os.path.join(d, "summary.json")
    cnv.main(["--sample", s, "--sample-name", "S", "--panel", *panel, "--calls", calls, "--summary", summ,
              "--min-depth", "10", *extra])
    rows = []
    for line in open(calls):
        if line.startswith("#"):
            continue
        t = line.rstrip("\n").split("\t")
        rows.append({"chrom": t[1], "start": int(t[2]), "type": t[5], "cn": int(t[6]), "targets": int(t[8]),
                     "confidence": t[13]})
    return rows, json.load(open(summ))


def test_deletions_and_duplication():
    ev = [("2", 100, 8, 0.5), ("5", 200, 10, 1.5), ("7", 50, 3, 0.0)]
    rows, summ = run(ev)
    found = {(r["chrom"], r["type"], r["cn"]) for r in rows}
    assert ("2", "DEL", 1) in found, rows
    assert ("5", "DUP", 3) in found, rows
    assert ("7", "DEL", 0) in found, rows
    other = [r for r in rows if r["chrom"] not in ("2", "5", "7")]
    assert len(other) <= 1, f"false calls on a clean simulation: {other}"
    assert summ["sex_complement"] == "XX" and summ["aneuploidy"] == []
    het = next(r for r in rows if r["chrom"] == "2")
    assert 6 <= het["targets"] <= 9 and het["confidence"] == "high"


def test_male_x_and_aneuploidy():
    rows, summ = run([("21", 0, PER_CHROM, 1.5), ("X", 30, 6, 0.0)], sex="male")
    assert summ["sex_complement"] == "XY", summ
    assert any(a.startswith("chromosome 21: gain of a whole copy") for a in summ["aneuploidy"]), summ
    x = [r for r in rows if r["chrom"] == "X"]
    assert x and x[0]["type"] == "DEL" and x[0]["cn"] == 0, "hemizygous loss on a male X"
    rows, summ = run([("X", 0, PER_CHROM, 2.0)], sex="male")
    assert summ["sex_complement"] == "XXY" and any("XXY" in a for a in summ["aneuploidy"])


def test_small_panel_not_analysed():
    rows, summ = run([], n_panel=2)
    assert rows == [] and summ["status"].startswith("not analysed")


def test_fixed_width_gc_bins():
    gc = np.concatenate([np.full(5000, 0.45), np.full(200, 0.82)])
    ratio = np.concatenate([np.ones(5000), np.full(200, 0.7)])
    corr = cnv.fixed_width_correction(ratio, gc)
    assert abs(corr[-1] - 0.7) < 1e-6 and abs(corr[0] - 1.0) < 1e-6, "GC extremes keep their own bin"
