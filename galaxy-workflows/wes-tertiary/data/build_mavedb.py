#!/usr/bin/env python3
"""Calibrated functional-assay evidence (ACMG PS3/BS3) from MaveDB for the clinical report.

Input: the MaveDB public bulk download (CC0), e.g. https://zenodo.org/records/20840937 (mavedb-dump.*.tar.gz). Only
score sets with a score calibration that assigns ACMG criteria (PS3/BS3 at a strength) to score ranges are used; their
variants are matched to the report by MaveDB's own mapping to GRCh38 (post-mapped genomic HGVS) and to the RefSeq
protein (post-mapped p. HGVS, which the report reaches through the MANE NP_ name).

Output columns (tab-separated):
    key  gene  urn  score  code  calibration  research_only  title
key            NC_...:g.X (nucleotide-level assays) or NP_...:p.X (protein-level assays)
code           e.g. PS3_Moderate or BS3_Supporting: the calibration's class for this score
research_only  1 when MaveDB marks the calibration research-use-only (e.g. ExCALIBR gene-level calibrations): the
               report shows such evidence but does not score it
Standard library only. Usage: build_mavedb.py mavedb-dump.tar.gz mave_evidence.tsv
"""
import csv
import hashlib
import io
import json
import os
import re
import sys
import tarfile

STRENGTH = {"SUPPORTING": "Supporting", "MODERATE": "Moderate", "MODERATE_PLUS": "Moderate", "STRONG": "Strong",
            "VERY_STRONG": "VeryStrong"}


def urn_file(urn):
    return urn.replace(":", "-")


def calibrations(main):
    """urn -> (genes, title, [(calibration title, research_only, [(lo, hi, lo_incl, hi_incl, code)])])"""
    out = {}
    for es in main["experimentSets"]:
        for e in es.get("experiments", []):
            for ss in e.get("scoreSets", []):
                cals = []
                for c in ss.get("scoreCalibrations") or []:
                    ranges = []
                    for f in c.get("functionalClassifications") or []:
                        a = f.get("acmgClassification")
                        if not a or not a.get("criterion") or not a.get("evidenceStrength") or not f.get("range"):
                            continue
                        lo, hi = f["range"]
                        code = f"{a['criterion']}_{STRENGTH.get(a['evidenceStrength'], a['evidenceStrength'].title())}"
                        ranges.append((lo, hi, f.get("inclusiveLowerBound", True), f.get("inclusiveUpperBound", False),
                                       code))
                    if ranges:
                        cals.append((c.get("title") or "calibration", bool(c.get("researchUseOnly")), ranges))
                if cals:
                    genes = [t.get("mappedHgncName") or t.get("name") or "" for t in ss.get("targetGenes", [])]
                    out[ss["urn"]] = (genes, ss.get("title", ""), cals)
    return out


def classify(score, ranges):
    for lo, hi, li, hi_i, code in ranges:
        if lo is not None and (score < lo or (score == lo and not li)):
            continue
        if hi is not None and (score > hi or (score == hi and not hi_i)):
            continue
        return code
    return None


def fnum(x):
    try:
        v = float(x)
        return v if v == v else None
    except (TypeError, ValueError):
        return None


def main(dump, out):
    want, scores, annots, cals = {}, {}, {}, None
    with tarfile.open(dump, "r:gz") as tf:
        for m in tf:
            name = os.path.basename(m.name)
            if name.startswith("._") or not m.isfile():
                continue
            if name == "main.json":
                cals = calibrations(json.load(tf.extractfile(m)))
                want = {urn_file(u): u for u in cals}
            elif name.endswith(".scores.csv") or name.endswith(".annotations.csv"):
                stem = name.rsplit(".", 2)[0]
                if cals is None:  # main.json comes first in the dump; keep everything until it is read
                    pass
                elif stem not in want:
                    continue
                rows = list(csv.DictReader(io.TextIOWrapper(tf.extractfile(m), encoding="utf-8")))
                (scores if name.endswith(".scores.csv") else annots)[stem] = rows
    n = 0
    with open(out, "w", encoding="utf-8", newline="\n") as f:
        f.write(f"# source\t{os.path.basename(dump)}\tsha256:{sha256(dump)}\n")
        f.write("key\tgene\turn\tscore\tcode\tcalibration\tresearch_only\ttitle\n")
        for stem, urn in sorted(want.items()):
            genes, title, cl = cals[urn]
            sc = {r.get("accession"): fnum(r.get("scores.score", r.get("score"))) for r in scores.get(stem, [])}
            for a in annots.get(stem, []):
                s = sc.get(a.get("accession"))
                if s is None:
                    continue
                keys = [k for k in (a.get("mavedb.post_mapped_hgvs_g"), a.get("mavedb.post_mapped_hgvs_p"))
                        if k and re.match(r"N[CP]_\d+\.\d+:[gp]\.", k)]
                if not keys:
                    continue
                key = keys[0]
                for ctitle, research, ranges in cl:
                    code = classify(s, ranges)
                    if code:
                        f.write("\t".join([key, ",".join(genes), urn, f"{s:.4g}", code, ctitle.replace("\t", " "),
                                           "1" if research else "0", title.replace("\t", " ")[:120]]) + "\n")
                        n += 1
    print(f"{len(want)} calibrated score sets, {n} variant-calibration rows", file=sys.stderr)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    main(*sys.argv[1:])
