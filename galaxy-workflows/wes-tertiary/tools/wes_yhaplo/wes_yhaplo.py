#!/usr/bin/env python3
"""Y-chromosome haplogroup from the reads of an exome (or genome), with Yleaf.

Yleaf (Ralf et al., Mol Biol Evol 2018; Yleaf 4) piles up the reads at the SNPs of the YFull tree and walks the tree
to the most derived haplogroup with support; it also gives the legacy ISOGG name. Exome captures cover the Y genes and
some off-target sequence, enough for a few hundred to a few thousand informative markers in a male.

The tool first checks that there is a Y chromosome: the Y/X ratio of mapped reads (samtools idxstats) is about 0.01 in
female exomes (reads mismapped to Y) and 0.2-0.3 in male exomes; below --min-y-ratio the sample is not typed.

Outputs: a JSON summary (haplogroup, YFull and ISOGG names, defining markers, path of derived nodes, Yleaf QC scores,
confidence) and a TSV of the derived nodes on the path.

    wes_yhaplo.py --bam reads.bam --reference ref.fa --json y.json --tsv y.tsv [--sample S] [--threads 4]
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile

VERSION = "0.1.0"
Y_NAMES = ("Y", "chrY")
X_NAMES = ("X", "chrX")


def idxstats(bam):
    out = subprocess.run(["samtools", "idxstats", bam], capture_output=True, text=True, check=True).stdout
    return {t[0]: int(t[2]) for t in (line.split("\t") for line in out.splitlines()) if len(t) > 2}


def clean(name):
    """'I1*(xI-Y33791,...)' -> 'I1*' (the excluded subclades are reported separately)."""
    return name.split("(x", 1)[0] if name else name


def yleaf(bam, ref, sample, threads, args):
    """Run Yleaf with the YFull and ISOGG trees -> {tree: sample report}, path of the YFull tree."""
    exe = shutil.which("Yleaf")
    if not exe:
        raise SystemExit("Yleaf not found")
    d = tempfile.mkdtemp()
    out = os.path.join(d, "out")
    rep = os.path.join(d, "report.json")
    cmd = [exe, "-bam", bam, "-rg", "hg38", "--ref-fasta", ref, "-o", out, "-t", str(threads), "-tree", "yfull",
           "isogg", "-force", "--no-update-check", "--report-json", rep,
           "-r", str(args.min_reads), "-q", str(args.min_mapq), "-b", str(args.base_majority),
           "-pq", str(args.prediction_quality)]
    p = subprocess.run(cmd, capture_output=True, text=True, env=dict(os.environ, YLEAF_NO_UPDATE_CHECK="1"))
    if p.returncode:
        raise SystemExit("Yleaf failed: " + (p.stderr or p.stdout)[-800:])
    res = {}
    for tree in ("yfull", "isogg"):
        f = os.path.join(d, f"report.{tree}.json")
        if os.path.exists(f):
            s = json.load(open(f))["samples"]
            res[tree] = s[0] if s else {}
    path = []
    for root, _, files in os.walk(out):
        if "haplogroup_path_yfull.json" in files:
            path = json.load(open(os.path.join(root, "haplogroup_path_yfull.json")))
    shutil.rmtree(d, ignore_errors=True)
    return res, path


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--bam", required=True)
    ap.add_argument("--reference", required=True, help="the FASTA the reads were aligned to (GRCh38)")
    ap.add_argument("--sample", default="sample")
    ap.add_argument("--json", required=True)
    ap.add_argument("--tsv", required=True)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--min-y-ratio", type=float, default=0.05,
                    help="lowest Y/X ratio of mapped reads for a Y chromosome (default 0.05; 0: no check)")
    ap.add_argument("--min-reads", type=int, default=10, help="Yleaf: reads needed at a marker (default 10)")
    ap.add_argument("--min-mapq", type=int, default=20, help="Yleaf: lowest mapping quality (default 20)")
    ap.add_argument("--base-majority", type=int, default=90, help="Yleaf: %% of reads for a base call (default 90)")
    ap.add_argument("--prediction-quality", type=float, default=0.95,
                    help="Yleaf: QC score for an accepted prediction (default 0.95)")
    a = ap.parse_args(argv)

    st = idxstats(a.bam)
    x = sum(st.get(c, 0) for c in X_NAMES)
    y = sum(st.get(c, 0) for c in Y_NAMES)
    ratio = round(y / x, 4) if x else None
    summary = {"sample": a.sample, "tool": f"wes_yhaplo {VERSION}", "program": "Yleaf 4.1.7",
               "y_reads": y, "x_reads": x, "y_x_ratio": ratio, "min_y_ratio": a.min_y_ratio}
    rows = []
    if a.min_y_ratio > 0 and (ratio is None or ratio < a.min_y_ratio):  # 0 turns the check off
        summary.update(status="not typed", haplogroup="", confidence="none",
                       note=f"no Y chromosome: Y/X read ratio {ratio} < {a.min_y_ratio} (female, or loss of Y)")
    else:
        res, path = yleaf(a.bam, a.reference, a.sample, a.threads, a)
        yf, ig = res.get("yfull", {}), res.get("isogg", {})
        hg = clean(yf.get("haplogroup", ""))
        if not hg or hg == "NA":
            summary.update(status="not determined", haplogroup="", confidence="none",
                           valid_markers=yf.get("valid_markers"),
                           note="Yleaf found too few informative markers for a haplogroup")
        else:
            qc = yf.get("qc_score")
            derived = [n for n in path if n.get("state") == "D"]
            rows = [(n["haplogroup"], n.get("marker") or "") for n in derived]
            summary.update(
                status="ok", haplogroup=hg, haplogroup_base=yf.get("haplogroup_base", hg),
                markers=yf.get("markers", []), excluded_subclades=len(yf.get("excluded_subclades", [])),
                isogg=clean(ig.get("haplogroup", "")), isogg_markers=ig.get("markers", []),
                valid_markers=yf.get("valid_markers"), qc_score=qc,
                qc=[yf.get("qc1"), yf.get("qc2"), yf.get("qc3")],
                confidence="high" if qc is not None and qc >= a.prediction_quality else "low",
                lineage=" > ".join(h for h, _ in rows),
                note=("paragroup: none of the " + str(len(yf.get("excluded_subclades", []))) +
                      " known subclades is supported (often lack of coverage in an exome)")
                if yf.get("excluded_subclades") else "")
    with open(a.json, "w") as f:
        json.dump(summary, f, indent=1)
    with open(a.tsv, "w") as f:
        f.write("haplogroup\tmarker\n")
        for r in rows:
            f.write("\t".join(r) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "lineage"}), file=sys.stderr)


if __name__ == "__main__":
    main()
