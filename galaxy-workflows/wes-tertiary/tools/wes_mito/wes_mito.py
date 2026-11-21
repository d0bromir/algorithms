#!/usr/bin/env python3
"""Mitochondrial DNA variants and heteroplasmy from the reads of an exome (or genome).

Exome captures pull in many mtDNA reads off target (about 50x in the GIAB IDT exomes), enough to call homoplasmic
variants and heteroplasmy down to ~10 % with read support. This tool piles up the reads on the mitochondrial
chromosome (MAPQ >= 20 removes most NUMT reads, base quality >= 20), calls every allele above the heteroplasmy cut-off
with enough supporting reads on both strands, and annotates it with:
  ClinVar (CLNSIG, review status, condition), gnomAD v3.1 chrM (homoplasmic and heteroplasmic frequency among 56,434
  genomes, maximum heteroplasmy seen, haplogroup-defining flag, common low-heteroplasmy flag, MitoTIP and PON-mt-tRNA
  predictions for tRNA variants) and the rCRS gene.

Haplogroup: with HaploGrep3 on the PATH (--haplogroup), the majority allele of every position covered by at least
--min-depth-hom reads is classified against PhyloTree 17 Forensic Update 1.2. Exomes cover the mtDNA unevenly, so
HaploGrep3 gets the covered ranges (HSD format): an uncovered position counts as unknown, not as the rCRS base.

Outputs: a TSV of calls (one row per allele) and a JSON summary (coverage, detection limit, counts, haplogroup).
Standard library plus samtools (and HaploGrep3 for the haplogroup).

    wes_mito.py --bam reads.bam --reference ref.fa [--gnomad-mt gnomad.chrM.vcf.bgz] [--clinvar clinvar.vcf.gz]
                --tsv calls.tsv --json summary.json [--min-hl 0.10] [--min-alt 3] [--min-depth 20] [--haplogroup]
"""
import argparse
import csv
import gzip
import os
import json
import re
import statistics
import shutil
import subprocess
import sys
import tempfile

VERSION = "0.2.0"
TAB = "\t"
TREE = "phylotree-fu-rcrs@1.2"  # PhyloTree 17 Forensic Update 1.2 (HaploGrep3 default)
# Below this fraction of the mtDNA covered the haplogroup is low confidence. Set on the tuning half of the 1000 Genomes
# exome panel (tests/haplogroup_eval.py): the first 0.1 step above every tuning sample whose call was less resolved
# than, or different from, the truth.
MIN_COVERED = 0.4
# rCRS (NC_012920.1) genes: start, end (1-based, inclusive), symbol
MT_GENES = [(577, 647, "MT-TF"), (648, 1601, "MT-RNR1"), (1602, 1670, "MT-TV"), (1671, 3229, "MT-RNR2"),
            (3230, 3304, "MT-TL1"), (3307, 4262, "MT-ND1"), (4263, 4331, "MT-TI"), (4329, 4400, "MT-TQ"),
            (4402, 4469, "MT-TM"), (4470, 5511, "MT-ND2"), (5512, 5579, "MT-TW"), (5587, 5655, "MT-TA"),
            (5657, 5729, "MT-TN"), (5761, 5826, "MT-TC"), (5826, 5891, "MT-TY"), (5904, 7445, "MT-CO1"),
            (7446, 7514, "MT-TS1"), (7518, 7585, "MT-TD"), (7586, 8269, "MT-CO2"), (8295, 8364, "MT-TK"),
            (8366, 8572, "MT-ATP8"), (8527, 9207, "MT-ATP6"), (9207, 9990, "MT-CO3"), (9991, 10058, "MT-TG"),
            (10059, 10404, "MT-ND3"), (10405, 10469, "MT-TR"), (10470, 10766, "MT-ND4L"), (10760, 12137, "MT-ND4"),
            (12138, 12206, "MT-TH"), (12207, 12265, "MT-TS2"), (12266, 12336, "MT-TL2"), (12337, 14148, "MT-ND5"),
            (14149, 14673, "MT-ND6"), (14674, 14742, "MT-TE"), (14747, 15887, "MT-CYB"), (15888, 15953, "MT-TT"),
            (15956, 16023, "MT-TP")]
COLS = ["pos", "ref", "alt", "gene", "depth", "alt_reads", "alt_fwd", "alt_rev", "heteroplasmy", "call",
        "gnomad_af_hom", "gnomad_af_het", "gnomad_max_hl", "haplogroup_defining", "common_low_heteroplasmy",
        "trna_prediction", "clinvar", "clinvar_review", "clinvar_condition", "flags"]
INDEL = re.compile(r"[+-](\d+)")
# MT-RNR1 variants of the CPIC aminoglycoside guideline (McDermott et al., Clin Pharmacol Ther 2022)
PGX_SITES = {1095: "T>C", 1494: "C>T", 1555: "A>G"}
MIXTURE_SITES = 3  # heteroplasmic calls at common haplogroup variants that flag a second lineage


def gene_of(pos):
    hits = [g for s, e, g in MT_GENES if s <= pos <= e]
    return ",".join(hits) if hits else ("control region" if pos >= 16024 or pos <= 576 else "intergenic")


def opener(path):
    with open(path, "rb") as f:
        gz = f.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if gz else open(path)


def mt_contig(bam):
    out = subprocess.run(["samtools", "idxstats", bam], capture_output=True, text=True, check=True).stdout
    for line in out.splitlines():
        c = line.split("\t")[0]
        if c in ("MT", "chrM", "M", "chrMT", "NC_012920.1"):
            return c
    raise SystemExit("no mitochondrial contig (MT/chrM) in the BAM")


def pileup(bam, ref, contig, min_mapq, min_bq):
    """Per position: depth, base counts by strand, indel counts by strand (allele as VCF ref/alt)."""
    cmd = ["samtools", "mpileup", "-r", contig, "-f", ref, "-q", str(min_mapq), "-Q", str(min_bq), "-B", "-d", "0",
           "-a", bam]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    for line in p.stdout:
        t = line.rstrip("\n").split("\t")
        pos, rb, depth, s = int(t[1]), t[2].upper(), int(t[3]), t[4] if len(t) > 4 else ""
        counts, indels, i = {}, {}, 0
        while i < len(s):
            c = s[i]
            if c == "^":
                i += 2
                continue
            if c == "$":
                i += 1
                continue
            if c in "+-":
                m = INDEL.match(s, i)
                n = int(m.group(1))
                seq = s[m.end():m.end() + n]
                i = m.end() + n
                key = (rb + seq.upper(), rb) if c == "-" else (rb, rb + seq.upper())
                fwd = seq.isupper()
                d = indels.setdefault(key, [0, 0])
                d[0 if fwd else 1] += 1
                continue
            if c in ".,":
                b, fwd = rb, c == "."
            elif c in "ACGTNacgtn":
                b, fwd = c.upper(), c.isupper()
            else:  # '*' (deleted base) and others
                i += 1
                continue
            d = counts.setdefault(b, [0, 0])
            d[0 if fwd else 1] += 1
            i += 1
        yield pos, rb, depth, counts, indels
    p.wait()
    if p.returncode:
        raise SystemExit("samtools mpileup failed: " + p.stderr.read()[:500])


def ranges(covered):
    """Sorted 1-based positions -> HaploGrep range string ("1-300;320-16569")."""
    out, start, prev = [], None, None
    for p in covered:
        if start is None:
            start = prev = p
        elif p == prev + 1:
            prev = p
        else:
            out.append(f"{start}-{prev}")
            start = prev = p
    if start is not None:
        out.append(f"{start}-{prev}")
    return ";".join(out)


def hsd_variants(pos, ref, alt):
    """VCF-style allele -> HaploGrep/PhyloTree notation (73G, 523d, 315.1C)."""
    if len(ref) == len(alt) == 1:
        return [f"{pos}{alt}"]
    if len(ref) > len(alt):  # deletion of ref[1:] after the anchor base
        return [f"{pos + k}d" for k in range(1, len(ref))]
    return [f"{pos}.{k}{b}" for k, b in enumerate(alt[1:], 1)]


def haplogroup(sample, covered, majority, tree, hits=3):
    """Classify the consensus profile with HaploGrep3 -> dict (best hit, alternatives, coverage of the tree)."""
    exe = shutil.which("haplogrep3")
    if not exe:
        return {"status": "not run: haplogrep3 not found"}
    if not covered:
        return {"status": "not determined: no mtDNA position covered"}
    d = tempfile.mkdtemp()
    hsd, out = os.path.join(d, "in.hsd"), os.path.join(d, "out.txt")
    # 3107 is the rCRS placeholder base (N): reads always "delete" it
    profile = [v for pos, ref, alt in majority for v in hsd_variants(pos, ref, alt) if v != "3107d"]
    with open(hsd, "w") as f:
        f.write(TAB.join([sample, ranges(covered), "?"] + profile) + "\n")
    p = subprocess.run([exe, "classify", "--tree", tree, "--in", hsd, "--out", out, "--extend-report",
                        "--hits", str(hits)], capture_output=True, text=True)
    if p.returncode or not os.path.exists(out):
        return {"status": "failed: " + (p.stderr or p.stdout)[-300:]}
    with open(out) as f:
        rows = list(csv.DictReader(f, delimiter=TAB))
    shutil.rmtree(d, ignore_errors=True)
    if not rows:
        return {"status": "failed: no HaploGrep3 result"}

    def split(x):
        """'152C (localPrivateMutation) 3107d (hotspot)' -> ['152C local private mutation', '3107d hotspot']"""
        out = []
        for v in (x or "").split():
            if v.startswith("(") and out:
                label = v.strip("()")
                out[-1] += " " + "".join(" " + c.lower() if c.isupper() else c for c in label).strip()
            else:
                out.append(v)
        return out
    best = rows[0]
    q = float(best["Quality"])
    alts = [{"haplogroup": r["Haplogroup"], "quality": round(float(r["Quality"]), 4)} for r in rows[1:]]
    return {"status": "ok", "tree": tree, "haplogroup": best["Haplogroup"], "quality": round(q, 4),
            "confidence": "low" if len(covered) / 16569 < MIN_COVERED else
                          "high" if q >= 0.9 else "medium" if q >= 0.8 else "low",
            "note": f"only {len(covered) / 16569:.0%} of the mtDNA covered: the haplogroup may be less resolved or "
                    "wrong (HaploGrep3's quality counts the covered positions only)"
                    if len(covered) / 16569 < MIN_COVERED else "",
            "covered_positions": len(covered), "covered_fraction": round(len(covered) / 16569, 4),
            "profile_variants": len(profile),
            "expected_not_found": split(best.get("Not_Found_Polys")),
            "extra_variants": split(best.get("Remaining_Polys")),
            "alternatives": alts}


def load_vcf_mt(path, fields):
    """(pos, ref, alt) -> {field: value} for the mitochondrial records of a VCF."""
    out = {}
    if not path:
        return out
    with opener(path) as f:
        for line in f:
            if line[0] == "#":
                continue
            t = line.rstrip("\r\n").split("\t", 8)
            if t[0] not in ("MT", "chrM", "M", "NC_012920.1"):
                continue
            info = {}
            for kv in t[7].split(";"):
                k, _, v = kv.partition("=")
                if k in fields:
                    info[k] = v if v else True
            for alt in t[4].split(","):
                out[(int(t[1]), t[3], alt)] = info
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--bam", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--gnomad-mt")
    ap.add_argument("--clinvar")
    ap.add_argument("--min-hl", type=float, default=0.10, help="lowest heteroplasmy reported (default 0.10)")
    ap.add_argument("--min-alt", type=int, default=3, help="supporting reads needed (default 3)")
    ap.add_argument("--min-depth", type=int, default=20, help="depth needed to call heteroplasmy (default 20)")
    ap.add_argument("--min-depth-hom", type=int, default=5,
                    help="depth needed to call a (near-)homoplasmic variant, >= 90 %% of reads (default 5)")
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--min-bq", type=int, default=20)
    ap.add_argument("--tsv", required=True)
    ap.add_argument("--json", required=True)
    ap.add_argument("--haplogroup", action="store_true", help="classify the mtDNA haplogroup with HaploGrep3")
    ap.add_argument("--tree", default=TREE, help=f"HaploGrep3 tree (default {TREE})")
    ap.add_argument("--sample", default="sample")
    a = ap.parse_args(argv)
    contig = mt_contig(a.bam)
    gn = load_vcf_mt(a.gnomad_mt, {"AF_hom", "AF_het", "max_hl", "hap_defining_variant", "common_low_heteroplasmy",
                                  "mitotip_trna_prediction", "pon_mt_trna_prediction", "filters"})
    cv = load_vcf_mt(a.clinvar, {"CLNSIG", "CLNREVSTAT", "CLNDN"})
    depths, rows, pgx_depth, covered, majority = [], [], {}, [], []
    for pos, rb, depth, counts, indels in pileup(a.bam, a.reference, contig, a.min_mapq, a.min_bq):
        depths.append(depth)
        if depth >= a.min_depth_hom and "N" not in rb:  # consensus for the haplogroup: majority allele (> 50 %)
            covered.append(pos)
            for (ref, alt), (f, r) in [((rb, b), v) for b, v in counts.items() if b not in (rb, "N")] + list(indels.items()):
                if (f + r) * 2 > depth and "N" not in alt:
                    majority.append((pos, ref, alt))
        if pos in PGX_SITES:
            pgx_depth[str(pos)] = depth
        alleles = [((rb, b), v) for b, v in counts.items() if b != rb and b != "N"] + list(indels.items())
        for (ref, alt), (f, r) in alleles:
            n = f + r
            if "N" in ref or "N" in alt or not depth:
                continue  # rCRS 3107N placeholder
            hl = n / depth
            homo = hl >= 0.9 and depth >= a.min_depth_hom
            if not homo and (depth < a.min_depth or n < a.min_alt or hl < a.min_hl):
                continue
            flags = []
            if min(f, r) == 0 and n >= 6:
                flags.append("one strand only")
            if pos in range(302, 317) or pos in range(513, 526) or pos in range(16182, 16195) or pos == 3107:
                flags.append("homopolymer/artefact-prone (302-316, 513-525, 16182-16194, 3107N)")
            if hl < 0.5 and len(ref) == len(alt):
                flags.append("heteroplasmic: confirm (NUMT reads or contamination can mimic low levels)")
            g = gn.get((pos, ref, alt), {})
            c = cv.get((pos, ref, alt), {})
            trna = "; ".join(x for x in (g.get("mitotip_trna_prediction") and f"MitoTIP {g['mitotip_trna_prediction']}",
                                         g.get("pon_mt_trna_prediction") and f"PON-mt-tRNA {g['pon_mt_trna_prediction']}")
                             if x and "None" not in x)
            rows.append([pos, ref, alt, gene_of(pos), depth, n, f, r, f"{hl:.3f}",
                         "homoplasmic" if hl >= 0.9 else "heteroplasmic", g.get("AF_hom", ""), g.get("AF_het", ""),
                         g.get("max_hl", ""), "yes" if g.get("hap_defining_variant") else "",
                         "yes" if g.get("common_low_heteroplasmy") else "", trna.replace("_", " "),
                         str(c.get("CLNSIG", "")), str(c.get("CLNREVSTAT", "")).replace("_", " "),
                         str(c.get("CLNDN", "")).replace("_", " "), "; ".join(flags)])
    with open(a.tsv, "w") as f:
        f.write("\t".join(COLS) + "\n")
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")
    cov = depths or [0]
    # A second mtDNA lineage (sample contamination, or NUMT reads in a low-coverage exome) shows as heteroplasmy at
    # positions that are common homoplasmic variants of other haplogroups (gnomAD AF_hom >= 1 %), like haplocheck.
    mix = [r for r in rows if r[9] == "heteroplasmic" and r[10] not in ("", None) and float(r[10]) >= 0.01]
    for r in mix:
        r[19] = "; ".join(x for x in (r[19], "common variant of another haplogroup at heteroplasmic level") if x)
    summary = {"tool": f"wes_mito {VERSION}", "contig": contig, "length": len(depths),
               "mean_depth": round(statistics.mean(cov), 1), "median_depth": statistics.median(cov),
               "fraction_ge_5x": round(sum(d >= 5 for d in cov) / max(len(cov), 1), 4),
               "fraction_ge_20x": round(sum(d >= 20 for d in cov) / max(len(cov), 1), 4),
               "fraction_ge_100x": round(sum(d >= 100 for d in cov) / max(len(cov), 1), 4),
               "min_hl": a.min_hl, "min_alt": a.min_alt, "min_depth": a.min_depth,
               "detection_limit_at_median": round(max(a.min_hl, a.min_alt / max(statistics.median(cov), 1)), 3),
               "homoplasmic": sum(r[9] == "homoplasmic" for r in rows),
               "heteroplasmic": sum(r[9] == "heteroplasmic" for r in rows),
               "mixture_sites": len(mix),
               "mixture": "possible second mtDNA lineage (contamination or NUMT reads): heteroplasmy at "
                          f"{len(mix)} common haplogroup variant sites" if len(mix) >= MIXTURE_SITES else "",
               "pgx_depth": pgx_depth, "gnomad_mt": bool(gn), "clinvar": bool(cv),
               "haplogroup": haplogroup(a.sample, covered, majority, a.tree) if a.haplogroup
               else {"status": "not run"}}
    with open(a.json, "w") as f:
        json.dump(summary, f, indent=1)
    print(json.dumps(summary), file=sys.stderr)


if __name__ == "__main__":
    main()
