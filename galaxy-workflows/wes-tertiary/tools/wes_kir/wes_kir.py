#!/usr/bin/env python3
"""KIR gene content, alleles, A/B genotype and B-content group from exome reads with T1K.

    wes_kir.py --bam sample.bam --kir-dat kir.dat --sample NAME --tsv kir.tsv --json kir.json [--threads 8]
        [--hla wes_hla_summary.json]

1. Builds T1K's allele reference from the IPD-KIR database file (kir.dat; t1k-build.pl, about 2 s).
2. Takes the reads of the KIR region (the leukocyte receptor complex on chromosome 19, named like the BAM: 19 or
   chr19) and the unmapped reads, and runs T1K (Song et al., Genome Res 2023; --preset kir-wes) on them. A KIR
   haplotype carries 4-14 genes and the reference genome only one haplotype, so reads of the other genes map to
   their paralogues in the region or not at all; T1K assigns every read to the alleles of all KIR genes.
3. A gene is **present** when its best allele has T1K quality >= --min-quality and an abundance of at least
   --min-abundance of the framework genes' mean (KIR3DL3, KIR3DP1, KIR2DL4, KIR3DL2: present on nearly every
   haplotype); otherwise absent. Alleles are given at three digits (protein level).
4. From the gene content: the A/B genotype (AA: only genes of the A haplotype), the centromeric and telomeric motifs
   (Cen-A/B, Tel-A/B) and the B-content group used in donor selection for haematopoietic cell transplantation
   (Cooley et al., Blood 2010: Best = Cen-B/B, Better = two or more B motifs, Neutral = otherwise).
5. With the wes_hla summary: the HLA ligands of the inhibitory KIR (HLA-C group C1/C2 for KIR2DL2/3 and KIR2DL1,
   HLA-A/B Bw4 for KIR3DL1) as given by wes_hla.
Standard library plus t1k-build.pl, run-t1k and samtools on the PATH.
"""
import argparse
import json
import os
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
import time

VERSION = "0.1.0"
# Leukocyte receptor complex, KIR part (GRCh38; KIR3DL3 starts at 54.72 Mb, KIR3DL2 ends at 54.86 Mb), with margin
KIR_REGION = (54_600_000, 55_000_000)
GENES = ["KIR3DL3", "KIR2DS2", "KIR2DL2", "KIR2DL3", "KIR2DL5B", "KIR2DS3", "KIR2DP1", "KIR2DL1", "KIR3DP1",
         "KIR2DL4", "KIR3DL1", "KIR3DS1", "KIR2DL5A", "KIR2DS5", "KIR2DS1", "KIR2DS4", "KIR3DL2"]
FRAMEWORK = ["KIR3DL3", "KIR3DP1", "KIR2DL4", "KIR3DL2"]
PSEUDO = {"KIR2DP1", "KIR3DP1"}
A_GENES = {"KIR3DL3", "KIR2DL3", "KIR2DP1", "KIR2DL1", "KIR3DP1", "KIR2DL4", "KIR3DL1", "KIR2DS4", "KIR3DL2"}
B_GENES = {"KIR2DS2", "KIR2DL2", "KIR2DL5A", "KIR2DL5B", "KIR2DS3", "KIR2DS5", "KIR3DS1", "KIR2DS1"}
MIN_QUALITY = 20       # set on the tuning half of the 1000 Genomes panel (tests/kir_eval.py)
MIN_ABUNDANCE = 0.15   # fraction of the framework genes' mean abundance (same)


def bam_chrom(bam):
    out = subprocess.run(["samtools", "view", "-H", bam], capture_output=True, text=True, check=True).stdout
    names = re.findall(r"^@SQ\tSN:(\S+)", out, re.M)
    return "chr19" if "chr19" in names else "19"


def build_reference(kir_dat, d):
    subprocess.run(["t1k-build.pl", "-d", kir_dat, "-o", d, "--prefix", "kir"], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return os.path.join(d, "kir_dna_seq.fa")


def db_release(kir_dat):
    """'CC   IPD-KIR Release Version 2.15.0' in the first record -> '2.15.0'."""
    with open(kir_dat) as f:
        for line in f:
            m = re.search(r"IPD-KIR Release Version (\d+\.\d+\.\d+)", line)
            if m:
                return m.group(1)
            if line.startswith("SQ"):
                break
    return "release not stated"


def extract(bam, regions, unmapped, work, threads):
    """Reads of the regions (and the unmapped reads) -> paired FASTQ files."""
    r1, r2 = os.path.join(work, "r1.fq"), os.path.join(work, "r2.fq")
    view = ["samtools", "view", "-u", "-F", "0x900", bam] + regions + (["*"] if unmapped else [])
    p1 = subprocess.Popen(view, stdout=subprocess.PIPE)
    p2 = subprocess.Popen(["samtools", "collate", "-u", "-O", "-", os.path.join(work, "collate")],
                          stdin=p1.stdout, stdout=subprocess.PIPE)
    p1.stdout.close()
    subprocess.run(["samtools", "fastq", "-1", r1, "-2", r2, "-s", "/dev/null", "-0", "/dev/null", "-"],
                   stdin=p2.stdout, check=True, stderr=subprocess.DEVNULL)
    p2.wait()
    if p1.wait() or p2.returncode:
        raise SystemExit("samtools could not read the KIR region of " + bam)
    with open(r1) as f:
        pairs = sum(1 for _ in f) // 4
    return r1, r2, pairs


def parse_genotype(path):
    calls = {}
    for line in open(path):
        t = line.rstrip("\n").split("\t")
        if len(t) < 8:
            continue
        calls[t[0]] = {"n": int(t[1]), "a1": t[2], "ab1": float(t[3]), "q1": int(t[4]), "a2": t[5],
                       "ab2": float(t[6]), "q2": int(t[7])}
    return calls


def genotypes(calls, min_q=MIN_QUALITY, min_ab=MIN_ABUNDANCE):
    """T1K calls -> one row per gene: presence, alleles (three digits), qualities, abundance relative to framework."""
    fw = [calls[g]["ab1"] + max(calls[g]["ab2"], 0) for g in FRAMEWORK if g in calls and calls[g]["q1"] >= 0]
    fw_mean = statistics.mean(fw) if fw else 0.0
    rows = []
    for g in GENES:
        c = calls.get(g)
        if not c or c["a1"] == ".":
            rows.append({"gene": g, "present": False, "alleles": "", "quality": -1, "relative_abundance": 0.0,
                         "note": "no reads assigned"})
            continue
        total = c["ab1"] + max(c["ab2"], 0)
        rel = round(total / fw_mean, 3) if fw_mean else 0.0
        best_q = max(c["q1"], c["q2"])
        present = best_q >= min_q and rel >= min_ab
        alleles = [a.split(",")[0] for a, q in ((c["a1"], c["q1"]), (c["a2"], c["q2"])) if a != "." and q >= min_q]
        note = "" if present else (f"absent: quality {best_q} < {min_q}" if best_q < min_q
                                   else f"absent: abundance {rel:.2f} of the framework genes < {min_ab}")
        rows.append({"gene": g, "present": present, "alleles": "/".join(alleles) if present else "",
                     "quality": best_q, "relative_abundance": rel, "note": note})
    return rows, fw_mean


def haplotypes(rows):
    """Gene content -> A/B genotype, centromeric/telomeric motifs and B-content group (Cooley et al. 2010)."""
    p = {r["gene"] for r in rows if r["present"]}
    out = {"genes_present": len(p), "framework_complete": all(g in p for g in FRAMEWORK)}
    out["genotype"] = "AA" if not (p & B_GENES) else "Bx"
    # centromeric: KIR2DL3 (A) / KIR2DL2 or KIR2DS2 (B); telomeric: KIR3DL1 + KIR2DS4 (A) / KIR3DS1 or KIR2DS1 (B)
    cen_a, cen_b = "KIR2DL3" in p, bool(p & {"KIR2DL2", "KIR2DS2"})
    tel_a, tel_b = bool(p & {"KIR3DL1", "KIR2DS4"}), bool(p & {"KIR3DS1", "KIR2DS1"})

    def motif(a, b):
        return "A/A" if a and not b else "B/B" if b and not a else "A/B" if a and b else "undetermined"
    out["centromeric"], out["telomeric"] = motif(cen_a, cen_b), motif(tel_a, tel_b)
    if "undetermined" in (out["centromeric"], out["telomeric"]):
        out["b_content_score"], out["b_content_group"] = None, "undetermined"
    else:
        score = out["centromeric"].count("B") + out["telomeric"].count("B")
        out["b_content_score"] = score
        out["b_content_group"] = ("Best" if out["centromeric"] == "B/B" else "Better" if score >= 2 else "Neutral")
    return out


def ligands(rows, hla):
    """KIR-HLA ligand pairs from the KIR gene content and the KIR ligand groups wes_hla gives per HLA allele."""
    p = {r["gene"] for r in rows if r["present"]}
    groups = {}
    for g in (hla or {}).get("genotypes", []):
        for k in ("kir_ligand1", "kir_ligand2"):
            if g.get(k) and g.get("confidence") == "high":
                groups.setdefault(g[k], []).append(g["allele1" if k.endswith("1") else "allele2"])
    pairs = []
    for kir, lig, kind in (("KIR2DL1", "C2", "inhibitory"), ("KIR2DL2", "C1", "inhibitory"),
                           ("KIR2DL3", "C1", "inhibitory"), ("KIR3DL1", "Bw4", "inhibitory"),
                           ("KIR2DS1", "C2", "activating"), ("KIR3DS1", "Bw4-80I", "activating (putative)")):
        have = sorted(set(a for grp, al in groups.items() if grp == lig or grp.startswith(lig + "-") for a in al))
        if lig == "Bw4-80I":
            have = sorted(set(groups.get("Bw4-80I", [])))
        pairs.append({"kir": kir, "kir_present": kir in p, "ligand": lig, "kind": kind,
                      "ligand_alleles": have, "status": ("ligand not typed" if not hla else
                                                         "pair present" if kir in p and have else
                                                         "KIR without its ligand" if kir in p else
                                                         "KIR absent")})
    return pairs


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--bam", required=True)
    ap.add_argument("--kir-dat", help="IPD-KIR kir.dat (without it nothing is typed)")
    ap.add_argument("--hla", help="wes_hla summary (JSON) for the KIR ligands")
    ap.add_argument("--sample", default="sample")
    ap.add_argument("--tsv", required=True)
    ap.add_argument("--json", required=True)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--min-quality", type=int, default=MIN_QUALITY)
    ap.add_argument("--min-abundance", type=float, default=MIN_ABUNDANCE)
    ap.add_argument("--regions", nargs="*", help="regions to read (default: the KIR region on 19/chr19)")
    ap.add_argument("--no-unmapped", action="store_true", help="do not add the unmapped reads")
    ap.add_argument("--preset", default="kir-wes", help="T1K preset (kir-wes; kir-wgs for genomes)")
    ap.add_argument("--keep", help="copy T1K's raw output files into this directory")
    a = ap.parse_args(argv)
    t0 = time.time()
    work = tempfile.mkdtemp(prefix="wes_kir_", dir=os.getcwd())
    hla = json.load(open(a.hla)) if a.hla else None
    pairs_n, fw_mean, calls = 0, 0.0, {}
    if a.kir_dat:
        seq = build_reference(a.kir_dat, os.path.join(work, "ref"))
        regions = a.regions or [f"{bam_chrom(a.bam)}:{KIR_REGION[0]}-{KIR_REGION[1]}"]
        r1, r2, pairs_n = extract(a.bam, regions, not a.no_unmapped, work, a.threads)
        if pairs_n:
            subprocess.run(["run-t1k", "-1", r1, "-2", r2, "-f", seq, "--preset", a.preset, "-t", str(a.threads),
                            "-o", "t1k", "--od", work], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
            calls = parse_genotype(os.path.join(work, "t1k_genotype.tsv"))
    rows, fw_mean = genotypes(calls, a.min_quality, a.min_abundance)
    if not a.kir_dat:
        for r in rows:
            r["note"] = "not typed: no KIR allele database (kir.dat) given"
    hap = haplotypes(rows) if calls else {"genotype": "not typed", "b_content_group": "undetermined"}
    if calls and not hap["framework_complete"]:
        hap["note"] = ("framework gene(s) " + ", ".join(g for g in FRAMEWORK if not any(
            r["gene"] == g and r["present"] for r in rows)) + " not detected: low coverage of the KIR region, "
            "or a rare haplotype; gene absences are uncertain")
    lig = ligands(rows, hla) if calls else []
    cols = ["gene", "present", "alleles", "quality", "relative_abundance", "note"]
    with open(a.tsv, "w", newline="\n") as o:
        o.write("sample\t" + "\t".join(cols) + "\n")
        for r in rows:
            o.write(a.sample + "\t" + "\t".join(str(r[c]) for c in cols) + "\n")
    json.dump({"sample": a.sample, "tool": f"wes_kir {VERSION}, T1K (--preset {a.preset})",
               "database": db_release(a.kir_dat) if a.kir_dat else "not given",
               "min_quality": a.min_quality, "min_abundance": a.min_abundance, "read_pairs": pairs_n,
               "framework_mean_abundance": round(fw_mean, 2), "genes": rows, "haplotypes": hap, "ligands": lig,
               "seconds": round(time.time() - t0, 1)}, open(a.json, "w"), indent=1)
    if a.keep:
        os.makedirs(a.keep, exist_ok=True)
        for f in os.listdir(work):
            if f.startswith("t1k_") and not f.endswith((".fq", ".fa")):
                shutil.copy(os.path.join(work, f), a.keep)
    shutil.rmtree(work, ignore_errors=True)
    print(f"{a.sample}: {sum(r['present'] for r in rows)} KIR genes present, genotype {hap.get('genotype')}, "
          f"{pairs_n} read pairs, {time.time() - t0:.0f} s", file=sys.stderr)


if __name__ == "__main__":
    main()
