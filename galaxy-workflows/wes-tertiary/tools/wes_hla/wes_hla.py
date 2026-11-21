#!/usr/bin/env python3
"""HLA typing from exome reads with T1K, pharmacogenomic HLA alleles for PharmCAT, and HLA disease associations.

    wes_hla.py --bam sample.bam --hla-dat hla.dat --sample NAME --tsv types.tsv --json types.json \\
        --pharmcat outside_calls.tsv [--threads 8] [--min-quality 20]

1. Builds T1K's allele reference from the IPD-IMGT/HLA database file (hla.dat; t1k-build.pl, about 15 s) and the
   gene-position file T1K needs for BAM input from the GRCh38 positions below (Ensembl 115, as t1k-build.pl -g
   would write them), named like the BAM's chromosomes (6 or chr6).
2. Runs T1K (Song et al., Genome Res 2023; --preset hla) on the reads of the HLA genes and the unmapped reads.
3. Writes per classical gene the two alleles at two-field resolution (protein level, e.g. A*01:01), T1K's quality
   (0-60) and every equally good candidate; a call is "low" confidence below --min-quality or when the candidates
   differ at two fields.
4. Lists the pharmacogenomic alleles of CPIC guidelines and the best-known disease associations (risk alleles, not
   diagnoses), and writes PharmCAT's outside-call file for HLA-A and HLA-B (only confident calls: PharmCAT then
   gives the CPIC/DPWG/FDA recommendations for abacavir, allopurinol, carbamazepine, oxcarbazepine, phenytoin).
Standard library plus t1k-build.pl, run-t1k and samtools on the PATH.
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time

VERSION = "0.1.0"
# GRCh38 gene positions (Ensembl 115 GTF, 1-based) used by T1K to extract the reads of each gene from a BAM.
# Genes absent from the primary assembly (DRB3/DRB4 haplotypes, pseudogenes R/Y) get no position: T1K then uses
# the unmapped reads only, as t1k-build.pl does.
GENE_POS = {
    "HLA-A": ("6", 29941260, 29949572, "+"), "HLA-B": ("6", 31353867, 31367067, "-"),
    "HLA-C": ("6", 31268749, 31272130, "-"), "HLA-DMA": ("6", 32948613, 32969094, "-"),
    "HLA-DMB": ("6", 32934629, 32941054, "-"), "HLA-DOA": ("6", 33004178, 33009594, "-"),
    "HLA-DOB": ("6", 32812763, 32820466, "-"), "HLA-DPA1": ("6", 33064569, 33080781, "-"),
    "HLA-DPA2": ("6", 33091753, 33097295, "-"), "HLA-DPB1": ("6", 33073786, 33089696, "+"),
    "HLA-DPB2": ("6", 33112516, 33129084, "+"), "HLA-DQA1": ("6", 32628179, 32647062, "+"),
    "HLA-DQA2": ("6", 32741391, 32747198, "+"), "HLA-DQB1": ("6", 32659467, 32668383, "-"),
    "HLA-DQB2": ("6", 32756094, 32763543, "-"), "HLA-DRA": ("6", 32439872, 32445047, "+"),
    "HLA-DRB1": ("6", 32577902, 32589848, "-"), "HLA-DRB5": ("6", 32517345, 32530369, "-"),
    "HLA-E": ("6", 30489503, 30494205, "+"), "HLA-F": ("6", 29722738, 29738528, "+"),
    "HLA-G": ("6", 29826967, 29831125, "+"), "HLA-H": ("6", 29887752, 29890482, "+"),
    "HLA-HFE": ("6", 26087226, 26098343, "+"), "HLA-J": ("6", 30006606, 30009539, "+"),
    "HLA-K": ("6", 29926459, 29929232, "+"), "HLA-L": ("6", 30259625, 30261703, "+"),
    "HLA-N": ("6", 30351416, 30351550, "+"), "HLA-P": ("6", 29800415, 29802425, "+"),
    "HLA-S": ("6", 31382074, 31382288, "-"), "HLA-T": ("6", 29896654, 29897786, "+"),
    "HLA-U": ("6", 29934101, 29934286, "+"), "HLA-V": ("6", 29792234, 29793136, "+"),
    "HLA-W": ("6", 29956596, 29958570, "+"),
}
CLASSICAL = ["HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5", "HLA-DQA1", "HLA-DQB1",
             "HLA-DPA1", "HLA-DPB1"]
PGX = [  # allele, drugs, effect, guideline
    ("B*57:01", "abacavir", "hypersensitivity reaction", "CPIC 2014 (Martin et al.)"),
    ("B*15:02", "carbamazepine, oxcarbazepine, phenytoin, fosphenytoin", "Stevens-Johnson syndrome / TEN",
     "CPIC 2017 (Phillips et al.)"),
    ("A*31:01", "carbamazepine", "SJS/TEN, DRESS and maculopapular exanthema", "CPIC 2017 (Phillips et al.)"),
    ("B*58:01", "allopurinol", "severe cutaneous adverse reactions (SCAR)", "CPIC 2015 (Saito et al.)"),
]
DISEASE = [  # allele (one or two fields), condition, note
    ("B*27", "ankylosing spondylitis / axial spondyloarthritis",
     "carried by about 90 % of patients of European ancestry; most B*27 carriers never develop the disease"),
    ("B*51", "Behcet disease", "the strongest genetic risk factor (odds ratio about 6)"),
    ("DQB1*06:02", "narcolepsy type 1", "carried by about 98 % of patients with hypocretin deficiency and by about "
     "25 % of the general population of European ancestry"),
]
CELIAC = [  # name, DQA1 prefix, DQB1 prefix
    ("DQ2.5", "DQA1*05", "DQB1*02:01"), ("DQ2.2", "DQA1*02:01", "DQB1*02:02"), ("DQ8", "DQA1*03", "DQB1*03:02"),
]


def two_field(a):
    """HLA-A*01:01:01:01 -> A*01:01 (keeps a trailing expression letter, e.g. DQB1*02:163N)."""
    g, _, f = a.partition("*")
    p = f.split(":")
    s = ":".join(p[:2])
    if len(p) > 2 and re.search(r"[A-Z]$", p[-1]) and not re.search(r"[A-Z]$", s):
        s += p[-1][-1]
    return g.replace("HLA-", "") + "*" + s


def allele_key(a):
    """A*02:01 < A*11:01 < A*11:01N: numeric order of the fields."""
    return [(int(x), "") if x.isdigit() else (-1, x) for x in re.findall(r"\d+|\D+", a.split("*")[1])]


def bam_chrom(bam):
    out = subprocess.run(["samtools", "view", "-H", bam], capture_output=True, text=True, check=True).stdout
    names = re.findall(r"^@SQ\tSN:(\S+)", out, re.M)
    return "chr6" if "chr6" in names else "6"


def build_reference(hla_dat, d, chrom6):
    subprocess.run(["t1k-build.pl", "-d", hla_dat, "-o", d, "--target", "HLA", "--prefix", "hla"], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    seq = os.path.join(d, "hla_dna_seq.fa")
    coord = os.path.join(d, "hla_dna_coord.fa")
    other = "chr19" if chrom6.startswith("chr") else "19"
    with open(seq) as f, open(coord, "w") as o:
        for line in f:
            if line.startswith(">"):
                allele = line[1:].split()[0]
                c, s, e, strand = GENE_POS.get(allele.split("*")[0], (None, -1, -1, "+"))
                o.write(f">{allele} {chrom6 if c else other} {s} {e} {strand}\n")
            else:
                o.write(line)
    return seq, coord


def db_release(hla_dat):
    """'IPD-IMGT/HLA Release Version 3.66.0' in the first record's comments -> '3.66.0'."""
    with open(hla_dat) as f:
        for line in f:
            m = re.search(r"Release Version (\d+\.\d+\.\d+)", line)
            if m:
                return m.group(1)
            if line.startswith("SQ"):
                break
    return "release not stated"


def parse_genotype(path):
    calls = {}
    for line in open(path):
        t = line.rstrip("\n").split("\t")
        if len(t) < 8:
            continue
        calls[t[0]] = {"n": int(t[1]), "a1": t[2], "ab1": float(t[3]), "q1": int(t[4]), "a2": t[5],
                       "ab2": float(t[6]), "q2": int(t[7])}
    return calls


CLASS_I = {"HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G"}


def ars_sequences(seq_fa, names):
    """Antigen-recognition-site exons (class I exons 2+3, class II exon 2) of T1K's candidate alleles.

    T1K names candidates by three fields (HLA-DQB1*02:01:01); the first reference allele with that prefix stands for
    it. The reference headers give the exon positions ('>allele n s1 e1 s2 e2 ...', t1k-build.pl)."""
    want = set(names)
    out, cur, keep = {}, None, False
    for line in open(seq_fa):
        if line[0] == ">":
            if cur and keep:
                out.setdefault(cur[0], (cur[1], "".join(cur[2])))
            t = line[1:].split()
            f = t[0].split(":")
            key = next((":".join(f[:k]) for k in range(len(f), 0, -1) if ":".join(f[:k]) in want), None)
            keep = key is not None and key not in out
            cur = (key, [int(x) for x in t[2:]], []) if keep else None
        elif keep:
            cur[2].append(line.strip())
    if cur and keep:
        out.setdefault(cur[0], (cur[1], "".join(cur[2])))
    ars = {}
    for name, (c, seq) in out.items():
        exons = (2, 3) if name.split("*")[0] in CLASS_I else (2,)
        ars[name] = "|".join(seq[c[2 * (e - 1)]:c[2 * (e - 1) + 1] + 1] for e in exons if 2 * e <= len(c))
    return ars


def resolve(field, ars):
    """T1K's candidates of one allele -> (call, two-field candidates, number of ARS groups, equivalent alleles).

    Candidates with identical antigen-recognition exons (one G group) cannot be told apart from exon reads: the call
    is then the expressed allele with the lowest number in that group (null and other expression-variant alleles,
    N/L/S/Q..., mostly differ outside these exons and are rare). Different ARS exons = a real ambiguity."""
    if field in (".", ""):
        return None, [], 0, []
    names = field.split(",")
    groups = []
    for n in names:
        seq = ars.get(n, n)  # an allele without a reference sequence is its own group
        for g in groups:
            if g[0] == seq:
                g[1].append(n)
                break
        else:
            groups.append((seq, [n]))
    first = [two_field(n) for n in groups[0][1]]
    expressed = sorted({x for x in first if not re.search(r"[A-Z]$", x)}, key=allele_key)
    call = (expressed or first)[0]
    cands = []
    for n in names:
        x = two_field(n)
        if x not in cands:
            cands.append(x)
    return call, cands, len(groups), sorted(set(first) - {call}, key=allele_key)


def genotypes(calls, min_q, ars=None):
    ars = ars or {}
    out = []
    for g in CLASSICAL:
        c = calls.get(g)
        row = {"gene": g, "allele1": "", "allele2": "", "quality1": "", "quality2": "", "confidence": "",
               "candidates1": "", "candidates2": "", "note": ""}
        if not c or c["n"] == 0:
            row["note"] = "not typed (no reads, or the gene is absent from this haplotype)" \
                if g in ("HLA-DRB3", "HLA-DRB4", "HLA-DRB5") else "not typed (too few reads)"
            row["confidence"] = "none"
            out.append(row)
            continue
        a1, c1, g1, e1 = resolve(c["a1"], ars)
        a2, c2, g2, e2 = resolve(c["a2"], ars) if c["n"] > 1 else (a1, c1, g1, e1)
        q1, q2 = c["q1"], (c["q2"] if c["n"] > 1 else c["q1"])
        if allele_key(a2) < allele_key(a1):
            a1, a2, c1, c2, q1, q2, g1, g2, e1, e2 = a2, a1, c2, c1, q2, q1, g2, g1, e2, e1
        row.update(allele1=a1, allele2=a2, quality1=q1, quality2=q2, candidates1=",".join(c1),
                   candidates2=",".join(c2))
        notes = []
        if c["n"] == 1 and g in ("HLA-DRB3", "HLA-DRB4", "HLA-DRB5"):
            notes.append("one allele found: one copy, or two of the same (each DR haplotype carries at most one of "
                         "DRB3, DRB4, DRB5)")
        elif c["n"] == 1:
            notes.append("one allele found: homozygous, or the second allele was missed")
        amb = g1 > 1 or g2 > 1
        if amb:
            notes.append("candidates differ in the antigen-recognition exons")
        eq = [x for x in dict.fromkeys(e1 + e2)]
        if eq:
            notes.append("same antigen-recognition exons (G group) as " + ", ".join(eq))
        low = min(q1, q2) < min_q or amb
        if min(q1, q2) < min_q:
            notes.append(f"quality below {min_q}")
        row["confidence"] = "low" if low else "high"
        row["note"] = "; ".join(notes)
        out.append(row)
    dr_haplotypes(out)
    return out


# Each DR haplotype carries at most one of the secondary DRB genes, fixed by its DRB1 allele group (Marsh et al.,
# The HLA FactsBook): DRB3 with DR3/11/12/13/14 (DR52), DRB4 with DR4/7/9 (DR53), DRB5 with DR15/16 (DR51), none
# with DR1/8/10.
DRB_PARTNER = {"03": "HLA-DRB3", "11": "HLA-DRB3", "12": "HLA-DRB3", "13": "HLA-DRB3", "14": "HLA-DRB3",
               "04": "HLA-DRB4", "07": "HLA-DRB4", "09": "HLA-DRB4", "15": "HLA-DRB5", "16": "HLA-DRB5"}


def dr_haplotypes(rows):
    """Mark DRB3/4/5 calls that the DRB1 alleles do not allow (reads of a related DRB gene) low confidence."""
    by = {r["gene"]: r for r in rows}
    d = by.get("HLA-DRB1")
    if not d or d["confidence"] != "high":
        return
    groups = [a.split("*")[1].split(":")[0] for a in (d["allele1"], d["allele2"])]
    expected = {DRB_PARTNER[x] for x in groups if x in DRB_PARTNER}
    for g in ("HLA-DRB3", "HLA-DRB4", "HLA-DRB5"):
        r = by[g]
        if r["confidence"] != "none" and g not in expected:
            r["confidence"] = "low"
            r["note"] = "; ".join(x for x in (r["note"], f"not expected with {d['allele1']}/{d['allele2']} (reads of "
                                              "a related DRB gene)") if x)
        elif r["confidence"] == "none" and g in expected:
            r["note"] = f"not typed, although {d['allele1']}/{d['allele2']} carry it (too few reads)"


def has(rows, allele):
    """Alleles of the genotype matching 'B*57:01' (two fields) or 'B*27' (allele group)."""
    gene = "HLA-" + allele.split("*")[0]
    n = 0
    for r in rows:
        if r["gene"] != gene or r["confidence"] == "none":
            continue
        for a in (r["allele1"], r["allele2"]):
            if a == allele or (":" not in allele and a.startswith(allele + ":")):
                n += 1
    return n


def associations(rows):
    conf = {r["gene"]: r["confidence"] for r in rows}
    out = []
    for allele, drugs, effect, src in PGX:
        n = has(rows, allele)
        gene = "HLA-" + allele.split("*")[0]
        out.append({"kind": "pharmacogenomic", "allele": allele, "copies": n, "condition": f"{drugs}: {effect}",
                    "status": ("present" if n else "absent") if conf.get(gene) == "high" else
                    ("present (low confidence)" if n else "not determined"), "source": src})
    for allele, cond, note in DISEASE:
        n = has(rows, allele)
        gene = "HLA-" + allele.split("*")[0]
        if n or conf.get(gene) != "high":
            out.append({"kind": "disease association", "allele": allele, "copies": n, "condition": cond,
                        "status": ("present" if n else "not determined") + (" (low confidence)" if n and conf.get(gene) != "high" else ""),
                        "source": note})
    dq = [name for name, a, b in CELIAC if has(rows, a) and has(rows, b)]
    sure = conf.get("HLA-DQA1") == "high" and conf.get("HLA-DQB1") == "high"
    # DQ2/DQ8 need DQB1*02:01, *02:02 or *03:02: a confident DQB1 without them excludes DQ2/DQ8 whatever DQA1 is
    if not sure and conf.get("HLA-DQB1") == "high" and not any(has(rows, b) for _n, _a, b in CELIAC):
        dq, sure = [], True
    out.append({"kind": "disease association", "allele": "+".join(dq) if dq else ("no DQ2/DQ8" if sure else "DQ2/DQ8"),
                "copies": len(dq), "condition": "celiac disease",
                "status": ("present" if dq else "absent") if sure else "not determined",
                "source": "DQ2.5/DQ2.2/DQ8 (DQA1 and DQB1 alleles, cis or trans) are carried by almost all patients; "
                          "without them celiac disease is very unlikely (negative predictive value > 99 %)"
                          if not dq else "present in about 30-40 % of the population; most carriers never develop "
                                         "celiac disease"})
    return out


def pharmcat_lines(rows):
    out = []
    for r in rows:
        if r["gene"] in ("HLA-A", "HLA-B") and r["confidence"] == "high":
            out.append(f"{r['gene']}\t*{r['allele1'].split('*')[1]}/*{r['allele2'].split('*')[1]}\n")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--bam", required=True)
    ap.add_argument("--hla-dat", help="IPD-IMGT/HLA hla.dat (without it nothing is typed)")
    ap.add_argument("--sample", default="sample")
    ap.add_argument("--tsv", required=True)
    ap.add_argument("--json", required=True)
    ap.add_argument("--pharmcat", required=True, help="PharmCAT outside-call file (HLA-A, HLA-B)")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--min-quality", type=int, default=20)
    ap.add_argument("--keep", help="copy T1K's raw output files into this directory")
    a = ap.parse_args(argv)
    t0 = time.time()
    work = tempfile.mkdtemp(prefix="wes_hla_", dir=os.getcwd())
    t_ref = 0.0
    if a.hla_dat:
        chrom6 = bam_chrom(a.bam)
        seq, coord = build_reference(a.hla_dat, os.path.join(work, "ref"), chrom6)
        t_ref = time.time() - t0
        subprocess.run(["run-t1k", "-b", a.bam, "-f", seq, "-c", coord, "--preset", "hla", "-t", str(a.threads),
                        "-o", "t1k", "--od", work], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        calls = parse_genotype(os.path.join(work, "t1k_genotype.tsv"))
        ars = ars_sequences(seq, [n for c in calls.values() for n in (c["a1"] + "," + c["a2"]).split(",") if n != "."])
    else:
        calls, ars = {}, {}
    rows = genotypes(calls, a.min_quality, ars)
    if not a.hla_dat:
        for r in rows:
            r["note"] = "not typed: no HLA allele database (hla.dat) given"
    assoc = associations(rows)
    cols = ["gene", "allele1", "allele2", "quality1", "quality2", "confidence", "candidates1", "candidates2", "note"]
    with open(a.tsv, "w", newline="\n") as o:
        o.write("sample\t" + "\t".join(cols) + "\n")
        for r in rows:
            o.write(a.sample + "\t" + "\t".join(str(r[c]) for c in cols) + "\n")
        o.write("\nsample\tkind\tallele\tcopies\tcondition\tstatus\tsource\n")
        for x in assoc:
            o.write("\t".join([a.sample, x["kind"], x["allele"], str(x["copies"]), x["condition"], x["status"],
                               x["source"]]) + "\n")
    with open(a.pharmcat, "w", newline="\n") as o:
        o.writelines(pharmcat_lines(rows))
    db = db_release(a.hla_dat) if a.hla_dat else "not given"
    json.dump({"sample": a.sample, "tool": f"wes_hla {VERSION}, T1K (--preset hla)", "database": db or "IPD-IMGT/HLA",
               "min_quality": a.min_quality, "genotypes": rows, "associations": assoc,
               "pharmcat_outside_calls": [x.strip() for x in pharmcat_lines(rows)],
               "seconds": {"reference": round(t_ref, 1), "total": round(time.time() - t0, 1)}},
              open(a.json, "w"), indent=1)
    if a.keep:
        os.makedirs(a.keep, exist_ok=True)
        for f in os.listdir(work):
            if f.startswith("t1k_") and not f.endswith((".fq", ".fa")):
                shutil.copy(os.path.join(work, f), a.keep)
    shutil.rmtree(work, ignore_errors=True)
    typed = sum(r["confidence"] == "high" for r in rows)
    print(f"{a.sample}: {typed} of {len(rows)} classical HLA genes typed with high confidence in "
          f"{time.time() - t0:.0f} s", file=sys.stderr)


if __name__ == "__main__":
    main()
