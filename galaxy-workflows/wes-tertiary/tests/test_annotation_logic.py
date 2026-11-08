#!/usr/bin/env python3
"""Functional test of the variant filtering/annotation logic of "WES tertiary analysis from VCF (v2)".

The bcftools and SnpSift command lines are rebuilt from the parameters stored in wes-tertiary-from-vcf.ga
(so the test exercises exactly what is shipped) and run on a small synthetic data set:

  bcftools norm     --check-ref <wf> --fasta-ref ref.fa --multiallelics <wf>
  bcftools filter   --regions-file padded.bed --regions-overlap <wf> --include <wf default quality filter>
  bcftools annotate --columns <wf ClinVar columns>    --annotations clinvar.vcf.gz
  bcftools annotate --columns <wf gnomAD columns>     --annotations gnomad.vcf.gz
  SnpSift filter    -e <wf default prioritization expression>       (on SnpEff-style ANN records)
  SnpSift extractFields <wf field list>

Requirements: pysam (bundles bcftools) and, for the SnpSift part, Java >= 21 plus SnpSift.jar
(e.g. from the bioconda "snpsift" package):

    python test_annotation_logic.py --snpsift-jar /path/to/SnpSift.jar
"""
import argparse
import json
import os
import random
import subprocess
import sys
import tempfile

import pysam
from pysam import bcftools

HERE = os.path.dirname(os.path.abspath(__file__))
WORKFLOW = os.path.join(HERE, "..", "wes-tertiary-from-vcf.ga")



def load_state(tool_state):
    """tool_state as a nested dict, whether plain JSON (gxformat2) or Galaxy's own encoding.

    Galaxy (and build.py, which stores Galaxy's complete states) JSON-encodes every top-level value
    a second time and always adds __rerun_remap_job_id__; gxformat2 writes one plain JSON object.
    """
    state = json.loads(tool_state) if isinstance(tool_state, str) else tool_state
    if "__rerun_remap_job_id__" in state:
        state = {k: (json.loads(v) if isinstance(v, str) and not k.startswith("__") else v)
                 for k, v in state.items()}
    return state

def workflow_settings():
    wf = json.load(open(WORKFLOW))
    steps = {s["label"]: s for s in wf["steps"].values()}
    state = {k: load_state(s["tool_state"]) for k, s in steps.items()}
    return {
        "check_ref": state["normalize"]["check_ref"],
        "multiallelics": state["normalize"]["multiallelics"]["mode"] + state["normalize"]["multiallelics"]["multiallelic_types"],
        "regions_overlap": state["target_filter"]["sec_restrict"]["regions"]["regions_overlap"],
        "quality": state["Quality filter (bcftools expression)"]["default"],
        "clinvar_columns": state["clinvar_annotation"]["sec_annofile"]["columns"],
        "gnomad_columns": state["population_annotation"]["sec_annofile"]["columns"],
        "prioritization": state["Prioritization filter (SnpSift expression)"]["default"],
        "fields": state["candidate_table"]["extract"],
    }


def alt_base(b):
    return {"A": "C", "C": "G", "G": "T", "T": "A"}[b]


def test_bcftools(cfg, work):
    random.seed(7)
    seqs = {"chr1": "".join(random.choice("ACGT") for _ in range(20000)),
            "chr2": "".join(random.choice("ACGT") for _ in range(8000))}
    with open("ref.fa", "w") as f:
        for name, s in seqs.items():
            f.write(f">{name}\n" + "\n".join(s[i:i + 60] for i in range(0, len(s), 60)) + "\n")
    pysam.faidx("ref.fa")
    ref = lambda c, p: seqs[c][p - 1]

    # targets padded by 50 bp, as bedtools slop -b 50 -g ref.fa.fai does
    targets = [("chr1", 1000, 1200), ("chr1", 5000, 5300), ("chr2", 100, 400)]
    with open("padded.bed", "w") as f:
        for c, s, e in targets:
            f.write(f"{c}\t{max(0, s - 50)}\t{min(len(seqs[c]), e + 50)}\n")

    # (chrom, pos, ref, alt, qual, filter, gt, dp, ad, keep, reason)
    recs = []
    def add(c, p, qual, flt, gt, dp, ad, keep, why, alt=None, ref_override=None):
        r = ref_override or ref(c, p)
        recs.append((c, p, r, alt or alt_base(ref(c, p)), qual, flt, gt, dp, ad, keep, why))
    add("chr1", 1050, 200, ".", "0/1", 40, "20,20", True, "on target, ClinVar P + common in gnomAD")
    add("chr1", 1100, 150, ".", "0/1", 35, "18,17", True, "on target; ClinVar/gnomAD only have another ALT here")
    add("chr1", 1150, 120, ".", "1/2", 30, "1,15,14", True, "multiallelic, split into two records",
        alt=alt_base(ref("chr1", 1150)) + "," + alt_base(alt_base(ref("chr1", 1150))))
    add("chr1", 1180, 10, ".", "0/1", 30, "15,15", False, "QUAL 10 < 20")
    add("chr1", 1190, 90, "LowQual", "0/1", 30, "15,15", False, "FILTER LowQual")
    add("chr1", 5010, 90, "PASS", "1/1", 6, "0,6", False, "FORMAT/DP 6 < 10")
    add("chr1", 5320, 90, "PASS", "0/1", 30, "15,15", True, "20 bp outside target, inside padding")
    add("chr1", 5500, 90, ".", "0/1", 30, "15,15", False, "200 bp outside target")
    add("chr2", 200, 90, ".", "0/1", 30, "15,15", False, "REF does not match reference",
        alt=ref("chr2", 200), ref_override=alt_base(ref("chr2", 200)))
    add("chr2", 250, 90, ".", "0/1", 30, "15,15", True, "on target, gnomAD AF 0.001")

    header = ["##fileformat=VCFv4.2"] + [f"##contig=<ID={n},length={len(s)}>" for n, s in seqs.items()] + [
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">',
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Caller allele frequency">',
        '##FILTER=<ID=LowQual,Description="Low quality">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tproband"]
    with open("input.vcf", "w") as f:
        f.write("\n".join(header) + "\n")
        for c, p, r, a, q, flt, gt, dp, ad, _k, _w in sorted(recs, key=lambda x: (x[0], x[1])):
            af = ",".join(["0.5"] * len(a.split(",")))
            f.write(f"{c}\t{p}\t.\t{r}\t{a}\t{q}\t{flt}\tDP={dp};AF={af}\tGT:DP:AD\t{gt}:{dp}:{ad}\n")

    r1050, r1100, r250 = ref("chr1", 1050), ref("chr1", 1100), ref("chr2", 250)
    with open("clinvar.vcf", "w") as f:
        f.write("##fileformat=VCFv4.1\n##contig=<ID=chr1>\n##contig=<ID=chr2>\n")
        for tag, typ, num in [("ALLELEID", "Integer", "1"), ("CLNDN", "String", "."), ("CLNREVSTAT", "String", "."),
                              ("CLNSIG", "String", "."), ("CLNSIGCONF", "String", ".")]:
            f.write(f'##INFO=<ID={tag},Number={num},Type={typ},Description="{tag}">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        f.write(f"chr1\t1050\t1\t{r1050}\t{alt_base(r1050)}\t.\t.\tALLELEID=1;CLNDN=Disease_A;"
                "CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts;CLNSIG=Pathogenic\n")
        f.write(f"chr1\t1100\t2\t{r1100}\t{alt_base(alt_base(r1100))}\t.\t.\tALLELEID=2;CLNDN=Disease_B;"
                "CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Likely_pathogenic\n")
    with open("gnomad.vcf", "w") as f:
        f.write("##fileformat=VCFv4.2\n##contig=<ID=chr1>\n##contig=<ID=chr2>\n")
        f.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Alternate allele frequency">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        f.write(f"chr1\t1050\t.\t{r1050}\t{alt_base(r1050)}\t.\tPASS\tAF=0.2\n")
        f.write(f"chr1\t1100\t.\t{r1100}\t{alt_base(alt_base(r1100))}\t.\tPASS\tAF=0.3\n")
        f.write(f"chr2\t250\t.\t{r250}\t{alt_base(r250)}\t.\tPASS\tAF=0.001\n")
    for name in ("input.vcf", "clinvar.vcf", "gnomad.vcf"):  # the wrappers bgzip + index their inputs
        pysam.tabix_compress(name, name + ".gz", force=True)
        bcftools.index(name + ".gz", catch_stdout=False)

    def run(fn, *args, out):
        fn(*args, "--output-type", "z", "-o", out, catch_stdout=False)
        bcftools.index(out, catch_stdout=False)

    run(bcftools.norm, "--check-ref", cfg["check_ref"], "--fasta-ref", "ref.fa",
        "--multiallelics", cfg["multiallelics"], "input.vcf.gz", out="norm.vcf.gz")
    run(bcftools.filter, "--regions-file", "padded.bed", "--regions-overlap", cfg["regions_overlap"],
        "--include", cfg["quality"], "norm.vcf.gz", out="filtered.vcf.gz")
    run(bcftools.annotate, "--columns", cfg["clinvar_columns"], "--annotations", "clinvar.vcf.gz",
        "filtered.vcf.gz", out="clinvar.annotated.vcf.gz")
    run(bcftools.annotate, "--columns", cfg["gnomad_columns"], "--annotations", "gnomad.vcf.gz",
        "clinvar.annotated.vcf.gz", out="annotated.vcf.gz")

    body = [l.split("\t") for l in pysam.BGZFile("annotated.vcf.gz").read().decode().splitlines()
            if l and not l.startswith("#")]
    kept = {(f[0], int(f[1]), f[4]): f[7] for f in body}
    failures = []
    for c, p, r, a, *_rest, keep, why in recs:
        got = any(k[0] == c and k[1] == p for k in kept)
        status = "ok " if got == keep else "BAD"
        if got != keep:
            failures.append(f"{c}:{p}")
        print(f"  [{status}] {c}:{p} {r}>{a:<5} expected {'KEEP' if keep else 'drop'}, got {'KEEP' if got else 'drop'}  ({why})")
    checks = [
        (sum(1 for k in kept if k[:2] == ("chr1", 1150)) == 2, "multiallelic record split into 2 biallelic records"),
        ("CLNSIG=Pathogenic" in kept.get(("chr1", 1050, alt_base(r1050)), ""), "ClinVar CLNSIG transferred (allele match)"),
        ("gnomAD_AF=0.2" in kept.get(("chr1", 1050, alt_base(r1050)), ""), "gnomAD AF transferred as gnomAD_AF"),
        ("AF=0.5" in kept.get(("chr1", 1050, alt_base(r1050)), ""), "caller INFO/AF left untouched"),
        ("CLNSIG" not in kept.get(("chr1", 1100, alt_base(r1100)), "CLNSIG"), "no ClinVar transfer from a different ALT"),
        ("gnomAD_AF" not in kept.get(("chr1", 1100, alt_base(r1100)), "gnomAD_AF"), "no gnomAD transfer from a different ALT"),
    ]
    for ok, what in checks:
        print(f"  [{'ok ' if ok else 'BAD'}] {what}")
        if not ok:
            failures.append(what)
    return failures


SNPEFF_HEADER = [
    '##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: \'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO\' ">',
    '##INFO=<ID=LOF,Number=.,Type=String,Description="Predicted loss of function effects for this variant. Format: \'Gene_Name | Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected\'">',
    '##INFO=<ID=NMD,Number=.,Type=String,Description="Predicted nonsense mediated decays effects for this variant. Format: \'Gene_Name | Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected\'">',
]


def ann(alt, effect, impact, gene, tx, hgvs_c, hgvs_p=""):
    return f"{alt}|{effect}|{impact}|{gene}|{gene}|transcript|{tx}|protein_coding|2/5|{hgvs_c}|{hgvs_p}|100/900|90/600|30/199||"


def test_snpsift(cfg, jar):
    rows = [  # chrom, pos, ref, alt, INFO, expected keep, reason
        ("chr1", 1050, "A", "C", "CLNSIG=Pathogenic;gnomAD_AF=0.2;ANN=" + ann("C", "missense_variant", "MODERATE", "GENEA", "TX1", "c.90A>C", "p.Lys30Asn"), True, "ClinVar pathogenic, common"),
        ("chr1", 1100, "G", "T", "LOF=(GENEA|GENEA|1|1.00);NMD=(GENEA|GENEA|1|1.00);ANN=" + ann("T", "stop_gained", "HIGH", "GENEA", "TX1", "c.140G>T", "p.Glu47*"), True, "HIGH impact, absent from gnomAD"),
        ("chr1", 1150, "T", "A", "ANN=" + ann("A", "synonymous_variant", "LOW", "GENEA", "TX1", "c.190T>A", "p.Leu63Leu"), False, "LOW impact, rare"),
        ("chr1", 1150, "T", "C", "ANN=" + ann("C", "missense_variant", "MODERATE", "GENEA", "TX1", "c.190T>C", "p.Leu63Pro") + "," + ann("C", "intron_variant", "MODIFIER", "GENEA", "TX2", "c.50+12T>C"), True, "MODERATE on one of two transcripts, rare"),
        ("chr1", 5320, "C", "G", "ANN=" + ann("G", "splice_region_variant&intron_variant", "LOW", "GENEB", "TX3", "c.300+5C>G"), False, "LOW impact splice region"),
        ("chr2", 250, "T", "A", "gnomAD_AF=0.001;ANN=" + ann("A", "missense_variant", "MODERATE", "GENEC", "TX4", "c.10T>A", "p.Ser4Thr"), True, "MODERATE, gnomAD 0.1%"),
        ("chr2", 300, "C", "T", "gnomAD_AF=0.05;ANN=" + ann("T", "missense_variant", "MODERATE", "GENEC", "TX4", "c.60C>T", "p.Ala20Val"), False, "MODERATE but gnomAD 5%"),
        ("chr2", 350, "G", "A", "CLNSIG=Benign;gnomAD_AF=0.3;ANN=" + ann("A", "missense_variant", "MODERATE", "GENEC", "TX4", "c.110G>A", "p.Arg37His"), False, "ClinVar benign, common"),
        ("chr2", 360, "A", "G", "CLNSIG=Conflicting_classifications_of_pathogenicity;gnomAD_AF=0.02;ANN=" + ann("G", "missense_variant", "MODERATE", "GENEC", "TX4", "c.120A>G", "p.Lys40Arg"), True, "ClinVar conflicting, 2%"),
        ("chr2", 370, "C", "A", "CLNSIG=Pathogenic/Likely_pathogenic;gnomAD_AF=0.04;ANN=" + ann("A", "synonymous_variant", "LOW", "GENEC", "TX4", "c.130C>A", "p.Ala44Ala"), True, "ClinVar P/LP, LOW impact, 4%"),
        ("chr2", 380, "T", "C", "gnomAD_AF=0.0099;ANN=" + ann("C", "frameshift_variant", "HIGH", "GENEC", "TX4", "c.140del", "p.Leu47fs"), True, "HIGH, gnomAD 0.99%"),
        ("chr2", 390, "T", "G", "gnomAD_AF=0.01;ANN=" + ann("G", "stop_gained", "HIGH", "GENEC", "TX4", "c.150T>G", "p.Tyr50*"), False, "HIGH, gnomAD exactly 1%"),
    ]
    header = ["##fileformat=VCFv4.2", "##contig=<ID=chr1>", "##contig=<ID=chr2>",
              '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">',
              '##INFO=<ID=AF,Number=A,Type=Float,Description="Caller allele frequency">',
              '##INFO=<ID=gnomAD_AF,Number=A,Type=Float,Description="gnomAD allele frequency">'] + [
              f'##INFO=<ID={t},Number=.,Type=String,Description="{t}">' for t in ("CLNSIG", "CLNREVSTAT", "CLNDN", "CLNSIGCONF")] + [
              '##INFO=<ID=ALLELEID,Number=1,Type=Integer,Description="ALLELEID">'] + SNPEFF_HEADER + [
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
              '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
              '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
              "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tproband"]
    with open("snpeff_annotated.vcf", "w") as f:
        f.write("\n".join(header) + "\n")
        for c, p, r, a, info, _k, _w in rows:
            f.write(f"{c}\t{p}\t.\t{r}\t{a}\t100\tPASS\tDP=30;AF=0.5;{info}\tGT:DP:AD\t0/1:30:15,15\n")
    # The Galaxy wrapper writes the expression to a file: SnpSift filter -f input -e exprFile
    open("expression.txt", "w").write(cfg["prioritization"])
    res = subprocess.run(["java", "-Xmx1g", "-jar", jar, "filter", "-f", "snpeff_annotated.vcf", "-e", "expression.txt"],
                         capture_output=True, text=True)
    if res.returncode:
        return ["SnpSift filter failed: " + res.stderr[-500:]]
    open("candidates.vcf", "w").write(res.stdout)
    kept = {(l.split("\t")[0], int(l.split("\t")[1]), l.split("\t")[4])
            for l in res.stdout.splitlines() if l and not l.startswith("#")}
    failures = []
    for c, p, r, a, _info, keep, why in rows:
        got = (c, p, a) in kept
        if got != keep:
            failures.append(f"{c}:{p}{r}>{a}")
        print(f"  [{'ok ' if got == keep else 'BAD'}] {c}:{p} {r}>{a}  expected {'KEEP' if keep else 'drop'}, "
              f"got {'KEEP' if got else 'drop'}  ({why})")
    # the wrapper runs: SnpSift extractFields -s ',' -e '.' - <fields>
    res = subprocess.run(["java", "-Xmx1g", "-jar", jar, "extractFields", "-s", ",", "-e", ".", "candidates.vcf"]
                         + cfg["fields"].split(), capture_output=True, text=True)
    if res.returncode:
        return failures + ["SnpSift extractFields failed: " + res.stderr[-500:]]
    lines = res.stdout.splitlines()
    print("  candidate table:")
    for line in lines:
        print("    " + line)
    if {len(l.split("\t")) for l in lines} != {len(cfg["fields"].split())}:
        failures.append("candidate table is ragged")
    if len(lines) - 1 != sum(1 for r in rows if r[5]):
        failures.append("candidate table row count")
    return failures


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--snpsift-jar", default=os.environ.get("SNPSIFT_JAR"), help="path to SnpSift.jar")
    args = parser.parse_args()
    cfg = workflow_settings()
    print("Settings read from", os.path.relpath(WORKFLOW))
    for k, v in cfg.items():
        print(f"  {k:16s} {v}")
    failures = []
    with tempfile.TemporaryDirectory() as work:
        os.chdir(work)
        print("\nbcftools: normalise -> restrict to padded targets + quality filter -> ClinVar -> gnomAD")
        failures += test_bcftools(cfg, work)
        if args.snpsift_jar:
            print("\nSnpSift: prioritisation expression and candidate table")
            failures += test_snpsift(cfg, args.snpsift_jar)
        else:
            print("\nSnpSift part skipped (no --snpsift-jar / SNPSIFT_JAR)")
    print("\nRESULT:", "PASS" if not failures else f"FAIL {failures}")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
