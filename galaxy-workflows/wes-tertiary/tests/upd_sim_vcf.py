"""Impose a maternal or paternal isodisomy on one chromosome of a proband VCF from a parent's calls (README
section 7, end-to-end UPD check): python3 tests/upd_sim_vcf.py child.vcf parent.vcf CHROM out.vcf. At each proband record of CHROM the genotype becomes homozygous for a
random allele of the parent (0/0 records are dropped); other chromosomes are unchanged."""
import importlib.util
import os
import random
import sys

spec = importlib.util.spec_from_file_location("r", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tools",
                                                                "wes_clinical_report", "wes_clinical_report.py"))
R = importlib.util.module_from_spec(spec)
spec.loader.exec_module(R)
child, parent, chrom, out = sys.argv[1:5]
par = R.load_parent_vcf(parent)
rng = random.Random(11)
kept = dropped = 0
with R.opener(child) as f, open(out, "w") as o:
    for line in f:
        if line[0] == "#":
            o.write(line)
            continue
        t = line.rstrip("\n").split("\t")
        if R.norm_chrom(t[0]) != chrom or "," in t[4]:
            o.write(line)
            continue
        call = par.get((chrom, int(t[1]), t[3], t[4]))
        m = 2 if call and call[0] == "hom" else (1 if call else 0)
        if not rng.random() < m / 2:
            dropped += 1
            continue
        fmt = t[8].split(":")
        val = t[9].split(":")
        d = dict(zip(fmt, val))
        dp = int(d["DP"]) if d.get("DP", ".").isdigit() else 40
        d["GT"] = "1/1"
        if "AD" in d:
            d["AD"] = f"0,{dp}"
        t[9] = ":".join(d.get(k, ".") for k in fmt)
        o.write("\t".join(t) + "\n")
        kept += 1
print(f"chr{chrom}: {kept} records made homozygous, {dropped} dropped")
