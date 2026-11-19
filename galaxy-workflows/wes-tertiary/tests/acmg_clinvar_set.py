"""ClinVar test set for the report's automated ACMG classification (README section 7): >= 2-star P/LP and B/LB
SNVs/indels (<= 50 bp) on 1-22 and X, 15,000 per side with up to 5,000 expert-panel/practice-guideline (3-4 star)
variants each, as one heterozygous pseudo-sample VCF (clinvar_set.vcf) plus the truth table (clinvar_truth.tsv).

    python3 tests/acmg_clinvar_set.py clinvar.GRCh38.vcf.gz GRCh38.ensembl115.primary_assembly.fa.fai OUTDIR

Run clinvar_set.vcf through the VCF workflow, then tests/acmg_clinvar_eval.py on its annotated VCF."""
import gzip, random, sys, collections
cv, fai, out = sys.argv[1:4]
TWO = {"criteria_provided,_multiple_submitters,_no_conflicts": 2, "reviewed_by_expert_panel": 3, "practice_guideline": 4}
CLS = {"Pathogenic": "P", "Likely_pathogenic": "LP", "Pathogenic/Likely_pathogenic": "P",
       "Benign": "B", "Likely_benign": "LB", "Benign/Likely_benign": "B"}
CHR = {str(i) for i in range(1, 23)} | {"X"}
pool = collections.defaultdict(list)
for line in gzip.open(cv, "rt"):
    if line[0] == "#":
        continue
    t = line.split("\t", 8)
    if t[0] not in CHR or len(t[3]) > 51 or len(t[4]) > 51 or t[4] in (".", "") or "," in t[4]:
        continue
    info = dict(kv.split("=", 1) for kv in t[7].split(";") if "=" in kv)
    st, sig = TWO.get(info.get("CLNREVSTAT")), CLS.get(info.get("CLNSIG"))
    if not st or not sig or "GENEINFO" not in info:
        continue
    side = "P" if sig in ("P", "LP") else "B"
    pool[(side, st >= 3)].append((t[0], int(t[1]), t[2], t[3], t[4], sig, st, info["GENEINFO"].split(":")[0]))
random.seed(1)
keep = []
for side in "PB":
    ex = pool[(side, True)]
    keep += random.sample(ex, min(5000, len(ex)))
    rest = pool[(side, False)]
    keep += random.sample(rest, min(15000 - min(5000, len(ex)), len(rest)))
    print(side, "pool: expert", len(ex), "multi-submitter", len(rest))
order = {c: i for i, c in enumerate([l.split("\t")[0] for l in open(fai)])}
keep.sort(key=lambda x: (order[x[0]], x[1]))
with open(f"{out}/clinvar_set.vcf", "w") as f, open(f"{out}/clinvar_truth.tsv", "w") as g:
    f.write("##fileformat=VCFv4.2\n")
    for l in open(fai):
        c, n = l.split("\t")[:2]
        f.write(f"##contig=<ID={c},length={n}>\n")
    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
            '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">\n'
            '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCLINVAR_SET\n")
    g.write("chrom\tpos\tref\talt\tclinvar_id\tclass\tstars\tgene\n")
    seen = set()
    for c, p, i, r, a, sig, st, gene in keep:
        if (c, p, r, a) in seen:
            continue
        seen.add((c, p, r, a))
        f.write(f"{c}\t{p}\t{i}\t{r}\t{a}\t100\tPASS\t.\tGT:DP:AD:GQ\t0/1:40:20,20:99\n")
        g.write(f"{c}\t{p}\t{r}\t{a}\t{i}\t{sig}\t{st}\t{gene}\n")
print("written", len(seen))
