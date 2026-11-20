"""The report's automated ACMG classification against ClinVar >= 2-star assertions (README section 7).

    python3 tests/acmg_clinvar_eval.py ANNOTATED.vcf clinvar_truth.tsv REF_DIR [--split tune|holdout|all] [--tsv OUT]
        [--mane MANE.summary.txt.gz] [--gene-disease gene_disease.tsv] [--vcep-rules vcep_rules.tsv]

ANNOTATED.vcf is the VCF workflow's annotated_vcf for tests/acmg_clinvar_set.py's pseudo-sample; REF_DIR holds the
reference files install.sh downloads (GenCC, gnomAD constraint, ClinGen dosage, ClinVar residue table, FASTA).
Odd chromosomes (and X) are the tuning half, even chromosomes the held-out half. The classifier never sees CLNSIG
(classification without phenotype terms, so PP4, PM3 and PS2/PM6 are not used). Standard library plus samtools."""
import argparse, collections, os, resource, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tools", "wes_clinical_report"))
import wes_clinical_report as R

ap = argparse.ArgumentParser()
ap.add_argument("vcf"); ap.add_argument("truth"); ap.add_argument("ref")
ap.add_argument("--split", default="tune"); ap.add_argument("--tsv")
ap.add_argument("--mane", help="MANE summary: classify on MANE transcripts (report >= 0.11)")
ap.add_argument("--gene-disease", help="data/build_gene_disease.py table: curated inheritance and mechanism")
ap.add_argument("--missense-policy", choices=("alphamissense", "revel", "concordant"),
                help="how AlphaMissense and REVEL are combined (default: the report's MISSENSE_POLICY)")
ap.add_argument("--vcep-rules", help="data/build_vcep.py rules: ClinGen expert-panel cut-offs, codes and strengths")
a = ap.parse_args()
t0 = time.time()
truth = {}
for l in open(a.truth):
    t = l.rstrip("\n").split("\t")
    if t[0] != "chrom":
        truth[(t[0], int(t[1]), t[2], t[3])] = (t[5], int(t[6]), t[7])
gencc = R.load_gencc(f"{a.ref}/gencc-submissions.tsv")
constraint = R.load_constraint(f"{a.ref}/gnomad.v4.1.constraint_metrics.tsv")
dosage = R.load_dosage(f"{a.ref}/ClinGen_gene_curation_list_GRCh38.tsv")
cvp = R.load_clinvar_protein(f"{a.ref}/clinvar_protein.tsv.gz")
R.MANE.update(R.load_mane(a.mane))
R.CURATED.update(R.load_gene_disease(a.gene_disease)[0])
R.VCEP.update(R.load_vcep_rules(a.vcep_rules)[0])
if a.missense_policy:
    R.MISSENSE_POLICY = a.missense_policy
_s, _c, variants, prefix = R.parse_vcf(a.vcf)


def half(c):
    odd = c == "X" or (c.isdigit() and int(c) % 2 == 1)
    return "tune" if odd else "holdout"


variants = [v for v in variants if (v.chrom, v.pos, v.ref, v.alt) in truth and (a.split == "all" or half(v.chrom) == a.split)]
inframe = [v for v in variants if any(e in (v.anns[0][1] if v.anns else "") for e in R.INFRAME_EFFECTS)]
ctx = R.reference_contexts(f"{a.ref}/GRCh38.ensembl115.primary_assembly.fa", prefix, inframe)
for v in inframe:
    if v.key() in ctx:
        v.repeat = R.repeat_flag(ctx[v.key()][0], ctx[v.key()][1], len(v.ref))


def moi_of(g):
    m = R.gene_moi(gencc.get(g, [])) | R.curated_moi(g)
    if not m and g in R.ACMG_SF:
        m = {"AR"} if R.ACMG_SF[g] in ("biallelic", "hfe") else ({"XL"} if R.ACMG_SF[g] == "xl" else {"AD"})
    return m


def kind(v, g):
    e = (R.ann_for_gene(v, g) or v.anns[0])[1] if v.anns else ""
    if R.is_lof(v, g):
        return "LoF"
    if "missense" in e:
        return "missense"
    if any(x in e for x in R.INFRAME_EFFECTS):
        return "in-frame/stop-lost"
    if "synonymous" in e:
        return "synonymous"
    if "splice_region" in e or "intron" in e:
        return "intronic/splice region"
    return "other (UTR, non-coding)"


def fbin(v):
    p = R.freq(v)
    return "absent" if p is None else ("< 0.01 %" if p < 1e-4 else ("< 1 %" if p < 0.01 else ">= 1 %"))


out = open(a.tsv, "w") if a.tsv else None
rows = []
for v in variants:
    cls_true, stars, tg = truth[(v.chrom, v.pos, v.ref, v.alt)]
    gs = R.genes_of(v)
    g = tg if tg in gs else next((x for x in gs[:3] if x in gencc), v.gene)
    moi = moi_of(g)
    pc, _pn = R.protein_evidence(v, g, cvp)
    ev = R.acmg_evidence(v, g, moi, constraint, False, 0.01, dosage, pc)
    cls, pts = R.acmg_classify(ev)
    auto = R.ACMG_ABBR[cls]
    side = "P/LP" if cls_true in ("P", "LP") else "B/LB"
    rows.append((side, cls_true, stars, auto, pts, ev, kind(v, g), fbin(v), bool(moi)))
    if out:
        out.write("\t".join(map(str, (v.chrom, v.pos, v.ref, v.alt, g, cls_true, stars, auto, pts, ",".join(ev),
                                      kind(v, g), fbin(v)))) + "\n")
secs = time.time() - t0
mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 ** 2
AUTO = ["P", "LP", "VUS", "LB", "B"]


def table(title, sel):
    print(f"\n{title}")
    print(f"  {'ClinVar':6} {'n':>6}  " + "  ".join(f"{x:>6}" for x in AUTO) + "   correct side  wrong side")
    for side in ("P/LP", "B/LB"):
        rs = [r for r in rows if r[0] == side and sel(r)]
        if not rs:
            continue
        c = collections.Counter(r[3] for r in rs)
        right = sum(c[x] for x in (("P", "LP") if side == "P/LP" else ("LB", "B"))) / len(rs)
        wrong = sum(c[x] for x in (("LB", "B") if side == "P/LP" else ("P", "LP"))) / len(rs)
        print(f"  {side:6} {len(rs):6}  " + "  ".join(f"{100 * c[x] / len(rs):5.1f}%" for x in AUTO)
              + f"   {100 * right:9.1f}%  {100 * wrong:8.1f}%")


print(f"split {a.split}: {len(rows)} ClinVar variants after the workflow's target and quality filters; "
      f"{secs:.0f} s, peak {mem:.2f} GiB")
table("All", lambda r: True)
table("Expert panel / practice guideline (3-4 stars)", lambda r: r[2] >= 3)
table("Rare (gnomAD < 1 % or absent)", lambda r: r[7] != ">= 1 %")
for k in ("LoF", "missense", "in-frame/stop-lost", "synonymous", "intronic/splice region", "other (UTR, non-coding)"):
    table(f"Rare, {k}", lambda r, k=k: r[6] == k and r[7] != ">= 1 %")
print("\nEvidence codes: share of ClinVar P/LP and B/LB variants with the code, and the likelihood ratio")
print("(P/LP share / B/LB share; ACMG point system expects ~2.1 Supporting, 4.3 Moderate, 19 Strong, 350 Very strong)")
nP, nB = sum(r[0] == "P/LP" for r in rows), sum(r[0] == "B/LB" for r in rows)
cc = collections.defaultdict(lambda: [0, 0])
for r in rows:
    for c in set(r[5]):
        cc[c][r[0] == "B/LB"] += 1
for c, (p, b) in sorted(cc.items(), key=lambda x: -sum(x[1])):
    lr = (p / nP) / (b / nB) if b else float("inf")
    print(f"  {c:16} P/LP {100 * p / nP:5.1f}%  B/LB {100 * b / nB:5.1f}%  LR {lr:8.2f}")
