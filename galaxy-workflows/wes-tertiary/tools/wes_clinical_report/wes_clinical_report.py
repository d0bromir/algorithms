#!/usr/bin/env python3
"""Human-readable clinical interpretation report for one germline exome.

Input is the annotated VCF of the "WES tertiary analysis from VCF" workflow (SnpEff ANN/LOF/NMD,
ClinVar CLNSIG/CLNREVSTAT/CLNDN/..., gnomAD v4.1 frequencies INFO/gnomAD_AF, _FAF95, _grpmax_AF, _nhomalt, ...),
plus optional knowledge
bases and QC outputs. Every statement in the report is derived by explicit, versioned rules from those
inputs; nothing is inferred by a statistical or language model. The report is decision support for a
qualified reviewer (clinical geneticist / laboratory scientist), not a diagnosis.

Sections: sample QC and sex check; diagnostic findings (ClinVar P/LP with inheritance-consistent
genotype, compound heterozygotes phased from the reads where possible); candidate variants not classified in
ClinVar (rare predicted loss-of-function, rare missense at pathogenic ClinVar residues, with high AlphaMissense
scores or in phenotype-matched genes, predicted splice-altering) with automated ACMG/AMP evidence codes and the
diseases that best explain the patient's phenotype; ACMG SF v3.3 secondary findings; carrier status;
pharmacogenomics (CPIC genes, single-nucleotide star alleles, PharmCAT); risk alleles and other ClinVar
associations; runs of homozygosity; coverage of clinically relevant genes; limits. Gene cells carry GenCC,
ClinGen dosage sensitivity and HPO disease features (onset, course, inheritance).

Standard library only (Python >= 3.8).
"""
import argparse
import bisect
import collections
import datetime
import gzip
import html
import json
import math
import re
import statistics
import subprocess
import sys

VERSION = "0.9.0"

# --------------------------------------------------------------------------------------------- rules
# ACMG secondary findings v3.3 (ACMG policy statement, Genet Med 2025;27:101454): the 81 genes of v3.2 (Miller et al.,
# Genet Med 2023;25:100866) plus ABCD1, CYP27A1 and PLN. Value: reporting rule.
#   "het" = one P/LP variant is reportable; "biallelic" = two P/LP alleles (hom or two het);
#   "xl" = hemizygous or homozygous, or two heterozygous P/LP variants (ABCD1; a single heterozygous variant in a
#   female is not reported); "hfe" = HFE p.Cys282Tyr homozygotes only; "ttn" = truncating variants only (TTN).
ACMG_SF = {g: "het" for g in (
    "APC RET BRCA1 BRCA2 PALB2 SDHD SDHAF2 SDHC SDHB MAX TMEM127 BMPR1A SMAD4 TP53 MLH1 MSH2 MSH6 PMS2 "
    "MEN1 NF2 STK11 PTEN RB1 TSC1 TSC2 VHL WT1 FBN1 TGFBR1 TGFBR2 SMAD3 ACTA2 MYH11 PKP2 DSP DSC2 TMEM43 "
    "DSG2 RYR2 CASQ2 TRDN CALM1 CALM2 CALM3 TNNT2 LMNA FLNC BAG3 DES RBM20 TNNC1 KCNQ1 KCNH2 SCN5A MYBPC3 "
    "MYH7 TNNI3 TPM1 MYL3 ACTC1 PRKAG2 MYL2 GLA COL3A1 LDLR APOB PCSK9 OTC RYR1 CACNA1S HNF1A TTR "
    "ACVRL1 ENG PLN").split()}
ACMG_SF.update({"MUTYH": "biallelic", "ATP7B": "biallelic", "BTD": "biallelic", "GAA": "biallelic",
                "RPE65": "biallelic", "CYP27A1": "biallelic", "ABCD1": "xl", "HFE": "hfe", "TTN": "ttn"})
ACMG_SF_VERSION = "ACMG SF v3.3"
ACMG_SF_CATEGORY = {}
for _cat, _genes in (
        ("Cancer predisposition", "APC RET BRCA1 BRCA2 PALB2 SDHD SDHAF2 SDHC SDHB MAX TMEM127 BMPR1A SMAD4 TP53 "
                                  "MLH1 MSH2 MSH6 PMS2 MEN1 MUTYH NF2 STK11 PTEN RB1 TSC1 TSC2 VHL WT1"),
        ("Cardiovascular", "FBN1 TGFBR1 TGFBR2 SMAD3 ACTA2 MYH11 PKP2 DSP DSC2 TMEM43 DSG2 RYR2 CASQ2 TRDN CALM1 "
                           "CALM2 CALM3 TNNT2 LMNA FLNC TTN BAG3 DES RBM20 TNNC1 KCNQ1 KCNH2 SCN5A MYBPC3 MYH7 "
                           "TNNI3 TPM1 MYL3 ACTC1 PRKAG2 MYL2 GLA COL3A1 LDLR APOB PCSK9 TTR PLN"),
        ("Inborn errors of metabolism", "OTC GAA BTD CYP27A1 ABCD1"),
        ("Miscellaneous", "ACVRL1 ENG ATP7B HFE RPE65 RYR1 CACNA1S HNF1A")):
    for _g in _genes.split():
        ACMG_SF_CATEGORY[_g] = _cat
assert len(ACMG_SF) == 84, len(ACMG_SF)

LOF_EFFECTS = ("stop_gained", "frameshift_variant", "splice_acceptor_variant", "splice_donor_variant",
               "start_lost", "exon_loss_variant", "transcript_ablation")
MISSENSE_EFFECTS = ("missense_variant", "inframe_insertion", "inframe_deletion", "disruptive_inframe_insertion",
                    "disruptive_inframe_deletion", "stop_lost")
REVIEW_STARS = {
    "practice_guideline": 4, "reviewed_by_expert_panel": 3,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2, "criteria_provided,_multiple_submitters": 2,
    "criteria_provided,_conflicting_classifications": 1, "criteria_provided,_conflicting_interpretations": 1,
    "criteria_provided,_single_submitter": 1,
}
GENCC_VALID = ("Definitive", "Strong", "Moderate", "Supportive")
GENCC_RANK = {"Definitive": 4, "Strong": 3, "Moderate": 2, "Supportive": 1}

# Pharmacogenomic positions (GRCh38; forward-strand alleles, checked against ClinVar and the reference).
# allele: (gene, star/name, rsID, chrom, pos, ref, alt, function)
PGX_ALLELES = [
    ("CYP2C19", "*2", "rs4244285", "10", 94781859, "G", "A", "no function"),
    ("CYP2C19", "*3", "rs4986893", "10", 94780653, "G", "A", "no function"),
    ("CYP2C19", "*17", "rs12248560", "10", 94761900, "C", "T", "increased function"),
    ("CYP2C9", "*2", "rs1799853", "10", 94942290, "C", "T", "decreased function"),
    ("CYP2C9", "*3", "rs1057910", "10", 94981296, "A", "C", "no function"),
    ("VKORC1", "c.-1639G>A", "rs9923231", "16", 31096368, "C", "T", "increased warfarin sensitivity"),
    ("CYP4F2", "*3", "rs2108622", "19", 15879621, "C", "T", "decreased function"),
    ("SLCO1B1", "c.521T>C", "rs4149056", "12", 21178615, "T", "C", "decreased function"),
    ("DPYD", "*2A (c.1905+1G>A)", "rs3918290", "1", 97450058, "C", "T", "no function"),
    ("DPYD", "*13 (c.1679T>G)", "rs55886062", "1", 97515787, "A", "C", "no function"),
    ("DPYD", "c.2846A>T", "rs67376798", "1", 97082391, "T", "A", "decreased function"),
    ("DPYD", "HapB3 (c.1236G>A)", "rs56038477", "1", 97573863, "C", "T", "decreased function"),
    ("TPMT", "*2 (c.238G>C)", "rs1800462", "6", 18143724, "C", "G", "no function"),
    ("TPMT", "*3B (c.460G>A)", "rs1800460", "6", 18138997, "C", "T", "no function"),
    ("TPMT", "*3C (c.719A>G)", "rs1142345", "6", 18130687, "T", "C", "no function"),
    ("NUDT15", "*3 (c.415C>T)", "rs116855232", "13", 48045719, "C", "T", "no function"),
    ("CYP3A5", "*3", "rs776746", "7", 99672916, "T", "C", "no function"),
    ("CYP2B6", "c.516G>T (*6/*9)", "rs3745274", "19", 41006936, "G", "T", "decreased function"),
    ("UGT1A1", "*6 (c.211G>A)", "rs4148323", "2", 233760498, "G", "A", "decreased function"),
    ("UGT1A1", "*80 (c.-364C>T, *28 proxy)", "rs887829", "2", 233759924, "C", "T", "decreased function"),
]
RISK_ALLELES = [
    ("F5", "Factor V Leiden (p.Arg534Gln)", "rs6025", "1", 169549811, "C", "T",
     "Venous thromboembolism risk (about 3-8 fold heterozygous, much higher homozygous); relevant for "
     "oestrogen-containing contraception, pregnancy, surgery."),
    ("F2", "Prothrombin c.*97G>A (G20210A)", "rs1799963", "11", 46739505, "G", "A",
     "Venous thromboembolism risk (about 2-3 fold heterozygous)."),
    ("HFE", "p.Cys282Tyr", "rs1800562", "6", 26092913, "G", "A",
     "Hereditary haemochromatosis: homozygotes are at risk of iron overload (incomplete penetrance); "
     "heterozygotes are carriers."),
    ("HFE", "p.His63Asp", "rs1799945", "6", 26090951, "C", "G",
     "Low-penetrance haemochromatosis allele; clinically relevant mainly as compound heterozygote with p.Cys282Tyr."),
    ("APOE", "rs429358 (APOE e4 tag)", "rs429358", "19", 44908684, "T", "C", None),
    ("APOE", "rs7412 (APOE e2 tag)", "rs7412", "19", 44908822, "C", "T", None),
]
PGX_DRUGS = {
    "CYP2C19": "clopidogrel, voriconazole, PPIs (omeprazole, pantoprazole...), SSRIs (citalopram, escitalopram, "
               "sertraline), tricyclic antidepressants",
    "CYP2C9": "warfarin (with VKORC1/CYP4F2), phenytoin/fosphenytoin, NSAIDs (celecoxib, ibuprofen, meloxicam...), "
              "siponimod",
    "VKORC1": "warfarin",
    "CYP4F2": "warfarin",
    "SLCO1B1": "statins (simvastatin, atorvastatin, rosuvastatin...)",
    "DPYD": "fluoropyrimidines (5-fluorouracil, capecitabine, tegafur)",
    "TPMT": "thiopurines (azathioprine, mercaptopurine, thioguanine)",
    "NUDT15": "thiopurines (azathioprine, mercaptopurine, thioguanine)",
    "CYP3A5": "tacrolimus",
    "CYP2B6": "efavirenz",
    "UGT1A1": "atazanavir, irinotecan",
}
NOT_ASSESSED_PGX = ("CYP2D6 (gene conversions/duplications with CYP2D7; codeine, tramadol, tamoxifen, "
                    "ondansetron...), HLA-A/HLA-B alleles (abacavir, carbamazepine, allopurinol...), "
                    "MT-RNR1 (aminoglycosides; mitochondrial genome), UGT1A1*28 TA-repeat (only the *80 proxy), "
                    "G6PD haplotypes beyond single ClinVar variants, IFNL3, NAT2")

PAR_GRCH38 = [(10001, 2781479), (155701383, 156030895)]

# Hypervariable genes, and clinically relevant genes with highly similar paralogs or pseudogenes where short-read
# alignment is unreliable (GIAB CMRG / "dark" medically relevant genes; Wagner et al. 2022, Ebbert et al. 2019).
POLYMORPHIC_GENES = {"MUC4", "MUC16", "MUC3A", "MUC5B", "MUC6", "MUC12", "MUC19", "MUC17", "OR4F5"}
PARALOG_GENES = {"PRSS1", "PRSS2", "PMS2", "SMN1", "SMN2", "CYP21A2", "TNXB", "GBA1", "GBA", "NEB", "STRC", "OTOA",
                 "HBA1", "HBA2", "IKBKG", "CHEK2", "PKD1", "NCF1", "CFC1", "SBDS", "HYDIN", "VWF", "TTN-AS1",
                 "CYP2D6", "CYP2A6", "RHD", "RHCE", "C4A", "C4B", "NPHP1", "FCGR2A", "FCGR3A", "FCGR3B", "LPA",
                 "KCNJ18", "MAP2K3", "OPN1LW", "OPN1MW", "CORO1A", "SPATA31A3", "ACTG1P"}


# ----------------------------------------------------------------------------------------- helpers
def opener(path):
    with open(path, "rb") as f:
        magic = f.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", errors="replace")
    return open(path, "rt", errors="replace")


def esc(x):
    return html.escape("" if x is None else str(x))


def norm_chrom(c):
    c = c[3:] if c.lower().startswith("chr") else c
    return "MT" if c == "M" else c


def fnum(x):
    try:
        v = float(x)
        return None if math.isnan(v) else v
    except (TypeError, ValueError):
        return None


def fmt_af(af):
    if af is None:
        return "absent"
    if af == 0:
        return "0"
    return f"{af:.2g}" if af < 0.001 else f"{af:.4f}".rstrip("0")


def freq(v):
    """Population frequency used by the frequency rules: the higher of the overall gnomAD AF and the highest
    filtering allele frequency of a genetic ancestry group (FAF95; Whiffin et al. 2017, the statistic of the ClinGen
    BA1/BS1 recommendations). FAF excludes bottlenecked groups (Ashkenazi, Finnish, Amish, remaining), so the overall
    AF still catches founder variants. None if the variant is absent from gnomAD (or no frequency file was given)."""
    xs = [x for x in (v.af, v.faf) if x is not None]
    return max(xs) if xs else None


def gnomad_caveats(v):
    out = []
    if v.gfilter and v.gfilter not in ("PASS", "."):
        out.append(f"gnomAD site filtered ({v.gfilter})")
    if v.gdisc is not None and v.gdisc < 1e-4:
        out.append("gnomAD exome/genome frequencies discordant")
    return out


def freq_cell(v):
    """HTML: overall AF, then grpmax AF, FAF95 and homozygote count where the v4.1 fields are present."""
    overall = v.af
    s = esc(fmt_af(overall))
    det = []
    if v.grpmax_af is not None and v.grpmax:
        det.append(f"max {esc(v.grpmax)} {fmt_af(v.grpmax_af)}")
    if v.faf is not None:
        det.append(f"FAF95 {fmt_af(v.faf)}" + (f" ({esc(v.faf_grp)})" if v.faf_grp else ""))
    if v.nhomalt:
        det.append(f"{int(v.nhomalt)} hom")
    if det:
        s += "<br><span class='small'>" + " &middot; ".join(det) + "</span>"
    return s


class Variant:
    __slots__ = ("chrom", "pos", "ref", "alt", "qual", "filter", "info", "gt", "dp", "ad", "anns", "lof",
                 "zyg", "gene", "effect", "impact", "hgvs_c", "hgvs_p", "feature", "af", "clnsig", "revstat",
                 "clndn", "stars", "rs", "ab", "flags",
                 "an", "nhomalt", "grpmax", "grpmax_af", "faf", "faf_grp", "gfilter", "gdisc", "artefact",
                 "am", "sai", "sai_gene", "repeat")

    def __init__(self):
        self.repeat = None  # repeat context of an in-frame change (repeat_flag), set when the reference is given

    def key(self):
        return f"{self.chrom}:{self.pos}:{self.ref}>{self.alt}"

    def hgvs(self):
        p = self.hgvs_p if self.hgvs_p else ""
        return f"{self.feature}:{self.hgvs_c}" + (f" ({p})" if p else "")


def parse_info(s):
    out = {}
    if s == ".":
        return out
    for kv in s.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
        else:
            out[kv] = True
    return out


def parse_vcf(path, sample_index=0):
    """Variant objects (one per ALT allele; the workflow already split multi-allelics), sample name,
    contig lengths and the chromosome prefix used by the file ("chr" or "")."""
    sample = None
    contigs = {}
    variants = []
    prefix = None
    with opener(path) as f:
        for line in f:
            if line.startswith("##contig="):
                m = re.search(r"ID=([^,>]+).*?length=(\d+)", line)
                if m:
                    contigs[norm_chrom(m.group(1))] = int(m.group(2))
                    if prefix is None:
                        prefix = "chr" if m.group(1).lower().startswith("chr") else ""
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                sample = cols[9 + sample_index] if len(cols) > 9 + sample_index else "sample"
                continue
            t = line.rstrip("\n").split("\t")
            if len(t) < 10:
                continue
            if prefix is None:
                prefix = "chr" if t[0].lower().startswith("chr") else ""
            fmt = t[8].split(":")
            vals = t[9 + sample_index].split(":")
            sd = dict(zip(fmt, vals))
            gt = sd.get("GT", "./.")
            alleles = [a for a in re.split(r"[/|]", gt)]
            if all(a in (".", "0") for a in alleles):
                continue  # reference or no-call
            info = parse_info(t[7])
            for ai, alt in enumerate(t[4].split(","), start=1):
                if alt in ("*", ".", "<*>", "<NON_REF>"):
                    continue
                n_alt = sum(1 for a in alleles if a == str(ai))
                if n_alt == 0:
                    continue
                v = Variant()
                v.chrom, v.pos, v.ref, v.alt = norm_chrom(t[0]), int(t[1]), t[3], alt
                v.qual = fnum(t[5])
                v.filter = t[6]
                v.info = info
                v.gt = gt
                v.zyg = "hom" if n_alt >= 2 or (len(alleles) == 1) else "het"
                v.dp = int(sd["DP"]) if sd.get("DP", ".").isdigit() else None
                ad = [int(x) if x.isdigit() else 0 for x in sd.get("AD", ".").split(",")] if sd.get("AD", ".") != "." else []
                v.ad = ad
                tot = sum(ad)
                v.ab = (ad[ai] / tot) if ad and len(ad) > ai and tot > 0 else None
                anns = []
                for a in info.get("ANN", "").split(","):
                    p = a.split("|")
                    if len(p) >= 11 and p[0] == alt:
                        anns.append(p)
                if not anns:
                    anns = [a.split("|") for a in info.get("ANN", "").split(",") if len(a.split("|")) >= 11]
                v.anns = anns
                first = anns[0] if anns else [""] * 16
                v.effect, v.impact, v.gene, v.feature = first[1], first[2], first[3], first[6]
                v.hgvs_c, v.hgvs_p = first[9], first[10]
                v.lof = info.get("LOF")
                gn = {k: str(x).split(",")[0] for k, x in info.items() if k.startswith("gnomAD_")}
                v.af = fnum(gn.get("gnomAD_AF"))
                v.an = fnum(gn.get("gnomAD_AN"))
                v.nhomalt = fnum(gn.get("gnomAD_nhomalt"))
                v.grpmax, v.grpmax_af = gn.get("gnomAD_grpmax"), fnum(gn.get("gnomAD_grpmax_AF"))
                v.faf, v.faf_grp = fnum(gn.get("gnomAD_FAF95")), gn.get("gnomAD_FAF95_grp")
                v.gfilter, v.gdisc = gn.get("gnomAD_filter"), fnum(gn.get("gnomAD_discordance_p"))
                v.am = fnum(str(info["AM_score"]).split(",")[0]) if info.get("AM_score") not in (None, True) else None
                v.sai, v.sai_gene = spliceai_max(info.get("SpliceAI"), alt)
                v.clnsig = str(info.get("CLNSIG", "")) if info.get("CLNSIG") not in (None, True, ".") else ""
                v.revstat = str(info.get("CLNREVSTAT", "")) if info.get("CLNREVSTAT") not in (None, True) else ""
                v.clndn = str(info.get("CLNDN", "")).replace("_", " ") if info.get("CLNDN") not in (None, True) else ""
                v.stars = REVIEW_STARS.get(v.revstat.split("|")[0], 0) if v.revstat else None
                rs = info.get("RS")
                v.rs = f"rs{rs}" if rs not in (None, True, ".") and not str(rs).startswith("rs") else rs
                v.flags = []
                variants.append(v)
    return sample, contigs, variants, prefix or ""


def spliceai_max(field, alt):
    """Highest SpliceAI delta score (DS_AG/AL/DG/DL) for this ALT over all its genes, and that gene."""
    best, gene = None, None
    if field in (None, True, "."):
        return best, gene
    for entry in str(field).split(","):
        p = entry.split("|")
        if len(p) < 6 or p[0] != alt:
            continue
        ds = [fnum(x) for x in p[2:6]]
        ds = [x for x in ds if x is not None]
        if ds and (best is None or max(ds) > best):
            best, gene = max(ds), p[1]
    return best, gene


# ClinGen-calibrated strengths (Bergquist et al., Genet Med 2025, Table 1) for AlphaMissense; the (+3)/(-3) points
# intervals are reported at Moderate, as the 2015 combining rules have no level between Moderate and Strong.
AM_PP3 = ((0.990, "PP3_Strong"), (0.906, "PP3_Moderate"), (0.792, "PP3_Supporting"))
AM_BP4 = ((0.099, "BP4_Moderate"), (0.169, "BP4_Supporting"))
CODE_RANK = {"PP3_Supporting": 1, "PP3_Moderate": 2, "PP3_Strong": 3}


def am_code(score):
    if score is None:
        return None
    s = round(score, 3)
    for t, c in AM_PP3:
        if s >= t:
            return c
    for t, c in AM_BP4:
        if s <= t:
            return c
    return None


def is_canonical_splice(v, gene):
    a = ann_for_gene(v, gene) if gene else (v.anns[0] if v.anns else None)
    return bool(a) and ("splice_acceptor_variant" in a[1] or "splice_donor_variant" in a[1])


def predictor_codes(v, gene):
    """PP3/BP4 from one calibrated missense predictor (AlphaMissense) and SpliceAI (Walker et al. 2023: max delta
    >= 0.2 PP3, <= 0.1 BP4, Supporting, outside the canonical +/-1,2 positions). The two PP3s are never summed: the
    stronger one is used. BP4 for a missense variant needs SpliceAI <= 0.1 (or no splice score)."""
    am = am_code(v.am) if is_missense(v, gene) and v.am is not None else None
    sp = None
    if v.sai is not None and not is_canonical_splice(v, gene):
        sp = "PP3_Supporting" if v.sai >= 0.2 else ("BP4_Supporting" if v.sai <= 0.1 else None)
    pp3 = [c for c in (am, sp) if c and c.startswith("PP3")]
    if pp3:
        return [max(pp3, key=CODE_RANK.get)]
    if am and am.startswith("BP4") and (v.sai is None or v.sai <= 0.1):
        return [am]
    if sp and not is_missense(v, gene) and not is_lof(v, gene):
        return [sp]  # BP4 from SpliceAI alone for a variant whose only predicted effect would be on splicing
    return []


def splice_candidate(v, gene, phenotype_match):
    """Predicted splice-altering outside the canonical sites: SpliceAI >= 0.5 (the developers' recommended cut-off),
    or >= 0.2 (the ClinGen PP3 cut-off) in a gene matching the patient's phenotype."""
    if v.sai is None or is_lof(v, gene) or (v.sai_gene and v.sai_gene != gene):
        return False
    return v.sai >= (0.2 if phenotype_match else 0.5)


def predictor_text(v):
    parts = []
    if v.am is not None:
        parts.append(f"AlphaMissense {v.am:.3f}")
    if v.sai is not None:
        parts.append(f"SpliceAI {v.sai:.2f}")
    return ", ".join(parts)


def genes_of(v):
    """All genes the variant is annotated in (SnpEff), first = most severe."""
    seen = []
    for a in v.anns:
        if a[3] and a[3] not in seen:
            seen.append(a[3])
    return seen


def ann_for_gene(v, gene):
    for a in v.anns:
        if a[3] == gene:
            return a
    return None


def main_ann(v, gene):
    """The gene's annotation on its principal transcript, approximated as the longest protein among the transcripts the
    variant hits (SnpEff field AA.pos/AA.length; Ensembl canonical/MANE Select is the longest isoform for most genes).
    SnpEff lists the most severe consequence first, often on a minor isoform; ACMG codes use this one instead.
    Without protein lengths (intronic and splice-site hits) the most severe annotation is kept."""
    best, best_len = None, 0
    for a in v.anns:
        if a[3] != gene:
            continue
        if best is None:
            best = a
        m = re.match(r"\d+/(\d+)", a[13]) if len(a) > 13 else None
        if m and int(m.group(1)) > best_len:
            best, best_len = a, int(m.group(1))
    return best


def is_lof(v, gene=None):
    a = ann_for_gene(v, gene) if gene else (v.anns[0] if v.anns else None)
    return bool(a) and any(e in a[1] for e in LOF_EFFECTS)


def is_missense(v, gene=None):
    a = ann_for_gene(v, gene) if gene else (v.anns[0] if v.anns else None)
    return bool(a) and any(e in a[1] for e in MISSENSE_EFFECTS)


def clin_terms(v):
    return [t.strip().lower() for t in re.split(r"[|/,]", v.clnsig.replace("_", " ")) if t.strip()]


def is_plp(v):
    """ClinVar pathogenic or likely pathogenic (any review status), not a conflicting assertion."""
    s = v.clnsig.lower()
    return "pathogenic" in s and "conflicting" not in s


def is_low_penetrance_only(v):
    terms = clin_terms(v)
    return terms and all(("low penetrance" in t) or ("risk" in t) or (t in ("other", "association")) for t in terms)


def stars_txt(n):
    return "" if n is None else ("&#9733;" * n + "&#9734;" * (4 - n))


# ------------------------------------------------------------------------------- knowledge bases
def load_gencc(path):
    """gene -> list of (disease, moi, classification, submitter)."""
    genes = collections.defaultdict(list)
    if not path:
        return genes
    with opener(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        ix = {h.strip('"'): i for i, h in enumerate(header)}
        need = ("gene_symbol", "disease_title", "classification_title", "moi_title", "submitter_title")
        if not all(k in ix for k in need):
            raise SystemExit(f"GenCC file lacks columns {need}")
        for line in f:
            t = [x.strip('"') for x in line.rstrip("\n").split("\t")]
            if len(t) < len(header) - 3:
                continue
            cls = t[ix["classification_title"]]
            if cls not in GENCC_VALID:
                continue
            genes[t[ix["gene_symbol"]]].append((t[ix["disease_title"]], t[ix["moi_title"]], cls,
                                                t[ix["submitter_title"]]))
    return genes


def gene_moi(entries):
    """Summarise modes of inheritance: set of 'AD', 'AR', 'XL', 'MT', 'other'."""
    out = set()
    for _d, moi, _c, _s in entries:
        m = moi.lower()
        if "x-linked" in m:
            out.add("XL")
        elif "recessive" in m:
            out.add("AR")
        elif "dominant" in m or "semidominant" in m:
            out.add("AD")
        elif "mitochondrial" in m:
            out.add("MT")
    return out


def gene_diseases(entries, limit=3):
    best = {}
    for d, moi, c, _s in entries:
        if d not in best or GENCC_RANK[c] > GENCC_RANK[best[d][1]]:
            best[d] = (moi, c)
    ranked = sorted(best.items(), key=lambda kv: -GENCC_RANK[kv[1][1]])
    return [(d, m, c) for d, (m, c) in ranked[:limit]], len(ranked)


def load_constraint(path):
    """gene -> (LOEUF, pLI, missense z) from gnomAD constraint (canonical/MANE transcript rows)."""
    out = {}
    if not path:
        return out
    with opener(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        ix = {h: i for i, h in enumerate(header)}
        gcol = "gene" if "gene" in ix else "gene_symbol"
        loeuf = next((c for c in ("lof.oe_ci.upper", "oe_lof_upper", "LOEUF") if c in ix), None)
        pli = next((c for c in ("lof.pLI", "pLI") if c in ix), None)
        misz = next((c for c in ("mis.z_score", "mis_z") if c in ix), None)
        canon = next((c for c in ("mane_select", "canonical") if c in ix), None)
        for line in f:
            t = line.rstrip("\n").split("\t")
            if canon and t[ix[canon]].lower() not in ("true", "1") and not (
                    canon == "mane_select" and "canonical" in ix and t[ix["canonical"]].lower() == "true"):
                continue
            g = t[ix[gcol]]
            rec = (fnum(t[ix[loeuf]]) if loeuf else None, fnum(t[ix[pli]]) if pli else None,
                   fnum(t[ix[misz]]) if misz else None)
            if g not in out or (rec[0] is not None and out[g][0] is None):
                out[g] = rec
    return out


def load_artefacts(path):
    """Recurrent-artefact list (data/build_artefact_list.sh): (chrom, pos, ref, alt) -> (samples with the allele,
    samples examined). Alleles called in many of the laboratory's own samples but rare or absent in gnomAD."""
    out = {}
    if not path:
        return out
    with opener(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            t = line.rstrip("\n").split("\t")
            if len(t) < 8:
                continue
            info = parse_info(t[7])
            n, total = info.get("ARTEFACT_N"), info.get("ARTEFACT_TOTAL")
            for alt in t[4].split(","):
                if alt not in ("*", "."):
                    out[(norm_chrom(t[0]), int(t[1]), t[3], alt)] = (
                        int(n) if str(n).isdigit() else None, int(total) if str(total).isdigit() else None)
    return out


def load_pharmcat(path):
    """PharmCAT calls-only TSV (-reporterCallsOnlyTsv): a 'PharmCAT <version>' line, a header, one row per gene."""
    out = {"version": "", "calls": []}
    if not path:
        return out
    with opener(path) as f:
        lines = [x.rstrip("\n") for x in f if x.strip()]
    if lines and lines[0].startswith("PharmCAT"):
        out["version"] = lines.pop(0).split(" ", 1)[-1].strip()
    if not lines:
        return out
    header = lines[0].split("\t")
    for line in lines[1:]:
        t = line.split("\t")
        out["calls"].append(dict(zip(header, t + [""] * (len(header) - len(t)))))
    return out


class Phenotype:
    """The patient's HPO terms matched against genes and diseases.

    gene_terms: gene -> matched patient terms; names: term -> name; ic: patient term -> information content
    (-ln of the fraction of annotated diseases that have the term or a descendant); disease_terms: disease ->
    matched patient terms; gene_diseases: gene -> diseases linked to it that match at least one term."""

    def __init__(self):
        self.gene_terms = collections.defaultdict(set)
        self.names = {}
        self.ic = {}
        self.disease_terms = collections.defaultdict(set)
        self.gene_diseases = collections.defaultdict(set)
        self.gene_all = collections.defaultdict(set)  # gene -> every disease linked to it
        self.terms = set()

    def score(self, disease):
        """Fraction of the patient's phenotype information explained by the disease (0-1)."""
        tot = sum(self.ic.get(t, 0.0) for t in self.terms)
        if not tot:
            return 0.0
        return sum(self.ic.get(t, 0.0) for t in self.disease_terms.get(disease, ())) / tot

    def best_diseases(self, gene, limit=2):
        ds = sorted(self.gene_diseases.get(gene, ()), key=lambda d: (-self.score(d), d))
        return [(self.score(d), d, sorted(self.disease_terms[d])) for d in ds[:limit]]


def load_hpo(path, terms):
    """Match the patient's HPO terms with the propagated phenotype_to_genes.txt (a general term matches genes and
    diseases annotated to any of its descendants). Columns: hpo_id, hpo_name, ncbi_gene_id, gene_symbol[, disease_id]."""
    ph = Phenotype()
    ph.terms = set(terms)
    if not path:
        return ph
    term_diseases = collections.defaultdict(set)
    all_diseases = set()
    with opener(path) as f:
        f.readline()
        for line in f:
            t = line.rstrip("\n").split("\t")
            if len(t) < 4:
                continue
            ph.names.setdefault(t[0], t[1])
            d = t[4] if len(t) > 4 and t[4] else None
            if d:
                all_diseases.add(d)
                ph.gene_all[t[3]].add(d)
            if t[0] in terms:
                ph.gene_terms[t[3]].add(t[0])
                if d:
                    term_diseases[t[0]].add(d)
                    ph.disease_terms[d].add(t[0])
                    ph.gene_diseases[t[3]].add(d)
    n = len(all_diseases)
    for term in terms:
        k = len(term_diseases.get(term, ()))
        ph.ic[term] = math.log(n / k) if n and k else 0.0
    return ph


# HPO clinical-course terms (aspect C) and inheritance terms (aspect I) used in phenotype.hpoa.
HPO_ONSET = {"HP:0030674": "antenatal", "HP:0011460": "embryonal", "HP:0011461": "fetal", "HP:0034199": "first trimester",
             "HP:0034198": "second trimester", "HP:0034197": "third trimester", "HP:0003577": "congenital",
             "HP:0003623": "neonatal", "HP:0003593": "infantile", "HP:0011463": "childhood", "HP:0003621": "juvenile",
             "HP:0011462": "young adult", "HP:0025708": "early young adult", "HP:0025709": "intermediate young adult",
             "HP:0025710": "late young adult", "HP:0003581": "adult", "HP:0003596": "middle age", "HP:0003584": "late"}
HPO_COURSE = {"HP:0003676": "progressive", "HP:0003677": "slowly progressive", "HP:0003811": "neonatal death",
              "HP:0001522": "death in infancy", "HP:0003819": "death in childhood", "HP:0003826": "stillbirth"}
HPO_INHERITANCE = {"HP:0000006": "AD", "HP:0000007": "AR", "HP:0001417": "XL", "HP:0001419": "XLR",
                   "HP:0001423": "XLD", "HP:0001427": "mitochondrial", "HP:0001450": "Y-linked",
                   "HP:0001426": "multifactorial", "HP:0010982": "polygenic", "HP:0010984": "digenic",
                   "HP:0001442": "somatic mosaicism", "HP:0003745": "sporadic"}


def load_hpoa(path):
    """HPO disease annotations (phenotype.hpoa): disease -> {name, onset, course, inheritance, n_features}."""
    out = {}
    if not path:
        return out
    with opener(path) as f:
        ix = None
        for line in f:
            if line.startswith("#"):
                continue
            t = line.rstrip("\n").split("\t")
            if ix is None:
                ix = {h: i for i, h in enumerate(t)}
                if "database_id" not in ix:
                    raise SystemExit("phenotype.hpoa lacks the database_id column")
                continue
            if len(t) < len(ix) - 1:
                continue
            d = out.setdefault(t[ix["database_id"]], {"name": t[ix["disease_name"]], "onset": [], "course": [],
                                                      "inheritance": [], "n_features": 0})
            hp, aspect = t[ix["hpo_id"]], t[ix["aspect"]]
            if t[ix["qualifier"]] == "NOT":
                continue
            if aspect == "P":
                d["n_features"] += 1
            elif aspect == "C" and hp in HPO_ONSET and HPO_ONSET[hp] not in d["onset"]:
                d["onset"].append(HPO_ONSET[hp])
            elif aspect == "C" and hp in HPO_COURSE and HPO_COURSE[hp] not in d["course"]:
                d["course"].append(HPO_COURSE[hp])
            elif aspect == "I" and hp in HPO_INHERITANCE and HPO_INHERITANCE[hp] not in d["inheritance"]:
                d["inheritance"].append(HPO_INHERITANCE[hp])
    return out


def load_gencc_disputed(path):
    """gene -> [(disease, classification, submitter)] for Disputed/Refuted/No Known Disease Relationship assertions."""
    out = collections.defaultdict(list)
    if not path:
        return out
    with opener(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        ix = {h.strip('"'): i for i, h in enumerate(header)}
        for line in f:
            t = [x.strip('"') for x in line.rstrip("\n").split("\t")]
            if len(t) < len(header) - 3:
                continue
            cls = t[ix["classification_title"]]
            if cls in ("Disputed Evidence", "Refuted Evidence", "No Known Disease Relationship"):
                out[t[ix["gene_symbol"]]].append((t[ix["disease_title"]], cls, t[ix["submitter_title"]]))
    return out


# ClinGen dosage sensitivity scores (Riggs et al. 2020).
DOSAGE_TEXT = {3: "sufficient evidence", 2: "some evidence", 1: "little evidence", 0: "no evidence",
               30: "autosomal recessive gene", 40: "dosage sensitivity unlikely"}


def load_dosage(path):
    """ClinGen gene curation list (ClinGen_gene_curation_list_GRCh38.tsv): gene -> (HI score, TS score)."""
    out = {}
    if not path:
        return out
    with opener(path) as f:
        ix = None
        for line in f:
            t = line.rstrip("\n").split("\t")
            if line.startswith("#Gene Symbol") or (ix is None and t[0] == "Gene Symbol"):
                ix = {h.lstrip("#"): i for i, h in enumerate(t)}
                continue
            if line.startswith("#") or ix is None or len(t) < len(ix) - 5:
                continue
            hi = t[ix["Haploinsufficiency Score"]]
            ts = t[ix["Triplosensitivity Score"]]
            out[t[ix["Gene Symbol"]]] = (int(hi) if hi.lstrip("-").isdigit() else None,
                                         int(ts) if ts.lstrip("-").isdigit() else None)
    return out


def load_clinvar_protein(path):
    """Protein-level ClinVar table (data/build_clinvar_protein.py): (gene, residue) -> [(ref_aa, alt_aa, class, stars,
    variation_id, chrom, pos, ref, alt)], class P (pathogenic/likely pathogenic) or B (benign/likely benign)."""
    out = collections.defaultdict(list)
    if not path:
        return out
    with opener(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            t = line.rstrip("\n").split("\t")
            if len(t) < 11 or not t[1].isdigit():
                continue
            out[(t[0], int(t[1]))].append((t[2], t[3], t[4], int(t[5]) if t[5].isdigit() else 0, t[6],
                                           norm_chrom(t[7]), int(t[8]) if t[8].isdigit() else 0, t[9], t[10]))
    return out


PM1_MIN_PATHOGENIC = 3  # pathogenic ClinVar missense variants at other residues within the window for PM1
PROTEIN_CHANGE = re.compile(r"^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$")


def protein_evidence(v, gene, cvp, window=5):
    """ACMG PS1 / PM5 from ClinVar missense variants at the same residue, and the ClinVar missense density around it.

    PS1: the same amino-acid change, by a different nucleotide change, is pathogenic/likely pathogenic (>= 1 star) and
    neither variant is predicted to affect splicing (SpliceAI < 0.2 for this one). PM5: a different missense change at
    the same residue is pathogenic/likely pathogenic (>= 1 star). Residues are compared only where the reference amino
    acid matches, i.e. where SnpEff's transcript and ClinVar's RefSeq transcript are numbered alike.
    Returns (codes, notes)."""
    if not cvp or not is_missense(v, gene):
        return [], []
    a = ann_for_gene(v, gene)
    m = PROTEIN_CHANGE.match(a[10]) if a and len(a) > 10 else None
    if not m:
        return [], []
    ref_aa, res, alt_aa = m.group(1), int(m.group(2)), m.group(3)
    here = cvp.get((gene, res), [])
    if here and not any(x[0] == ref_aa for x in here):
        return [], [f"ClinVar residue {res} is {here[0][0]}, not {ref_aa}: transcripts numbered differently"]
    codes, notes = [], []
    same = [x for x in here if x[1] == alt_aa and x[2] == "P" and x[3] >= 1 and (x[5], x[6], x[7], x[8]) != (v.chrom, v.pos, v.ref, v.alt)]
    if same and (v.sai is None or v.sai < 0.2):
        codes.append("PS1")
        notes.append(f"same change p.{ref_aa}{res}{alt_aa} pathogenic in ClinVar (VCV{int(same[0][4]):09d})")
    other = sorted({x[1] for x in here if x[1] != alt_aa and x[2] == "P" and x[3] >= 1})
    if other:
        codes.append("PM5")
        notes.append(f"other pathogenic missense at residue {res}: " + ", ".join(f"p.{ref_aa}{res}{x}" for x in other))
    benign_same = [x for x in here if x[1] == alt_aa and x[2] == "B" and x[3] >= 1]
    if benign_same:
        notes.append(f"same amino-acid change benign in ClinVar (VCV{int(benign_same[0][4]):09d})")
    n_p = n_b = n_p_near = 0
    for r in range(res - window, res + window + 1):
        for x in cvp.get((gene, r), []):
            if x[3] >= 1:
                n_p += x[2] == "P"
                n_b += x[2] == "B"
                n_p_near += x[2] == "P" and r != res
    if n_p_near >= PM1_MIN_PATHOGENIC and n_b == 0:
        # PM1 (mutational hot spot), Supporting: pathogenic missense at other residues nearby and no benign one; the
        # variant's own residue is left to PS1/PM5 so that the same ClinVar entries are not counted twice
        codes.append("PM1_Supporting")
    if n_p >= 3:
        notes.append(f"ClinVar missense within {window} residues: {n_p} pathogenic, {n_b} benign"
                     + (" (pathogenic cluster)" if n_b == 0 else ""))
    return codes, notes


# ----------------------------------------------------------------------------------- read-backed phase
CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def read_alleles(start, cigar, seq, variants):
    """For one aligned read: variant key -> 'alt' / 'ref' / None (not covered or unclear). SNVs are read from the
    aligned base; indels from an insertion/deletion of the right length right after the VCF anchor base."""
    aligned = {}   # ref pos -> read base
    ins = {}       # ref pos of the base before the insertion -> inserted length
    dels = {}      # first deleted ref pos -> deleted length
    rpos, qpos = start, 0
    for n, op in CIGAR.findall(cigar):
        n = int(n)
        if op in "M=X":
            for i in range(n):
                aligned[rpos + i] = seq[qpos + i]
            rpos += n
            qpos += n
        elif op == "I":
            ins[rpos - 1] = n
            qpos += n
        elif op in "DN":
            dels[rpos] = n
            rpos += n
        elif op == "S":
            qpos += n
    out = {}
    for v in variants:
        if len(v.ref) == 1 and len(v.alt) == 1:
            b = aligned.get(v.pos)
            out[v.key()] = None if b is None else ("alt" if b == v.alt else ("ref" if b == v.ref else None))
            continue
        span = range(v.pos, v.pos + max(len(v.ref), 2))
        if not all(p in aligned or p in dels or any(d <= p < d + dels[d] for d in dels) for p in span):
            out[v.key()] = None
            continue
        if len(v.ref) > len(v.alt):
            hit = dels.get(v.pos + 1) == len(v.ref) - len(v.alt)
        else:
            hit = ins.get(v.pos) == len(v.alt) - len(v.ref)
        clean = not any(v.pos <= p <= v.pos + len(v.ref) for p in list(ins) + list(dels))
        out[v.key()] = "alt" if hit else ("ref" if clean else None)
    return out


def phase_pair(bam, prefix, v1, v2, max_distance=1000):
    """Read-backed phase of two heterozygous variants: ('trans'|'cis'|'unclear'|'too far apart', reads cis, reads trans).
    A read or read pair that covers both positions supports cis (alt-alt or ref-ref) or trans (alt-ref or ref-alt)."""
    if v1.chrom != v2.chrom or abs(v1.pos - v2.pos) > max_distance:
        return "too far apart", 0, 0
    lo, hi = sorted((v1.pos, v2.pos))
    try:
        sam = run_samtools(["view", "-F", "3852", "-q", "20", bam, f"{prefix}{v1.chrom}:{lo}-{hi + 1}"])
    except (Exception, SystemExit):
        return "unclear", 0, 0
    per = collections.defaultdict(dict)
    for line in sam.splitlines():
        t = line.split("\t")
        if len(t) < 10 or t[5] == "*":
            continue
        for k, s in read_alleles(int(t[3]), t[5], t[9], (v1, v2)).items():
            if s is None:
                continue
            if per[t[0]].get(k, s) != s:
                per[t[0]][k] = "conflict"
            else:
                per[t[0]][k] = s
    cis = trans = 0
    for d in per.values():
        a, b = d.get(v1.key()), d.get(v2.key())
        if a in ("alt", "ref") and b in ("alt", "ref"):
            if a == b:
                cis += 1
            else:
                trans += 1
    if cis + trans >= 4 and trans >= 0.9 * (cis + trans):
        return "trans", cis, trans
    if cis + trans >= 4 and cis >= 0.9 * (cis + trans):
        return "cis", cis, trans
    return "unclear", cis, trans


def phase_text(ph):
    state, cis, trans = ph
    if state == "trans":
        return f"in trans (read-backed: {trans} of {cis + trans} read pairs)"
    if state == "cis":
        return f"in cis on the same allele (read-backed: {cis} of {cis + trans} read pairs)"
    if state == "too far apart":
        return "phase unknown (too far apart for read-backed phasing; test the parents)"
    return f"phase unknown (read-backed phasing unclear: {cis} cis, {trans} trans read pairs; test the parents)"


def load_roh(path):
    segs = []
    if not path:
        return segs
    with opener(path) as f:
        for line in f:
            t = line.rstrip("\n").split("\t")
            if t and t[0] == "RG" and len(t) >= 8:
                segs.append((norm_chrom(t[2]), int(t[3]), int(t[4]), int(t[5]), int(t[6]), fnum(t[7])))
    return segs


def load_thresholds(path):
    """mosdepth thresholds BED -> list of (chrom, start, end, name, {threshold: bases})."""
    rows = []
    if not path:
        return rows
    with opener(path) as f:
        header = None
        for line in f:
            t = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                header = t
                continue
            if len(t) < 5:
                continue
            cols = header[4:] if header else [f"{i}X" for i in range(len(t) - 4)]
            rows.append((norm_chrom(t[0]), int(t[1]), int(t[2]), t[3],
                         {c.upper(): int(x) for c, x in zip(cols, t[4:])}))
    return rows


def load_mosdepth_summary(path):
    out = {}
    if not path:
        return out
    with opener(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            t = line.rstrip("\n").split("\t")
            if len(t) >= 4:
                rec = dict(zip(header, t))
                out[rec.get("chrom", t[0])] = rec
    return out


def run_samtools(args):
    r = subprocess.run(["samtools"] + args, capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"samtools {' '.join(args[:2])} failed: {r.stderr.strip()[:500]}")
    return r.stdout


def site_depths(bam, sites, prefix):
    """Read depth (MAPQ >= 20, base quality >= 20, duplicates excluded) at single positions."""
    out = {}
    for chrom, pos in sites:
        d = 0
        for line in run_samtools(["depth", "-a", "-Q", "20", "-q", "20", "-r", f"{prefix}{chrom}:{pos}-{pos}",
                                  bam]).splitlines():
            t = line.split("\t")
            if len(t) >= 3 and int(t[1]) == pos:
                d = int(t[2])
        out[(chrom, pos)] = d
    return out


def reference_contexts(fasta, prefix, variants, flank=20):
    """Reference sequence around each variant: {key: (seq, offset of POS in seq)}."""
    regions, keys = [], []
    for v in variants:
        start = max(1, v.pos - flank)
        regions.append(f"{prefix}{v.chrom}:{start}-{v.pos + len(v.ref) + flank}")
        keys.append((v.key(), v.pos - start))
    out = {}
    for i in range(0, len(regions), 200):
        seqs, cur = [], None
        for line in run_samtools(["faidx", fasta] + regions[i:i + 200]).splitlines():
            if line.startswith(">"):
                if cur is not None:
                    seqs.append("".join(cur).upper())
                cur = []
            else:
                cur.append(line.strip())
        if cur is not None:
            seqs.append("".join(cur).upper())
        for (k, off), s in zip(keys[i:i + 200], seqs):
            out[k] = (s, off)
    return out


def repeat_flag(seq, off, ref_len):
    """Homopolymer (>= 6) or short tandem repeat (unit 2-4, >= 4 copies) at or next to the variant."""
    best_h, base = 0, ""
    for i in range(max(0, off - 1), min(len(seq), off + max(ref_len, 1) + 1)):
        j, k = i, i
        while j > 0 and seq[j - 1] == seq[i]:
            j -= 1
        while k + 1 < len(seq) and seq[k + 1] == seq[i]:
            k += 1
        if k - j + 1 > best_h:
            best_h, base = k - j + 1, seq[i]
    if best_h >= 6:
        return f"homopolymer ({base} x{best_h})"
    for u in (2, 3, 4):
        for s in range(max(0, off - 3 * u), off + 1):
            unit = seq[s:s + u]
            if len(unit) < u or len(set(unit)) == 1:
                continue
            n = 1
            while seq[s + n * u:s + (n + 1) * u] == unit:
                n += 1
            if n >= 4 and s <= off <= s + n * u:
                return f"tandem repeat (({unit}) x{n})"
    return None


# ------------------------------------------------------------------------------------- analyses
def sample_qc(variants, contigs, msum, target_cov):
    snv = [v for v in variants if len(v.ref) == 1 and len(v.alt) == 1]
    indel = [v for v in variants if not (len(v.ref) == 1 and len(v.alt) == 1)]
    ts = sum(1 for v in snv if {v.ref, v.alt} in ({"A", "G"}, {"C", "T"}))
    tv = len(snv) - ts
    auto = [v for v in snv if v.chrom not in ("X", "Y", "MT")]
    het = sum(1 for v in auto if v.zyg == "het")
    hom = len(auto) - het
    abs_ = sorted(v.ab for v in auto if v.zyg == "het" and v.ab is not None and (v.dp or 0) >= 20)
    qc = {
        "variants": len(variants), "snv": len(snv), "indel": len(indel),
        "tstv": ts / tv if tv else None, "het_hom": het / hom if hom else None,
        "het_ab_median": abs_[len(abs_) // 2] if abs_ else None,
        "het_ab_low": (sum(1 for a in abs_ if a < 0.3) / len(abs_)) if abs_ else None,
    }
    # sex: chrX non-PAR heterozygosity and (if available) X/Y coverage relative to autosomes
    xs = [v for v in snv if v.chrom == "X" and not any(a <= v.pos <= b for a, b in PAR_GRCH38)
          and (v.qual or 0) >= 30 and (v.dp or 0) >= 15]
    xhet = sum(1 for v in xs if v.zyg == "het")
    qc["x_snvs"] = len(xs)
    qc["x_het_frac"] = xhet / len(xs) if xs else None
    auto_cov = x_cov = y_cov = None
    if msum:
        def mean_of(keys):
            for k in keys:
                if k in msum:
                    return fnum(msum[k].get("mean"))
            return None
        auto_cov = mean_of(["total_region", "total"])
        x_cov = mean_of(["X_region", "chrX_region"])
        y_cov = mean_of(["Y_region", "chrY_region"])
        qc["mean_target_cov"] = auto_cov
    qc["x_ratio"] = (x_cov / auto_cov) if (x_cov and auto_cov) else None
    qc["y_ratio"] = (y_cov / auto_cov) if (y_cov is not None and auto_cov) else None
    votes = []
    if qc["x_het_frac"] is not None and len(xs) >= 20:
        votes.append("male" if qc["x_het_frac"] < 0.1 else ("female" if qc["x_het_frac"] > 0.2 else "?"))
    if qc["y_ratio"] is not None:
        votes.append("male" if qc["y_ratio"] > 0.1 else ("female" if qc["y_ratio"] < 0.02 else "?"))
    if votes and all(x == votes[0] for x in votes) and votes[0] != "?":
        qc["sex"] = votes[0]
    elif votes:
        qc["sex"] = "inconclusive (" + ", ".join(votes) + ")"
    else:
        qc["sex"] = "not determined"
    if target_cov:
        tot = sum(r[2] - r[1] for r in target_cov)
        for thr in ("10X", "20X", "30X"):
            if thr in target_cov[0][4]:
                qc[f"pct_{thr}"] = 100.0 * sum(r[4][thr] for r in target_cov) / tot if tot else None
    return qc


def variant_flags(v):
    f = []
    if v.dp is not None and v.dp < 20:
        f.append(f"low depth ({v.dp}x)")
    if v.zyg == "het" and v.ab is not None and (v.ab < 0.25 or v.ab > 0.75):
        f.append(f"allele balance {v.ab:.2f}")
    if v.zyg == "hom" and v.ab is not None and v.ab < 0.85:
        f.append(f"hom call with alt fraction {v.ab:.2f}")
    if v.qual is not None and v.qual < 50:
        f.append(f"QUAL {v.qual:.0f}")
    if len(v.ref) != len(v.alt):
        f.append("indel")
    g = v.gene or ""
    if g.startswith("HLA-") or g.startswith("KIR") or g in POLYMORPHIC_GENES:
        f.append("difficult/polymorphic region")
    elif g in PARALOG_GENES:
        f.append("segmental duplication / pseudogene: short-read calls unreliable")
    return f + gnomad_caveats(v)


def lof_mechanism(gene, moi, constraint, dosage):
    """Is loss of function a known disease mechanism? Recessive or X-linked inheritance, ClinGen haploinsufficiency
    score 3 or 2, or a LoF-intolerant gene (LOEUF < 0.6) unless ClinGen scores haploinsufficiency as unlikely (40)."""
    hi = dosage.get(gene, (None, None))[0]
    if "AR" in moi or "XL" in moi or hi in (2, 3):
        return True
    loeuf = constraint.get(gene, (None, None, None))[0]
    return hi != 40 and loeuf is not None and loeuf < 0.6


def nmd_predicted(v, gene):
    """SnpEff's NMD prediction (INFO/NMD lists the genes whose transcripts are predicted to undergo NMD)."""
    nmd = v.info.get("NMD")
    if nmd in (None, True):
        return nmd is True
    return any(e.strip("()").split("|")[0] == gene for e in str(nmd).split(","))


def nmd_by_exon(a):
    """NMD expected from the exon of the stop (SnpEff field Rank/Total): a premature stop before the penultimate exon
    is > 50 nt upstream of the last exon junction. Stops in the last two exons and single-exon genes are left to the
    protein-fraction rule (the 3' end of the penultimate exon escapes NMD, and the exon boundaries are not in the VCF)."""
    m = re.match(r"(\d+)/(\d+)", a[8]) if len(a) > 8 else None
    return bool(m) and int(m.group(1)) < int(m.group(2)) - 1


def protein_fraction_lost(a):
    """Share of the protein after a truncating variant, from SnpEff's AA.pos / AA.length, or None."""
    m = re.match(r"(\d+)/(\d+)", a[13]) if len(a) > 13 else None
    if not m or int(m.group(2)) == 0:
        return None
    return 1 - int(m.group(1)) / int(m.group(2))


def pvs1_code(v, gene):
    """PVS1 strength after the ClinGen SVI decision tree (Abou Tayoun et al., Hum Mutat 2018), the parts that a VCF
    supports: nonsense/frameshift predicted to undergo NMD, canonical splice-site, exon or transcript loss: Very strong;
    nonsense/frameshift escaping NMD (last exon, end of the penultimate exon) that removes > 10 % of the protein:
    Strong, otherwise Moderate; start loss: Moderate (another in-frame start codon may be used)."""
    a = main_ann(v, gene) or v.anns[0]
    eff = a[1]
    if "start_lost" in eff and not any(e in eff for e in ("stop_gained", "frameshift", "splice_")):
        return "PVS1_Moderate"
    if ("stop_gained" in eff or "frameshift" in eff) and "splice_" not in eff and not nmd_predicted(v, gene) \
            and not nmd_by_exon(a):
        lost = protein_fraction_lost(a)
        return "PVS1_Strong" if lost is None or lost > 0.10 else "PVS1_Moderate"
    return "PVS1"


# Variants excluded from BA1 by the ClinGen SVI (Ghosh et al., Hum Mutat 2018): common, but pathogenic
BA1_EXCEPTIONS = {("ACAD9", "c.-44_-41dupTAAG"), ("BTD", "c.1330G>C"), ("GJB2", "c.-23+1G>A"), ("GJB2", "c.109G>A"),
                  ("HFE", "c.187C>G"), ("HFE", "c.845G>A"), ("MEFV", "c.1105C>T"), ("MEFV", "c.1223G>A"),
                  ("PIBF1", "c.1214G>A")}
INFRAME_EFFECTS = ("inframe_insertion", "inframe_deletion", "disruptive_inframe_insertion",
                   "disruptive_inframe_deletion", "stop_lost")
BS1_DOMINANT = 1e-3  # BS1 cut-off for genes with only dominant disorders (the report's 'too common for dominant')


def acmg_evidence(v, gene, moi, constraint, hpo_hit, max_af, dosage=None, protein_codes=()):
    """Automated ACMG/AMP 2015 evidence codes with the ClinGen SVI strength modifications, the subset a single VCF
    with population, predictor and ClinVar-residue data supports (no functional, segregation or case-count evidence).
    v.repeat (set by the caller from the reference sequence) separates PM4 from BP3 for in-frame changes."""
    codes = []
    a = main_ann(v, gene) or (v.anns[0] if v.anns else [""] * 16)
    utr_splice = "splice_" in a[1] and a[9].startswith(("c.-", "c.*"))  # splice site of an intron in a UTR
    if any(e in a[1] for e in LOF_EFFECTS) and not utr_splice and lof_mechanism(gene, moi, constraint, dosage or {}):
        codes.append(pvs1_code(v, gene))
    p = freq(v)
    dominant_only = bool(moi) and not ({"AR", "XL"} & moi)
    if p is None or p < (1e-4 if dominant_only else 1e-3):
        codes.append("PM2_supporting")
    if p is not None and p > 0.05 and (gene, a[9]) not in BA1_EXCEPTIONS:
        codes.append("BA1")
    elif p is not None and (p > max_af or (dominant_only and p > BS1_DOMINANT)) and (gene, a[9]) not in BA1_EXCEPTIONS:
        codes.append("BS1")
    if v.nhomalt and v.nhomalt >= 2 and moi == {"AR"} and v.chrom != "X":
        codes.append("BS2_Supporting")  # homozygotes in gnomAD for a recessive disorder (penetrance unknown)
    codes += list(protein_codes)
    if any(e in a[1] for e in INFRAME_EFFECTS) and not is_lof(v, gene):
        codes.append("BP3" if getattr(v, "repeat", None) and "stop_lost" not in a[1] else "PM4")
    if is_missense(v, gene):
        z = constraint.get(gene, (None, None, None))[2]
        if z is not None and z > 3.09:
            codes.append("PP2")
    if not any(c.startswith("PVS1") for c in codes):  # PP3 is not combined with PVS1
        pc = predictor_codes(v, gene)
        if "synonymous_variant" in a[1] and pc == ["BP4_Supporting"]:
            pc = ["BP7"]  # synonymous, no predicted splice effect (Walker et al. 2023: BP7 instead of BP4)
        codes += pc
    if hpo_hit:
        codes.append("PP4")
    return codes


# Points of the ACMG/AMP criteria (Tavtigian et al., Genet Med 2018; Hum Mutat 2020): Supporting 1, Moderate 2,
# Strong 4, Very strong 8, benign criteria negative; BA1 is stand-alone benign
STRENGTH_POINTS = {"supporting": 1, "moderate": 2, "strong": 4, "verystrong": 8}
DEFAULT_STRENGTH = {"PVS": "verystrong", "PS": "strong", "PM": "moderate", "PP": "supporting", "BS": "strong",
                    "BP": "supporting"}
ACMG_CLASSES = [(10, "pathogenic"), (6, "likely pathogenic"), (0, "uncertain significance"), (-6, "likely benign")]


def code_points(code):
    """Points of one evidence code such as PVS1, PM2_supporting or PP3_Strong (None for BA1)."""
    base, _, strength = code.partition("_")
    if base == "BA1":
        return None
    pts = STRENGTH_POINTS[strength.lower().replace("_", "") or DEFAULT_STRENGTH[base[:-1] if base[:3] != "PVS" else "PVS"]]
    return -pts if base.startswith("B") else pts


def acmg_classify(codes):
    """(class, points) from evidence codes by the point system. PM1 and PP3 together count at most 4 points (Strong;
    Pejaver et al. 2022). BA1 alone makes a variant benign; a variant with BA1 and pathogenic codes is not benign."""
    codes = list(dict.fromkeys(codes))
    pts = {c: code_points(c) for c in codes if c != "BA1"}
    pm1 = sum(p for c, p in pts.items() if c.startswith("PM1"))
    pp3 = sum(p for c, p in pts.items() if c.startswith("PP3"))
    total = sum(pts.values()) - max(0, pm1 + pp3 - 4) if pm1 and pp3 else sum(pts.values())
    if "BA1" in codes:
        return ("benign", total) if total <= 0 or not any(p > 0 for p in pts.values()) else ("conflicting (BA1)", total)
    for t, name in ACMG_CLASSES:
        if total >= t:
            return name, total
    return "benign", total


ACMG_ABBR = {"pathogenic": "P", "likely pathogenic": "LP", "uncertain significance": "VUS", "likely benign": "LB",
             "benign": "B", "conflicting (BA1)": "VUS"}


def acmg_tag(cls, pts):
    kind = "t-red" if cls in ("pathogenic", "likely pathogenic") else ("t-grn" if "benign" in cls and "BA1" not in cls else "t-amb")
    return tag(f"{cls[0].upper() + cls[1:]} ({pts:+d})", kind)


# --------------------------------------------------------------------------------------- report
class Report:
    def __init__(self):
        self.parts = []
        self.rows = []  # machine-readable findings

    def h(self, s):
        self.parts.append(s)

    def table(self, header, rows, cls=""):
        if not rows:
            self.h('<p class="none">None.</p>')
            return
        self.h(f'<table class="{cls}"><thead><tr>' + "".join(f"<th>{esc(x)}</th>" for x in header)
               + "</tr></thead><tbody>")
        for r in rows:
            self.h("<tr>" + "".join(f"<td>{x}</td>" for x in r) + "</tr>")
        self.h("</tbody></table>")


CSS = """
body{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;margin:24px auto;max-width:1180px;
color:#1d2430;line-height:1.45;padding:0 16px}
h1{font-size:1.6em;margin-bottom:0}h2{font-size:1.25em;border-bottom:2px solid #d8dee8;padding-bottom:4px;
margin-top:2em}h3{font-size:1.05em;margin-bottom:.3em}
.sub{color:#5b6576;margin-top:2px}.box{border:1px solid #d8dee8;border-radius:8px;padding:12px 16px;
background:#f7f9fc;margin:12px 0}.warn{border-color:#e7b85c;background:#fff8e8}.alert{border-color:#d9534f;
background:#fdf1f0}.ok{border-color:#6fb26f;background:#f1f9f1}
table{border-collapse:collapse;width:100%;font-size:.86em;margin:8px 0 14px}th,td{border:1px solid #d8dee8;
padding:4px 6px;vertical-align:top;text-align:left}th{background:#eef2f7}tr:nth-child(even) td{background:#fbfcfe}
.tag{display:inline-block;padding:0 6px;border-radius:9px;font-size:.8em;background:#e6ebf2;margin:1px}
.t-red{background:#f6d3d1}.t-amb{background:#fbe7bf}.t-grn{background:#d7eed7}.t-blu{background:#d6e6f7}
.none{color:#5b6576;font-style:italic}.small{font-size:.84em;color:#5b6576}code{font-size:.92em}
ul.kv{list-style:none;padding:0;columns:2}ul.kv li{margin:2px 0}
"""


def tag(text, kind=""):
    return f'<span class="tag {kind}">{esc(text)}</span>'


def var_cell(v):
    s = f"<b>{esc(v.gene)}</b><br><code>{esc(v.hgvs_c)}</code>"
    if v.hgvs_p:
        s += f"<br><code>{esc(v.hgvs_p)}</code>"
    s += f'<br><span class="small">{esc(v.chrom)}:{v.pos} {esc(v.ref)}&gt;{esc(v.alt)}'
    if v.rs:
        s += f" {esc(v.rs)}"
    s += f"<br>{esc(v.feature)}</span>"
    return s


def low_confidence(v):
    """Reasons the call itself is unreliable (validated on GIAB HG001/HG002/HG005: every false-positive finding had
    one of these, no confirmed finding did): het alt fraction < 0.25, QUAL < 30, an indel in a homopolymer or
    tandem repeat of >= 8, or a paralogous/hypervariable gene); and an allele on the laboratory's recurrent-artefact
    list."""
    why = []
    if v.zyg == "het" and v.ab is not None and v.ab < 0.25:
        why.append("alt fraction < 0.25")
    if v.qual is not None and v.qual < 30:
        why.append("QUAL < 30")
    if len(v.ref) != len(v.alt):
        for f in v.flags:
            m = re.search(r"x(\d+)\)", f)
            if (f.startswith("homopolymer") or f.startswith("tandem")) and m and int(m.group(1)) >= 8 - (
                    3 if f.startswith("tandem") else 0):
                why.append("indel in a repeat")
    if any(f.startswith(("segmental duplication", "difficult/polymorphic", "REF does not match")) for f in v.flags):
        why.append("paralogous / hypervariable region")
    if any(f.startswith("recurrent artefact") for f in v.flags):
        why.append("recurrent assay artefact")
    return why


def conf_sorted(vs, key=lambda x: x):
    """Confident calls first, keeping the original order within each group."""
    return sorted(vs, key=lambda x: bool(low_confidence(key(x))))


def call_cell(v, sex):
    z = v.zyg
    if v.chrom == "X" and sex == "male" and z == "hom" and not any(a <= v.pos <= b for a, b in PAR_GRCH38):
        z = "hemizygous"
    ad = ",".join(map(str, v.ad)) if v.ad else "."
    s = f"{z}<br><span class='small'>GT {esc(v.gt)} DP {v.dp if v.dp is not None else '.'} AD {ad}"
    if v.ab is not None:
        s += f" (alt {v.ab:.2f})"
    s += "</span>"
    lc = low_confidence(v)
    if lc:
        s += "<br>" + tag("LOW CONFIDENCE: " + ", ".join(lc), "t-red")
    if v.flags:
        s += "<br>" + "".join(tag(f, "t-amb") for f in v.flags)
    return s


def clinvar_cell(v):
    if not v.clnsig:
        return '<span class="small">not in ClinVar</span>'
    s = esc(v.clnsig.replace("_", " ").replace("|", "; "))
    if v.revstat:
        s += f"<br><span class='small'>{stars_txt(v.stars)} {esc(v.revstat.replace('_', ' ').replace('|', ', '))}</span>"
    if v.clndn:
        dn = [d for d in v.clndn.split("|") if d and d not in ("not provided", "not specified")]
        if dn:
            s += f"<br><span class='small'>{esc('; '.join(dn[:4]))}{' ...' if len(dn) > 4 else ''}</span>"
    return s


class Knowledge:
    """Gene- and disease-level knowledge bases given to the report."""

    def __init__(self, gencc=None, disputed=None, pheno=None, hpoa=None, dosage=None):
        self.gencc = gencc or {}
        self.disputed = disputed or {}
        self.pheno = pheno or Phenotype()
        self.hpoa = hpoa or {}
        self.dosage = dosage or {}


def disease_features(d, info):
    bits = []
    if info.get("onset"):
        order = list(HPO_ONSET.values())
        bits.append("onset " + "/".join(sorted(info["onset"], key=order.index)[:3]))
    if info.get("inheritance"):
        bits.append(", ".join(info["inheritance"]))
    if info.get("course"):
        bits.append(", ".join(info["course"][:2]))
    return "; ".join(bits)


def disease_cell(gene, K, phenotype=True):
    entries = K.gencc.get(gene, [])
    s = ""
    if entries:
        ds, n = gene_diseases(entries)
        s = "<br>".join(f"{esc(d)} <span class='small'>({esc(m)}; {esc(c)})</span>" for d, m, c in ds)
        if n > len(ds):
            s += f"<br><span class='small'>+{n - len(ds)} more</span>"
    else:
        s = '<span class="small">no established gene-disease association in GenCC</span>' if K.gencc else ""
    for d, cls, sub in K.disputed.get(gene, [])[:2]:
        s += "<br>" + tag(f"{cls.replace(' Evidence', '').lower()}: {d} ({sub})", "t-amb")
    hi, ts = K.dosage.get(gene, (None, None))
    dz = []
    if hi in (1, 2, 3, 30, 40):
        dz.append(f"haploinsufficiency: {DOSAGE_TEXT[hi]}")
    if ts in (2, 3):
        dz.append(f"triplosensitivity: {DOSAGE_TEXT[ts]}")
    if dz:
        s += f"<br><span class='small'>ClinGen {esc('; '.join(dz))}</span>"
    P = K.pheno
    if phenotype and P.gene_terms.get(gene):
        s += "<br>" + "".join(tag(f"phenotype match: {P.names.get(t, t)}", "t-blu") for t in sorted(P.gene_terms[gene]))
        for sc, d, terms in P.best_diseases(gene, 1):
            info = K.hpoa.get(d, {})
            s += (f"<br><span class='small'>best match: <b>{esc(info.get('name', d))}</b> ({esc(d)}): "
                  f"{len(terms)} of {len(P.terms)} terms, {sc:.0%} of the phenotype information"
                  + (f"; {esc(disease_features(d, info))}" if disease_features(d, info) else "") + "</span>")
    elif K.hpoa and P.gene_all.get(gene):
        feats = []
        for d in sorted(P.gene_all[gene], key=lambda d: -K.hpoa.get(d, {}).get("n_features", 0))[:2]:
            info = K.hpoa.get(d)
            if info and disease_features(d, info):
                feats.append(f"{info['name']}: {disease_features(d, info)}")
        if feats:
            s += "<br><span class='small'>HPO: " + esc(" | ".join(feats)) + "</span>"
    return s


def pgx_pheno(gene, alle):
    """CPIC-style genotype and phenotype from allele copy numbers: alle = [(name, rs, function, copies, ...)]."""
    cnt = {a[0]: a[3] for a in alle}
    if any(c is None for c in cnt.values()):
        unk = [k for k, c in cnt.items() if c is None]
        return "indeterminate", "not assessable: " + ", ".join(unk) + " (off-target or low coverage)"
    if gene == "CYP2C19":
        nf = cnt["*2"] + cnt["*3"]
        inc = cnt["*17"]
        dip = "/".join((["*2"] * cnt["*2"] + ["*3"] * cnt["*3"] + ["*17"] * inc + ["*1"] * 2)[:2])
        if nf >= 2:
            return dip, "Poor metabolizer"
        if nf == 1:
            return dip, "Intermediate metabolizer"
        return dip, {0: "Normal metabolizer", 1: "Rapid metabolizer"}.get(inc, "Ultrarapid metabolizer")
    if gene == "CYP2C9":
        score = 2 - 0.5 * cnt["*2"] - 1.0 * cnt["*3"]
        dip = "/".join((["*2"] * cnt["*2"] + ["*3"] * cnt["*3"] + ["*1"] * 2)[:2])
        ph = "Normal metabolizer" if score >= 2 else ("Intermediate metabolizer" if score >= 1 else "Poor metabolizer")
        return dip, f"{ph} (activity score {score:g})"
    if gene == "DPYD":
        losses = sorted([0.0] * (cnt["*2A (c.1905+1G>A)"] + cnt["*13 (c.1679T>G)"]) +
                        [0.5] * (cnt["c.2846A>T"] + cnt["HapB3 (c.1236G>A)"]))
        alleles = [1.0, 1.0]
        for i, val in enumerate(losses[:2]):
            alleles[i] = val
        score = sum(alleles)
        ph = "Normal metabolizer" if score >= 2 else ("Intermediate metabolizer" if score >= 1 else "Poor metabolizer")
        note = " (two variants assumed in trans)" if len(losses) >= 2 else ""
        return f"activity score {score:g}", ph + note
    if gene == "TPMT":
        b, c, two = cnt["*3B (c.460G>A)"], cnt["*3C (c.719A>G)"], cnt["*2 (c.238G>C)"]
        a3 = min(b, c)
        nf = two + a3 + (b - a3) + (c - a3)
        dip = "/".join(((["*3A"] * a3) + ["*2"] * two + ["*3B"] * (b - a3) + ["*3C"] * (c - a3) + ["*1"] * 2)[:2])
        return dip, {0: "Normal metabolizer", 1: "Intermediate metabolizer"}.get(nf, "Possible/Poor metabolizer")
    if gene == "NUDT15":
        n = cnt["*3 (c.415C>T)"]
        return "/".join((["*3"] * n + ["*1"] * 2)[:2]), {0: "Normal metabolizer", 1: "Intermediate metabolizer"}.get(n, "Poor metabolizer")
    if gene == "CYP3A5":
        n = cnt["*3"]
        return "/".join((["*3"] * n + ["*1"] * 2)[:2]), {2: "Poor metabolizer (CYP3A5 non-expresser; the usual genotype outside African ancestry)",
                                                         1: "Intermediate metabolizer", 0: "Normal metabolizer (CYP3A5 expresser)"}[n]
    if gene == "SLCO1B1":
        n = cnt["c.521T>C"]
        return {0: "c.521TT", 1: "c.521TC", 2: "c.521CC"}[n], {0: "Normal function", 1: "Decreased function",
                                                              2: "Poor function"}[n]
    if gene == "VKORC1":
        n = cnt["c.-1639G>A"]
        return {0: "-1639GG", 1: "-1639GA", 2: "-1639AA"}[n], {0: "typical warfarin dose requirement",
                                                              1: "increased warfarin sensitivity",
                                                              2: "high warfarin sensitivity"}[n]
    if gene == "CYP4F2":
        n = cnt["*3"]
        return "/".join((["*3"] * n + ["*1"] * 2)[:2]), "decreased function (slightly higher warfarin dose)" if n else "normal"
    if gene == "CYP2B6":
        n = cnt["c.516G>T (*6/*9)"]
        return {0: "516GG", 1: "516GT", 2: "516TT"}[n], {0: "Normal metabolizer (for this allele)",
                                                        1: "Intermediate metabolizer (likely)",
                                                        2: "Poor metabolizer (likely)"}[n]
    if gene == "UGT1A1":
        n6, n80 = cnt["*6 (c.211G>A)"], cnt["*80 (c.-364C>T, *28 proxy)"]
        dec = n6 + n80
        return f"*6 x{n6}, *80 x{n80}", {0: "Normal metabolizer (likely)", 1: "Intermediate metabolizer (likely)"}.get(dec, "Poor metabolizer (likely; Gilbert genotype)")
    return "", ""


# ----------------------------------------------------------------------------------------------- trio
def decompose(pos, ref, alt):
    """Minimal representation(s) of one ALT, matching the workflow's bcftools norm --atomize where no reference is
    needed: common suffix and prefix trimmed (keeping one anchor base for indels), MNPs split into SNVs."""
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    if len(ref) == len(alt):
        return [(pos + i, r, a) for i, (r, a) in enumerate(zip(ref, alt)) if r != a]
    return [(pos, ref, alt)]


def load_parent_vcf(path):
    """A parent's calls (any single-sample VCF, raw or filtered): (chrom, pos, ref, alt) -> (zygosity, depth, alt
    fraction, QUAL). Multi-allelic and complex records are split and trimmed like the workflow's normalisation."""
    out = {}
    if not path:
        return out
    with opener(path) as f:
        for line in f:
            if line[0] == "#":
                continue
            t = line.rstrip("\n").split("\t")
            if len(t) < 10:
                continue
            sd = dict(zip(t[8].split(":"), t[9].split(":")))
            gt = re.split(r"[/|]", sd.get("GT", "./."))
            ad = [int(x) if x.isdigit() else 0 for x in sd.get("AD", "").split(",")] if sd.get("AD", ".") not in (".", "") else []
            dp = int(sd["DP"]) if sd.get("DP", ".").isdigit() else (sum(ad) if ad else None)
            for ai, alt in enumerate(t[4].split(","), start=1):
                n = sum(1 for a in gt if a == str(ai))
                if n == 0 or alt in ("*", ".", "<*>", "<NON_REF>"):
                    continue
                zyg = "hom" if n >= 2 or len(gt) == 1 else "het"
                ab = ad[ai] / sum(ad) if ad and len(ad) > ai and sum(ad) else None
                for p, r, a in decompose(int(t[1]), t[3], alt):
                    k = (norm_chrom(t[0]), p, r, a)
                    if k not in out or out[k][0] == "het":
                        out[k] = (zyg, dp, ab, fnum(t[5]))
    return out


def allele_counts(bam, prefix, v):
    """Reads (MAPQ >= 20, duplicates excluded) carrying the ALT and the REF allele of a variant."""
    try:
        sam = run_samtools(["view", "-F", "3844", "-q", "20", bam,
                            f"{prefix}{v.chrom}:{v.pos}-{v.pos + max(len(v.ref), 2)}"])
    except (Exception, SystemExit):
        return None
    alt = ref = 0
    for line in sam.splitlines():
        t = line.split("\t")
        if len(t) < 10 or t[5] == "*":
            continue
        s = read_alleles(int(t[3]), t[5], t[9], (v,))[v.key()]
        alt += s == "alt"
        ref += s == "ref"
    return alt, ref


class Trio:
    """Parental genotypes for the proband's variants (one or both parents; VCF and optionally BAM)."""

    MIN_READS = 15

    def __init__(self, mother_vcf=None, father_vcf=None, mother_bam=None, father_bam=None, prefix="",
                 affected="unknown"):
        self.vcf = {"mother": load_parent_vcf(mother_vcf) if mother_vcf else None,
                    "father": load_parent_vcf(father_vcf) if father_vcf else None}
        self.bam = {"mother": mother_bam, "father": father_bam}
        self.prefix = prefix
        self.affected = affected
        self._counts = {}
        self.qc = {}
        self.parentage_ok = True
        self.parentage_checked = False

    def has(self, who):
        return self.vcf[who] is not None or bool(self.bam[who])

    def present(self):
        return [w for w in ("mother", "father") if self.has(w)]

    def counts(self, who, v):
        k = (who, v.key())
        if k not in self._counts:
            self._counts[k] = allele_counts(self.bam[who], self.prefix, v) if self.bam[who] else None
        return self._counts[k]

    def status(self, who, v, use_reads=True):
        """(state, text): state het | hom | present (reads carry it, not called) | ref (reads) | low (too few reads) |
        absent (not called, coverage unknown) | None (no data for this parent)."""
        if not self.has(who):
            return None, ""
        call = self.vcf[who].get((v.chrom, v.pos, v.ref, v.alt)) if self.vcf[who] is not None else None
        if call:
            zyg, dp, ab, _q = call
            return zyg, f"{zyg}" + (f" (DP {dp}" + (f", alt {ab:.2f}" if ab is not None else "") + ")" if dp else "")
        c = self.counts(who, v) if use_reads else None
        if c:
            alt, ref = c
            n = alt + ref
            if alt >= 2 and alt >= 0.1 * n:
                return "present", f"not called, but {alt} of {n} reads carry it"
            if n < self.MIN_READS:
                return "low", f"not called, only {n} reads"
            return "ref", f"reference, {n} reads" + (f" ({alt} alt read)" if alt else "")
        return "absent", "not called (coverage unknown)"

    def origin(self, v, sex, use_reads=True):
        """(label, text) for a proband variant. Labels: maternal, paternal, both, de novo, possible de novo,
        one parent only, unknown."""
        st = {w: self.status(w, v, use_reads) for w in ("mother", "father")}
        carry = {w for w, (s, _t) in st.items() if s in ("het", "hom", "present")}
        txt = "; ".join(f"{w}: {t}" for w, (s, t) in st.items() if s)
        chrx = v.chrom == "X" and not any(a <= v.pos <= b for a, b in PAR_GRCH38)
        if chrx and sex == "male":  # hemizygous: from the mother
            s = st["mother"][0]
            if "mother" in carry:
                return "maternal", txt
            if s == "ref":
                return "de novo", txt
            return ("possible de novo" if s == "absent" else "unknown"), txt
        if v.zyg == "hom":
            if carry == {"mother", "father"}:
                return "both", txt
            refs = {w for w, (s, _t) in st.items() if s == "ref"}
            if carry and refs:
                return "one parent only", txt
            return "unknown", txt  # a homozygote needs both parents' alleles: one carrier alone says nothing
        if carry == {"mother", "father"}:
            return "both", txt
        if carry == {"mother"}:
            return "maternal", txt
        if carry == {"father"}:
            return "paternal", txt
        states = {s for s, _t in st.values()}
        if len(self.present()) == 2 and states == {"ref"}:
            return "de novo", txt
        if len(self.present()) == 2 and states <= {"ref", "absent", "low"}:
            return "possible de novo", txt
        return "unknown", txt

    def check(self, variants):
        """Parentage / sample-identity check on common autosomal SNVs (gnomAD AF >= 5 %, QUAL >= 30, DP >= 20): the
        share of the proband's homozygous sites where each parent carries the allele, and of heterozygous sites where
        at least one parent does. Calls only (no reads), so missed parental calls lower the rates slightly."""
        hom = {"mother": [0, 0], "father": [0, 0]}
        het = [0, 0]
        for v in variants:
            if v.chrom in ("X", "Y", "MT") or len(v.ref) != 1 or len(v.alt) != 1:
                continue
            if (v.af or 0) < 0.05 or (v.qual or 0) < 30 or (v.dp or 0) < 20:
                continue
            k = (v.chrom, v.pos, v.ref, v.alt)
            has = {w: self.vcf[w] is not None and k in self.vcf[w] for w in ("mother", "father")}
            if v.zyg == "hom":
                for w in ("mother", "father"):
                    if self.vcf[w] is not None:
                        hom[w][0] += has[w]
                        hom[w][1] += 1
            elif all(self.vcf[w] is not None for w in ("mother", "father")):
                het[0] += has["mother"] or has["father"]
                het[1] += 1
        self.qc = {w: (hom[w][0] / hom[w][1] if hom[w][1] else None, hom[w][1]) for w in ("mother", "father")}
        self.qc["het"] = (het[0] / het[1] if het[1] else None, het[1])
        rates = [r for r, n in self.qc.values() if r is not None and n >= 100]
        self.parentage_checked = bool(rates)
        self.parentage_ok = not rates or min(rates) >= 0.9
        return self.qc


def trio_phase(trio, a, b, sex):
    """Phase of two heterozygous proband variants from the parents: 'trans', 'cis' or None."""
    if trio is None or a.zyg != "het" or b.zyg != "het":
        return None
    oa, ob = trio.origin(a, sex)[0], trio.origin(b, sex)[0]
    if {oa, ob} == {"maternal", "paternal"}:
        return "trans"
    if oa == ob and oa in ("maternal", "paternal"):
        return "cis"
    return None


def de_novo_codes(o, pheno_match, pheno_given, parentage_ok=True):
    """ClinGen SVI de novo points for one proband: confirmed de novo with a phenotype consistent with the gene = 1 point
    (PS2_Moderate); assumed de novo = 0.5 (PM6_Supporting). Without phenotype terms the strength cannot be judged."""
    if not pheno_given or not pheno_match or not parentage_ok:
        return []
    return ["PS2_Moderate"] if o == "de novo" else (["PM6_Supporting"] if o == "possible de novo" else [])


# ---------------------------------------------------------------- uniparental disomy and parental origin (trio)
# (chromosome, parent who gave both copies) -> imprinting disorder (Eggermann et al., Clin Epigenetics 2015; Nakka et
# al., Am J Hum Genet 2019)
IMPRINTING_UPD = {
    ("6", "father"): "transient neonatal diabetes mellitus (6q24)",
    ("7", "mother"): "Silver-Russell syndrome",
    ("11", "father"): "Beckwith-Wiedemann syndrome (11p15)",
    ("11", "mother"): "Silver-Russell syndrome (11p15, rare)",
    ("14", "mother"): "Temple syndrome (14q32)",
    ("14", "father"): "Kagami-Ogata syndrome (14q32)",
    ("15", "mother"): "Prader-Willi syndrome",
    ("15", "father"): "Angelman syndrome",
    ("16", "mother"): "growth restriction (often with confined placental trisomy 16)",
    ("20", "father"): "pseudohypoparathyroidism type 1B (20q13)",
    ("20", "mother"): "Mulchandani-Bhoj-Conlin syndrome",
}
UPD_P = (0.3, 0.1)  # exclusion rates of a parent's alleles inside an isodisomy and a heterodisomy (common SNVs)
UPD_MIN_LOD = 10.0  # log10 likelihood ratio of the segment against the genome-wide error rate
UPD_MIN_EXCL = 5   # excluding sites in the segment
UPD_MIN_SPAN = 2_000_000  # UPD spans many genes; shorter runs are inherited deletions or paralogous loci
UPD_CLUSTER_BP = 100_000  # excluding sites within this distance of the previous one count once (HLA, IGH, deletions)
OTHER = {"mother": "father", "father": "mother"}


class Covered:
    """Targets with >= 80 % of bases at >= 20x in the proband (mosdepth thresholds): where a missing call means the
    proband carries no ALT allele."""

    def __init__(self, rows):
        by = collections.defaultdict(list)
        for c, s, e, _n, th in rows:
            if th.get("20X", 0) >= 0.8 * (e - s):
                by[c].append((s, e))
        self.iv = {c: sorted(v) for c, v in by.items()}
        self.starts = {c: [s for s, _e in v] for c, v in self.iv.items()}

    def __contains__(self, cp):
        c, p = cp
        st = self.starts.get(c)
        if not st:
            return False
        i = bisect.bisect_right(st, p - 1) - 1  # 0-based start < p (1-based position)
        return i >= 0 and self.iv[c][i][1] >= p

    def __bool__(self):
        return bool(self.iv)


def upd_chromosomes(sex):
    return [str(i) for i in range(1, 23)] + (["X"] if sex == "female" else [])


def upd_sites(variants, trio, sex, covered=None):
    """Informative common SNVs for parent-of-origin, per chromosome and parent: [(pos, excluded)] where excluded means
    the proband did not get that parent's allele (further excluding sites within UPD_CLUSTER_BP count once): (a) the proband is homozygous for an ALT allele (gnomAD AF >= 1 %)
    that the parent's calls lack; (b) the parent is homozygous for an ALT allele that the proband, covered at >= 20x
    there, does not carry. Calls with QUAL >= 30 and depth >= 20 only; X only for a female proband (outside PAR)."""
    chroms = set(upd_chromosomes(sex))
    sites = collections.defaultdict(list)
    child = {}
    for v in variants:
        if v.chrom in chroms and len(v.ref) == 1 and len(v.alt) == 1:
            child[(v.chrom, v.pos, v.ref, v.alt)] = v
    par = [(10001, 2781479), (155701383, 156030895)]
    for who in ("mother", "father"):
        calls = trio.vcf.get(who)
        if calls is None:
            continue
        for k, v in child.items():  # (a)
            if v.zyg == "hom" and (v.af or 0) >= 0.01 and (v.qual or 0) >= 30 and (v.dp or 0) >= 20:
                if k[0] == "X" and any(a <= k[1] <= b for a, b in par):
                    continue
                sites[(k[0], who)].append((k[1], k not in calls))
        if covered:
            for k, (zyg, dp, _ab, q) in calls.items():  # (b)
                if zyg != "hom" or k[0] not in chroms or len(k[2]) != 1 or len(k[3]) != 1:
                    continue
                if (q or 0) < 30 or (dp or 0) < 20 or (k[0], k[1]) not in covered:
                    continue
                if k[0] == "X" and any(a <= k[1] <= b for a, b in par):
                    continue
                sites[(k[0], who)].append((k[1], k not in child))
    for key, s in sites.items():
        s.sort()
        kept, last = [], -UPD_CLUSTER_BP
        for p, x in s:
            if x and p - last < UPD_CLUSTER_BP:
                continue  # one locus with several Mendelian errors (a deletion, a paralog) counts once
            if x:
                last = p
            kept.append((p, x))
        sites[key] = kept
    return sites


def best_segment(xs, p0, p1s=UPD_P):
    """Maximum-scoring run (Kadane) of a 0/1 exclusion series: (log10 LR, i, j) for xs[i..j], the best over the
    isodisomy and heterodisomy exclusion rates."""
    return max(_kadane(xs, p0, p1) for p1 in p1s)


def _kadane(xs, p0, p1):
    hit, miss = math.log10(p1 / p0), math.log10((1 - p1) / (1 - p0))
    best, cur, start, bi, bj = 0.0, 0.0, 0, 0, -1
    for k, x in enumerate(xs):
        s = hit if x else miss
        if cur <= 0:
            cur, start = s, k
        else:
            cur += s
        if cur > best:
            best, bi, bj = cur, start, k
    return best, bi, bj


def upd_scan(variants, trio, sex, covered=None, cnv_calls=(), roh=()):
    """Uniparental disomy from a trio: per chromosome and parent, the run of informative sites where the proband
    lacks that parent's alleles most often, against that parent's genome-wide (median per-chromosome) rate of such
    Mendelian errors. Returns (segments, rates); each segment is a dict with chrom, start, end, from (the parent who
    gave both copies), missing (the other parent), excluded, sites, lod, het (proband heterozygosity in the segment
    relative to the rest of the genome), kind (isodisomy / heterodisomy / mixed), deletion (a CNV loss call
    covering it), roh_mb."""
    if trio is None or any(trio.vcf.get(w) is None for w in ("mother", "father")):
        return [], {}
    sites = upd_sites(variants, trio, sex, covered)
    rates = {}
    for who in ("mother", "father"):
        per = [sum(x for _p, x in s) / len(s) for (c, w), s in sites.items() if w == who and len(s) >= 50]
        rates[who] = max(statistics.median(per) if per else 0.01, 0.002)
    common = collections.defaultdict(lambda: [0, 0])  # chrom -> [het, all] common SNV calls of the proband
    het_pos = collections.defaultdict(list)
    for v in variants:
        if len(v.ref) == 1 and len(v.alt) == 1 and (v.af or 0) >= 0.01 and (v.qual or 0) >= 30 and (v.dp or 0) >= 20:
            het_pos[v.chrom].append((v.pos, v.zyg == "het"))
    for c in het_pos:
        het_pos[c].sort()
    allhet = [x for c in upd_chromosomes(sex) for _p, x in het_pos.get(c, [])]
    base_het = sum(allhet) / len(allhet) if allhet else None
    segs = []
    for (c, who), s in sites.items():
        lod, i, j = best_segment([x for _p, x in s], rates[who])
        n_ex = sum(x for _p, x in s[i:j + 1]) if j >= i else 0
        start, end = (s[i][0], s[j][0]) if j >= i else (0, 0)
        if lod < UPD_MIN_LOD or n_ex < UPD_MIN_EXCL or end - start < UPD_MIN_SPAN:
            continue
        inside = [x for p, x in het_pos.get(c, []) if start <= p <= end]
        het = (sum(inside) / len(inside) / base_het) if inside and base_het else None
        roh_bp = sum(max(0, min(e, end) - max(b, start)) for rc, b, e, *_ in roh if rc == c)
        iso = roh_bp / max(end - start, 1)
        kind = ("isodisomy" if (het is not None and het < 0.15) or iso > 0.8 else
                "heterodisomy" if (het is None or het > 0.6) and iso < 0.2 else "mixed iso- and heterodisomy")
        dele = [x for x in cnv_calls if x["chrom"] == c and x["type"] == "DEL" and x["confidence"] == "high"
                and min(x["end"], end) - max(x["start"], start) > 0.5 * (end - start)]
        segs.append({"chrom": c, "start": start, "end": end, "from": OTHER[who], "missing": who, "excluded": n_ex,
                     "sites": j - i + 1, "lod": lod, "het": het, "kind": kind, "deletion": dele,
                     "roh_mb": roh_bp / 1e6, "rate": rates[who],
                     "disorder": IMPRINTING_UPD.get((c, OTHER[who]), "")})
    order = {c: n for n, c in enumerate(upd_chromosomes(sex))}
    segs.sort(key=lambda d: (order.get(d["chrom"], 99), d["start"]))
    return segs, rates


def upd_label(d):
    """'maternal isodisomy of chromosome 15 (Prader-Willi syndrome)', with the region for a segment"""
    seg = f" {d['start'] / 1e6:.1f}-{d['end'] / 1e6:.1f} Mb"
    return (f"{'maternal' if d['from'] == 'mother' else 'paternal'} {d['kind']} of chromosome {d['chrom']}{seg}"
            + (f" ({d['disorder']})" if d["disorder"] else "")
            + (" or a deletion of the " + d["missing"] + "'s copy" if d["deletion"] else ""))


def in_segments(v, segs):
    return next((d for d in segs if d["chrom"] == v.chrom and d["start"] <= v.pos <= d["end"]), None)


def gain_origin(variants, trio, chrom):
    """Parent of origin of an extra copy of a chromosome: at SNVs where one parent is homozygous ALT and the other's
    calls lack the ALT, the share of the proband's reads from the mother's allele is ~0.5 with two copies, ~0.67 with
    an extra maternal copy and ~0.33 with an extra paternal one. Returns (median maternal share, sites) or None."""
    if trio is None or any(trio.vcf.get(w) is None for w in ("mother", "father")):
        return None
    fm = []
    for v in variants:
        if v.chrom != chrom or v.zyg != "het" or v.ab is None or (v.dp or 0) < 20 or len(v.ref) != 1 or len(v.alt) != 1:
            continue
        k = (v.chrom, v.pos, v.ref, v.alt)
        m, f = trio.vcf["mother"].get(k), trio.vcf["father"].get(k)
        if m and m[0] == "hom" and not f:
            fm.append(v.ab)
        elif f and f[0] == "hom" and not m:
            fm.append(1 - v.ab)
    return (statistics.median(fm), len(fm)) if len(fm) >= 10 else None


def gain_origin_text(o):
    if o is None:
        return "parent of origin not determined (fewer than 10 informative sites)"
    share, n = o
    who = "maternal" if share >= 0.6 else ("paternal" if share <= 0.4 else "unclear (allele balance near 0.5: "
                                                                            "mosaic, or not a whole-copy gain)")
    return f"extra copy {who} (proband reads from the mother's allele {share:.2f}, {n} sites)"


def x_parents(variants, trio, covered=None):
    """Which parents' X chromosomes the proband carries (outside PAR), from sites where only one parent has the ALT
    and is homozygous (hemizygous in the father): share of such sites where the proband carries it, per parent.
    A typical female carries both (~1, ~1); a male, 45,X or maternal XXY only the mother's; paternal XXY both."""
    if trio is None or any(trio.vcf.get(w) is None for w in ("mother", "father")):
        return {}
    child = {(v.chrom, v.pos, v.ref, v.alt) for v in variants if v.chrom == "X"}
    par = [(10001, 2781479), (155701383, 156030895)]
    out = {}
    for who in ("mother", "father"):
        got = n = 0
        for k, (zyg, dp, _ab, q) in trio.vcf[who].items():
            if k[0] != "X" or zyg != "hom" or len(k[2]) != 1 or len(k[3]) != 1 or (q or 0) < 30 or (dp or 0) < 10:
                continue
            if any(a <= k[1] <= b for a, b in par) or k in trio.vcf[OTHER[who]]:
                continue
            if k in child:
                got += 1
            elif covered and ("X", k[1]) in covered:
                pass
            else:
                continue
            n += 1
        out[who] = (got / n if n else None, n)
    return out


def x_parents_text(xp):
    parts = []
    for who in ("mother", "father"):
        r, n = xp.get(who, (None, 0))
        if r is not None and n >= 10:
            parts.append(f"{'an X from the ' + who if r >= 0.5 else 'no X from the ' + who} ({r:.0%} of {n} sites)")
    return "; ".join(parts) if parts else "not determined (too few informative X sites)"


CLUSTER_TEXT = "3 or more parent-absent changes within 10 bp, typical of reads from a paralogous copy"
NOREF_TEXT = ("no reference reads in a heterozygous call: the other allele is another variant at this locus (a "
              "repeat or complex site the parents' calls may write differently)")


def mosaic_p(v):
    """One-sided binomial P(alt reads <= observed | heterozygous, 50 %) for a heterozygous call, or None."""
    if v.zyg != "het" or v.ab is None or not v.dp:
        return None
    n, k = int(v.dp), int(round(v.ab * v.dp))
    return sum(math.comb(n, i) for i in range(k + 1)) / 2 ** n


def denovo_clusters(screened, span=10, n=3):
    """Variants of the de novo screen (position-sorted) with >= n distinct positions within +/- span bp."""
    out = []
    for i, v in enumerate(screened):
        near = {w.pos for w in screened[max(0, i - 2 * n):i + 2 * n + 1] if w.chrom == v.chrom and abs(w.pos - v.pos) <= span}
        if len(near) >= n:
            out.append(v)
    return out


ORIGIN_TEXT = {"maternal": "inherited from the mother", "paternal": "inherited from the father",
               "both": "both parents carry it", "de novo": "de novo (absent from both parents' reads)",
               "possible de novo": "possibly de novo (not called in either parent; reads missing, too few or not conclusive)",
               "one parent only": "only one parent carries it: deletion of the other allele or uniparental disomy",
               "unknown": "parental origin unknown"}


# ---- copy-number variants (wes_cnv calls from exome read depth)
# immunoglobulin and T-cell receptor loci (GRCh38): somatic V(D)J rearrangement in lymphocytes deletes parts of them,
# so blood- or cell-line-derived DNA can show deletions there that are not germline
SOMATIC_REARRANGING = [("14", 105_586_437, 106_879_844, "IGH"), ("2", 88_857_361, 90_235_368, "IGK"),
                       ("22", 22_026_076, 22_922_913, "IGL"), ("14", 21_621_904, 22_552_132, "TRA/TRD"),
                       ("7", 142_299_011, 142_813_287, "TRB"), ("7", 38_240_024, 38_368_055, "TRG")]
ANEUPLOIDY_NAMES = {"21": "Down syndrome", "18": "Edwards syndrome", "13": "Patau syndrome", "XXY": "Klinefelter syndrome",
                    "X": "Turner syndrome (45,X) or a mosaic", "XXX": "triple X syndrome", "XYY": "47,XYY",
                    "XXYY": "48,XXYY"}
CNV_CLASSES = [(0.99, "pathogenic"), (0.90, "likely pathogenic"), (-0.89, "uncertain significance"),
               (-0.98, "likely benign"), (-99, "benign")]


def load_cnv(path):
    """wes_cnv calls TSV -> (calls, status). status is None when analysed, else the reason it was not."""
    if not path:
        return [], None
    calls, status, cols = [], None, None
    with opener(path) as f:
        for line in f:
            if line.startswith("##"):
                if "not analysed" in line:
                    status = line.split("; ", 1)[1].strip()
                continue
            t = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                cols = [x.lstrip("#") for x in t]
                continue
            if not cols or len(t) < len(cols) - 1:
                continue
            c = dict(zip(cols, t + [""] * (len(cols) - len(t))))
            for k in ("start", "end", "length_bp", "copy_number", "normal_copies", "targets", "panel_carriers"):
                c[k] = int(c[k])
            for k in ("mean_ratio", "log10_bf", "segdup_fraction"):
                c[k] = float(c[k])
            c["chrom"] = norm_chrom(c["chrom"])
            c["gene_exons"] = []
            for g in filter(None, c.get("genes", "").split(",")):
                name, _, frac = g.rpartition(":")
                k, _, n = frac.partition("/")
                if name and k.isdigit() and n.isdigit():
                    c["gene_exons"].append((name, int(k), int(n)))
                else:
                    c["gene_exons"].append((g, 0, 0))
            calls.append(c)
    return calls, status


def load_cnv_summary(path):
    if not path:
        return {}
    try:
        return json.load(opener(path))
    except (OSError, ValueError):
        return {}


def load_clingen_regions(path):
    """ClinGen region curation list (ClinGen_region_curation_list_GRCh38.tsv): recurrent CNV regions with their
    haploinsufficiency and triplosensitivity scores, including population regions scored 40 (dosage sensitivity
    unlikely)."""
    out = []
    if not path:
        return out
    with opener(path) as f:
        ix = None
        for line in f:
            t = line.rstrip("\n").split("\t")
            if line.startswith("#ISCA ID") or (ix is None and t[0] == "ISCA ID"):
                ix = {h.lstrip("#"): i for i, h in enumerate(t)}
                continue
            if line.startswith("#") or ix is None or len(t) < 5:
                continue
            m = re.match(r"(?:chr)?(\w+):(\d+)-(\d+)", t[ix["Genomic Location"]].replace(",", "").strip())
            if not m:
                continue

            def score(col):
                x = t[ix[col]].strip() if col in ix and ix[col] < len(t) else ""
                return int(x) if x.lstrip("-").isdigit() else None
            out.append({"id": t[ix["ISCA ID"]], "name": t[ix["ISCA Region Name"]], "cyto": t[ix.get("cytoBand", 2)],
                        "chrom": norm_chrom(m.group(1)), "start": int(m.group(2)), "end": int(m.group(3)),
                        "hi": score("Haploinsufficiency Score"), "ts": score("Triplosensitivity Score")})
    return out


def cnv_class(total):
    for thr, name in CNV_CLASSES:
        if total >= thr:
            return name
    return "benign"


def cnv_loc(c):
    size = c["length_bp"]
    sz = f"{size / 1e6:.2f} Mb" if size >= 1e6 else (f"{size / 1e3:.1f} kb" if size >= 1e3 else f"{size} bp")
    return f"{c['chrom']}:{c['start'] + 1:,}-{c['end']:,} ({sz})"


def load_cnv_targets(path):
    """wes_cnv per-target TSV -> {chrom: [(start, end, ratio, call_copy), ...]} for the targets used."""
    out = collections.defaultdict(list)
    if not path:
        return None
    with opener(path) as f:
        head = f.readline().rstrip("\n").split("\t")
        ix = {h: i for i, h in enumerate(head)}
        if "ratio" not in ix:
            return None
        for line in f:
            t = line.rstrip("\n").split("\t")
            if t[ix["used"]] != "1" or t[ix["ratio"]] == "NA":
                continue
            cc = t[ix["call_copy"]]
            out[norm_chrom(t[0])].append((int(t[1]), int(t[2]), float(t[ix["ratio"]]),
                                          int(cc) if cc.lstrip("-").isdigit() else None))
    for c in out:
        out[c].sort()
    return out if out else None


def cnv_parent_depth(c, targets):
    """A parent's read depth over the targets of the child's CNV: 'present', 'absent' or 'unclear', with the mean
    ratio. Present: the parent's own call covers half of the targets, the mean ratio moves the same way by >= 0.25,
    or >= 2 consecutive targets each move the same way by >= 0.2 (a parent's CNV can be shorter than the child's
    call); absent: >= 3 targets, no such run and the mean ratio within 0.15 of normal."""
    rows = [t for t in targets.get(c["chrom"], ()) if t[0] < c["end"] and c["start"] < t[1]]
    if not rows:
        return "unclear", None
    r = sum(x[2] for x in rows) / len(rows)
    loss = c["type"] == "DEL"
    same = sum(1 for x in rows if x[3] is not None and ((x[3] < 2) if loss else (x[3] > 2)))
    if c["chrom"] == "X":  # copy numbers differ between the sexes on X; the ratios are on one scale: use them only
        same = 0
    run = best = 0
    for x in rows:
        run = run + 1 if ((x[2] <= 0.8) if loss else (x[2] >= 1.2)) else 0
        best = max(best, run)
    shift = (r <= 0.75) if loss else (r >= 1.25)
    if same >= 0.5 * len(rows) or shift or best >= 2:
        return "present", r
    if len(rows) >= 3 and abs(r - 1) < 0.15:
        return "absent", r
    return "unclear", r


def cnv_parent_status(c, parents):
    """parents: {'mother': per-target ratios or None, ...} -> {parent: (status, ratio)} for the parents given."""
    return {w: cnv_parent_depth(c, t) for w, t in parents.items() if t is not None}


def cnv_evidence(c, K, constraint, regions, origin):
    """ACMG/ClinGen CNV points (Riggs et al., Genet Med 2020), the subset that exome read depth can support.
    origin: 'de novo', 'inherited (unaffected)', 'inherited (affected)', other or None. Returns (codes, notes):
    codes = [(code, points, text)]."""
    loss = c["type"] == "DEL"
    codes, notes = [], []
    key = "hi" if loss else "ts"
    for r in regions:
        if r["chrom"] != c["chrom"] or r["end"] <= c["start"] or c["end"] <= r["start"]:
            continue
        ov = (min(r["end"], c["end"]) - max(r["start"], c["start"])) / max(r["end"] - r["start"], 1)
        sc = r[key]
        if sc == 3 and ov >= 0.8:
            codes.append(("2A", 1.0, f"covers the established {'haploinsufficient' if loss else 'triplosensitive'} "
                                     f"region {r['name']} ({r['cyto']})"))
        elif sc in (3, 2) and ov > 0:
            notes.append(f"overlaps {max(ov, 0.01):.0%} of the ClinGen region {r['name']} ({r['cyto']}; "
                         f"{'HI' if loss else 'TS'} score {sc}): breakpoints decide, review")
        elif sc == 40 and r["start"] <= c["start"] and c["end"] <= r["end"]:
            codes.append(("2F", -1.0, f"inside the ClinGen benign population region {r['name']} "
                                      "(dosage sensitivity unlikely)"))
    names = [g for g, _, _ in c["gene_exons"]]
    if not any(code == "2A" for code, _, _ in codes):
        for g, k, n in c["gene_exons"]:
            hi, ts = K.dosage.get(g, (None, None))
            sc = hi if loss else ts
            whole = n and k >= n
            if sc == 3 and whole:
                codes.append(("2A", 1.0, f"all {n} coding exons of {g}, established "
                                         f"{'haploinsufficient' if loss else 'triplosensitive'} gene (ClinGen 3)"))
                break
            if loss and sc == 3 and k:
                codes.append(("2E", 0.45, f"{k} of {n} coding exons of {g} (ClinGen haploinsufficient): an intragenic "
                                          "deletion; strength 0.45 applied, the exons and reading frame decide (PVS1)"))
                break
            if not loss and k and not whole and hi == 3:
                codes.append(("2I", 0.45, f"intragenic duplication of {k} of {n} exons of {g} (ClinGen haploinsufficient):"
                                          " a tandem duplication can disrupt the gene; position and orientation unknown"))
                break
        else:
            if loss and any((constraint.get(g, (None,))[0] or 9) < 0.35 for g in names):
                g = min(names, key=lambda g: constraint.get(g, (9,))[0] or 9)
                codes.append(("2H", 0.15, f"{g} is loss-of-function intolerant (gnomAD LOEUF "
                                          f"{constraint[g][0]:.2f})"))
    n_genes = len(names)
    lim = (25, 35) if loss else (35, 50)
    if n_genes >= lim[1]:
        codes.append(("3C", 0.9, f"{n_genes} protein-coding genes"))
    elif n_genes >= lim[0]:
        codes.append(("3B", 0.45, f"{n_genes} protein-coding genes"))
    if origin == "de novo":
        codes.append(("5A", 0.45 if any(K.pheno.gene_terms.get(g) for g in names) else 0.15,
                      "de novo (absent from both parents' read depth); points depend on the phenotype match"))
    elif origin == "inherited (unaffected)":
        codes.append(("5B", -0.3, "inherited from a parent reported unaffected"))
    elif origin == "inherited (affected)":
        notes.append("inherited from an affected parent: segregates with the condition (5D needs more family members)")
    return codes, notes


def build(args):
    sample, contigs, variants, prefix = parse_vcf(args.vcf)
    sample = args.sample_name or sample or "sample"
    grch38 = contigs.get("1") == 248956422 if contigs else None
    gencc = load_gencc(args.gencc)
    constraint = load_constraint(args.constraint)
    hpo_terms = {t.strip().upper() for t in re.split(r"[,;\s]+", args.hpo_terms or "") if t.strip()}
    hpo_terms = {t if t.startswith("HP:") else "HP:" + t for t in hpo_terms if re.match(r"(HP:)?\d{7}$", t)}
    pheno = load_hpo(args.hpo_genes, hpo_terms)
    hpo_match, hpo_names = pheno.gene_terms, pheno.names
    hpoa = load_hpoa(args.hpo_annotations)
    dosage = load_dosage(args.clingen_dosage)
    cvp = load_clinvar_protein(args.clinvar_protein)
    K = Knowledge(gencc, load_gencc_disputed(args.gencc), pheno, hpoa, dosage)
    roh = load_roh(args.roh)
    target_cov = load_thresholds(args.target_coverage)
    gene_cov_rows = load_thresholds(args.gene_coverage)
    msum = load_mosdepth_summary(args.coverage_summary)
    qc = sample_qc(variants, contigs, msum, target_cov)
    sex = qc["sex"]
    trio = None
    if args.mother_vcf or args.father_vcf or args.mother_bam or args.father_bam:
        trio = Trio(args.mother_vcf, args.father_vcf, args.mother_bam, args.father_bam, prefix, args.affected_parents)
        if trio.vcf["mother"] is not None or trio.vcf["father"] is not None:
            trio.check(variants)

    doubt = {}  # id -> why a parent-absent variant is not taken as de novo (filled after the artefact flags)

    def origin(v):
        if not trio:
            return None, ""
        o, txt = trio.origin(v, sex)
        if id(v) in doubt and o in ("de novo", "possible de novo"):
            return "possible de novo", doubt[id(v)] + "; " + txt
        return o, txt

    def dn_codes(v, o, g):
        if not trio or id(v) in doubt:
            return []
        return de_novo_codes(o, bool(hpo_match.get(g)), bool(hpo_terms), trio.parentage_ok)

    def unaffected(who):
        return args.affected_parents == "none" or (args.affected_parents in ("mother", "father")
                                                   and args.affected_parents != who)

    def affected(who):
        return args.affected_parents in (who, "both")

    def origin_note(v, dominant):
        """Text on the parental origin of a finding, judged against the parents' affection status."""
        o, _t = origin(v)
        if not o:
            return ""
        who = {"maternal": "mother", "paternal": "father"}.get(o)
        txt = ORIGIN_TEXT[o]
        if who and dominant and unaffected(who):
            txt += (f" (reported unaffected): reduced penetrance or variable expressivity, a milder phenotype in the "
                    f"{who}, or not the cause")
        elif who and affected(who):
            txt += " (affected): segregates with the condition in this meiosis"
        return txt

    def gt_cell(v):
        s = call_cell(v, sex)
        if trio:
            o, txt = origin(v)
            kind = "t-red" if o == "de novo" else ("t-amb" if o in ("possible de novo", "one parent only") else "t-blu")
            s += "<br>" + tag(ORIGIN_TEXT[o], kind) + (f"<br><span class='small'>{esc(txt)}</span>" if txt else "")
        return s

    def with_origin(v, note, dominant=False):
        t = origin_note(v, dominant)
        return f"{note}; {t}" if t and note else (t or note)
    artefacts = load_artefacts(args.artefacts)
    n_artefact = 0
    for v in variants:
        v.flags = variant_flags(v)
        v.artefact = artefacts.get((v.chrom, v.pos, v.ref, v.alt))
        if v.artefact:
            n_artefact += 1
            n, total = v.artefact
            seen = f"{n}/{total}" if n and total else "several"
            v.flags.append(f"recurrent artefact: called in {seen} of the laboratory's samples, rare or absent in gnomAD")
    by_pos = {(v.chrom, v.pos, v.alt): v for v in variants}

    # de novo screen (trio): good-quality rare variants absent from both parents' calls; >= 3 of them within 10 bp are
    # reads of a paralogous copy (or a CNV) mapped here, not a de novo event (de novo multinucleotide variants are
    # almost always 1-2 bp), so they lose the de novo label
    screened = []
    if trio and len(trio.present()) == 2:
        for v in variants:
            hemi = v.chrom == "X" and sex == "male" and not any(a <= v.pos <= b for a, b in PAR_GRCH38)
            if v.chrom in ("Y", "MT") or (v.zyg != "het" and not hemi):
                continue
            if (freq(v) is not None and freq(v) >= 1e-4) or (v.qual or 0) < 30 or (v.dp or 0) < 20 or v.artefact:
                continue
            if v.zyg == "het" and v.ab is not None and v.ab < 0.25:
                continue
            if trio.origin(v, sex, use_reads=False)[0] == "possible de novo":
                screened.append(v)
    doubt.update((id(v), CLUSTER_TEXT) for v in denovo_clusters(screened))
    doubt.update((id(v), NOREF_TEXT) for v in screened if v.zyg == "het" and v.ab is not None and v.ab >= 0.9
                 and id(v) not in doubt)

    v_chrom_of_gene = {v.gene: v.chrom for v in variants}

    phase_cache = {}

    def phase(a, b):
        """Read-backed phase of two heterozygous variants, or None without aligned reads."""
        if not args.bam or a.zyg != "het" or b.zyg != "het":
            return None
        k = tuple(sorted((a.key(), b.key())))
        if k not in phase_cache:
            phase_cache[k] = phase_pair(args.bam, prefix, a, b)
        return phase_cache[k]

    def second_allele(v, others):
        """(may a second allele be in trans?, text). others: the other qualifying variants in the gene."""
        others = [x for x in others if x is not v]
        if not others:
            return False, ""
        if v.zyg == "hom" or any(x.zyg == "hom" for x in others):
            return True, ""
        tp = [(x, trio_phase(trio, v, x, sex)) for x in others]
        if any(p == "trans" for _x, p in tp):
            x = next(x for x, p in tp if p == "trans")
            return True, f"with {x.hgvs_p or x.hgvs_c} in trans (one inherited from each parent)"
        if all(p == "cis" for _x, p in tp):
            x = tp[0][0]
            return False, f"with {x.hgvs_p or x.hgvs_c} in cis (both inherited from the same parent): one allele only"
        phases = [(x, phase(v, x)) for x in others]
        trans = [(x, p) for x, p in phases if p and p[0] == "trans"]
        if trans:
            x, p = trans[0]
            return True, f"with {x.hgvs_p or x.hgvs_c} {phase_text(p)}"
        if all(p and p[0] == "cis" for _x, p in phases):
            x, p = phases[0]
            return False, f"with {x.hgvs_p or x.hgvs_c} {phase_text(p)}: one allele only"
        p = next((p for _x, p in phases if p), None)
        return True, phase_text(p) if p else "phase unknown"

    def moi_of(g):
        m = gene_moi(gencc.get(g, []))
        if not m and g in ACMG_SF:
            m = {"AR"} if ACMG_SF[g] in ("biallelic", "hfe") else ({"XL"} if ACMG_SF[g] == "xl" else {"AD"})
        if not m and v_chrom_of_gene.get(g) == "X":
            m = {"XL"}
        return m

    def in_roh(v):
        return any(c == v.chrom and s <= v.pos <= e for c, s, e, *_ in roh)

    R = Report()
    today = datetime.date.today().isoformat()

    # ---------------------------------------------------------------- diagnostic (ClinVar P/LP)
    plp = [v for v in variants if is_plp(v)]
    diag, carrier, lowpen, review = [], [], [], []
    carrier_note = {}
    by_gene_plp = collections.defaultdict(list)
    for v in plp:
        by_gene_plp[v.gene].append(v)
    for v in plp:
        g = v.gene
        moi = moi_of(g)
        nonpar_x = v.chrom == "X" and not any(a <= v.pos <= b for a, b in PAR_GRCH38)
        second, phase_note = second_allele(v, by_gene_plp[g])
        biallelic = v.zyg == "hom" or second
        two = ("homozygous" if v.zyg == "hom" else
               "two P/LP alleles in this gene (possible compound heterozygote" + (f", {phase_note}" if phase_note else "")
               + ")")
        if v.zyg == "hom" and v.nhomalt and v.nhomalt >= 2:
            two += (f" ({int(v.nhomalt)} homozygotes in gnomAD: the condition is probably mild, late-onset or of "
                    "reduced penetrance)")
        dominant_only = bool(moi) and "AR" not in moi
        if is_low_penetrance_only(v):
            lowpen.append(v)
        elif freq(v) is not None and (freq(v) > 0.05 or (dominant_only and freq(v) > 0.001)):
            review.append((v, f"ClinVar P/LP, but gnomAD frequency {fmt_af(freq(v))} is too common for a highly penetrant "
                              + ("dominant disorder" if dominant_only else "disorder")
                              + ": a low-penetrance/risk allele, an outdated assertion or a mapping artefact"))
        elif not moi:
            review.append((v, "no established gene-disease association or mode of inheritance (GenCC)"))
        elif "AR" in moi and not ({"AD", "XL"} & moi):
            if biallelic:
                diag.append((v, "autosomal recessive: " + two))
            else:
                carrier.append(v)
                if phase_note:
                    carrier_note[v.key()] = phase_note
        elif "XL" in moi and nonpar_x:
            if sex == "male" or v.zyg == "hom":
                diag.append((v, "X-linked: " + ("hemizygous" if sex == "male" else "homozygous")))
            else:
                review.append((v, "X-linked gene, heterozygous female: usually a carrier of an X-linked recessive "
                                  "condition; affected if the condition is X-linked dominant or with skewed X-inactivation"))
        elif "AR" not in moi:
            diag.append((v, "autosomal dominant: one allele is sufficient (penetrance may be incomplete)"))
        elif biallelic:
            diag.append((v, "gene with dominant and recessive disorders: " + two))
        else:
            review.append((v, "gene with both dominant and recessive disorders: a heterozygote is affected only by "
                              "the dominant one; check which condition the ClinVar assertion refers to"))

    # ---------------------------------------------------------------- automated ACMG classification
    # repeat context of the rare in-frame changes (PM4 outside repeats, BP3 inside)
    inframe = [v for v in variants if (freq(v) is None or freq(v) <= 0.05)
               and any(e in (v.anns[0][1] if v.anns else "") for e in INFRAME_EFFECTS)]
    if args.reference and inframe:
        ictx = reference_contexts(args.reference, prefix, inframe)
        for v in inframe:
            if v.key() in ictx:
                v.repeat = repeat_flag(ictx[v.key()][0], ictx[v.key()][1], len(v.ref))

    def pm3_codes(v, g, moi):
        """ClinGen SVI PM3 points (v1.0) for a recessive gene in a patient whose phenotype matches it: in trans with a
        ClinVar pathogenic variant 1.0 (likely pathogenic 0.5), phase unknown 0.5 (0.25), homozygous 0.5; 1 point is
        PM3 (Moderate), 0.5 PM3_Supporting."""
        if "AR" not in moi or not hpo_match.get(g):
            return []
        best = 0.5 if v.zyg == "hom" and not (v.nhomalt or 0) >= 2 else 0.0
        for x in by_gene_plp.get(g, []):
            if x is v or x.key() == v.key() or x.zyg == "hom":
                continue
            full = "pathogenic" in clin_terms(x)
            ph = trio_phase(trio, v, x, sex)
            if ph is None and phase(v, x):
                ph = phase(v, x)[0]
            pts = (1.0 if full else 0.5) if ph == "trans" else (0.0 if ph == "cis" else (0.5 if full else 0.25))
            best = max(best, pts)
        return ["PM3"] if best >= 1 else (["PM3_Supporting"] if best >= 0.5 else [])

    R.acmg = {}

    def classify(v, g, extra=()):
        """Evidence codes and class of variant v in gene g; kept for the findings table."""
        moi = moi_of(g)
        pc, pn = protein_evidence(v, g, cvp)
        ev = acmg_evidence(v, g, moi, constraint, bool(hpo_match.get(g)), args.max_af, dosage,
                           pc + list(extra) + pm3_codes(v, g, moi))
        cls, pts = acmg_classify(ev)
        R.acmg[id(v)] = (cls, pts, ev)
        return ev, pn, cls, pts

    def auto_cell(v):
        """ClinVar cell with the automated classification under it (a cross-check of the ClinVar assertion)."""
        g = v.gene
        o = origin(v)[0] if trio else None
        ev, _pn, cls, pts = classify(v, g, dn_codes(v, o, g))
        s = clinvar_cell(v) + "<br><span class='small'>automated: </span>" + acmg_tag(cls, pts)
        if "benign" in cls and "BA1" not in cls:
            s += "<br><span class='small'>disagrees with ClinVar: " + esc(",".join(ev)) + "</span>"
        return s

    # ---------------------------------------------------------------- candidates not in ClinVar
    cands = []
    rare_by_gene = collections.defaultdict(list)
    for v in variants:
        if v.clnsig and not any(t in v.clnsig.lower() for t in ("uncertain", "not_provided", "conflicting")):
            continue  # classified in ClinVar (handled elsewhere)
        if freq(v) is not None and freq(v) >= args.max_af:
            continue
        for g in genes_of(v)[:3]:
            if g in gencc or g in ACMG_SF:
                if (is_lof(v, g) or (is_missense(v, g) and hpo_match.get(g))
                        or (is_missense(v, g) and CODE_RANK.get(am_code(v.am), 0) >= 2)
                        or splice_candidate(v, g, bool(hpo_match.get(g)))
                        or {"PS1", "PM5"} & set(protein_evidence(v, g, cvp)[0])):
                    rare_by_gene[g].append(v)
                break
    for g, vs in rare_by_gene.items():
        moi = moi_of(g)
        loeuf = constraint.get(g, (None, None, None))[0]
        for v in vs:
            lof = is_lof(v, g)
            reasons, consistent = [], False
            if "AR" in moi:  # recessive model: two alleles needed
                others = [x for x in vs + by_gene_plp.get(g, []) if x is not v]
                ok, ph_note = second_allele(v, others)
                ok = v.zyg == "hom" or ok
                hom_pop = v.zyg == "hom" and (v.nhomalt or 0) >= 2
                consistent |= ok and not hom_pop
                reasons.append("recessive: " + (f"homozygous, but {int(v.nhomalt)} homozygotes in gnomAD" if hom_pop else
                                                "homozygous" if v.zyg == "hom" else
                                                f"second rare allele in the gene ({ph_note})" if ok else
                                                f"second rare allele {ph_note}" if others else
                                                "single heterozygous allele (carrier at most)"))
                if trio and v.zyg == "hom":
                    o = origin(v)[0]
                    if o == "both":
                        reasons[-1] += ", both parents heterozygous carriers"
                    elif o == "one parent only":
                        reasons[-1] += f", but {ORIGIN_TEXT[o]}"
            if not moi or {"AD", "XL", "MT"} & moi:  # dominant / X-linked model: one rare allele
                strict = freq(v) is None or freq(v) < 1e-4
                hi = dosage.get(g, (None, None))[0]
                tolerant = lof and "XL" not in moi and hi not in (2, 3) and (
                    hi == 40 or (loeuf is not None and loeuf >= 0.6))
                o = origin(v)[0] if trio else None
                who = {"maternal": "mother", "paternal": "father"}.get(o)
                from_unaffected = bool(who) and unaffected(who) and not (v.chrom == "X" and sex == "male")
                consistent |= strict and not tolerant and not from_unaffected
                r = ("dominant/X-linked" if moi else "inheritance unknown")
                if not strict:
                    r += ": too frequent for a dominant disorder"
                elif tolerant:
                    r += (": ClinGen: haploinsufficiency unlikely" if dosage.get(g, (None,))[0] == 40
                          else f": LoF-tolerant gene (LOEUF {loeuf:.2f})")
                elif from_unaffected:
                    r += f": inherited from the unaffected {who}"
                elif o in ("de novo", "possible de novo"):
                    r += f": {o}"
                elif who and affected(who):
                    r += f": inherited from the affected {who}"
                reasons.append(r)
            why = "; ".join(reasons)
            o = origin(v)[0] if trio else None
            ev, pnotes, cls, pts = classify(v, g, dn_codes(v, o, g))
            if id(v) in doubt:
                o = None  # no de novo score
            best = pheno.best_diseases(g, 1)
            score = (2 if hpo_match.get(g) else 0) + (1 if best and best[0][0] >= 0.5 else 0) + \
                (2 if consistent else 0) + (1 if lof else 0) + \
                (1 if "PVS1" in ev else 0) + (1 if CODE_RANK.get(next((c for c in ev if c.startswith("PP3")), ""), 0) >= 2 else 0) + \
                (2 if "PS1" in ev else (1 if "PM5" in ev else 0)) + \
                (3 if o == "de novo" and (not moi or {"AD", "XL"} & moi) else (1 if o == "possible de novo" else 0)) - \
                (1 if v.flags and any("depth" in f or "balance" in f for f in v.flags) else 0)
            cands.append((score, v, g, consistent, why, ev, pnotes, cls, pts))
    cands.sort(key=lambda x: (-x[0], -x[8], x[2], x[1].pos))

    # ---------------------------------------------------------------- de novo variants (trio)
    denovo, n_denovo_any, n_denovo_cluster = [], len(screened), 0
    for v in screened:
        if not (v.impact in ("HIGH", "MODERATE") or (v.sai or 0) >= 0.2):
            continue
        o, txt = origin(v)
        if o not in ("de novo", "possible de novo"):
            continue
        clustered = id(v) in doubt
        mp = mosaic_p(v)
        if not clustered and mp is not None and mp < 1e-3:
            txt = (f"alt fraction {v.ab:.2f} of {v.dp} reads, below a heterozygote's 0.5 (p = {mp:.0e}): possibly "
                   "post-zygotic (mosaic), or acquired in blood or cell culture; " + txt)
        n_denovo_cluster += clustered
        g = next((x for x in genes_of(v)[:3] if x in gencc), v.gene)
        moi = moi_of(g)
        ev, pn, cls, pts = classify(v, g, dn_codes(v, o, g))
        rank = -10 if clustered else (2 if o == "de novo" else 0) + (
            2 if g in gencc and (not moi or {"AD", "XL"} & moi) else (1 if g in gencc else 0)) +             (2 if hpo_match.get(g) else 0) + (1 if v.impact == "HIGH" else 0)
        denovo.append((rank, v, g, o, txt, ev, pn, cls, pts))
    denovo.sort(key=lambda x: (-x[0], x[2], x[1].pos))

    # ---------------------------------------------------------------- ACMG SF
    sf = []
    for v in variants:
        for g in genes_of(v)[:2]:
            if g not in ACMG_SF:
                continue
            rule = ACMG_SF[g]
            plp_v = is_plp(v) and not is_low_penetrance_only(v)
            if rule == "hfe":
                if v.pos == 26092913 and v.alt == "A" and v.zyg == "hom" and grch38:
                    sf.append((v, g, "reportable", "HFE p.Cys282Tyr homozygous"))
                continue
            if rule == "ttn" and plp_v and not is_lof(v, g):
                continue
            if plp_v:
                if rule in ("biallelic", "xl"):
                    hemi = rule == "xl" and sex == "male"
                    n = sum(2 if (x.zyg == "hom" or hemi) else 1 for x in variants if is_plp(x) and g in genes_of(x)[:2])
                    sf.append((v, g, "reportable" if n >= 2 else "not reportable (single allele: carrier)",
                               "hemizygous, homozygous or two heterozygous P/LP variants only" if rule == "xl"
                               else "biallelic only"))
                else:
                    sf.append((v, g, "reportable", "P/LP in ClinVar"))
            elif (not v.clnsig or "uncertain" in v.clnsig.lower()) and is_lof(v, g) and (freq(v) is None or freq(v) < 1e-4):
                sf.append((v, g, "review", "rare predicted loss-of-function, not classified in ClinVar"))
            break

    # ---------------------------------------------------------------- pharmacogenomics
    tc_index = collections.defaultdict(list)
    for r in target_cov:
        tc_index[r[0]].append(r)

    def target_fraction_20x(chrom, pos):
        for c, s, e, _n, th in tc_index.get(chrom, []):
            if s < pos <= e:
                return th.get("20X", 0) / (e - s) if e > s else 0.0
        return None  # off-target

    fixed_sites = [(a[3], a[4]) for a in PGX_ALLELES] + [(a[3], a[4]) for a in RISK_ALLELES]
    depth_at = site_depths(args.bam, fixed_sites, prefix) if (args.bam and grch38) else {}

    def lookup(chrom, pos, ref, alt):
        """(alt allele copies or None if not assessable, variant, evidence text)."""
        v = by_pos.get((chrom, pos, alt))
        if depth_at:
            dp = depth_at.get((chrom, pos), 0)
            if v:
                return (2 if v.zyg == "hom" else 1), v, f"called, {dp}x"
            if dp >= args.min_ref_depth:
                return 0, None, f"reference, {dp}x"
            return None, None, f"only {dp}x"
        if target_cov:
            frac = target_fraction_20x(chrom, pos)
            txt = "off-target" if frac is None else f"target {frac:.0%} at 20x"
            if v:
                return (2 if v.zyg == "hom" else 1), v, f"called, DP {v.dp}"
            return (0 if frac is not None and frac >= 0.9 else None), None, txt
        if v:
            return (2 if v.zyg == "hom" else 1), v, f"called, DP {v.dp}"
        return 0, None, "no call; coverage unknown (assumed reference)"

    pgx = collections.OrderedDict()
    for gene, name, rs, chrom, pos, ref, alt, func in PGX_ALLELES:
        n, v, cov = lookup(chrom, pos, ref, alt) if grch38 is not False else (None, None, None)
        pgx.setdefault(gene, []).append((name, rs, func, n, v, cov))
    risk = []
    for gene, name, rs, chrom, pos, ref, alt, text in RISK_ALLELES:
        n, v, cov = lookup(chrom, pos, ref, alt) if grch38 is not False else (None, None, None)
        risk.append((gene, name, rs, n, v, cov, text))

    clinvar_drug = [v for v in variants if "drug_response" in v.clnsig.lower()]
    other_assoc = [v for v in variants if v.clnsig and not is_plp(v) and any(
        t in v.clnsig.lower() for t in ("risk_factor", "association", "protective", "affects", "risk_allele"))
        and "drug_response" not in v.clnsig.lower()]

    # ---------------------------------------------------------------- sequence context of reported variants
    if args.reference:
        shown = {v.key(): v for v in [x[0] for x in diag] + [x[0] for x in review] + carrier + lowpen
                 + [c[1] for c in cands] + [x[0] for x in sf] + [d[1] for d in denovo]}
        ctx = reference_contexts(args.reference, prefix, list(shown.values()))
        for k, v in shown.items():
            if k in ctx:
                seq, off = ctx[k]
                if seq[off:off + len(v.ref)] != v.ref.upper():
                    v.flags.append("REF does not match reference")
                f = repeat_flag(seq, off, len(v.ref))
                if f and f not in v.flags:
                    v.flags.append(f)

    # confident calls first in every table (low-confidence calls stay listed, marked)
    diag, review = conf_sorted(diag, key=lambda x: x[0]), conf_sorted(review, key=lambda x: x[0])
    carrier, lowpen = conf_sorted(carrier), conf_sorted(lowpen)
    sf = conf_sorted(sf, key=lambda x: x[0])
    cands = conf_sorted(cands, key=lambda x: x[1])
    denovo = conf_sorted(denovo, key=lambda x: x[1])

    def short(v, gene=True):
        txt = f"{v.gene} {v.hgvs_p or v.hgvs_c}" if gene else v.gene
        return esc(txt) + (" <i>(low confidence)</i>" if low_confidence(v) else "")

    # ---------------------------------------------------------------- copy-number variants
    cnv_calls, cnv_status = load_cnv(args.cnv_calls)
    cnv_summ = load_cnv_summary(args.cnv_summary)
    cnv_regions = load_clingen_regions(args.clingen_regions)
    cnv_parents = {"mother": load_cnv_targets(args.mother_cnv), "father": load_cnv_targets(args.father_cnv)}
    small_by_gene = collections.defaultdict(list)  # P/LP and candidate small variants per gene (compound het)
    for v, _w in diag + review:
        small_by_gene[v.gene].append(v)
    for v in carrier:
        small_by_gene[v.gene].append(v)
    for _s, v, g, *_r in cands:
        small_by_gene[g].append(v)
    cnv_rows = []  # (rank, call, codes, notes, total, class, fit, origin)
    for c in cnv_calls:
        loss, cn, normal = c["type"] == "DEL", c["copy_number"], c["normal_copies"]
        par = cnv_parent_status(c, cnv_parents)
        c["parents"] = par
        origin_c = None
        if par:
            present = [w for w, (st, _r) in par.items() if st == "present"]
            if present:
                origin_c = ("inherited (unaffected)" if all(unaffected(w) for w in present) else
                            ("inherited (affected)" if any(affected(w) for w in present) else "inherited"))
            elif len(par) == 2 and all(st == "absent" for st, _r in par.values()):
                # de novo CNVs are rare: a weak call absent from both parents is more often noise than a new event
                if c["log10_bf"] < 10 or c["targets"] < 5:
                    origin_c = "absent from both parents (weak call: review)"
                elif trio and not trio.parentage_ok:
                    origin_c = "absent from both parents (parentage failed)"
                else:
                    origin_c = "de novo"
            elif any(st == "unclear" for st, _r in par.values()):
                origin_c = "unclear in a parent"
            else:
                origin_c = "absent from the parent given"
        codes, notes = cnv_evidence(c, K, constraint, cnv_regions, origin_c)
        for ch, s_, e_, locus in SOMATIC_REARRANGING:
            if c["chrom"] == ch and c["start"] < e_ and s_ < c["end"]:
                notes.append(f"{locus} locus: deletions here arise somatically by V(D)J rearrangement in lymphocytes "
                             "and cell lines; not evidence of a germline CNV")
                c["confidence"] = (c["confidence"] + ", " if c["confidence"] != "high" else "low: ") + f"{locus} locus"
        total = round(sum(p for _, p, _ in codes), 2)
        fit = []  # clinical fit with the genes' inheritance
        for g, k, n in c["gene_exons"]:
            moi = gene_moi(K.gencc.get(g, []))
            if not moi:
                continue
            hemi_x = c["chrom"] == "X" and normal == 1
            if loss and cn == 0 and ("AR" in moi or (hemi_x and "XL" in moi)):
                fit.append(f"{g}: {'hemizygous' if hemi_x else 'homozygous'} loss of {k}/{n} exons in a "
                           f"{'recessive' if 'AR' in moi else 'X-linked'} gene")
            elif loss and cn == 1 and "AD" in moi and lof_mechanism(g, moi, constraint, K.dosage):
                fit.append(f"{g}: heterozygous loss in a dominant gene where loss of function causes disease")
            elif loss and cn == 1 and "AR" in moi:
                other = list({id(v): v for v in small_by_gene.get(g, [])}.values())
                inside_hom = [v for v in other if v.chrom == c["chrom"] and c["start"] < v.pos <= c["end"] and v.zyg == "hom"]
                if inside_hom:
                    fit.append(f"{g}: {', '.join(short(v) for v in inside_hom)} appears homozygous inside the "
                               "heterozygous deletion, so it is hemizygous: in trans with the deletion (compound "
                               "heterozygous)")
                elif other:
                    fit.append(f"{g}: heterozygous deletion plus {', '.join(short(v) for v in other[:3])}: possible "
                               "compound heterozygote (phase unknown; parents' testing decides)")
                else:
                    fit.append(f"{g}: heterozygous deletion in a recessive gene: carrier")
            elif not loss and "AD" in moi and K.dosage.get(g, (None, None))[1] in (2, 3):
                fit.append(f"{g}: duplication of a triplosensitive dominant gene")
        cls = cnv_class(total)
        rank = (c["confidence"] != "high", -total, -len(fit), -c["log10_bf"])
        cnv_rows.append((rank, c, codes, notes, total, cls, fit, origin_c))
    cnv_rows.sort(key=lambda x: x[0])
    cnv_high = [x for x in cnv_rows if x[1]["confidence"] == "high"]
    cnv_report = [x for x in cnv_high if x[5] in ("pathogenic", "likely pathogenic") or x[6]]
    cnv_sex = cnv_summ.get("sex_complement")
    cnv_aneu = cnv_summ.get("aneuploidy", [])
    covered = Covered(target_cov) if target_cov else None
    both_parents = bool(trio) and all(trio.vcf.get(w) is not None for w in ("mother", "father"))
    upd_segs, upd_rates = upd_scan(variants, trio, sex, covered, cnv_calls, roh) if both_parents else ([], {})

    def aneu_origin(text):
        """parent of origin of a reported chromosome gain or sex-chromosome aneuploidy (trio)"""
        if not both_parents:
            return ""
        m = re.match(r"chromosome (\w+): (gain|partial gain)", text)
        if m:
            return gain_origin_text(gain_origin(variants, trio, m.group(1)))
        if text.startswith("sex chromosomes") or re.match(r"chromosome X: ", text):
            return "X chromosomes: " + x_parents_text(x_parents(variants, trio, covered))
        return ""

    def aneu_name(text):
        m = re.match(r"chromosome (\w+): gain of a whole copy", text)
        if m and m.group(1) in ANEUPLOIDY_NAMES:
            return ANEUPLOIDY_NAMES[m.group(1)]
        m = re.match(r"sex chromosomes (\w+)", text)
        return ANEUPLOIDY_NAMES.get(m.group(1), "") if m else ""

    # ---------------------------------------------------------------- coverage of clinically relevant genes
    gene_cov = collections.defaultdict(lambda: [0, 0, 0])
    for c, s, e, name, th in gene_cov_rows:
        g = gene_cov[name]
        g[0] += e - s
        g[1] += th.get("20X", 0)
        g[2] += th.get("10X", 0)

    # ================================================================ HTML
    R.h(f"<!doctype html><html><head><meta charset='utf-8'><title>Exome report: {esc(sample)}</title>"
        f"<style>{CSS}</style></head><body>")
    R.h(f"<h1>Exome interpretation report</h1><p class='sub'>Sample <b>{esc(sample)}</b> &middot; generated {today} "
        f"by wes_clinical_report {VERSION} &middot; reference {'GRCh38' if grch38 else ('unknown build' if grch38 is None else 'non-GRCh38')}</p>")
    R.h("<div class='box warn'><b>Research use / decision support only.</b> Findings come from automated rules applied "
        "to ClinVar, GenCC, gnomAD and HPO data. Every finding that may be acted on must be reviewed by a qualified "
        "clinical scientist and confirmed by an orthogonal method (e.g. Sanger sequencing) in an accredited "
        "laboratory before it is reported to a patient or used for treatment decisions.</div>")

    # summary
    sf_rep = [x for x in sf if x[2] == "reportable"]
    n_pgx_act = 0
    pgx_summary = []
    for gene, alle in pgx.items():
        dip, ph = pgx_pheno(gene, alle)
        pgx_summary.append((gene, dip, ph))
        if not re.match(r"^(Normal|normal|typical|indeterminate)", ph) and dip != "indeterminate":
            n_pgx_act += 1
    R.h("<h2>Summary</h2><div class='box'><ul>")
    R.h(f"<li><b>Diagnostic findings</b> (ClinVar pathogenic/likely pathogenic, genotype consistent with inheritance): "
        f"<b>{len(diag)}</b>" + (": " + ", ".join(short(v) for v, _ in diag) if diag else "") + "</li>")
    R.h(f"<li><b>Pathogenic variants with unclear inheritance fit</b> (review): <b>{len(review)}</b>"
        + (": " + ", ".join(short(v) for v, _ in review) if review else "") + "</li>")
    n_lp = [c for c in cands if c[7] in ("pathogenic", "likely pathogenic")]
    R.h(f"<li><b>Candidate variants needing review</b> (not classified in ClinVar): <b>{sum(1 for c in cands if c[3])}</b> "
        f"inheritance-consistent of {len(cands)}; automated ACMG class pathogenic or likely pathogenic: <b>{len(n_lp)}</b>"
        + (" (" + ", ".join(short(c[1]) for c in n_lp[:6]) + ")" if n_lp else "") + "</li>")
    R.h(f"<li><b>Secondary findings ({ACMG_SF_VERSION})</b>: <b>{len(sf_rep)}</b> reportable"
        + (": " + ", ".join(short(v) for v, *_ in sf_rep) if sf_rep else "") + "</li>")
    R.h(f"<li><b>Carrier status</b> (one P/LP allele in a recessive gene): <b>{len(carrier)}</b>"
        + (": " + ", ".join(short(v, gene=False) for v in carrier) if carrier else "") + "</li>")
    R.h(f"<li><b>Pharmacogenomics</b>: {n_pgx_act} of {len(pgx_summary)} genes with a non-normal phenotype; "
        f"{len(clinvar_drug)} ClinVar drug-response variants</li>")
    R.h(f"<li><b>Risk alleles / associations</b>: {sum(1 for r in risk if r[3])} of {len(risk)} tested alleles carried; "
        f"{len(other_assoc)} other ClinVar risk/association variants</li>")
    if roh:
        big = [s for s in roh if s[3] >= 1_000_000]
        R.h(f"<li><b>Runs of homozygosity</b>: {len(big)} segments &ge; 1 Mb, total {sum(s[3] for s in big) / 1e6:.1f} Mb</li>")
    R.h(f"<li><b>Sex (genetic)</b>: {esc(sex)}" + (f"; sex chromosomes from read depth: {esc(cnv_sex)}" if cnv_sex else "")
        + "</li>")
    if args.cnv_calls:
        if cnv_status:
            R.h(f"<li><b>Copy-number variants</b>: {esc(cnv_status)}</li>")
        else:
            pl = [x for x in cnv_report if x[5] in ("pathogenic", "likely pathogenic")]
            R.h(f"<li><b>Copy-number variants</b> (exome read depth): <b>{len(cnv_report)}</b> of clinical interest "
                f"({len(pl)} pathogenic/likely pathogenic by CNV points)"
                + (": " + ", ".join(f"{x[1]['type']} {esc(','.join(g for g, _, _ in x[1]['gene_exons'][:3]))}"
                                    for x in cnv_report[:5]) if cnv_report else "")
                + f"; {len(cnv_high)} high-confidence and {len(cnv_calls) - len(cnv_high)} low-confidence calls in all</li>")
        if cnv_aneu:
            R.h("<li><b>Chromosome dosage</b>: <b>" + "; ".join(esc(a + (f" ({aneu_name(a)})" if aneu_name(a) else "")
                                                                    + (f"; {aneu_origin(a)}" if aneu_origin(a) else ""))
                                                      for a in cnv_aneu) + "</b></li>")
    if trio:
        conf = sum(1 for d in denovo if d[3] == "de novo")
        R.h(f"<li><b>Trio</b> ({' and '.join(trio.present())} given; affected parents: {esc(args.affected_parents)}): "
            + (f"<b>{conf}</b> de novo and {len(denovo) - conf} possibly de novo protein-altering variants"
               + (": " + ", ".join(short(d[1]) for d in denovo[:6]) if denovo else "")
               if len(trio.present()) == 2 else "one parent only, no de novo analysis")
            + ("" if trio.parentage_ok or not trio.qc else "; <b>parentage check failed</b>") + "</li>")
        if both_parents:
            R.h("<li><b>Uniparental disomy</b> (trio, section 9): "
                + ("<b>" + "; ".join(esc(upd_label(d)) for d in upd_segs) + "</b>" if upd_segs else "none found")
                + "</li>")
    if hpo_terms:
        R.h(f"<li><b>Phenotype</b>: {', '.join(esc(t + ' ' + hpo_names.get(t, '(unknown term)')) for t in sorted(hpo_terms))}; "
            f"{len(hpo_match)} genes associated</li>")
    else:
        R.h("<li><b>Phenotype</b>: no HPO terms given; candidates are not phenotype-ranked</li>")
    R.h("</ul></div>")

    # sample QC
    R.h("<h2>1. Sample and data quality</h2><ul class='kv'>")
    for k, label, f in (("variants", "QC-passed on-target variants", "{:,}"), ("snv", "SNVs", "{:,}"),
                        ("indel", "Indels / other", "{:,}"), ("tstv", "Ts/Tv (SNVs; ~2.4-3.3 for exomes, lower with UTR/flanking sites)", "{:.2f}"),
                        ("het_hom", "Het/hom ratio (autosomal SNVs; expectation 1.3-2.0)", "{:.2f}"),
                        ("het_ab_median", "Median alt allele fraction of het SNVs (expect ~0.5)", "{:.2f}"),
                        ("het_ab_low", "Het SNVs with alt fraction < 0.3 (contamination hint if high)", "{:.1%}"),
                        ("x_het_frac", "chrX non-PAR SNV heterozygosity", "{:.2f}"),
                        ("y_ratio", "chrY / autosomal target coverage", "{:.3f}"),
                        ("mean_target_cov", "Mean target coverage", "{:.1f}x"),
                        ("pct_10X", "Target bases >= 10x", "{:.1f}%"), ("pct_20X", "Target bases >= 20x", "{:.1f}%"),
                        ("pct_30X", "Target bases >= 30x", "{:.1f}%")):
        if qc.get(k) is not None:
            R.h(f"<li>{esc(label)}: <b>{f.format(qc[k])}</b></li>")
    R.h(f"<li>Genetic sex: <b>{esc(sex)}</b></li></ul>")
    if trio and trio.qc:
        parts = []
        for w in ("mother", "father"):
            r, n = trio.qc[w]
            if r is not None:
                parts.append(f"{w} carries the allele at {r:.1%} of {n:,} common homozygous sites")
        r, n = trio.qc["het"]
        if r is not None:
            parts.append(f"a parent carries it at {r:.1%} of {n:,} common heterozygous sites")
        R.h("<p>Trio check (common SNVs, parents' calls): " + esc("; ".join(parts)) + "."
            + ("" if trio.parentage_checked else " Too few common sites (&lt; 100) to check parentage.") + "</p>")
        if not trio.parentage_ok:
            R.h("<div class='box alert'>The parents' genotypes do not fit the proband's (expected &ge; 90%, usually "
                "&ge; 97%): sample mix-up, non-paternity/non-maternity or a parental sample of poor quality. De novo "
                "status is not used as evidence (no PS2).</div>")
    warn = []
    if qc.get("tstv") is not None and not (2.3 <= qc["tstv"] <= 3.5):
        warn.append("Ts/Tv outside the usual exome range: check for false positives or off-target calls.")
    if qc.get("het_ab_low") is not None and qc["het_ab_low"] > 0.15:
        warn.append("Many heterozygous calls with low alt fraction: possible sample contamination or mosaicism.")
    if qc.get("pct_20X") is not None and qc["pct_20X"] < 90:
        warn.append(f"Only {qc['pct_20X']:.0f}% of target bases reach 20x (clinical exomes usually require >= 90-95%): "
                    "absence of findings is not reliable in poorly covered genes (see section 10).")
    if warn:
        R.h("<div class='box warn'>" + "<br>".join(esc(w) for w in warn) + "</div>")

    # diagnostic
    R.h("<h2>2. Diagnostic findings</h2><p class='small'>ClinVar pathogenic or likely pathogenic variants (not "
        "conflicting) whose genotype fits the gene's mode of inheritance (GenCC): one allele in dominant or X-linked "
        "genes, two alleles in recessive genes. Whether a finding explains the patient's condition depends on the "
        "clinical picture.</p>")
    R.table(["Variant", "Genotype", "ClinVar", "gnomAD frequency", "Gene-disease (GenCC)", "Interpretation"],
            [[var_cell(v), gt_cell(v), auto_cell(v), freq_cell(v),
              disease_cell(v.gene, K), esc(with_origin(v, why, "autosomal dominant" in why or "X-linked" in why))]
             for v, why in diag])
    for v, why in diag:
        R.rows.append(("diagnostic", v, with_origin(v, why, "autosomal dominant" in why or "X-linked" in why)))
    R.h("<h3>Pathogenic variants whose relevance depends on the condition (review)</h3><p class='small'>Genes with "
        "both dominant and recessive disorders, X-linked genes in heterozygous females, and genes without an "
        "established disease association in GenCC.</p>")
    R.table(["Variant", "Genotype", "ClinVar", "gnomAD frequency", "Gene-disease (GenCC)", "Why review"],
            [[var_cell(v), gt_cell(v), auto_cell(v), freq_cell(v),
              disease_cell(v.gene, K), esc(with_origin(v, why))] for v, why in review])
    for v, why in review:
        R.rows.append(("pathogenic_review", v, with_origin(v, why)))
    if lowpen:
        R.h("<h3>Low-penetrance pathogenic alleles</h3>")
        R.table(["Variant", "Genotype", "ClinVar", "gnomAD frequency"],
                [[var_cell(v), gt_cell(v), clinvar_cell(v), freq_cell(v)] for v in lowpen])
        for v in lowpen:
            R.rows.append(("low_penetrance", v, ""))

    # candidates
    R.h("<h2>3. Candidate variants not classified in ClinVar</h2><p class='small'>Rare (gnomAD frequency &lt; "
        f"{args.max_af:g}) predicted loss-of-function variants in genes with an established disease association "
        "(GenCC definitive/strong/moderate/supportive)" + (", rare missense/in-frame variants in genes matching "
                                                          "the patient's phenotype" if hpo_terms else "")
        + ", rare missense variants with AlphaMissense at PP3_Moderate or stronger (>= 0.906) and rare variants "
        "predicted to alter splicing outside the canonical sites (SpliceAI >= 0.5; >= 0.2 in phenotype-matched genes)"
        + (" and rare missense variants at a residue where ClinVar has a pathogenic missense variant (PS1/PM5)" if cvp else "")
        + ". Each variant gets an automated ACMG/AMP class from the evidence codes a VCF supports, combined by the "
        "point system (pathogenic &ge; 10 points, likely pathogenic 6-9, uncertain 0-5, likely benign -1 to -6, benign "
        "&le; -7; methods in section 11). It is a draft for the analyst: segregation, functional, case-count and "
        "literature evidence are missing, and PM3/PP4/PS2 depend on the phenotype terms given. Rows are ranked by "
        "phenotype match, inheritance fit and variant consequence, then by points.</p>")
    rows = []
    for score, v, g, consistent, why, ev, pnotes, cls, pts in cands[:args.max_candidates]:
        a = ann_for_gene(v, g) or v.anns[0]
        loeuf = constraint.get(g, (None, None, None))[0]
        rows.append([var_cell(v).replace(f"<b>{esc(v.gene)}</b>", f"<b>{esc(g)}</b>"), gt_cell(v),
                     esc(a[1].replace("&", " & ")) + (f"<br><span class='small'>LOEUF {loeuf:.2f}</span>" if loeuf is not None else ""),
                     freq_cell(v), disease_cell(g, K),
                     ("" if consistent else tag("inheritance not met", "t-grn") + "<br>") + esc(why),
                     acmg_tag(cls, pts) + "<br>"
                     + " ".join(tag(c, "t-red" if c.startswith(("PVS", "PS")) else ("t-grn" if c.startswith("B") else "")) for c in ev)
                     + (f"<br><span class='small'>{esc(predictor_text(v))}</span>" if predictor_text(v) else "")
                     + "".join(f"<br><span class='small'>{esc(n)}</span>" for n in pnotes)])
        R.rows.append(("candidate" if consistent else "candidate_inheritance_not_met", v,
                       with_origin(v, "; ".join([why, ",".join(ev)] + pnotes))))
    R.table(["Variant", "Genotype", "Consequence", "gnomAD frequency", "Gene-disease (GenCC)", "Inheritance fit",
             "ACMG class (automated) and evidence"], rows)
    if len(cands) > args.max_candidates:
        R.h(f"<p class='small'>{len(cands) - args.max_candidates} more candidates in the TSV output.</p>")
        for score, v, g, consistent, why, ev, pnotes, cls, pts in cands[args.max_candidates:]:
            R.rows.append(("candidate" if consistent else "candidate_inheritance_not_met", v,
                           with_origin(v, "; ".join([why, ",".join(ev)] + pnotes))))

    # phenotype match by disease
    if hpo_terms and pheno.disease_terms:
        R.h("<h3>Diseases matching the phenotype in genes with a finding</h3><p class='small'>For every gene with a "
            "diagnostic, review or candidate variant: the linked disease (HPO annotations) that explains the largest share "
            "of the patient's phenotype, weighting each term by its information content (rare features count more than "
            "common ones). Terms the disease does not explain argue against it, or for a second diagnosis.</p>")
        fgenes = collections.OrderedDict()
        for v, _w in diag + review:
            fgenes.setdefault(v.gene, []).append(v)
        for c in cands:
            fgenes.setdefault(c[2], []).append(c[1])
        rows = []
        for g, vs in fgenes.items():
            for sc, d, terms in pheno.best_diseases(g, 1):
                if sc <= 0:
                    continue
                info = hpoa.get(d, {})
                miss = sorted(pheno.terms - set(terms))
                rows.append((sc, [f"<b>{esc(info.get('name', d))}</b><br><span class='small'>{esc(d)}"
                                  + (f"; {esc(disease_features(d, info))}" if disease_features(d, info) else "") + "</span>",
                                  f"<b>{esc(g)}</b><br><span class='small'>"
                                  + esc(", ".join(dict.fromkeys(x.hgvs_p or x.hgvs_c for x in vs))) + "</span>",
                                  f"{sc:.0%}", esc(", ".join(f"{t} {hpo_names.get(t, '')}" for t in terms)),
                                  esc(", ".join(f"{t} {hpo_names.get(t, '')}" for t in miss)) or "none"]))
                R.rows.append(("phenotype_match", None, (g, d, f"{sc:.2f} {info.get('name', '')}")))
        rows.sort(key=lambda r: -r[0])
        R.table(["Disease", "Gene and variants", "Phenotype explained", "Patient terms explained", "Not explained"],
                [r[1] for r in rows[:30]])

    # de novo variants
    if trio:
        R.h("<h3>De novo variants (trio)</h3>")
        if not trio.parentage_ok:
            R.h("<div class='box alert'><b>Parentage check failed</b> (section 1): the parents' genotypes do not fit "
                "the proband's, so most variants below are inherited from the biological parents, not de novo. Confirm "
                "the sample identities before using this table; no de novo evidence (PS2/PM6) is applied.</div>")
        if len(trio.present()) < 2:
            R.h("<p class='none'>Only one parent given: de novo variants cannot be identified.</p>")
        else:
            R.h("<p class='small'>Rare (gnomAD &lt; 0.01 %), protein-altering or predicted splice-altering variants of "
                "good quality (QUAL &ge; 30, depth &ge; 20, alt fraction &ge; 0.25) that neither parent's calls contain"
                + (", checked in the parents' reads: <b>de novo</b> when both parents have &ge; 15 reads and at most one "
                   "alt read, <b>possibly de novo</b> when a parent's reads are missing or too few"
                   if trio.bam["mother"] and trio.bam["father"] else
                   "; without the parents' reads, absence from their calls may be missing coverage (<b>possibly de "
                   "novo</b>)")
                + f". {n_denovo_any} variants of any consequence passed the quality and frequency filters and are absent "
                "from both parents' calls (an exome usually has 1-3 true coding de novo variants; many more point to "
                "missed parental calls or a sample problem). PS2 (confirmed de novo) is applied at Moderate strength "
                "when the gene matches the phenotype (ClinGen SVI de novo points), PM6_Supporting for a possible de "
                "novo variant in a phenotype-matched gene."
                + (f" {n_denovo_cluster} listed variants lie in a cluster of 3 or more parent-absent changes within "
                   "10 bp (the pattern of reads from a paralogous copy) or have no reference reads in a heterozygous "
                   "call (a repeat or complex locus the parents' calls may write differently): they are listed last, "
                   "as possibly de novo, without de novo points." if n_denovo_cluster else "") + "</p>")
            rows = []
            for rank, v, g, o, txt, ev, pn, cls, pts in denovo:
                a = ann_for_gene(v, g) or (v.anns[0] if v.anns else [""] * 11)
                rows.append([var_cell(v).replace(f"<b>{esc(v.gene)}</b>", f"<b>{esc(g)}</b>"), call_cell(v, sex),
                             esc(a[1].replace("&", " & ")), freq_cell(v), disease_cell(g, K),
                             tag(ORIGIN_TEXT[o], "t-red" if o == "de novo" else "t-amb") + f"<br><span class='small'>{esc(txt)}</span>",
                             acmg_tag(cls, pts) + "<br>"
                             + " ".join(tag(c, "t-red" if c.startswith(("PVS", "PS")) else ("t-grn" if c.startswith("B") else ""))
                                        for c in ev) + "".join(f"<br><span class='small'>{esc(n)}</span>" for n in pn)])
                R.rows.append(("de_novo" if o == "de novo" else "de_novo_possible", v,
                               "; ".join([ORIGIN_TEXT[o], txt, ",".join(ev)] + pn)))
            R.table(["Variant", "Genotype", "Consequence", "gnomAD frequency", "Gene-disease (GenCC)", "Parents",
                     "ACMG class (automated) and evidence"], rows)

    # copy-number variants and chromosome dosage
    R.h("<h2>4. Copy-number variants and chromosome dosage</h2>")
    if not args.cnv_calls:
        R.h("<p class='none'>Not analysed: no copy-number calls given (the FASTQ workflow calls them from the "
            "targets' read depth against a panel of at least 3 other samples of the same kit).</p>")
    elif cnv_status:
        R.h(f"<p class='none'>{esc(cnv_status[:1].upper() + cnv_status[1:])}.</p>")
    else:
        sm = cnv_summ
        R.h("<p class='small'>Deletions and duplications called from the read depth of the capture targets against a "
            f"panel of {len(sm.get('panel', []))} samples of the same kit"
            + (f" ({esc(', '.join(sm.get('panel', [])))})" if sm.get("panel") else "")
            + f"; {sm.get('targets_used', '?'):,} targets used, ratio noise (robust SD) {sm.get('ratio_sd_autosomes', 0):.3f}. "
            "A call is <b>high confidence</b> with at least 3 targets, mostly outside segmental duplications, not seen "
            "in the panel and log10 Bayes factor &ge; 5. In validation at 50x mean target depth with a 6-sample panel, "
            "heterozygous deletions of 5 or more targets were found 83 % of the time (10 targets: all), of 3 targets "
            "28 %; single-exon deletions are not reliably detected. A negative result does not exclude a CNV. "
            "CNV points: the ACMG/ClinGen CNV scoring (Riggs et al., Genet Med 2020) for the evidence read depth and "
            "ClinGen dosage can give (sections 2, 3 and 5); &ge; 0.99 pathogenic, 0.90-0.98 likely pathogenic. "
            "Confirm by an orthogonal method (MLPA, qPCR, chromosomal microarray) before clinical use.</p>")
        dos = sm.get("chromosome_dosage", {})
        if cnv_aneu:
            R.h("<div class='box alert'><b>Chromosome dosage</b>: " + "<br>".join(
                esc(a) + (f" &rarr; <b>{esc(aneu_name(a))}</b>" if aneu_name(a) else "")
                + (f"; {esc(aneu_origin(a))}" if aneu_origin(a) else "") for a in cnv_aneu)
                + "<br><span class='small'>From exome read depth: confirm by karyotype or chromosomal microarray. A "
                  "partial change can be mosaicism or a large CNV.</span></div>")
        elif dos:
            lo, hi = min(dos.values()), max(dos.values())
            R.h(f"<p>Chromosome dosage: no whole-chromosome gain or loss (median ratio per chromosome {lo:.2f}-{hi:.2f}); "
                f"sex chromosomes {esc(cnv_sex or '?')}.</p>")
        if cnv_sex and sex in ("male", "female") and cnv_sex in ("XX", "XY") and \
                (cnv_sex == "XY") != (sex == "male"):
            R.h(f"<div class='box warn'>Read depth gives {esc(cnv_sex)} but chrX heterozygosity gives {esc(sex)}: "
                "check for a sample swap or contamination.</div>")

        def cnv_row(x):
            _r, c, codes, notes, total, cls, fit, origin_c = x
            genes = c["gene_exons"]
            gtxt = ", ".join(f"<b>{esc(g)}</b> <span class='small'>{k}/{n}</span>" for g, k, n in genes[:12])
            if len(genes) > 12:
                gtxt += f" <span class='small'>+{len(genes) - 12} more</span>"
            dz = [g for g, _, _ in genes if K.gencc.get(g) or K.dosage.get(g, (None, None)) != (None, None)]
            know = "<br>".join(f"<b>{esc(g)}</b>: {disease_cell(g, K)}" for g in dz[:3])
            if len(dz) > 3:
                know += f"<br><span class='small'>+{len(dz) - 3} more genes with disease or dosage data</span>"
            kind = "t-red" if c["type"] == "DEL" else "t-blu"
            cell = (tag(f"{c['type']} copy number {c['copy_number']}", kind) + f"<br>{esc(cnv_loc(c))}<br>"
                    f"<span class='small'>{c['targets']} targets, ratio {c['mean_ratio']:.2f}, log10 BF "
                    f"{c['log10_bf']:.0f}</span>")
            if c["confidence"] != "high":
                cell += "<br>" + tag(c["confidence"], "t-amb")
            ev = " ".join(tag(f"{code} {p:+.2f}", "t-red" if p > 0 else "t-grn") for code, p, _ in codes)
            ev += "".join(f"<br><span class='small'>{code}: {esc(t)}</span>" for code, _, t in codes)
            ev += "".join(f"<br><span class='small'>{esc(n)}</span>" for n in notes)
            klass = tag(f"{cls} ({total:+.2f})", "t-red" if total >= 0.9 else ("t-grn" if total <= -0.9 else ""))
            par = (tag(origin_c, "t-red" if origin_c == "de novo" else "t-blu") if origin_c else
                   "<span class='small'>parents not given</span>")
            for w, (st, r) in (c.get("parents") or {}).items():
                par += f"<br><span class='small'>{w}: {st}" + (f" (ratio {r:.2f})" if r is not None else "") + "</span>"
            return [cell, gtxt or "<span class='small'>no coding exons</span>", know, "<br>".join(fit) or "",
                    par, klass + "<br>" + ev]
        hdr = ["CNV", "Genes (coding exons in the CNV / of the gene)", "Gene-disease, dosage", "Clinical fit",
               "Parents", "CNV points"]
        R.h(f"<h3>CNVs of clinical interest ({len(cnv_report)})</h3><p class='small'>High-confidence calls that are "
            "pathogenic/likely pathogenic by CNV points or fit the inheritance of a disease gene.</p>")
        R.table(hdr, [cnv_row(x) for x in cnv_report])
        other_high = [x for x in cnv_high if x not in cnv_report]
        R.h(f"<h3>Other high-confidence CNVs ({len(other_high)})</h3>")
        R.table(hdr, [cnv_row(x) for x in other_high[:40]])
        low_dz = [x for x in cnv_rows if x[1]["confidence"] != "high"
                  and any(K.gencc.get(g) for g, _, _ in x[1]["gene_exons"])]
        R.h(f"<h3>Low-confidence calls in disease genes ({len(low_dz)} of {len(cnv_rows) - len(cnv_high)} "
            "low-confidence calls)</h3><p class='small'>Segmental duplications, polymorphic regions seen in the panel, "
            "1-2 targets or weak evidence: often artefacts or common CNVs. Listed for review, not as findings.</p>")
        R.table(hdr, [cnv_row(x) for x in low_dz[:30]])
        for x in cnv_rows:
            _r, c, codes, notes, total, cls, fit, origin_c = x
            cat = ("cnv_" + cls.replace(" ", "_") if x in cnv_report else
                   ("cnv_high_confidence" if c["confidence"] == "high" else "cnv_low_confidence"))
            note = "; ".join([cnv_loc(c), f"{cls} ({total:+.2f}): " + ",".join(f"{k} {p:+.2f}" for k, p, _ in codes),
                              c["confidence"]] + [re.sub("<[^>]+>", "", f) for f in fit] + notes
                             + ([origin_c] if origin_c else []))
            R.rows.append((cat, None, (",".join(g for g, _, _ in c["gene_exons"][:20]),
                                         f"{c['type']} CN{c['copy_number']}", note,
                                         (c["chrom"], c["start"] + 1, c["end"]))))

    # ACMG SF
    R.h(f"<h2>5. Secondary findings: {ACMG_SF_VERSION} ({len(ACMG_SF)} genes)</h2><p class='small'>Medically "
        "actionable genes that the ACMG recommends reporting regardless of the indication for testing, if the patient "
        "has opted in. Gene-specific rules applied: biallelic only for MUTYH, ATP7B, BTD, GAA, RPE65, CYP27A1; ABCD1 "
        "hemizygous, homozygous or two heterozygous variants (not a single heterozygous variant in a female); HFE "
        "p.Cys282Tyr homozygotes only; TTN truncating variants only.</p>")
    if sf:
        R.table(["Variant", "Genotype", "ClinVar", "gnomAD frequency", "Category", "Status"],
                [[var_cell(v), gt_cell(v), auto_cell(v), freq_cell(v), esc(ACMG_SF_CATEGORY.get(g, "")),
                  tag(st, "t-red" if st == "reportable" else "t-amb") + "<br>" + esc(note)] for v, g, st, note in sf])
        for v, g, st, note in sf:
            R.rows.append(("acmg_sf_" + st.split()[0], v, with_origin(v, note)))
    else:
        R.h("<div class='box ok'>No reportable secondary findings: no ClinVar pathogenic/likely pathogenic variants "
            f"and no rare predicted loss-of-function variants in the {len(ACMG_SF)} {ACMG_SF_VERSION} genes.</div>")
    if gene_cov:
        low = sorted((g for g in ACMG_SF if g in gene_cov and gene_cov[g][0] and gene_cov[g][1] / gene_cov[g][0] < 0.95),
                     key=lambda g: gene_cov[g][1] / gene_cov[g][0])
        if low:
            R.h("<div class='box warn'><b>ACMG SF genes with &lt; 95% of coding bases at &ge; 20x</b> (a negative result is "
                "less reliable here): " + ", ".join(f"{esc(g)} {100 * gene_cov[g][1] / gene_cov[g][0]:.0f}%" for g in low)
                + "</div>")

    # carriers
    R.h("<h2>6. Carrier status for recessive conditions</h2><p class='small'>One pathogenic/likely pathogenic "
        "allele in a gene associated with autosomal recessive disease. Carriers are usually unaffected; the result "
        "matters for reproductive risk if a partner carries a pathogenic variant in the same gene.</p>")
    R.table(["Variant", "Genotype", "ClinVar", "gnomAD frequency", "Condition (GenCC)"],
            [[var_cell(v), gt_cell(v) + (f"<br>{tag(carrier_note[v.key()], 't-blu')}" if v.key() in carrier_note else ""),
              auto_cell(v), freq_cell(v), disease_cell(v.gene, K, phenotype=False)]
             for v in carrier])
    for v in carrier:
        R.rows.append(("carrier", v, with_origin(v, carrier_note.get(v.key(), ""))))

    # PGx
    R.h("<h2>7. Pharmacogenomics</h2><p class='small'>Single-nucleotide star alleles of CPIC genes, read at fixed "
        "GRCh38 positions. "
        + (f"Absence of a variant is read as the reference (*1) allele only where the aligned reads cover the position "
           f"with &ge; {args.min_ref_depth} reads (MAPQ &ge; 20, base quality &ge; 20)" if depth_at else
           "Absence of a variant is read as the reference (*1) allele only where the position lies in a target with "
           "&ge; 90% of bases at 20x" if target_cov else
           "<b>No alignment or coverage data given: absence of a variant is assumed to be the reference allele, which "
           "is wrong where the position was not sequenced</b>")
        + ". Phases are unknown; two variants are assumed to be on different alleles. Phenotypes follow the CPIC "
        "allele-function and phenotype tables; dosing must follow the current CPIC/DPWG guideline.</p>")
    rows = []
    for (gene, dip, ph), (_g, alle) in zip(pgx_summary, pgx.items()):
        det = []
        for name, rs, func, n, v, cov in alle:
            s = "not assessable" if n is None else {0: "not carried", 1: "heterozygous", 2: "homozygous"}[n]
            det.append(f"{esc(name)} <span class='small'>{esc(rs)}, {esc(func)}</span>: {esc(s)} "
                       f"<span class='small'>({esc(cov)})</span>")
        kind = "t-grn" if re.match(r"^(Normal|normal|typical)", ph) else ("" if dip == "indeterminate" else "t-amb")
        rows.append([f"<b>{esc(gene)}</b>", "<br>".join(det), esc(dip), tag(ph, kind), esc(PGX_DRUGS.get(gene, ""))])
        R.rows.append(("pgx", None, (gene, dip, ph)))
    R.table(["Gene", "Alleles tested", "Genotype", "Predicted phenotype", "Affected drugs"], rows)
    R.h(f"<p class='small'><b>Not assessed by this exome pipeline:</b> {esc(NOT_ASSESSED_PGX)}.</p>")
    pharmcat = load_pharmcat(args.pharmcat)
    if args.pharmcat:
        R.h(f"<h3>PharmCAT {esc(pharmcat['version'])}</h3><p class='small'>Diplotypes from PharmCAT's named-allele "
            "matcher (all defining positions of each star allele, including multi-variant haplotypes), with its "
            "CPIC/DPWG phenotypes; recommendations are in the PharmCAT report. <b>Missing positions</b> were not in the "
            "VCF (often off-target in an exome): an allele defined there cannot be excluded. Where PharmCAT and the "
            "table above disagree, the PharmCAT call uses more positions but assumes nothing about missing ones.</p>")
        rows = []
        for r in pharmcat["calls"]:
            ph = r.get("Phenotype", "")
            kind = "t-grn" if re.match(r"^(Normal|normal)", ph) else ("" if not ph or "no call" in r.get("Source Diplotype", "") else "t-amb")
            miss = r.get("Missing positions", "")
            rows.append([f"<b>{esc(r.get('Gene', ''))}</b>", esc(r.get("Source Diplotype", "")),
                         tag(ph, kind) if ph else "", esc(r.get("Activity Score", "")),
                         f"<span class='small'>{esc(miss)}</span>" if miss else "none"])
            R.rows.append(("pgx_pharmcat", None, (r.get("Gene", ""), r.get("Source Diplotype", ""), ph)))
        R.table(["Gene", "Diplotype (PharmCAT)", "Phenotype", "Activity score", "Missing positions"], rows)
    if clinvar_drug:
        shown = sorted([v for v in clinvar_drug if (v.stars or 0) >= 1], key=lambda v: (-(v.stars or 0), v.gene))
        R.h(f"<h3>ClinVar drug-response variants</h3><p class='small'>{len(shown)} with review status &ge; 1 star shown "
            f"(most reliable first); {len(clinvar_drug) - len(shown)} without assertion criteria are listed in the TSV only.</p>")
        R.table(["Variant", "Genotype", "ClinVar"], [[var_cell(v), call_cell(v, sex), clinvar_cell(v)] for v in shown])
        for v in clinvar_drug:
            R.rows.append(("clinvar_drug_response", v, ""))

    # risk alleles
    R.h("<h2>8. Risk alleles and other associations</h2>")
    rows = []
    e4 = e2 = None
    for gene, name, rs, n, v, cov, text in risk:
        if gene == "APOE":
            if "e4" in name:
                e4 = n
            else:
                e2 = n
            continue
        s = "not assessable" if n is None else {0: "not carried", 1: "heterozygous", 2: "homozygous"}[n]
        rows.append([f"<b>{esc(gene)}</b>", f"{esc(name)} <span class='small'>{esc(rs)}</span>",
                     tag(s, "t-amb" if n else ("t-grn" if n == 0 else "")) + f"<br><span class='small'>{esc(cov)}</span>",
                     esc(text)])
        if n:
            R.rows.append(("risk_allele", v, f"{gene} {name}: {s}"))
    if e4 is None or e2 is None:
        rows.append(["<b>APOE</b>", "rs429358 / rs7412", tag("not assessable"), ""])
    else:
        apoe = "/".join((["e4"] * e4 + ["e2"] * e2 + ["e3"] * 2)[:2])
        txt = ("APOE genotype is associated with late-onset Alzheimer disease risk (e4) and lipid levels. It is not "
               "predictive, is not a standard return-of-results finding, and should only be disclosed after genetic "
               "counselling and with the person's explicit consent.")
        rows.append(["<b>APOE</b>", "rs429358 / rs7412 <span class='small'>(e4 / e2 tags)</span>",
                     tag(apoe, "t-amb" if e4 else "t-grn"), esc(txt)])
        R.rows.append(("apoe", None, ("APOE", apoe, "")))
    R.table(["Gene", "Allele", "Genotype", "Meaning"], rows)
    if other_assoc:
        shown = sorted([v for v in other_assoc if (v.stars or 0) >= 1], key=lambda v: (-(v.stars or 0), v.gene))
        R.h(f"<h3>Other ClinVar risk factor / association / protective variants</h3><p class='small'>{len(shown)} with "
            f"review status &ge; 1 star shown; {len(other_assoc) - len(shown)} without assertion criteria (mostly common "
            "polymorphisms from early association studies) are listed in the TSV only.</p>")
        R.table(["Variant", "Genotype", "ClinVar", "gnomAD frequency"],
                [[var_cell(v), call_cell(v, sex), clinvar_cell(v), freq_cell(v)] for v in shown[:100]])
        for v in other_assoc:
            R.rows.append(("clinvar_association", v, ""))

    # ROH
    R.h("<h2>9. Runs of homozygosity and uniparental disomy</h2>")
    if roh:
        big = sorted([s for s in roh if s[3] >= 1_000_000], key=lambda s: -s[3])
        total = sum(s[3] for s in big)
        per_chrom = collections.Counter()
        for s in big:
            per_chrom[s[0]] += s[3]
        R.h(f"<p>{len(big)} segments &ge; 1 Mb, total <b>{total / 1e6:.1f} Mb</b> "
            f"(about {100 * total / 2.875e9:.1f}% of the autosomal genome; estimated from exome variant sites with "
            "bcftools roh).</p><p class='small'>Interpretation guide: a total above ~70-100 Mb suggests parental "
            "relatedness (first cousins: ~6%, 180 Mb), increasing the prior for recessive disease; a single very long "
            "run confined to one chromosome can indicate uniparental isodisomy. Exome-based estimates are approximate.</p>")
        if big and big[0][3] > 20_000_000 and per_chrom.most_common(1)[0][1] > 0.8 * total and total > 25_000_000:
            R.h("<div class='box warn'>Homozygosity is concentrated on chromosome "
                f"{esc(per_chrom.most_common(1)[0][0])}: consider uniparental disomy if clinically relevant.</div>")
        R.table(["Chromosome", "Start", "End", "Length (Mb)", "Sites", "Quality"],
                [[esc(c), f"{s:,}", f"{e:,}", f"{ln / 1e6:.2f}", str(n), f"{q:.0f}" if q is not None else ""]
                 for c, s, e, ln, n, q in big[:40]])
        rhom = [v for v in variants if v.zyg == "hom" and in_roh(v) and (is_plp(v) or any(c[1] is v for c in cands))]
        if rhom:
            R.h("<p>Diagnostic or candidate homozygous variants inside a run of homozygosity: "
                + ", ".join(esc(f"{v.gene} {v.hgvs_p or v.hgvs_c}") for v in rhom) + "</p>")
    else:
        R.h("<p class='none'>Not computed (no bcftools roh output given).</p>")
    R.h("<h3>Uniparental disomy (trio)</h3>")
    if not both_parents:
        R.h("<p class='none'>Needs both parents' calls.</p>")
    else:
        R.h("<p class='small'>Common SNVs (gnomAD &ge; 1 %, QUAL &ge; 30, depth &ge; 20) where the proband lacks a "
            "parent's allele: homozygous for an allele that parent lacks, or (in targets the proband covers at "
            "&ge; 20x) without an allele the parent is homozygous for. Each chromosome is scanned for the run where "
            "such sites are most frequent, against the trio's genome-wide rate (mother "
            f"{upd_rates.get('mother', 0):.2%}, father {upd_rates.get('father', 0):.2%}); a run of &ge; 2 Mb with "
            f"&ge; {UPD_MIN_EXCL} such sites and log10 likelihood ratio &ge; {UPD_MIN_LOD:.0f} is reported. "
            "Isodisomy (two copies of one parental homologue) leaves the run homozygous; heterodisomy (both "
            "homologues of one parent) does not. In validation with simulated UPD in three GIAB trios, whole-chromosome "
            "iso- and heterodisomy and isodisomic segments of 8-30 Mb were all found; heterodisomic segments under ~30 Mb "
            "mostly were not. A deletion of one parent's copy looks like isodisomy: copy number tells them apart. "
            "Confirm by methylation or microsatellite analysis.</p>")
        if not upd_segs:
            R.h("<p class='none'>No uniparental disomy found.</p>")
        else:
            rows = []
            for d in upd_segs:
                notes = []
                if d["deletion"]:
                    notes.append("a high-confidence CNV deletion covers it: more likely a deletion of the "
                                 f"{d['missing']}'s copy than UPD")
                elif not args.cnv_calls:
                    notes.append("no copy-number data: a deletion of the " + d["missing"] + "'s copy would look "
                                 "the same")
                if d["roh_mb"]:
                    notes.append(f"{d['roh_mb']:.1f} Mb in runs of homozygosity")
                inside = [v for v in variants if v.zyg == "hom" and in_segments(v, [d]) and
                          (is_plp(v) or any(c[1] is v for c in cands))]
                if inside:
                    notes.append("homozygous diagnostic or candidate variants here (only the " + d["from"]
                                 + " carries them, low recurrence risk): " + ", ".join(
                                     f"{v.gene} {v.hgvs_p or v.hgvs_c}" for v in inside[:6]))
                rows.append([f"<b>{esc(d['chrom'])}</b>", f"{d['start']:,}-{d['end']:,} "
                             f"({(d['end'] - d['start']) / 1e6:.1f} Mb)", f"<b>{esc(d['from'])}</b>", esc(d["kind"]),
                             f"{d['excluded']} of {d['sites']} sites lack the {esc(d['missing'])}'s allele; "
                             f"log10 LR {d['lod']:.0f}"
                             + (f"; heterozygosity {d['het']:.0%} of the genome's" if d["het"] is not None else ""),
                             f"<b>{esc(d['disorder'])}</b>" if d["disorder"] else "no imprinting disorder known",
                             esc("; ".join(notes))])
                R.rows.append(("upd", None, ("",
                                             f"UPD({'mat' if d['from'] == 'mother' else 'pat'}) {d['kind']}",
                                             f"both copies from the {d['from']}, {d['excluded']}/{d['sites']} sites, "
                                             f"log10 LR {d['lod']:.0f}"
                                             + (f", {d['disorder']}" if d["disorder"] else "")
                                             + (", " + ", ".join(notes) if notes else ""),
                                             (d["chrom"], d["start"], d["end"]))))
            R.table(["Chromosome", "Region", "Both copies from", "Type", "Evidence", "Imprinting disorder", "Notes"],
                    rows)

    # coverage
    R.h("<h2>10. Coverage of clinically relevant genes</h2>")
    if gene_cov:
        rel = set(ACMG_SF) | set(pgx) | set(g for g in hpo_match) | {v.gene for v, _ in diag}
        rows = []
        for g in sorted(rel):
            if g in gene_cov and gene_cov[g][0]:
                tot, b20, b10 = gene_cov[g]
                if b20 / tot < 0.95 or g in hpo_match:
                    rows.append([esc(g), f"{tot:,}", f"{100 * b10 / tot:.1f}%", f"{100 * b20 / tot:.1f}%",
                                 esc(ACMG_SF_CATEGORY.get(g, "pharmacogene" if g in pgx else ("phenotype gene" if g in hpo_match else "")))])
        allg = [(g, x) for g, x in gene_cov.items() if x[0]]
        well = sum(1 for g, x in allg if x[1] / x[0] >= 0.95)
        R.h(f"<p>{well:,} of {len(allg):,} protein-coding genes have &ge; 95% of their coding bases at &ge; 20x. "
            "Clinically relevant genes (ACMG SF, pharmacogenes, phenotype genes, genes with findings) below 95%:</p>")
        R.table(["Gene", "Coding bases", ">= 10x", ">= 20x", "Group"], rows)
    else:
        R.h("<p class='none'>Per-gene coverage not given (VCF-only mode, or no coding-region BED). A negative result "
            "cannot be distinguished from missing coverage.</p>")

    # methods & limits
    R.h("<h2>11. Methods, data sources and limitations</h2><ul>")
    R.h("<li>Variant annotation: SnpEff (consequence, HGVS), ClinVar (clinical significance and review status), "
        "gnomAD (population allele frequency) as provided to the workflow.</li>")
    R.h("<li>Frequency rules use the higher of the overall gnomAD allele frequency and the highest filtering allele "
        "frequency (FAF95) of a genetic ancestry group, where the frequency file provides it (gnomAD v4.1 joint "
        "exome + genome extract): BA1 above 5%, BS1 above the rare-variant cut-off, PM2_supporting below 0.01% "
        "(dominant) or 0.1% (recessive). Frequencies from gnomAD sites that failed gnomAD's own filters or whose exome "
        "and genome frequencies disagree are flagged. A variant absent from gnomAD is not proof of rarity where gnomAD "
        "coverage is low.</li>")
    R.h(f"<li>Gene-disease validity and inheritance: GenCC ({'loaded, ' + str(len(gencc)) + ' genes' if gencc else 'not provided'}"
        f"{'; disputed, refuted and no-known-relationship assertions are shown next to the gene' if K.disputed else ''}); "
        f"gene constraint: gnomAD LOEUF ({'loaded, ' + str(len(constraint)) + ' genes' if constraint else 'not provided'}); "
        f"dosage sensitivity: ClinGen ({'loaded, ' + str(len(dosage)) + ' genes' if dosage else 'not provided'}); "
        f"phenotype: HPO phenotype_to_genes ({'loaded' if args.hpo_genes else 'not provided'}) and disease annotations "
        f"({'loaded, ' + str(len(hpoa)) + ' diseases' if hpoa else 'not provided'}).</li>")
    R.h("<li>Loss of function counts as a disease mechanism (PVS1, and a LoF candidate in a dominant gene) for recessive "
        "or X-linked genes, ClinGen haploinsufficiency score 3 or 2, or LOEUF &lt; 0.6 unless ClinGen scores "
        "haploinsufficiency as unlikely (Abou Tayoun et al. 2018; Riggs et al. 2020).</li>")
    if cvp:
        R.h(f"<li>Protein-level ClinVar evidence ({len(cvp):,} residues with pathogenic or benign missense variants, "
            "&ge; 1 star): PS1 when the same amino-acid change by a different nucleotide change is pathogenic and "
            "SpliceAI is below 0.2; PM5 when another missense change at the residue is pathogenic. Residues are compared "
            "only where the reference amino acid agrees between SnpEff's and ClinVar's transcripts. The density of "
            "pathogenic and benign missense variants within 5 residues is shown for orientation, not scored "
            "(a PM1 hotspot needs a gene-specific definition).</li>")
    if args.bam:
        R.h("<li>Two heterozygous variants in a recessive gene are phased from the reads where a read or read pair "
            "covers both (within 1 kb; MAPQ &ge; 20; &ge; 4 informative reads, &ge; 90% agreeing). Variants in cis "
            "are one allele: the patient is then a carrier, not affected. Otherwise the phase stays unknown and parental "
            "testing decides.</li>")
    if hpo_terms and pheno.disease_terms:
        R.h("<li>Disease-level phenotype match: for each gene, the linked HPO-annotated disease that explains the largest "
            "share of the patient's terms, each weighted by its information content (-ln of the fraction of diseases "
            "annotated with the term or a descendant). A match of &ge; 50% adds one rank point to a candidate.</li>")
    if args.artefacts:
        R.h(f"<li>Recurrent-artefact list: {len(artefacts)} alleles called in many of the laboratory's own samples but "
            f"rare or absent in gnomAD (systematic errors of the assay); {n_artefact} calls of this sample match it and "
            "are marked low confidence, not removed.</li>")
    else:
        R.h("<li>No recurrent-artefact list was provided: recurrent false calls of the assay (e.g. in GC-rich targets) "
            "are not recognised.</li>")
    R.h("<li>Rules: ACMG SF v3.3 gene list and reporting rules (Genet Med 2025;27:101454); automated ACMG/AMP 2015 "
        "classification (Richards et al. 2015) with the ClinGen SVI modifications, combined by points (Tavtigian et "
        "al. 2018, 2020): PVS1 by the Abou Tayoun et al. 2018 decision tree as far as SnpEff supports it, on the principal "
        "transcript (longest protein among those hit; NMD when the stop is before the penultimate exon or SnpEff "
        "predicts it): nonsense/frameshift with NMD, canonical splice sites outside UTR introns and exon loss Very "
        "strong; NMD-escaping truncation Strong if it removes &gt; 10 % of the protein, else Moderate; start loss "
        "Moderate); PM2_Supporting (gnomAD &lt; 0.01 % for "
        "dominant-only genes, &lt; 0.1 % otherwise); BA1 &gt; 5 % (ClinGen exception list excepted), BS1 above the "
        "rare-variant cut-off (0.1 % for dominant-only genes); BS2_Supporting for &ge; 2 gnomAD homozygotes in a "
        "recessive-only autosomal gene; PS1, PM5 and PM1_Supporting (&ge; 3 pathogenic and no benign ClinVar missense "
        "variants at other residues within 5 residues; with PP3 at most 4 points) from the ClinVar residue table; PM4 "
        "for in-frame changes and stop loss outside repeats, BP3 inside repeats; PP2 (missense Z &gt; 3.09); PP3/BP4 "
        "and BP7 (synonymous, SpliceAI &le; 0.1) from the calibrated predictors; PM3 (ClinGen SVI points; in trans, "
        "phase unknown or homozygous with a ClinVar pathogenic variant), PS2/PM6 and PP4 only when phenotype terms "
        "match the gene; CPIC allele-function and "
        "phenotype assignments for single-nucleotide alleles.</li>")
    n_am = sum(1 for v in variants if v.am is not None)
    n_sai = sum(1 for v in variants if v.sai is not None)
    R.h("<li>In-silico predictors (one calibrated tool per mechanism, ClinGen SVI): missense PP3/BP4 from AlphaMissense "
        "(Cheng et al. 2023) at the calibrated strengths of Bergquist et al. 2025 (PP3 Supporting &ge; 0.792, Moderate "
        "&ge; 0.906, Strong &ge; 0.990; BP4 Supporting &le; 0.169, Moderate &le; 0.099); splicing PP3 (max delta &ge; 0.2) "
        "and BP4 (&le; 0.1) at Supporting strength from SpliceAI outside the canonical &plusmn;1,2 positions (Walker et al. "
        "2023). The two PP3s are not summed, PP3 is not combined with PVS1, and missense BP4 needs SpliceAI &le; 0.1. "
        f"Scores present for this sample: AlphaMissense {n_am}, SpliceAI {n_sai} variants"
        + ("" if n_am and n_sai else " (a missing predictor gives no PP3/BP4)") + ".</li>")
    if trio:
        R.h("<li>Trio: parental genotypes from the parents' calls (multi-allelic and complex records split and trimmed "
            "to match the proband's normalised alleles)" + (" and, where a parent has no call, from the parent's reads "
                                                             "(MAPQ &ge; 20, duplicates excluded)"
                                                             if trio.bam["mother"] or trio.bam["father"] else "")
            + ". Compound heterozygotes are phased by parental origin first, then from the reads. A dominant finding "
            "inherited from a parent reported unaffected is kept but marked (reduced penetrance or not causal); such "
            "a candidate does not count as fitting the inheritance. De novo evidence follows the ClinGen SVI "
            "recommendation (confirmed de novo, phenotype consistent with the gene: PS2_Moderate; assumed: "
            "PM6_Supporting) and is withheld when the parentage check fails."
            + (" Uniparental disomy (section 9): for each chromosome and parent, the maximum-scoring run (log10 "
               "likelihood ratio, exclusion rate 30 % for isodisomy or 10 % for heterodisomy against the trio's "
               "genome-wide Mendelian-error rate) of common SNVs where the proband lacks that parent's allele; errors "
               "within 100 kb of each other count once. Parent of origin of a chromosome gain: the proband's read "
               "fraction of the mother's allele where one parent is homozygous ALT and the other lacks it (0.5 two "
               "copies, ~0.67 extra maternal, ~0.33 extra paternal copy); of the X chromosomes: whether the proband "
               "carries alleles only the mother or only the father has." if both_parents else "") + "</li>")
    if args.cnv_calls and not cnv_status:
        R.h("<li>Copy-number variants (section 4): wes_cnv read-depth calls (GC- and length-corrected target depth "
            "against a same-kit panel, one principal component of the panel removed, Student-t hidden Markov model). "
            f"ClinGen dosage regions: {'loaded, ' + str(len(cnv_regions)) + ' regions' if cnv_regions else 'not provided'}."
            " Small (1-2 exon) and heterozygous single-exon CNVs, balanced rearrangements, CNVs in segmental "
            "duplications and the breakpoints of any CNV are beyond exome read depth.</li>")
    R.h("<li><b>Not detected by this pipeline:</b> "
        + ("single- and two-exon copy-number changes (CNVs of several exons: section 4), "
           if args.cnv_calls and not cnv_status else "copy-number variants and exon deletions/duplications, ")
        + "structural variants, repeat expansions (e.g. FMR1, HTT, C9orf72), mitochondrial DNA variants, deep intronic and most "
        "regulatory variants, low-level mosaicism, methylation/imprinting defects, variants in segmental duplications "
        "and pseudogene regions (e.g. SMN1/SMN2, CYP21A2, PMS2 exons 11-15, GBA1), and genes or exons with low coverage. "
        "A negative report does not exclude a genetic diagnosis.</li>")
    R.h("<li>The automated ACMG/AMP class uses only the evidence a VCF and the knowledge files support. Functional "
        "studies (PS3/BS3), case counts (PS4), segregation (PP1/BS4), an alternative cause (BP5), gene-specific "
        "ClinGen VCEP rules and the literature are not used, so a variant can move by several points on manual "
        "review; a likely pathogenic or pathogenic class needs that review before it is reported.</li>")
    R.h("<li>ClinVar classifications change: findings should be re-evaluated against the current ClinVar release "
        "before reporting, and the data periodically reanalysed.</li></ul>")
    R.h("</body></html>")
    return R, sample


def write_tsv(R, path, sample):
    with open(path, "w") as f:
        f.write("sample\tcategory\tgene\tchrom\tpos\tref\talt\tgenotype\tdepth\talt_fraction\thgvs_c\thgvs_p\teffect\t"
                "gnomad_af\tgnomad_grpmax_af\tgnomad_faf95\tgnomad_nhomalt\talphamissense\tspliceai\tclinvar\t"
                "clinvar_review\tclinvar_condition\tqc_flags\tconfidence\tnote\tacmg_class\tacmg_points\tacmg_codes\n")
        acmg = getattr(R, "acmg", {})
        for cat, v, note in R.rows:
            if v is None and len(note) == 4:  # CNV: (genes, type and copy number, note, (chrom, start, end))
                genes, genotype, txt, (chrom, start, end) = note
                f.write("\t".join([sample, cat, genes, chrom, str(start), "", f"<{genotype.split()[0]}>", genotype]
                                  + [""] * 15 + [txt.split("; ")[2] if txt.count("; ") >= 2 else "",
                                                 txt.replace("\t", " ") + f"; end {end}", "", "", ""]) + "\n")
                continue
            if v is None:  # gene-level result: (gene, genotype, phenotype)
                gene, genotype, pheno = note
                f.write("\t".join([sample, cat, gene, "", "", "", "", genotype] + [""] * 16 + [pheno, "", "", ""]) + "\n")
                continue
            lc = low_confidence(v)
            f.write("\t".join(str(x) for x in (
                sample, cat, v.gene, v.chrom, v.pos, v.ref, v.alt, v.gt, v.dp if v.dp is not None else "",
                f"{v.ab:.3f}" if v.ab is not None else "", v.hgvs_c, v.hgvs_p, v.effect, fmt_af(v.af),
                fmt_af(v.grpmax_af) if v.grpmax_af is not None else "", fmt_af(v.faf) if v.faf is not None else "",
                int(v.nhomalt) if v.nhomalt is not None else "",
                f"{v.am:.3f}" if v.am is not None else "", f"{v.sai:.2f}" if v.sai is not None else "", v.clnsig,
                v.revstat, v.clndn, ";".join(v.flags), "low: " + ", ".join(lc) if lc else "high",
                note.replace("\t", " "))) + "\t"
                + ("\t".join([acmg[id(v)][0], str(acmg[id(v)][1]), ",".join(acmg[id(v)][2])]) if id(v) in acmg else "\t\t")
                + "\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--vcf", required=True, help="annotated VCF (SnpEff ANN, ClinVar CLN*, INFO/gnomAD_AF and the other gnomAD_* fields)")
    ap.add_argument("--sample-name")
    ap.add_argument("--gencc", help="GenCC submissions TSV")
    ap.add_argument("--constraint", help="gnomAD gene constraint TSV")
    ap.add_argument("--hpo-genes", help="HPO phenotype_to_genes.txt")
    ap.add_argument("--hpo-terms", help="patient HPO terms, comma/space separated (HP:0001250 ...)")
    ap.add_argument("--hpo-annotations", help="HPO disease annotations (phenotype.hpoa): disease names, onset, course")
    ap.add_argument("--clingen-dosage", help="ClinGen gene curation list (ClinGen_gene_curation_list_GRCh38.tsv)")
    ap.add_argument("--clinvar-protein", help="protein-level ClinVar table (data/build_clinvar_protein.py): PS1/PM5")
    ap.add_argument("--mother-vcf", help="the mother's variant calls (any single-sample VCF of the same build)")
    ap.add_argument("--father-vcf", help="the father's variant calls")
    ap.add_argument("--mother-bam", help="the mother's indexed BAM: reads at variants absent from her calls")
    ap.add_argument("--father-bam", help="the father's indexed BAM")
    ap.add_argument("--affected-parents", choices=("unknown", "none", "mother", "father", "both"), default="unknown",
                    help="which parents have the proband's condition")
    ap.add_argument("--cnv-calls", help="copy-number calls of the sample (wes_cnv calls TSV)")
    ap.add_argument("--cnv-summary", help="wes_cnv summary JSON: panel, noise, chromosome dosage, sex chromosomes")
    ap.add_argument("--clingen-regions", help="ClinGen region curation list (ClinGen_region_curation_list_GRCh38.tsv)")
    ap.add_argument("--mother-cnv", help="the mother's wes_cnv per-target ratios (--targets output): her read depth "
                                         "at the proband's CNVs, inherited or de novo")
    ap.add_argument("--father-cnv", help="the father's wes_cnv per-target ratios")
    ap.add_argument("--roh", help="bcftools roh output (RG lines)")
    ap.add_argument("--artefacts", help="recurrent-artefact VCF of the laboratory (data/build_artefact_list.sh)")
    ap.add_argument("--pharmcat", help="PharmCAT calls-only TSV (pharmcat_pipeline -reporterCallsOnlyTsv)")
    ap.add_argument("--target-coverage", help="mosdepth thresholds BED over the exome targets")
    ap.add_argument("--gene-coverage", help="mosdepth thresholds BED over coding exons named by gene")
    ap.add_argument("--coverage-summary", help="mosdepth summary")
    ap.add_argument("--bam", help="indexed BAM: read depth at pharmacogenomic/risk positions (needs samtools)")
    ap.add_argument("--reference", help="indexed reference FASTA: repeat context of reported variants (needs samtools)")
    ap.add_argument("--min-ref-depth", type=int, default=15,
                    help="reads needed to read an uncalled pharmacogenomic/risk position as reference")
    ap.add_argument("--max-af", type=float, default=0.01, help="rare-variant cut-off, applied to the higher of gnomAD AF and grpmax FAF95")
    ap.add_argument("--max-candidates", type=int, default=60)
    ap.add_argument("--html", required=True)
    ap.add_argument("--tsv", required=True)
    args = ap.parse_args(argv)
    R, sample = build(args)
    with open(args.html, "w") as f:
        f.write("\n".join(R.parts))
    write_tsv(R, args.tsv, sample)


if __name__ == "__main__":
    sys.exit(main())
