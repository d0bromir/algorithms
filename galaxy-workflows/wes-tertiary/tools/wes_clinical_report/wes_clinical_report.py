#!/usr/bin/env python3
"""Human-readable clinical interpretation report for one germline exome.

Input is the annotated VCF of the "WES tertiary analysis from VCF" workflow (SnpEff ANN/LOF/NMD,
ClinVar CLNSIG/CLNREVSTAT/CLNDN/..., gnomAD v4.1 frequencies INFO/gnomAD_AF, _FAF95, _grpmax_AF, _nhomalt, ...),
plus optional knowledge
bases and QC outputs. Every statement in the report is derived by explicit, versioned rules from those
inputs; nothing is inferred by a statistical or language model. The report is decision support for a
qualified reviewer (clinical geneticist / laboratory scientist), not a diagnosis.

Sections: sample QC and sex check; diagnostic findings (ClinVar P/LP with inheritance-consistent
genotype); candidate variants not classified in ClinVar (rare predicted loss-of-function, and rare
missense in phenotype-matched genes) with automated ACMG/AMP evidence codes; ACMG SF v3.3 secondary
findings; carrier status; pharmacogenomics (CPIC genes, single-nucleotide star alleles); risk alleles
and other ClinVar associations; runs of homozygosity; coverage of clinically relevant genes; limits.

Standard library only (Python >= 3.8).
"""
import argparse
import collections
import datetime
import gzip
import html
import math
import re
import subprocess
import sys

VERSION = "0.2.0"

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
                 "an", "nhomalt", "grpmax", "grpmax_af", "faf", "faf_grp", "gfilter", "gdisc")

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
                v.clnsig = str(info.get("CLNSIG", "")) if info.get("CLNSIG") not in (None, True, ".") else ""
                v.revstat = str(info.get("CLNREVSTAT", "")) if info.get("CLNREVSTAT") not in (None, True) else ""
                v.clndn = str(info.get("CLNDN", "")).replace("_", " ") if info.get("CLNDN") not in (None, True) else ""
                v.stars = REVIEW_STARS.get(v.revstat.split("|")[0], 0) if v.revstat else None
                rs = info.get("RS")
                v.rs = f"rs{rs}" if rs not in (None, True, ".") and not str(rs).startswith("rs") else rs
                v.flags = []
                variants.append(v)
    return sample, contigs, variants, prefix or ""


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


def load_hpo(path, terms):
    """For the patient's HPO terms: gene -> set of matched (term id, term name). Uses the propagated
    phenotype_to_genes.txt, so a general term matches genes annotated to any of its descendants."""
    match = collections.defaultdict(set)
    names = {}
    if not path or not terms:
        return match, names
    with opener(path) as f:
        f.readline()
        for line in f:
            t = line.rstrip("\n").split("\t")
            if len(t) < 4:
                continue
            names.setdefault(t[0], t[1])
            if t[0] in terms:
                match[t[3]].add(t[0])
    return match, names


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


def acmg_evidence(v, gene, moi, constraint, hpo_hit, max_af):
    """Automated subset of ACMG/AMP 2015 evidence codes (provisional; for triage, not classification)."""
    codes = []
    loeuf = constraint.get(gene, (None, None, None))[0]
    if is_lof(v, gene) and ("AR" in moi or "XL" in moi or (loeuf is not None and loeuf < 0.6)):
        ann = ann_for_gene(v, gene)
        nmd = "NMD" in v.info
        codes.append("PVS1" if (nmd or "splice" in ann[1] or "start_lost" in ann[1] or "exon_loss" in ann[1]
                                or "stop_gained" in ann[1] or "frameshift" in ann[1]) else "PVS1_moderate")
    p = freq(v)
    if p is None or p < (1e-4 if ("AD" in moi and "AR" not in moi) else 1e-3):
        codes.append("PM2_supporting")
    if p is not None and p > 0.05:
        codes.append("BA1")
    elif p is not None and p > max_af:
        codes.append("BS1")
    if is_missense(v, gene):
        z = constraint.get(gene, (None, None, None))[2]
        if z is not None and z > 3.09:
            codes.append("PP2")
    if hpo_hit:
        codes.append("PP4")
    return codes


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
    tandem repeat of >= 8, or a paralogous/hypervariable gene."""
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


def disease_cell(gene, gencc, hpo_match, hpo_names):
    entries = gencc.get(gene, [])
    s = ""
    if entries:
        ds, n = gene_diseases(entries)
        s = "<br>".join(f"{esc(d)} <span class='small'>({esc(m)}; {esc(c)})</span>" for d, m, c in ds)
        if n > len(ds):
            s += f"<br><span class='small'>+{n - len(ds)} more</span>"
    else:
        s = '<span class="small">no established gene-disease association in GenCC</span>' if gencc else ""
    if hpo_match.get(gene):
        s += "<br>" + "".join(tag(f"phenotype match: {hpo_names.get(t, t)}", "t-blu") for t in sorted(hpo_match[gene]))
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


def build(args):
    sample, contigs, variants, prefix = parse_vcf(args.vcf)
    sample = args.sample_name or sample or "sample"
    grch38 = contigs.get("1") == 248956422 if contigs else None
    gencc = load_gencc(args.gencc)
    constraint = load_constraint(args.constraint)
    hpo_terms = {t.strip().upper() for t in re.split(r"[,;\s]+", args.hpo_terms or "") if t.strip()}
    hpo_terms = {t if t.startswith("HP:") else "HP:" + t for t in hpo_terms if re.match(r"(HP:)?\d{7}$", t)}
    hpo_match, hpo_names = load_hpo(args.hpo_genes, hpo_terms)
    roh = load_roh(args.roh)
    target_cov = load_thresholds(args.target_coverage)
    gene_cov_rows = load_thresholds(args.gene_coverage)
    msum = load_mosdepth_summary(args.coverage_summary)
    qc = sample_qc(variants, contigs, msum, target_cov)
    sex = qc["sex"]
    for v in variants:
        v.flags = variant_flags(v)
    by_pos = {(v.chrom, v.pos, v.alt): v for v in variants}

    v_chrom_of_gene = {v.gene: v.chrom for v in variants}

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
    by_gene_plp = collections.defaultdict(list)
    for v in plp:
        by_gene_plp[v.gene].append(v)
    for v in plp:
        g = v.gene
        moi = moi_of(g)
        nonpar_x = v.chrom == "X" and not any(a <= v.pos <= b for a, b in PAR_GRCH38)
        biallelic = v.zyg == "hom" or len(by_gene_plp[g]) >= 2
        two = ("homozygous" if v.zyg == "hom" else
               "two P/LP alleles in this gene (possible compound heterozygote, phase unknown)")
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
                if is_lof(v, g) or (is_missense(v, g) and hpo_match.get(g)):
                    rare_by_gene[g].append(v)
                break
    for g, vs in rare_by_gene.items():
        moi = moi_of(g)
        loeuf = constraint.get(g, (None, None, None))[0]
        for v in vs:
            lof = is_lof(v, g)
            reasons, consistent = [], False
            if "AR" in moi:  # recessive model: two alleles needed
                ok = v.zyg == "hom" or len(vs) + len(by_gene_plp.get(g, [])) >= 2
                hom_pop = v.zyg == "hom" and (v.nhomalt or 0) >= 2
                consistent |= ok and not hom_pop
                reasons.append("recessive: " + (f"homozygous, but {int(v.nhomalt)} homozygotes in gnomAD" if hom_pop else
                                                "homozygous" if v.zyg == "hom" else
                                                "second rare allele in the gene (phase unknown)" if ok else
                                                "single heterozygous allele (carrier at most)"))
            if not moi or {"AD", "XL", "MT"} & moi:  # dominant / X-linked model: one rare allele
                strict = freq(v) is None or freq(v) < 1e-4
                tolerant = lof and loeuf is not None and loeuf >= 0.6 and "XL" not in moi
                consistent |= strict and not tolerant
                r = ("dominant/X-linked" if moi else "inheritance unknown")
                if not strict:
                    r += ": too frequent for a dominant disorder"
                elif tolerant:
                    r += f": LoF-tolerant gene (LOEUF {loeuf:.2f})"
                reasons.append(r)
            why = "; ".join(reasons)
            ev = acmg_evidence(v, g, moi, constraint, bool(hpo_match.get(g)), args.max_af)
            score = (2 if hpo_match.get(g) else 0) + (2 if consistent else 0) + (1 if lof else 0) + \
                (1 if "PVS1" in ev else 0) - (1 if v.flags and any("depth" in f or "balance" in f for f in v.flags) else 0)
            cands.append((score, v, g, consistent, why, ev))
    cands.sort(key=lambda x: (-x[0], x[2], x[1].pos))

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
                 + [c[1] for c in cands] + [x[0] for x in sf]}
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

    def short(v, gene=True):
        txt = f"{v.gene} {v.hgvs_p or v.hgvs_c}" if gene else v.gene
        return esc(txt) + (" <i>(low confidence)</i>" if low_confidence(v) else "")

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
    R.h(f"<li><b>Candidate variants needing review</b> (not classified in ClinVar): <b>{sum(1 for c in cands if c[3])}</b> "
        f"inheritance-consistent of {len(cands)}</li>")
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
    R.h(f"<li><b>Sex (genetic)</b>: {esc(sex)}</li>")
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
    warn = []
    if qc.get("tstv") is not None and not (2.3 <= qc["tstv"] <= 3.5):
        warn.append("Ts/Tv outside the usual exome range: check for false positives or off-target calls.")
    if qc.get("het_ab_low") is not None and qc["het_ab_low"] > 0.15:
        warn.append("Many heterozygous calls with low alt fraction: possible sample contamination or mosaicism.")
    if qc.get("pct_20X") is not None and qc["pct_20X"] < 90:
        warn.append(f"Only {qc['pct_20X']:.0f}% of target bases reach 20x (clinical exomes usually require >= 90-95%): "
                    "absence of findings is not reliable in poorly covered genes (see section 9).")
    if warn:
        R.h("<div class='box warn'>" + "<br>".join(esc(w) for w in warn) + "</div>")

    # diagnostic
    R.h("<h2>2. Diagnostic findings</h2><p class='small'>ClinVar pathogenic or likely pathogenic variants (not "
        "conflicting) whose genotype fits the gene's mode of inheritance (GenCC): one allele in dominant or X-linked "
        "genes, two alleles in recessive genes. Whether a finding explains the patient's condition depends on the "
        "clinical picture.</p>")
    R.table(["Variant", "Genotype", "ClinVar", "gnomAD frequency", "Gene-disease (GenCC)", "Interpretation"],
            [[var_cell(v), call_cell(v, sex), clinvar_cell(v), freq_cell(v),
              disease_cell(v.gene, gencc, hpo_match, hpo_names), esc(why)] for v, why in diag])
    for v, why in diag:
        R.rows.append(("diagnostic", v, why))
    R.h("<h3>Pathogenic variants whose relevance depends on the condition (review)</h3><p class='small'>Genes with "
        "both dominant and recessive disorders, X-linked genes in heterozygous females, and genes without an "
        "established disease association in GenCC.</p>")
    R.table(["Variant", "Genotype", "ClinVar", "gnomAD frequency", "Gene-disease (GenCC)", "Why review"],
            [[var_cell(v), call_cell(v, sex), clinvar_cell(v), freq_cell(v),
              disease_cell(v.gene, gencc, hpo_match, hpo_names), esc(why)] for v, why in review])
    for v, why in review:
        R.rows.append(("pathogenic_review", v, why))
    if lowpen:
        R.h("<h3>Low-penetrance pathogenic alleles</h3>")
        R.table(["Variant", "Genotype", "ClinVar", "gnomAD frequency"],
                [[var_cell(v), call_cell(v, sex), clinvar_cell(v), freq_cell(v)] for v in lowpen])
        for v in lowpen:
            R.rows.append(("low_penetrance", v, ""))

    # candidates
    R.h("<h2>3. Candidate variants not classified in ClinVar</h2><p class='small'>Rare (gnomAD frequency &lt; "
        f"{args.max_af:g}) predicted loss-of-function variants in genes with an established disease association "
        "(GenCC definitive/strong/moderate/supportive)" + (", and rare missense/in-frame variants in genes matching "
                                                          "the patient's phenotype" if hpo_terms else "")
        + ". Evidence codes are an automated subset of the ACMG/AMP 2015 criteria (PVS1, PM2, BA1, BS1, PP2, PP4) "
        "for triage; they are not a classification. Rows are ranked by phenotype match, inheritance fit and "
        "variant consequence.</p>")
    rows = []
    for score, v, g, consistent, why, ev in cands[:args.max_candidates]:
        a = ann_for_gene(v, g) or v.anns[0]
        loeuf = constraint.get(g, (None, None, None))[0]
        rows.append([var_cell(v).replace(f"<b>{esc(v.gene)}</b>", f"<b>{esc(g)}</b>"), call_cell(v, sex),
                     esc(a[1].replace("&", " & ")) + (f"<br><span class='small'>LOEUF {loeuf:.2f}</span>" if loeuf is not None else ""),
                     freq_cell(v), disease_cell(g, gencc, hpo_match, hpo_names),
                     ("" if consistent else tag("inheritance not met", "t-grn") + "<br>") + esc(why),
                     " ".join(tag(c, "t-red" if c.startswith("PVS") else ("t-grn" if c.startswith("B") else "")) for c in ev)])
        R.rows.append(("candidate" if consistent else "candidate_inheritance_not_met", v, why + "; " + ",".join(ev)))
    R.table(["Variant", "Genotype", "Consequence", "gnomAD frequency", "Gene-disease (GenCC)", "Inheritance fit", "Evidence"], rows)
    if len(cands) > args.max_candidates:
        R.h(f"<p class='small'>{len(cands) - args.max_candidates} more candidates in the TSV output.</p>")
        for score, v, g, consistent, why, ev in cands[args.max_candidates:]:
            R.rows.append(("candidate" if consistent else "candidate_inheritance_not_met", v, why + "; " + ",".join(ev)))

    # ACMG SF
    R.h(f"<h2>4. Secondary findings: {ACMG_SF_VERSION} ({len(ACMG_SF)} genes)</h2><p class='small'>Medically "
        "actionable genes that the ACMG recommends reporting regardless of the indication for testing, if the patient "
        "has opted in. Gene-specific rules applied: biallelic only for MUTYH, ATP7B, BTD, GAA, RPE65, CYP27A1; ABCD1 "
        "hemizygous, homozygous or two heterozygous variants (not a single heterozygous variant in a female); HFE "
        "p.Cys282Tyr homozygotes only; TTN truncating variants only.</p>")
    if sf:
        R.table(["Variant", "Genotype", "ClinVar", "gnomAD frequency", "Category", "Status"],
                [[var_cell(v), call_cell(v, sex), clinvar_cell(v), freq_cell(v), esc(ACMG_SF_CATEGORY.get(g, "")),
                  tag(st, "t-red" if st == "reportable" else "t-amb") + "<br>" + esc(note)] for v, g, st, note in sf])
        for v, g, st, note in sf:
            R.rows.append(("acmg_sf_" + st.split()[0], v, note))
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
    R.h("<h2>5. Carrier status for recessive conditions</h2><p class='small'>One pathogenic/likely pathogenic "
        "allele in a gene associated with autosomal recessive disease. Carriers are usually unaffected; the result "
        "matters for reproductive risk if a partner carries a pathogenic variant in the same gene.</p>")
    R.table(["Variant", "Genotype", "ClinVar", "gnomAD frequency", "Condition (GenCC)"],
            [[var_cell(v), call_cell(v, sex), clinvar_cell(v), freq_cell(v), disease_cell(v.gene, gencc, {}, {})]
             for v in carrier])
    for v in carrier:
        R.rows.append(("carrier", v, ""))

    # PGx
    R.h("<h2>6. Pharmacogenomics</h2><p class='small'>Single-nucleotide star alleles of CPIC genes, read at fixed "
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
    if clinvar_drug:
        shown = sorted([v for v in clinvar_drug if (v.stars or 0) >= 1], key=lambda v: (-(v.stars or 0), v.gene))
        R.h(f"<h3>ClinVar drug-response variants</h3><p class='small'>{len(shown)} with review status &ge; 1 star shown "
            f"(most reliable first); {len(clinvar_drug) - len(shown)} without assertion criteria are listed in the TSV only.</p>")
        R.table(["Variant", "Genotype", "ClinVar"], [[var_cell(v), call_cell(v, sex), clinvar_cell(v)] for v in shown])
        for v in clinvar_drug:
            R.rows.append(("clinvar_drug_response", v, ""))

    # risk alleles
    R.h("<h2>7. Risk alleles and other associations</h2>")
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
    R.h("<h2>8. Runs of homozygosity</h2>")
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

    # coverage
    R.h("<h2>9. Coverage of clinically relevant genes</h2>")
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
    R.h("<h2>10. Methods, data sources and limitations</h2><ul>")
    R.h("<li>Variant annotation: SnpEff (consequence, HGVS), ClinVar (clinical significance and review status), "
        "gnomAD (population allele frequency) as provided to the workflow.</li>")
    R.h("<li>Frequency rules use the higher of the overall gnomAD allele frequency and the highest filtering allele "
        "frequency (FAF95) of a genetic ancestry group, where the frequency file provides it (gnomAD v4.1 joint "
        "exome + genome extract): BA1 above 5%, BS1 above the rare-variant cut-off, PM2_supporting below 0.01% "
        "(dominant) or 0.1% (recessive). Frequencies from gnomAD sites that failed gnomAD's own filters or whose exome "
        "and genome frequencies disagree are flagged. A variant absent from gnomAD is not proof of rarity where gnomAD "
        "coverage is low.</li>")
    R.h(f"<li>Gene-disease validity and inheritance: GenCC ({'loaded, ' + str(len(gencc)) + ' genes' if gencc else 'not provided'}); "
        f"gene constraint: gnomAD LOEUF ({'loaded, ' + str(len(constraint)) + ' genes' if constraint else 'not provided'}); "
        f"phenotype: HPO phenotype_to_genes ({'loaded' if args.hpo_genes else 'not provided'}).</li>")
    R.h("<li>Rules: ACMG SF v3.3 gene list and reporting rules (Genet Med 2025;27:101454); automated ACMG/AMP 2015 evidence "
        "subset (Richards et al. 2015; PVS1 after Abou Tayoun et al. 2018, simplified); CPIC allele-function and "
        "phenotype assignments for single-nucleotide alleles.</li>")
    R.h("<li><b>Not detected by this pipeline:</b> copy-number variants and exon deletions/duplications, structural "
        "variants, repeat expansions (e.g. FMR1, HTT, C9orf72), mitochondrial DNA variants, deep intronic and most "
        "regulatory variants, low-level mosaicism, methylation/imprinting defects, variants in segmental duplications "
        "and pseudogene regions (e.g. SMN1/SMN2, CYP21A2, PMS2 exons 11-15, GBA1), and genes or exons with low coverage. "
        "A negative report does not exclude a genetic diagnosis.</li>")
    R.h("<li>Candidate variants without a ClinVar classification need manual ACMG/AMP classification with all evidence "
        "types (segregation, de novo status, functional data, computational predictions, literature).</li>")
    R.h("<li>ClinVar classifications change: findings should be re-evaluated against the current ClinVar release "
        "before reporting, and the data periodically reanalysed.</li></ul>")
    R.h("</body></html>")
    return R, sample


def write_tsv(R, path, sample):
    with open(path, "w") as f:
        f.write("sample\tcategory\tgene\tchrom\tpos\tref\talt\tgenotype\tdepth\talt_fraction\thgvs_c\thgvs_p\teffect\t"
                "gnomad_af\tgnomad_grpmax_af\tgnomad_faf95\tgnomad_nhomalt\tclinvar\tclinvar_review\tclinvar_condition\t"
                "qc_flags\tconfidence\tnote\n")
        for cat, v, note in R.rows:
            if v is None:  # gene-level result: (gene, genotype, phenotype)
                gene, genotype, pheno = note
                f.write("\t".join([sample, cat, gene, "", "", "", "", genotype] + [""] * 14 + [pheno]) + "\n")
                continue
            lc = low_confidence(v)
            f.write("\t".join(str(x) for x in (
                sample, cat, v.gene, v.chrom, v.pos, v.ref, v.alt, v.gt, v.dp if v.dp is not None else "",
                f"{v.ab:.3f}" if v.ab is not None else "", v.hgvs_c, v.hgvs_p, v.effect, fmt_af(v.af),
                fmt_af(v.grpmax_af) if v.grpmax_af is not None else "", fmt_af(v.faf) if v.faf is not None else "",
                int(v.nhomalt) if v.nhomalt is not None else "", v.clnsig,
                v.revstat, v.clndn, ";".join(v.flags), "low: " + ", ".join(lc) if lc else "high",
                note.replace("\t", " "))) + "\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--vcf", required=True, help="annotated VCF (SnpEff ANN, ClinVar CLN*, INFO/gnomAD_AF and the other gnomAD_* fields)")
    ap.add_argument("--sample-name")
    ap.add_argument("--gencc", help="GenCC submissions TSV")
    ap.add_argument("--constraint", help="gnomAD gene constraint TSV")
    ap.add_argument("--hpo-genes", help="HPO phenotype_to_genes.txt")
    ap.add_argument("--hpo-terms", help="patient HPO terms, comma/space separated (HP:0001250 ...)")
    ap.add_argument("--roh", help="bcftools roh output (RG lines)")
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
