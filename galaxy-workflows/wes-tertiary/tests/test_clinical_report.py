#!/usr/bin/env python3
"""Tests of the clinical report tool (tools/wes_clinical_report) on small synthetic inputs.

    python -m pytest tests/test_clinical_report.py      (or run this file directly)

Only the standard library is needed; the BAM/reference features (samtools) are not exercised here.
"""
import csv
import importlib.util
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
TOOL = os.path.join(HERE, "..", "tools", "wes_clinical_report")
DATA = os.path.join(TOOL, "test-data")

spec = importlib.util.spec_from_file_location("wes_clinical_report", os.path.join(TOOL, "wes_clinical_report.py"))
rep = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rep)


def run_report(extra=()):
    out = tempfile.mkdtemp()
    html, tsv = os.path.join(out, "r.html"), os.path.join(out, "r.tsv")
    rep.main(["--vcf", os.path.join(DATA, "annotated.vcf"), "--gencc", os.path.join(DATA, "gencc.tsv"),
              "--constraint", os.path.join(DATA, "constraint.tsv"), "--html", html, "--tsv", tsv, *extra])
    with open(tsv) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    with open(html) as f:
        return rows, f.read()


def categories(rows, gene):
    return {r["category"] for r in rows if r["gene"] == gene}


def test_acmg_sf_list():
    assert len(rep.ACMG_SF) == 84, "ACMG SF v3.3"
    assert rep.ACMG_SF["MUTYH"] == "biallelic" and rep.ACMG_SF["HFE"] == "hfe" and rep.ACMG_SF["TTN"] == "ttn"
    assert rep.ACMG_SF["PLN"] == "het" and rep.ACMG_SF["CYP27A1"] == "biallelic" and rep.ACMG_SF["ABCD1"] == "xl"
    assert set(rep.ACMG_SF) == set(rep.ACMG_SF_CATEGORY)


def test_sf_v33_genes():
    rows, html = run_report()
    assert categories(rows, "PLN") == {"diagnostic", "acmg_sf_reportable"}, "PLN P/LP het: reportable (dominant)"
    assert "acmg_sf_not" in categories(rows, "ABCD1"), "a single heterozygous ABCD1 variant is not an SF"
    assert "acmg_sf_reportable" not in categories(rows, "ABCD1")
    assert "ACMG SF v3.3 (84 genes)" in html


def test_gnomad_v41_frequency_rules():
    rows, html = run_report()
    assert categories(rows, "GENEF") == {"pathogenic_review"}, \
        "overall AF 0.0004 but grpmax FAF95 0.004: too common for a dominant disorder"
    f = next(r for r in rows if r["gene"] == "GENEF")
    assert f["gnomad_faf95"] == "0.004" and f["gnomad_grpmax_af"] == "0.006"
    assert "gnomAD site filtered (EXOMES_FILTERED)" in f["qc_flags"] and f["confidence"] == "high"
    assert categories(rows, "GENER") == {"candidate_inheritance_not_met"}, "6 homozygotes in gnomAD"
    assert "6 homozygotes in gnomAD" in next(r["note"] for r in rows if r["gene"] == "GENER")
    assert "FAF95 0.004 (afr)" in html


def test_freq():
    v = rep.Variant()
    v.af, v.faf = None, None
    assert rep.freq(v) is None
    v.af, v.faf = 0.02, 0.001
    assert rep.freq(v) == 0.02, "founder variants: FAF excludes bottlenecked groups"
    v.af, v.faf = 0.0004, 0.004
    assert rep.freq(v) == 0.004


def test_findings_by_inheritance():
    rows, html = run_report()
    assert categories(rows, "BRCA2") == {"diagnostic", "acmg_sf_reportable"}, "dominant P/LP: diagnostic + SF"
    assert categories(rows, "PAH") == {"carrier"}, "single P/LP allele in a recessive gene: carrier"
    assert "carrier" in categories(rows, "MUTYH") and "acmg_sf_reportable" not in categories(rows, "MUTYH"), \
        "MUTYH is reportable as a secondary finding only when biallelic"
    assert categories(rows, "GENEX") == {"candidate"}, "rare stop-gain in a LoF-intolerant dominant gene"
    assert not categories(rows, "TP53"), "benign common variant must not be reported"
    assert "Research use / decision support only" in html


def test_common_pathogenic_is_reviewed():
    rows, _ = run_report()
    assert categories(rows, "GENED") == {"pathogenic_review"}, "ClinVar P but gnomAD AF 0.2 must not be diagnostic"


def test_mixed_inheritance_homozygous_lof():
    rows, _ = run_report()
    assert categories(rows, "GENEM") == {"candidate"}, \
        "homozygous LoF in a gene with a recessive disorder fits, even if the gene is LoF-tolerant for dominant disease"


def test_artefact_list():
    rows, html = run_report()
    g = next(r for r in rows if r["gene"] == "GENEX")
    assert g["confidence"] == "high" and "recurrent artefact" not in g["qc_flags"]
    assert "No recurrent-artefact list was provided" in html
    rows, html = run_report(["--artefacts", os.path.join(DATA, "artefacts.vcf")])
    g = next(r for r in rows if r["gene"] == "GENEX")
    assert "recurrent artefact: called in 12/40" in g["qc_flags"], "matched on CHROM (chr prefix ignored)/POS/REF/ALT"
    assert g["confidence"].startswith("low") and "recurrent assay artefact" in g["confidence"]
    assert categories(rows, "GENEX") == {"candidate"}, "flagged, not removed"
    assert "2 alleles called in many" in html and "1 calls of this sample match it" in html


def test_gencc_formats(tmp_path=None):
    """The GenCC export keyed on sgc_id/version_number (current) and the legacy uuid one parse alike."""
    new = rep.load_gencc(os.path.join(DATA, "gencc.tsv"))
    d = tmp_path or tempfile.mkdtemp()
    legacy = os.path.join(str(d), "gencc_legacy.tsv")
    with open(os.path.join(DATA, "gencc.tsv")) as f, open(legacy, "w") as g:
        for i, line in enumerate(f):
            g.write(("uuid" if i == 0 else "x") + "\t" + line.split("\t", 2)[2])
    assert new and dict(new) == dict(rep.load_gencc(legacy))
    assert ("phenylketonuria", "Autosomal recessive", "Definitive", "Test submitter") in new["PAH"]


def test_predictor_candidates():
    rows, html = run_report()
    p = next(r for r in rows if r["gene"] == "GENEP")
    assert p["category"] == "candidate" and "PP3_Strong" in p["note"] and p["alphamissense"] == "0.995"
    s = next(r for r in rows if r["gene"] == "GENES")
    assert s["category"] == "candidate" and "PP3_Supporting" in s["note"] and s["spliceai"] == "0.71"
    assert "AlphaMissense 0.995" in html and "Scores present for this sample: AlphaMissense 1, SpliceAI 2" in html


def test_predictor_codes():
    assert [rep.am_code(x) for x in (0.990, 0.989, 0.906, 0.905, 0.792, 0.791, 0.170, 0.169, 0.100, 0.099)] == [
        "PP3_Strong", "PP3_Moderate", "PP3_Moderate", "PP3_Supporting", "PP3_Supporting", None, None,
        "BP4_Supporting", "BP4_Supporting", "BP4_Moderate"], "Bergquist 2025 intervals"
    assert rep.spliceai_max("A|G1|0.1|0.3|0|0|1|2|3|4,T|G1|0.9|0|0|0|1|1|1|1,A|G2|0.05|0|0|0.4|1|1|1|1", "A") == (0.4, "G2")

    def mis(am, sai, effect="missense_variant"):
        v = rep.Variant()
        v.anns = [["A", effect, "MODERATE", "G", "G", "transcript", "T", "protein_coding", "1/2", "c.1A>T", "p.X1Y"]]
        v.am, v.sai = am, sai
        return rep.predictor_codes(v, "G")
    assert mis(0.95, 0.3) == ["PP3_Moderate"], "missense and splice PP3 are not summed"
    assert mis(0.05, 0.5) == ["PP3_Supporting"], "a splice prediction overrides a benign missense score"
    assert mis(0.05, 0.15) == [], "missense BP4 needs SpliceAI <= 0.1"
    assert mis(0.05, None) == ["BP4_Moderate"]
    assert mis(None, 0.9, "splice_donor_variant") == [], "canonical +/-1,2: PVS1 decision tree, not PP3"


def test_candidate_evidence_codes():
    rows, _ = run_report()
    g = next(r for r in rows if r["gene"] == "GENEX")
    assert "PVS1" in g["note"] and "PM2_supporting" in g["note"]
    assert g["acmg_class"] in ("likely pathogenic", "pathogenic") and int(g["acmg_points"]) >= 6
    assert all(r["acmg_class"] == "" for r in rows if r["category"].startswith("pgx"))


def _var(effect, gene="G", aa="", info=None, af=None, nhomalt=None, hgvs_c="c.1A>T"):
    v = rep.Variant()
    v.anns = [["A", effect, "HIGH", gene, gene, "transcript", "T", "protein_coding", "1/2", hgvs_c, "p.X1Y", "", "", aa]]
    v.info, v.af, v.faf, v.nhomalt, v.am, v.sai, v.chrom = info or {}, af, None, nhomalt, None, None, "1"
    return v


def test_acmg_points_and_classes():
    cl = rep.acmg_classify
    assert cl(["PVS1", "PM2_supporting"]) == ("likely pathogenic", 9), "PVS1 + PM2_Supporting is LP (Tavtigian)"
    assert cl(["PVS1", "PS1"]) == ("pathogenic", 12)
    assert cl(["PM2_supporting", "PP3_Moderate", "PM5"]) == ("uncertain significance", 5)
    assert cl(["PM2_supporting", "PP3_Strong", "PM1_Supporting", "PM5"]) == ("likely pathogenic", 7), \
        "PM1 + PP3 count at most 4 points"
    assert cl(["BS1", "BP4_Moderate", "BP7"]) == ("benign", -7)
    assert cl(["BS2_Supporting", "BP4_Supporting"]) == ("likely benign", -2)
    assert cl(["BA1"]) == ("benign", 0) and cl(["BA1", "PP3_Strong"])[0] == "conflicting (BA1)"
    assert [rep.code_points(c) for c in ("PVS1_Strong", "PVS1_Moderate", "PM3_Supporting", "PS2_Moderate", "BP3")] == \
        [4, 2, 1, 2, -1]


def test_acmg_criteria():
    nmd = _var("stop_gained", aa="10/400", info={"NMD": "(G|G|1|1.00)"})
    assert rep.pvs1_code(nmd, "G") == "PVS1"
    assert rep.pvs1_code(_var("stop_gained", aa="100/400"), "G") == "PVS1_Strong", "escapes NMD, removes 75 %"
    assert rep.pvs1_code(_var("frameshift_variant", aa="390/400"), "G") == "PVS1_Moderate", "removes < 10 %"
    assert rep.pvs1_code(_var("start_lost"), "G") == "PVS1_Moderate"
    assert rep.pvs1_code(_var("splice_donor_variant&intron_variant"), "G") == "PVS1"
    ev = lambda v, moi={"AR"}: rep.acmg_evidence(v, "G", moi, {}, False, 0.01, {})
    assert "BA1" in ev(_var("missense_variant", af=0.2))
    hfe = rep.acmg_evidence(_var("missense_variant", gene="HFE", af=0.06, hgvs_c="c.845G>A"), "HFE", {"AR"}, {}, False,
                            0.01, {})
    assert "BA1" not in hfe and "BS1" not in hfe, "ClinGen BA1 exception list"
    assert "BS1" in ev(_var("missense_variant", af=0.002), {"AD"}) and "BS1" not in ev(_var("missense_variant", af=0.002))
    assert "BS2_Supporting" in ev(_var("missense_variant", af=0.004, nhomalt=3))
    assert "BS2_Supporting" not in ev(_var("missense_variant", af=0.004, nhomalt=3), {"AD", "AR"})
    inf = _var("inframe_deletion")
    assert "PM4" in ev(inf)
    inf.repeat = "tandem repeat ((GCA) x6)"
    assert "BP3" in ev(inf) and "PM4" not in ev(inf)
    syn = _var("synonymous_variant")
    syn.sai = 0.02
    assert "BP7" in ev(syn) and "BP4_Supporting" not in ev(syn)
    early = _var("frameshift_variant", aa="100/400")
    early.anns[0][8] = "3/10"
    assert rep.pvs1_code(early, "G") == "PVS1", "stop before the penultimate exon: NMD (50-nt rule)"
    minor = _var("stop_gained", aa="20/120")  # LoF on a short isoform only; the principal transcript is missense
    minor.anns.append(["A", "missense_variant", "MODERATE", "G", "G", "transcript", "T2", "protein_coding", "2/9",
                       "c.60A>T", "p.K20N", "", "", "20/800"])
    assert rep.main_ann(minor, "G")[6] == "T2"
    assert not [c for c in ev(minor) if c.startswith("PVS1")]
    assert "PVS1" in ev(_var("splice_donor_variant&intron_variant", hgvs_c="c.100+1G>A"))
    assert not [c for c in ev(_var("splice_donor_variant&intron_variant", hgvs_c="c.-187+2T>C")) if
                c.startswith("PVS1")], "splice site of a 5' UTR intron"


def test_pm1_from_clinvar_residues():
    v = _var("missense_variant")
    v.anns[0][10] = "p.Arg50Trp"
    v.pos, v.ref, v.alt = 100, "C", "T"
    cvp = {("G", r): [("Gly", "Asp", "P", 2, "1", "1", 1, "G", "A")] for r in (47, 48, 52)}
    assert "PM1_Supporting" in rep.protein_evidence(v, "G", cvp)[0]
    cvp[("G", 53)] = [("Ala", "Thr", "B", 2, "2", "1", 2, "G", "A")]
    assert "PM1_Supporting" not in rep.protein_evidence(v, "G", cvp)[0], "a benign missense nearby: no hot spot"
    cvp2 = {("G", 50): [("Arg", "Gln", "P", 2, "1", "1", 1, "G", "A")] * 3}
    codes = rep.protein_evidence(v, "G", cvp2)[0]
    assert "PM5" in codes and "PM1_Supporting" not in codes, "the variant's own residue counts for PM5 only"


def test_pharmacogenomics_vcf_only():
    rows, html = run_report()
    pgx = {r["gene"]: (r["genotype"], r["note"]) for r in rows if r["category"] == "pgx"}
    assert pgx["CYP2C19"] == ("*2/*2", "Poor metabolizer")
    assert pgx["CYP2C9"][1].startswith("Normal metabolizer")
    assert all(len(r) == 31 and None not in r.values() for r in rows)
    assert "absence of a variant is assumed to be the reference allele" in html


def test_pharmcat_calls():
    rows, html = run_report(["--pharmcat", os.path.join(DATA, "pharmcat.tsv")])
    pc = {r["gene"]: (r["genotype"], r["note"]) for r in rows if r["category"] == "pgx_pharmcat"}
    assert pc["CYP2C19"] == ("*1/*2", "Intermediate Metabolizer"), "PharmCAT 3.4.0 on HG001 (GeT-RM: *1/*2)"
    assert pc["CYP3A5"] == ("*3/*3", "Poor Metabolizer")
    assert "<h3>PharmCAT 3.4.0</h3>" in html and "Missing positions" in html
    assert all(len(r) == 31 for r in rows)
    _, html = run_report()
    assert "PharmCAT" not in html.split("<h2>7. Pharmacogenomics")[1].split("<h2>8.")[0]


def test_hpo_ranking(tmp_path=None):
    d = tmp_path or tempfile.mkdtemp()
    p = os.path.join(str(d), "p2g.txt")
    with open(p, "w") as f:
        f.write("hpo_id\thpo_name\tncbi_gene_id\tgene_symbol\tdisease_id\nHP:0000001\tExample\t1\tGENEX\tOMIM:1\n")
    rows, html = run_report(["--hpo-genes", p, "--hpo-terms", "HP:0000001"])
    assert "phenotype match: Example" in html
    note = next(r["note"] for r in rows if r["gene"] == "GENEX")
    assert "PP4" in note


def test_low_confidence():
    def var(ref, alt, zyg="het", ab=0.5, qual=200.0, flags=()):
        v = rep.Variant()
        v.ref, v.alt, v.zyg, v.ab, v.qual, v.flags = ref, alt, zyg, ab, qual, list(flags)
        return v
    assert rep.low_confidence(var("C", "T")) == []
    assert rep.low_confidence(var("C", "T", ab=0.30, qual=36, flags=["low depth (10x)"])) == [], "TPP1-like true call"
    assert "alt fraction < 0.25" in rep.low_confidence(var("C", "T", ab=0.17))
    assert "QUAL < 30" in rep.low_confidence(var("C", "T", qual=20))
    assert "indel in a repeat" in rep.low_confidence(var("CT", "C", ab=0.28, flags=["homopolymer (T x10)"]))
    assert rep.low_confidence(var("CT", "C", flags=["homopolymer (T x6)"])) == []
    assert rep.low_confidence(var("C", "T", flags=["segmental duplication / pseudogene: short-read calls unreliable"]))


def test_repeat_flag():
    assert rep.repeat_flag("ACGTTTTTTTTGCA", 3, 2) == "homopolymer (T x8)"
    assert rep.repeat_flag("GGACACACACACTT", 3, 1).startswith("tandem repeat")
    assert rep.repeat_flag("ACGTACGGTCAT", 5, 1) is None


def test_pgx_rules():
    # DPYD: one no-function and one decreased-function allele (assumed in trans) -> activity score 0.5
    dpyd = [("*2A (c.1905+1G>A)", "", "", 1, None, ""), ("*13 (c.1679T>G)", "", "", 0, None, ""),
            ("c.2846A>T", "", "", 1, None, ""), ("HapB3 (c.1236G>A)", "", "", 0, None, "")]
    assert rep.pgx_pheno("DPYD", dpyd) == ("activity score 0.5", "Poor metabolizer (two variants assumed in trans)")
    cyp2c19 = [("*2", "", "", 1, None, ""), ("*3", "", "", 0, None, ""), ("*17", "", "", 1, None, "")]
    assert rep.pgx_pheno("CYP2C19", cyp2c19) == ("*2/*17", "Intermediate metabolizer")
    cyp2c9 = [("*2", "", "", 1, None, ""), ("*3", "", "", 1, None, "")]
    assert rep.pgx_pheno("CYP2C9", cyp2c9) == ("*2/*3", "Poor metabolizer (activity score 0.5)")
    tpmt = [("*2 (c.238G>C)", "", "", 0, None, ""), ("*3B (c.460G>A)", "", "", 1, None, ""), ("*3C (c.719A>G)", "", "", 1, None, "")]
    assert rep.pgx_pheno("TPMT", tpmt) == ("*3A/*1", "Intermediate metabolizer")
    unknown = [("*2", "", "", None, None, ""), ("*3", "", "", 0, None, ""), ("*17", "", "", 0, None, "")]
    assert rep.pgx_pheno("CYP2C19", unknown)[0] == "indeterminate"


def test_protein_evidence():
    cvp = rep.load_clinvar_protein(os.path.join(DATA, "clinvar_protein.tsv"))
    _s, _c, variants, _p = rep.parse_vcf(os.path.join(DATA, "annotated.vcf"))
    v = next(x for x in variants if x.gene == "GENEP")
    codes, notes = rep.protein_evidence(v, "GENEP", cvp)
    assert codes == ["PS1", "PM5"], "p.Arg67His pathogenic by another nucleotide change; p.Arg67Cys pathogenic"
    assert "p.Arg67Cys" in notes[1] and "VCV000001001" in notes[0]
    assert any("3 pathogenic, 1 benign" in n for n in notes), "0-star p.Ala68Val is not counted"
    f = next(x for x in variants if x.gene == "GENEF")
    assert rep.protein_evidence(f, "GENEF", cvp) == ([], ["ClinVar residue 4 is Gly, not Arg: transcripts numbered differently"])
    v.sai = 0.3
    assert rep.protein_evidence(v, "GENEP", cvp)[0] == ["PM5"], "no PS1 if the variant may act through splicing"
    rows, html = run_report(["--clinvar-protein", os.path.join(DATA, "clinvar_protein.tsv")])
    note = next(r["note"] for r in rows if r["gene"] == "GENEP")
    assert "PS1" in note and "PM5" in note and ">PS1</span>" in html and "pathogenic cluster" not in html


def test_clingen_dosage():
    d = rep.load_dosage(os.path.join(DATA, "clingen_dosage.tsv"))
    assert d == {"GENEX": (40, 0), "BRCA2": (3, 1), "PAH": (30, None)}
    assert rep.lof_mechanism("BRCA2", {"AD"}, {}, d), "HI 3 without constraint data"
    assert not rep.lof_mechanism("GENEX", {"AD"}, {"GENEX": (0.2, 1.0, 0.0)}, d), "HI 40 overrides LOEUF"
    rows, html = run_report()
    assert "PVS1" in next(r["note"] for r in rows if r["gene"] == "GENEX"), "LOEUF < 0.6 in the test constraint file"
    rows, html = run_report(["--clingen-dosage", os.path.join(DATA, "clingen_dosage.tsv")])
    g = next(r for r in rows if r["gene"] == "GENEX")
    assert g["category"] == "candidate_inheritance_not_met" and "haploinsufficiency unlikely" in g["note"]
    assert "PVS1" not in g["note"]
    assert "ClinGen haploinsufficiency: sufficient evidence" in html


def test_disputed_gencc(tmp_path=None):
    d = tmp_path or tempfile.mkdtemp()
    p = os.path.join(str(d), "gencc.tsv")
    with open(os.path.join(DATA, "gencc.tsv")) as f, open(p, "w") as out:
        out.write(f.read())
        out.write("SGC-9\t1\tHGNC:9\tGENEP\tMONDO:9\tGene P cardiomyopathy\tRefuted Evidence\tAutosomal dominant\tLab X\n")
    assert rep.load_gencc_disputed(p) == {"GENEP": [("Gene P cardiomyopathy", "Refuted Evidence", "Lab X")]}
    out = tempfile.mkdtemp()
    rep.main(["--vcf", os.path.join(DATA, "annotated.vcf"), "--gencc", p, "--html", os.path.join(out, "r.html"),
              "--tsv", os.path.join(out, "r.tsv")])
    assert "refuted: Gene P cardiomyopathy (Lab X)" in open(os.path.join(out, "r.html")).read()


def test_disease_phenotype_match():
    terms = {"HP:0001250", "HP:0001263"}
    ph = rep.load_hpo(os.path.join(DATA, "p2g.txt"), terms)
    assert ph.gene_terms["GENEP"] == terms and ph.gene_all["PAH"] == {"OMIM:100003"}
    best = ph.best_diseases("GENEP", 2)
    assert best[0][1] == "OMIM:100001" and abs(best[0][0] - 1.0) < 1e-9, "both terms explained"
    assert best[1][1] == "OMIM:100002" and 0 < best[1][0] < 1
    assert ph.ic["HP:0001263"] > ph.ic["HP:0001250"], "the rarer term carries more information"
    hpoa = rep.load_hpoa(os.path.join(DATA, "phenotype.hpoa"))
    assert hpoa["OMIM:100001"] == {"name": "Gene P syndrome", "onset": ["infantile"], "course": [],
                                   "inheritance": ["AD"], "n_features": 2}
    rows, html = run_report(["--hpo-genes", os.path.join(DATA, "p2g.txt"), "--hpo-terms", "HP:0001250,HP:0001263",
                             "--hpo-annotations", os.path.join(DATA, "phenotype.hpoa")])
    assert "best match: <b>Gene P syndrome</b> (OMIM:100001): 2 of 2 terms, 100% of the phenotype information; " \
           "onset infantile; AD" in html
    pm = [r for r in rows if r["category"] == "phenotype_match"]
    assert pm and pm[0]["gene"] == "GENEP" and pm[0]["genotype"] == "OMIM:100001"
    assert "Diseases matching the phenotype in genes with a finding" in html
    _rows, html = run_report(["--hpo-genes", os.path.join(DATA, "p2g.txt"),
                              "--hpo-annotations", os.path.join(DATA, "phenotype.hpoa")])
    assert "HPO: Phenylketonuria: onset congenital; AR" in html, "disease features without patient terms"


def sam(name, pos, cigar, seq):
    return "\t".join([name, "99", "1", str(pos), "60", cigar, "=", str(pos + 100), "200", seq, "I" * len(seq)])


def test_read_alleles():
    def var(pos, ref, alt):
        v = rep.Variant()
        v.chrom, v.pos, v.ref, v.alt = "1", pos, ref, alt
        return v
    snv, dele, ins = var(103, "C", "T"), var(106, "GA", "G"), var(108, "A", "AT")
    # ref: 101 AAC GTG AAC 109; read with T at 103, deletion of A at 107 and insertion of T after 108
    r = rep.read_alleles(101, "6M1D1M1I2M", "AATGTGATCC", (snv, dele, ins))
    assert r == {snv.key(): "alt", dele.key(): "alt", ins.key(): "alt"}
    r = rep.read_alleles(101, "10M", "AACGTGAACC", (snv, dele, ins))
    assert r == {snv.key(): "ref", dele.key(): "ref", ins.key(): "ref"}
    assert rep.read_alleles(104, "5M", "GTGAA", (snv,)) == {snv.key(): None}, "not covered"


def test_phase_pair():
    def var(pos, ref, alt):
        v = rep.Variant()
        v.chrom, v.pos, v.ref, v.alt, v.zyg = "1", pos, ref, alt, "het"
        return v
    a, b = var(103, "C", "T"), var(110, "G", "A")
    trans = [sam(f"r{i}", 101, "12M", "AATGTGAACGTT") for i in range(3)] + \
            [sam(f"s{i}", 101, "12M", "AACGTGAACATT") for i in range(3)]
    old = rep.run_samtools
    try:
        rep.run_samtools = lambda args: "\n".join(trans)
        assert rep.phase_pair("x.bam", "", a, b) == ("trans", 0, 6)
        rep.run_samtools = lambda args: "\n".join([sam(f"r{i}", 101, "12M", "AATGTGAACATT") for i in range(5)])
        assert rep.phase_pair("x.bam", "", a, b)[0] == "cis"
        assert rep.phase_pair("x.bam", "", a, var(5000, "G", "A"))[0] == "too far apart"
        assert "in trans" in rep.phase_text(("trans", 0, 6)) and "test the parents" in rep.phase_text(("unclear", 1, 1))
    finally:
        rep.run_samtools = old


def test_decompose():
    assert rep.decompose(100, "TCC", "TTC") == [(101, "C", "T")], "MNP split into SNVs"
    assert rep.decompose(100, "ACGT", "AGGA") == [(101, "C", "G"), (103, "T", "A")]
    assert rep.decompose(100, "ACCT", "ACT") == [(100, "AC", "A")], "suffix then prefix trimmed, anchor kept"
    assert rep.decompose(100, "G", "GA") == [(100, "G", "GA")]
    p = rep.load_parent_vcf(os.path.join(DATA, "mother.vcf"))
    assert p[("1", 5000006, "C", "T")][0] == "het" and p[("1", 5000000, "ACGTAAC", "G")][0] == "het", "1/2 split"


def trio_report(extra=()):
    return run_report(["--mother-vcf", os.path.join(DATA, "mother.vcf"), "--father-vcf", os.path.join(DATA, "father.vcf"),
                       "--hpo-genes", os.path.join(DATA, "p2g.txt"), "--hpo-terms", "HP:0001250,HP:0001263", *extra])


def test_trio_origin():
    rows, html = trio_report(["--affected-parents", "none"])
    x = next(r for r in rows if r["gene"] == "GENEX")
    assert x["category"] == "candidate_inheritance_not_met" and "inherited from the unaffected father" in x["note"]
    assert "inherited from the mother" in next(r["note"] for r in rows if r["gene"] == "PAH")
    assert "inherited from the father" in next(r["note"] for r in rows if r["gene"] == "MUTYH")
    p = next(r for r in rows if r["gene"] == "GENEP" and r["category"].startswith("candidate"))
    assert "possible de novo" in p["note"] and "PM6_Supporting" in p["note"], "no parental reads: assumed de novo"
    assert any(r["category"] == "de_novo_possible" and r["gene"] == "GENEP" for r in rows)
    assert "De novo variants (trio)" in html and "Too few common sites" in html
    rows, _ = trio_report(["--affected-parents", "father"])
    x = next(r for r in rows if r["gene"] == "GENEX")
    assert x["category"] == "candidate" and "inherited from the affected father" in x["note"]
    rows, html = run_report(["--father-vcf", os.path.join(DATA, "father.vcf")])
    assert "Only one parent given" in html and not any(r["category"].startswith("de_novo") for r in rows)


def test_trio_denovo_cluster():
    """three parent-absent changes within 10 bp (reads of a paralogous copy): possibly de novo, no de novo points"""
    out = tempfile.mkdtemp()
    vcf = os.path.join(out, "cluster.vcf")
    with open(os.path.join(DATA, "annotated.vcf")) as f, open(vcf, "w") as w:
        for line in f:
            w.write(line)
            if line.startswith("2\t2000100\t"):
                for pos, ref, alt in ((2000104, "C", "T"), (2000108, "G", "T")):
                    t = line.split("\t")
                    t[1], t[3], t[4] = str(pos), ref, alt
                    t[7] = t[7].replace("p.Arg67His", f"p.Gly{pos - 2000100 + 67}Val").replace("SpliceAI=A|", f"SpliceAI={alt}|")
                    t[7] = t[7].replace("ANN=A|", f"ANN={alt}|")
                    w.write("\t".join(t))
                t = line.split("\t")  # a heterozygous call without reference reads (other allele: another variant)
                t[1], t[9] = "3000100", "0/1:40:0,40\n"
                w.write("\t".join(t))
    html, tsv = os.path.join(out, "r.html"), os.path.join(out, "r.tsv")
    rep.main(["--vcf", vcf, "--gencc", os.path.join(DATA, "gencc.tsv"), "--constraint", os.path.join(DATA, "constraint.tsv"),
              "--mother-vcf", os.path.join(DATA, "mother.vcf"), "--father-vcf", os.path.join(DATA, "father.vcf"),
              "--hpo-genes", os.path.join(DATA, "p2g.txt"), "--hpo-terms", "HP:0001250,HP:0001263",
              "--affected-parents", "none", "--html", html, "--tsv", tsv])
    rows = list(csv.DictReader(open(tsv), delimiter="\t"))
    dn = [r for r in rows if r["category"].startswith("de_novo") and r["gene"] == "GENEP"]
    assert len(dn) == 4 and all(r["category"] == "de_novo_possible" for r in dn)
    assert sum("paralogous copy" in r["note"] for r in dn) == 3
    assert "no reference reads" in next(r["note"] for r in dn if r["pos"] == "3000100")
    assert not any("PM6" in r["note"] or "PS2" in r["note"] for r in rows if r["gene"] == "GENEP")
    assert "4 listed variants lie in a cluster" in open(html).read()
    s = [rep.Variant() for _ in range(4)]
    for v, p in zip(s, (100, 104, 108, 200)):
        v.chrom, v.pos = "1", p
    assert rep.denovo_clusters(s) == s[:3] and rep.denovo_clusters(s[:2] + s[3:]) == []
    m = rep.Variant()
    m.zyg, m.dp, m.ab = "het", 120, 0.267
    assert rep.mosaic_p(m) < 1e-3, "32 of 120 reads: below a heterozygote"
    m.dp, m.ab = 40, 0.35
    assert rep.mosaic_p(m) > 1e-3, "14 of 40 reads: consistent with a heterozygote"


def test_trio_reads_and_phase():
    def var(chrom, pos, ref, alt, zyg="het"):
        v = rep.Variant()
        v.chrom, v.pos, v.ref, v.alt, v.zyg = chrom, pos, ref, alt, zyg
        return v
    t = rep.Trio(os.path.join(DATA, "mother.vcf"), os.path.join(DATA, "father.vcf"), "m.bam", "f.bam")
    counts = {("mother", "2:2000100:G>A"): (0, 30), ("father", "2:2000100:G>A"): (1, 25),
              ("mother", "3:5:G>A"): (0, 8), ("father", "3:5:G>A"): (0, 30),
              ("mother", "3:6:G>A"): (6, 20), ("father", "3:6:G>A"): (0, 30)}
    t.counts = lambda who, v: counts.get((who, v.key()))
    assert t.origin(var("2", 2000100, "G", "A"), "female")[0] == "de novo", "parents' reads: 0 and 1 alt"
    assert t.origin(var("3", 5, "G", "A"), "female")[0] == "possible de novo", "mother has too few reads"
    assert t.origin(var("3", 6, "G", "A"), "female")[0] == "maternal", "mother's reads carry it, not called"
    assert t.origin(var("1", 1000100, "C", "T", "hom"), "female")[0] == "unknown"
    m, f = var("12", 102840474, "T", "C"), var("1", 1000100, "C", "T")
    m.gene = f.gene = "G"
    assert rep.trio_phase(t, m, f, "female") == "trans", "one from each parent"
    m2 = var("1", 45331556, "C", "T")
    assert rep.trio_phase(t, f, m2, "female") == "cis", "both paternal"
    assert rep.de_novo_codes("de novo", True, True) == ["PS2_Moderate"]
    assert rep.de_novo_codes("de novo", True, True, parentage_ok=False) == []
    assert rep.de_novo_codes("possible de novo", True, True) == ["PM6_Supporting"]
    assert rep.de_novo_codes("de novo", False, True) == [] and rep.de_novo_codes("de novo", True, False) == []


if __name__ == "__main__":
    failed = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"ok   {name}")
            except AssertionError as e:
                failed += 1
                print(f"FAIL {name}: {e}")
    sys.exit(1 if failed else 0)


def _cnv_inputs(d):
    """Synthetic wes_cnv outputs: calls, summary, ClinGen regions and both parents' per-target ratios."""
    cols = ["sample", "chrom", "start", "end", "length_bp", "type", "copy_number", "normal_copies", "targets",
            "mean_ratio", "log10_bf", "segdup_fraction", "panel_carriers", "confidence", "genes"]
    calls = [
        ("13", 32315000, 32400000, "DEL", 1, 27, 0.5, 40.0, 0.0, 0, "high", "BRCA2:27/27"),  # whole HI-3 gene
        ("12", 102860000, 102900000, "DEL", 1, 3, 0.48, 9.0, 0.0, 0, "high", "PAH:3/13"),  # + PAH P het: comp. het?
        ("1", 1000, 5000, "DEL", 1, 4, 0.5, 8.0, 0.0, 0, "high", "GENEX:2/9"),  # inside a benign population region
        ("22", 22100000, 22200000, "DEL", 1, 6, 0.5, 12.0, 0.0, 0, "high", "IGLL5:1/1"),  # IGL locus
        ("16", 74410535, 74412246, "DEL", 1, 4, 0.52, 7.0, 0.9, 0, "low: segmental duplication", "CLEC18B:2/4"),
    ]
    p = os.path.join(d, "calls.tsv")
    with open(p, "w") as f:
        f.write("##wes_cnv 0.1.0; panel 6 samples\n#" + "\t".join(cols) + "\n")
        for c, s, e, t, cn, n, r, bf, sdf, pc, conf, genes in calls:
            f.write("\t".join(map(str, ["S", c, s, e, e - s, t, cn, 2, n, r, bf, sdf, pc, conf, genes])) + "\n")
    js = os.path.join(d, "summary.json")
    with open(js, "w") as f:
        f.write('{"status": "analysed", "panel": ["A", "B", "C", "D", "E", "F"], "targets_used": 190000, '
                '"ratio_sd_autosomes": 0.16, "sex_complement": "XY", "chromosome_dosage": {"21": 1.48, "1": 1.0}, '
                '"aneuploidy": ["chromosome 21: gain of a whole copy (trisomy), median ratio 1.48"]}')
    reg = os.path.join(d, "regions.tsv")
    with open(reg, "w") as f:
        f.write("#ClinGen Region Curation Results\n#ISCA ID\tISCA Region Name\tcytoBand\tGenomic Location\t"
                "Haploinsufficiency Score\tHaploinsufficiency Description\tTriplosensitivity Score\n")
        f.write("ISCA-1\t1p36 population region\t1p36.33\tchr1:500-9000\t40\tDosage sensitivity unlikely\t40\n")
        f.write("ISCA-2\t13q13 test region\t13q13.1\tchr13:32300000-32410000\t3\tSufficient evidence\t1\n")
    parents = []
    for who, r_brca2 in (("mother", 1.01), ("father", 0.98)):
        t = os.path.join(d, f"{who}.targets.tsv")
        with open(t, "w") as f:
            f.write("chrom\tstart\tend\tdepth\tpanel_depth\tratio\tsd\tused\tcall_copy\tgc\n")
            for i in range(27):  # BRCA2: normal in both parents -> de novo
                s = 32315000 + i * 3000
                f.write(f"13\t{s}\t{s + 150}\t50\t50\t{r_brca2}\t0.1\t1\t2\t0.4\n")
            for i in range(3):  # PAH: the mother carries the deletion
                s = 102860000 + i * 10000
                f.write(f"12\t{s}\t{s + 150}\t50\t50\t{0.5 if who == 'mother' else 1.0}\t0.1\t1\t"
                        f"{1 if who == 'mother' else 2}\t0.4\n")
        parents += [f"--{who}-cnv", t]
    return ["--cnv-calls", p, "--cnv-summary", js, "--clingen-regions", reg,
            "--clingen-dosage", os.path.join(DATA, "clingen_dosage.tsv"), *parents, "--affected-parents", "none"]


def test_cnv_interpretation():
    d = tempfile.mkdtemp()
    rows, html = run_report(_cnv_inputs(d))
    cnv = {r["gene"]: r for r in rows if r["category"].startswith("cnv")}
    b = cnv["BRCA2"]
    assert b["category"] == "cnv_pathogenic" and b["chrom"] == "13" and b["alt"] == "<DEL>", b
    assert "2A +1.00" in b["note"] and "5A" in b["note"] and "de novo" in b["note"], "whole HI-3 gene/region, de novo"
    assert "heterozygous loss in a dominant gene" in b["note"]
    p = cnv["PAH"]
    assert "possible compound heterozygote" in p["note"] and "inherited (unaffected)" in p["note"] \
        and "5B -0.30" in p["note"], "PAH deletion from the mother plus a pathogenic PAH variant"
    assert "2F -1.00" in cnv["GENEX"]["note"] and cnv["GENEX"]["category"] == "cnv_high_confidence", \
        "inside a ClinGen population region: benign evidence, not of interest"
    assert cnv["IGLL5"]["confidence"].startswith("low") and "IGL locus" in cnv["IGLL5"]["note"]
    assert cnv["CLEC18B"]["category"] == "cnv_low_confidence"
    assert "<h2>4. Copy-number variants and chromosome dosage</h2>" in html
    assert "Down syndrome" in html and "Chromosome dosage</b>: <b>chromosome 21" in html
    assert "<h2>5. Secondary findings" in html and "<h2>13. Methods" in html


def test_cnv_not_analysed():
    d = tempfile.mkdtemp()
    p = os.path.join(d, "calls.tsv")
    with open(p, "w") as f:
        f.write("##wes_cnv 0.1.0; not analysed: 1 panel samples other than S (at least 3 of the same kit are needed)\n"
                "#sample\tchrom\tstart\tend\n")
    rows, html = run_report(["--cnv-calls", p])
    assert not [r for r in rows if r["category"].startswith("cnv")]
    assert "Not analysed: not analysed: 1 panel samples" not in html and "at least 3 of the same kit" in html
    rows, html = run_report()
    assert "no copy-number calls given" in html


def _upd_trio(event=None, seed=3, n=300, step=200_000, chroms=("1", "2", "7", "15", "21")):
    """Synthetic trio of common SNVs (AF 0.3) spaced `step` apart. event: (chrom, kind, donor) with kind iso | het |
    tri (two donor alleles plus one of the other parent). Returns (child variants, Trio, Covered)."""
    import random
    rng = random.Random(seed)
    trio = rep.Trio()
    trio.vcf = {"mother": {}, "father": {}}
    kids, rows = [], []
    for c in chroms:
        for i in range(n):
            pos = 1_000_000 + i * step
            k = (c, pos, "A", "G")
            g = {w: (rng.random() < 0.3) + (rng.random() < 0.3) for w in ("mother", "father")}
            for w, a in g.items():
                if a:
                    trio.vcf[w][k] = ("hom" if a == 2 else "het", 40, a / 2, 60.0)
            take = lambda w: rng.random() < g[w] / 2
            if event and event[0] == c:
                donor, other = event[2], rep.OTHER[event[2]]
                alt, tot = {"iso": (2 * take(donor), 2), "het": (g[donor], 2),
                            "tri": (g[donor] + take(other), 3)}[event[1]]
            else:
                alt, tot = take("mother") + take("father"), 2
            rows.append((c, pos - 11, pos + 10, "t", {"20X": 21}))
            if alt:
                v = rep.Variant()
                v.chrom, v.pos, v.ref, v.alt, v.af, v.qual, v.dp = c, pos, "A", "G", 0.3, 60.0, 40
                v.zyg, v.ab, v.anns, v.gene = ("hom" if alt == tot else "het"), alt / tot, [], ""
                kids.append(v)
    return kids, trio, rep.Covered(rows)


def test_upd_isodisomy_and_heterodisomy():
    kids, trio, cov = _upd_trio()
    assert rep.upd_scan(kids, trio, "male", cov)[0] == [], "a normal trio: no UPD"
    kids, trio, cov = _upd_trio(("15", "iso", "mother"))
    segs, _ = rep.upd_scan(kids, trio, "male", cov)
    assert [(d["chrom"], d["from"], d["kind"]) for d in segs] == [("15", "mother", "isodisomy")]
    assert segs[0]["disorder"].startswith("Prader-Willi") and segs[0]["end"] - segs[0]["start"] > 40_000_000
    assert "maternal isodisomy of chromosome 15" in rep.upd_label(segs[0])
    kids, trio, cov = _upd_trio(("7", "het", "father"))
    segs, _ = rep.upd_scan(kids, trio, "male", cov)
    assert [(d["chrom"], d["from"], d["kind"]) for d in segs] == [("7", "father", "heterodisomy")]
    assert segs[0]["disorder"] == "", "paternal UPD7 has no known imprinting disorder"
    kids, trio, cov = _upd_trio(("15", "iso", "mother"))
    cnv = [{"chrom": "15", "start": 0, "end": 100_000_000, "type": "DEL", "confidence": "high"}]
    assert rep.upd_scan(kids, trio, "male", cov, cnv)[0][0]["deletion"], "a CNV deletion explains it"


def test_upd_ignores_error_clusters():
    kids, trio, cov = _upd_trio()
    for i in range(12):  # 12 Mendelian errors within 60 kb (a paralogous locus): one event, not UPD
        v = rep.Variant()
        v.chrom, v.pos, v.ref, v.alt, v.af, v.qual, v.dp = "2", 5_000_001 + 5000 * i, "C", "T", 0.3, 60.0, 40
        v.zyg, v.ab, v.anns, v.gene = "hom", 1.0, [], ""
        kids.append(v)
    assert rep.upd_scan(kids, trio, "male", cov)[0] == []


def test_aneuploidy_and_x_origin():
    kids, trio, _cov = _upd_trio(("21", "tri", "mother"), n=400)
    share, n = rep.gain_origin(kids, trio, "21")
    assert 0.6 < share < 0.75 and n >= 10 and "extra copy maternal" in rep.gain_origin_text((share, n))
    kids, trio, _cov = _upd_trio(("21", "tri", "father"), n=400)
    assert "extra copy paternal" in rep.gain_origin_text(rep.gain_origin(kids, trio, "21"))
    kids, trio, _cov = _upd_trio(n=400)
    assert "unclear" in rep.gain_origin_text(rep.gain_origin(kids, trio, "21")), "two copies: allele balance 0.5"
    assert "not determined" in rep.gain_origin_text(None)
    kx, tx, cx = _upd_trio(chroms=("X",))
    xp = rep.x_parents(kx, tx, cx)
    assert xp["mother"][0] > 0.9 and xp["father"][0] > 0.9, "a female: one X from each parent"
    male = [v for v in kx if not tx.vcf["father"].get((v.chrom, v.pos, v.ref, v.alt))
            or tx.vcf["mother"].get((v.chrom, v.pos, v.ref, v.alt))]
    assert "no X from the father" in rep.x_parents_text(rep.x_parents(male, tx, cx))


def _write_tsv(rows, path):
    cols = list(rows[0].keys())
    with open(path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(r[c] for c in cols) + "\n")


def test_reanalysis_changes():
    rows, _ = run_report()
    tmp = tempfile.mkdtemp()
    same = os.path.join(tmp, "same.tsv")
    _write_tsv(rows, same)
    again, html = run_report(["--previous-tsv", same])
    assert not [r for r in again if r["category"].startswith("change_")]  # same inputs: nothing changed
    assert "0</b> new, <b>0</b> changed and <b>0</b> no longer" in html

    prev = [dict(r) for r in rows if r["category"] != "report_meta"]
    prev = [r for r in prev if r["gene"] != "PLN"]  # PLN: new diagnostic and secondary finding
    for r in prev:
        if r["gene"] == "BRCA2":  # BRCA2: previously a VUS in ClinVar, listed as a candidate
            r.update(category="candidate", clinvar="Uncertain_significance")
        if r["gene"] == "CYP2C9":
            r["genotype"] = "*1/*3"
    gone = dict(prev[0], gene="GENEZ", chrom="3", pos="3000100", ref="A", alt="G", category="diagnostic")
    called = dict(prev[0], gene="GENEX", chrom="1", pos="1000100", ref="C", alt="T", category="diagnostic",
                  clinvar="Likely_pathogenic")
    prev = [r for r in prev if r["gene"] != "GENEX"] + [gone, called]
    old = os.path.join(tmp, "old.tsv")
    _write_tsv(prev, old)
    cur, html = run_report(["--previous-tsv", old, "--previous-name", "report of 2025"])
    ch = {(r["category"], r["gene"]): r["note"] for r in cur if r["category"].startswith("change_")}
    assert ("change_new", "PLN") in ch and "tier 3" in ch[("change_new", "PLN")]
    assert "ClinVar Pathogenic (was Uncertain significance)" in ch[("change_changed", "BRCA2")]
    assert "genotype *1/*1 (was *1/*3)" in ch[("change_changed", "CYP2C9")]
    assert "no longer in the variant calls" in ch[("change_removed", "GENEZ")]
    assert "ClinVar not listed (was Likely pathogenic)" in ch[("change_changed", "GENEX")]
    assert "Changes since the previous report" in html and "report of 2025" in html
    assert "a version before 0.10.0" in html  # no report_meta rows in the old TSV
    assert "another sample" not in html
    assert len([k for k in ch if k[1] not in ("PLN", "BRCA2", "CYP2C9", "GENEZ", "GENEX")]) == 0


def test_reanalysis_other_sample():
    rows, _ = run_report()
    tmp = tempfile.mkdtemp()
    meta = [r for r in rows if r["category"] == "report_meta"]
    first = next(r for r in rows if r["category"] != "report_meta" and r["pos"])
    other = [dict(first, pos=str(5_000_000 + i), gene=f"OTHER{i}", gnomad_af="absent") for i in range(24)]
    path = os.path.join(tmp, "other.tsv")
    _write_tsv(meta + other, path)
    _, html = run_report(["--previous-tsv", path])
    assert "the previous report is probably of another sample" in html


def test_hgvs_g():
    v = rep.Variant()
    v.chrom, v.pos, v.ref, v.alt = "1", 1000, "C", "T"
    assert rep.hgvs_g(v) == "NC_000001.11:g.1000C>T"
    # reference ...G[A]AAAT...: POS 1000 is the G (offset 2 in seq)
    seq, off = "TTGAAAATCC", 2
    v.ref, v.alt = "GA", "G"  # one A deleted: 3'-most A is at 1004
    assert rep.hgvs_g(v, (seq, off)) == "NC_000001.11:g.1004del"
    assert rep.hgvs_g(v).endswith("g.1001del (not 3'-shifted)")
    v.ref, v.alt = "G", "GA"  # one more A: a duplication of the last A
    assert rep.hgvs_g(v, (seq, off)) == "NC_000001.11:g.1004dup"
    v.ref, v.alt = "G", "GC"  # C inserted between G and A: not a copy of a neighbour
    assert rep.hgvs_g(v, (seq, off)) == "NC_000001.11:g.1000_1001insC"
    seq2 = "TTGCACAGG"  # POS 1000 = G; CA CA repeat; inserting CA after G -> dup of the last CA (1003_1004)
    v.ref, v.alt = "G", "GCA"
    assert rep.hgvs_g(v, (seq2, 2)) == "NC_000001.11:g.1003_1004dup"
    v.chrom, v.ref, v.alt = "X", "GAT", "GCC"
    assert rep.hgvs_g(v) == "NC_000023.11:g.1001_1002delinsCC"


def test_mane_and_disease_context():
    extra = ["--mane", os.path.join(DATA, "mane.tsv"), "--gene-disease", os.path.join(DATA, "gene_disease.tsv")]
    rows, html = run_report(extra)
    brca2 = next(r for r in rows if r["gene"] == "BRCA2")
    assert brca2["hgvs_refseq"].startswith("NM_000059.4:c.1234C>T") and brca2["transcript"] == "TXB"
    assert "TXB (MANE Select)" in html
    pah = [r for r in rows if r["gene"] == "PAH"]
    assert pah and all("check the numbering" in r["hgvs_refseq"] for r in pah), \
        "MANE pairs TXP.2, SnpEff annotated TXP: the RefSeq name carries a warning"
    assert brca2["hgvs_g"] == "NC_000013.11:g.32340301C>T"
    assert "<h2>11. Disease context of the reported genes</h2>" in html
    assert "BRCA2-related Fanconi anemia" in html and "MONDO:0013697" in html
    assert "autosomal recessive: one allele found" in html, "BRCA2 het: carrier for Fanconi anemia"
    gener = next(r for r in rows if r["gene"] == "GENER")
    assert "PVS1" not in gener["acmg_codes"], "GENER: curated mechanism is gain of function only"
    rows0, _ = run_report()
    assert "PVS1" in next(r for r in rows0 if r["gene"] == "GENER")["acmg_codes"]
    assert "loss of function is not a curated mechanism" in html


def test_curated_moi_and_fit():
    rep.CURATED.clear()
    rep.CURATED.update(rep.load_gene_disease(os.path.join(DATA, "gene_disease.tsv"))[0])
    try:
        assert rep.curated_moi("BRCA2") == {"AD", "AR"}
        assert rep.curated_moi("GENER") == {"AD"}, "Orphanet's disease-group inheritance (AD,AR,XLR) is not used"
        assert rep.lof_not_mechanism("GENER") and not rep.lof_not_mechanism("PAH")
        v = _var("missense_variant")
        v.zyg, v.pos = "het", 100
        assert rep.genotype_fit(["AR"], [v], "female")[0] == "carrier"
        assert rep.genotype_fit(["AD"], [v], "female")[0] == "fits"
        v.chrom = "X"
        assert rep.genotype_fit(["XLR"], [v], "male")[0] == "fits"
        assert rep.genotype_fit(["XLR"], [v], "female")[0] == "carrier"
        w = _var("missense_variant")
        w.zyg, w.pos = "het", 200
        assert rep.genotype_fit(["AR"], [v, w], "female", "in cis (reads)")[0] == "carrier"
    finally:
        rep.CURATED.clear()


def test_vcep_rules():
    rep.VCEP.clear()
    rep.VCEP.update(rep.load_vcep_rules(os.path.join(DATA, "vcep_rules.tsv"))[0])
    try:
        adj = rep.vcep_adjust
        assert adj("BRCA2", ["PM2_supporting", "PP2", "PP3_Strong"], None)[0] == ["PM2_supporting", "PP3_Strong"], \
            "PP2 not used by the panel; the calibrated PP3 strength is kept (VCEP_KEEP_STRENGTH)"
        keep = rep.VCEP_KEEP_STRENGTH
        rep.VCEP_KEEP_STRENGTH = set()
        try:
            assert adj("BRCA2", ["PP3_Strong"], None)[0] == ["PP3_Supporting"], "with caps applied: Supporting only"
        finally:
            rep.VCEP_KEEP_STRENGTH = keep
        assert adj("BRCA2", ["PM2_supporting"], 0.0005)[0] == ["BS1"], "panel BS1 cut-off 0.0001"
        assert adj("BRCA2", ["BS1"], 0.002)[0] == ["BA1"], "panel BA1 cut-off 0.001 (default 5 %)"
        assert adj("GENER", ["PVS1", "PM2_supporting"], None)[0] == ["PM2_supporting"]
        assert adj("OTHER", ["PVS1", "PP2"], None)[0] == ["PVS1", "PP2"], "genes without a specification unchanged"
        assert rep.vcep_frequency("BRCA2") == {"BA1": 0.001, "BS1": 0.0001}
    finally:
        rep.VCEP.clear()


def test_vcep_and_erepo_in_report():
    extra = ["--vcep-rules", os.path.join(DATA, "vcep_rules.tsv"), "--erepo", os.path.join(DATA, "erepo_variants.tsv")]
    rows, html = run_report(extra)
    assert "ClinGen ENIGMA BRCA1 and BRCA2 VCEP: Likely Pathogenic" in html, "ERepo matched by genomic HGVS"
    assert "ClinGen VCEP rules: ENIGMA BRCA1 and BRCA2 VCEP v1.1" in html
    gener = next(r for r in rows if r["gene"] == "GENER")
    assert "PVS1" not in gener["acmg_codes"]
    meta = {r["gene"]: r["note"] for r in rows if r["category"] == "report_meta"}
    assert "vcep_rules" in meta and "erepo" in meta


def test_missense_policy():
    v = _var("missense_variant")
    v.am, v.revel = 0.995, 0.80
    assert rep.revel_code(0.80) == "PP3_Moderate" and rep.revel_code(0.01) == "BP4_Strong"
    assert rep.revel_code(0.29) == "BP4_Supporting" and rep.revel_code(0.5) is None
    assert rep.missense_code(v, "alphamissense") == "PP3_Strong"
    assert rep.missense_code(v, "revel") == "PP3_Moderate"
    assert rep.missense_code(v, "concordant") == "PP3_Moderate", "both pathogenic: the weaker strength"
    v.revel = 0.1
    assert rep.missense_code(v, "concordant") is None, "disagreeing predictors give no PP3/BP4"
    v.revel = None
    assert rep.missense_code(v, "concordant") == "PP3_Strong", "one score alone counts as itself"
    v.revel = 0.80
    assert rep.missense_code(v) == "PP3_Moderate", "default policy: REVEL first"


def test_cnv_collins_dosage():
    d = tempfile.mkdtemp()
    rows, _ = run_report(_cnv_inputs(d) + ["--collins-dosage", os.path.join(DATA, "collins_dosage.tsv")])
    pah = next(r for r in rows if r["category"].startswith("cnv") and r["gene"] == "PAH")
    assert "Collins et al. 2022 pHaplo): PAH 0.91" in pah["note"]
    brca2 = next(r for r in rows if r["category"].startswith("cnv") and r["gene"] == "BRCA2")
    assert "Collins" not in brca2["note"], "pHaplo 0.42 is below the 0.86 cut-off"


def test_mavedb_functional_evidence():
    rows, html = run_report(["--mave", os.path.join(DATA, "mave_evidence.tsv")])
    brca2 = next(r for r in rows if r["gene"] == "BRCA2" and r["acmg_codes"])
    assert "PS3_Moderate" in brca2["acmg_codes"] and "PS3_Strong" not in brca2["acmg_codes"], \
        "the clinical calibration is scored; the research-only one is shown, not scored"
    assert "MaveDB PS3 Strong (research only)" in html


def test_mitochondrial_section():
    rows, html = run_report(["--mito-tsv", os.path.join(DATA, "mito_calls.tsv"),
                             "--mito-json", os.path.join(DATA, "mito_summary.json")])
    assert "<h2>10. Mitochondrial DNA</h2>" in html
    assert "m.3243A&gt;G" in html and "MELAS" in html, "ClinVar pathogenic heteroplasmic variant listed"
    assert "carried, 100%" in html and "not assessable (2 reads)" in html and "not carried (33 reads)" in html
    cats = {(r["category"], r["gene"]) for r in rows}
    assert ("mito_pathogenic", "MT-TL1") in cats and ("pgx", "MT-RNR1") in cats
    assert "MT-RNR1 (aminoglycosides" not in html.split("<h2>8.")[0], "MT-RNR1 no longer listed as not assessed"
    _, html0 = run_report()
    assert "Not analysed: no mitochondrial calls given" in html0


def test_cnv_gnomad_frequency():
    d = tempfile.mkdtemp()
    rows, _ = run_report(_cnv_inputs(d) + ["--gnomad-cnv", os.path.join(DATA, "gnomad_cnv.vcf")])
    cnv = {r["gene"]: r for r in rows if r["category"].startswith("cnv")}
    assert "4O -1.00" in cnv["PAH"]["note"], "40 of 3,000 NFE exomes (1.3 %): common, benign evidence"
    assert "2 of 464,000 carriers" in cnv["BRCA2"]["note"] and "4O" not in cnv["BRCA2"]["note"]
    assert "not among the rare gnomAD v4.1 exome CNVs" in cnv["GENEX"]["note"]


def test_uniprot_features():
    _, html = run_report(["--uniprot-features", os.path.join(DATA, "uniprot_features.tsv")])
    assert "UniProt P51587: domain BRC repeat test (410-420); region Interaction with NPM1 (300-500)" in html, \
        "BRCA2 p.Gln412*: the narrower domain first"


def test_gnomad_coverage_and_pext():
    extra = ["--gnomad-coverage", os.path.join(DATA, "gnomad_coverage.tsv"), "--pext", os.path.join(DATA, "gnomad_pext.tsv")]
    rows, html = run_report(extra)
    brca2 = next(r for r in rows if r["gene"] == "BRCA2" and r["category"] == "diagnostic")
    assert "PM2" not in brca2["acmg_codes"], "absent from gnomAD where only 10 % of gnomAD exomes reach 20x"
    assert "absence uninformative, no PM2" in html
    gener = next(r for r in rows if r["gene"] == "GENER")
    assert "PVS1" in gener["acmg_codes"], "pext is shown, not scored (PEXT_LOW = None)"
    assert "pext 0.05 (highest tissue 0.05; gene max 0.90): little-expressed exon in GTEx" in html
    old = rep.PEXT_LOW
    rep.PEXT_LOW = 0.1
    try:
        rows, _ = run_report(extra)
        assert "PVS1" not in next(r for r in rows if r["gene"] == "GENER")["acmg_codes"]
    finally:
        rep.PEXT_LOW = old
