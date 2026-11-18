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
    note = next(r["note"] for r in rows if r["gene"] == "GENEX")
    assert "PVS1" in note and "PM2_supporting" in note


def test_pharmacogenomics_vcf_only():
    rows, html = run_report()
    pgx = {r["gene"]: (r["genotype"], r["note"]) for r in rows if r["category"] == "pgx"}
    assert pgx["CYP2C19"] == ("*2/*2", "Poor metabolizer")
    assert pgx["CYP2C9"][1].startswith("Normal metabolizer")
    assert all(len(r) == 25 and None not in r.values() for r in rows)
    assert "absence of a variant is assumed to be the reference allele" in html


def test_pharmcat_calls():
    rows, html = run_report(["--pharmcat", os.path.join(DATA, "pharmcat.tsv")])
    pc = {r["gene"]: (r["genotype"], r["note"]) for r in rows if r["category"] == "pgx_pharmcat"}
    assert pc["CYP2C19"] == ("*1/*2", "Intermediate Metabolizer"), "PharmCAT 3.4.0 on HG001 (GeT-RM: *1/*2)"
    assert pc["CYP3A5"] == ("*3/*3", "Poor Metabolizer")
    assert "<h3>PharmCAT 3.4.0</h3>" in html and "Missing positions" in html
    assert all(len(r) == 25 for r in rows)
    _, html = run_report()
    assert "PharmCAT" not in html.split("<h2>6. Pharmacogenomics")[1].split("<h2>7.")[0]


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
