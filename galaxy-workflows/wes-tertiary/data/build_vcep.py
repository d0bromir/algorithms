#!/usr/bin/env python3
"""ClinGen expert-panel (VCEP) rules and classifications for the clinical report.

Inputs (both free, ClinGen terms of use):
  CSPEC_DIR  the released criteria specifications of the ClinGen CSpec registry, one JSON per specification:
             https://cspec.genome.network/cspec/api/svis lists them; each is fetched from its @id with ?detail=high
             (--fetch downloads them into CSPEC_DIR first)
  EREPO.tsv  the ClinGen Evidence Repository's expert-panel classifications:
             https://erepo.clinicalgenome.org/evrepo/api/classifications/all?format=tabbed

Outputs:
  vcep_rules.tsv     gene  spec  version  panel  code  applicable  value  text
      applicable  the strengths the panel allows for the code (Supporting, Moderate, Strong, Very Strong,
                  Stand Alone; ';'-separated), '' when the code is not applicable to the gene at any strength
      value       for BA1, BS1 and PM2: the allele-frequency cut-off parsed from the panel's text (fraction), else ''
      text        the panel's text for the code (first 400 characters; review it before relying on a value)
  erepo_variants.tsv hgvs_g  gene  disease  mondo  moi  assertion  met  not_met  panel  approved  link
      hgvs_g      the GRCh38 genomic HGVS (NC_0000xx.yy:g....) the repository gives; the report matches it with its own
                  3'-shifted genomic names

Only 'Released' specifications are used; a gene in several (e.g. RYR1, ACTA1) keeps all of them, and the report
combines them conservatively. Standard library only.

Usage: build_vcep.py [--fetch] CSPEC_DIR EREPO.tsv OUT_RULES.tsv OUT_EREPO.tsv
"""
import csv
import glob
import hashlib
import json
import os
import re
import sys
import urllib.request

STRENGTHS = ("Supporting", "Moderate", "Strong", "Very Strong", "Stand Alone")
FREQ_CODES = ("BA1", "BS1", "PM2")


def fetch(cspec_dir):
    os.makedirs(cspec_dir, exist_ok=True)
    with urllib.request.urlopen("https://cspec.genome.network/cspec/api/svis?detail=low") as r:
        for x in json.load(r)["data"]:
            sid = x["@id"].rsplit("/", 1)[1]
            p = os.path.join(cspec_dir, f"{sid}.json")
            if not os.path.exists(p):
                urllib.request.urlretrieve(x["@id"] + "?detail=high", p)


CMP = re.compile(r"(≥|>=|≤|<=|>|<|greater than(?: or equal to)?|less than|lower than|above|below|exceeds?|"
                 r"cut-?off (?:of )?|threshold (?:of )?|frequency (?:of )?)\s*(\d*\.?\d+)(?:\s*[x×]\s*10\s*<sup>\s*([-−]\d+)"
                 r"\s*</sup>|\s*[x×]\s*10\^?([-−]\d+)|e([-−]\d+))?\s*(%?)", re.I)


def freq_value(code, texts):
    """Allele-frequency cut-off from a panel's BA1/BS1/PM2 text: the first number after a comparison ('≥0.001',
    'greater than 0.0163', '1.11 x 10<sup>-5</sup>', '0.017%'), as a fraction. None when the text gives no such
    number (the report then keeps its default rule). Bare numbers are not used: they are prevalences, penetrances,
    confidence levels and allele counts."""
    for t in texts:
        t = t.replace("\\<", "<").replace("\\>", ">").replace(",", "").replace("\\_", "_")
        for m in CMP.finditer(t):
            x = float(m.group(2))
            exp = m.group(3) or m.group(4) or m.group(5)
            if exp:
                x *= 10 ** int(exp.replace("−", "-"))
            if m.group(6) == "%":
                x /= 100
            if 0 < x < 1:
                return x
    return None


def check_order(rows):
    """Drop parsed cut-offs that contradict each other within a specification (BA1 > BS1 > PM2 is required)."""
    by = {}
    for r in rows:
        if r[4] in FREQ_CODES and r[6]:
            by.setdefault((r[0], r[1]), {})[r[4]] = float(r[6])
    bad = set()
    for k, d in by.items():
        if "BA1" in d and "BS1" in d and d["BS1"] >= d["BA1"]:
            bad |= {(k, "BA1"), (k, "BS1")}
        if "PM2" in d and min(d.get("BS1", 1), d.get("BA1", 1)) <= d["PM2"]:
            bad.add((k, "PM2"))
    for r in rows:
        if ((r[0], r[1]), r[4]) in bad:
            r[6] = ""
    return len(bad)


def strength_label(s):
    s = s.strip().lower()
    return next((x for x in STRENGTHS if x.lower() == s), s.title())


def rules(cspec_dir):
    out = []
    for f in sorted(glob.glob(os.path.join(cspec_dir, "*.json"))):
        d = json.load(open(f, encoding="utf-8"))
        if d.get("currentStatus") != "Released":
            continue
        sid = d["@id"].rsplit("/", 1)[1]
        panel = d.get("affiliation", {}).get("label", "")
        ver = str(d.get("version", ""))
        for rs in d.get("ruleSets", []):
            genes = [g["label"] for g in rs.get("genes", [])]
            for c in rs.get("criteriaCodes", []):
                code = c["label"]
                ok = [s for s in c.get("evidenceStrengths", [])
                      if str(s.get("applicability", "")).strip().lower().startswith("applicable")]
                texts = [s.get("description") or "" for s in ok] + [c.get("description") or ""]
                val = freq_value(code, texts[:-1] or texts) if code in FREQ_CODES and ok else None
                text = " | ".join(f"{strength_label(s['label'])}: {s.get('description') or ''}" for s in ok)
                text = re.sub(r"\s+", " ", text)[:400]
                for g in genes:
                    out.append([g, sid, ver, panel, code, ";".join(strength_label(s["label"]) for s in ok),
                                f"{val:g}" if val is not None else "", text])
    return out


def erepo(path):
    out = []
    with open(path, encoding="utf-8", newline="") as f:
        r = csv.reader(f, delimiter="\t")
        head = next(r)
        ix = {h.lstrip("#"): i for i, h in enumerate(head)}
        for t in r:
            if len(t) < len(head) or t[ix["Retracted"]].strip().lower() == "true":
                continue
            hg = [h.strip() for h in t[ix["HGVS Expressions"]].split(",")]
            g38 = [h for h in hg if ":g." in h and h.split(":")[0] in GRCH38]
            for g in g38[:1]:
                out.append([g, t[ix["HGNC Gene Symbol"]], t[ix["Disease"]], t[ix["Mondo Id"]],
                            t[ix["Mode of Inheritance"]], t[ix["Assertion"]],
                            t[ix["Applied Evidence Codes (Met)"]].replace(" ", ""),
                            t[ix["Applied Evidence Codes (Not Met)"]].replace(" ", ""), t[ix["Expert Panel"]],
                            t[ix["Approval Date"]][:10], t[ix["Evidence Repo Link"]]])
    return out


# GRCh38 RefSeq chromosome accessions (the repository lists GRCh37 and NCBI36 names too)
GRCH38 = {f"NC_0000{i:02d}.{v}" for i, v in zip(range(1, 23), (11, 12, 12, 12, 10, 12, 14, 11, 12, 11, 10, 12, 11, 9,
                                                                 10, 10, 11, 10, 10, 11, 9, 11))} | {
    "NC_000023.11", "NC_000024.10", "NC_012920.1"}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main(argv):
    if argv and argv[0] == "--fetch":
        fetch(argv[1])
        argv = argv[1:]
    if len(argv) != 4:
        raise SystemExit(__doc__)
    cdir, ep, out_rules, out_erepo = argv
    rr = rules(cdir)
    n_bad = check_order(rr)
    digest = hashlib.sha256()
    for f in sorted(glob.glob(os.path.join(cdir, "*.json"))):
        digest.update(open(f, "rb").read())
    with open(out_rules, "w", encoding="utf-8", newline="\n") as f:
        f.write(f"# source\tCSpec registry ({len(glob.glob(os.path.join(cdir, '*.json')))} specifications)\t"
                f"sha256:{digest.hexdigest()}\n")
        f.write("gene\tspec\tversion\tpanel\tcode\tapplicable\tvalue\ttext\n")
        for r in rr:
            f.write("\t".join(x.replace("\t", " ") for x in r) + "\n")
    er = erepo(ep)
    with open(out_erepo, "w", encoding="utf-8", newline="\n") as f:
        f.write(f"# source\t{os.path.basename(ep)}\tsha256:{sha256(ep)}\n")
        f.write("hgvs_g\tgene\tdisease\tmondo\tmoi\tassertion\tmet\tnot_met\tpanel\tapproved\tlink\n")
        for r in er:
            f.write("\t".join(x.replace("\t", " ") for x in r) + "\n")
    print(f"{len(rr)} rule rows for {len({r[0] for r in rr})} genes ({n_bad} contradictory cut-offs dropped); "
          f"{len(er)} expert classifications", file=sys.stderr)


if __name__ == "__main__":
    main(sys.argv[1:])
