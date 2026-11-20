#!/usr/bin/env python3
"""Checks the report's reanalysis section against the inputs, independently of the tool's code.

    reanalysis_eval.py OUT_DIR SAMPLE [SAMPLE ...]

For each SAMPLE, OUT_DIR holds (made by rerunning the report, see the README):
  SAMPLE.cur.tsv    report on the current annotated VCF
  SAMPLE.old.tsv    report on the same VCF re-annotated with an older ClinVar release
  SAMPLE.diff.tsv   report on the current VCF with --previous-tsv SAMPLE.old.tsv
  SAMPLE.same.tsv   report on the current VCF with --previous-tsv SAMPLE.cur.tsv
  SAMPLE.oldcv.vcf  the re-annotated VCF; the current VCF path is read from SAMPLE.cur.vcfpath

Expected changes: every listed variant (in either report) whose ClinVar CLNSIG differs between the two VCFs, and
every finding whose category set differs between the two reports. Found changes: the change_* rows. Reasons:
each "ClinVar X (was Y)" must quote the two VCFs. Same-input diffs must be empty.
"""
import collections
import csv
import os
import re
import sys

LOW = {"cnv_low_confidence"}


def rows(path):
    with open(path) as f:
        return [r for r in csv.DictReader(f, delimiter="\t") if r["category"] != "report_meta"]


def key(r):
    c = r["category"]
    if c == "upd" or c.startswith("cnv_"):
        return ("region", r["chrom"], r["pos"], r["genotype"])
    if r["pos"]:
        return ("var", r["chrom"], r["pos"], r["ref"], r["alt"])
    return ("gene", c, r["gene"])


def findings(path):
    out = collections.defaultdict(set)
    for r in rows(path):
        if r["category"] not in LOW and not r["category"].startswith("change_"):
            out[key(r)].add(r["category"])
    return out


def clnsig(vcf, keys):
    want = {(k[1], k[2], k[3], k[4]) for k in keys if k[0] == "var"}
    got = {}
    with open(vcf) as f:
        for line in f:
            if line[0] == "#":
                continue
            t = line.split("\t", 8)
            k = (t[0], t[1], t[3], t[4])
            if k in want:
                m = re.search(r"(?:^|;)CLNSIG=([^;\t]*)", t[7])
                got[k] = m.group(1) if m else ""
    return got


def main(out_dir, samples):
    tot = collections.Counter()
    print("sample\texpected\tfound\tmissed\textra\tclinvar_reasons\treasons_wrong\tsame_input_changes\t"
          "tier3_new\ttier3_changed\ttier3_removed\ttier2\ttier1")
    for s in samples:
        p = os.path.join(out_dir, s)
        cur, old = findings(p + ".cur.tsv"), findings(p + ".old.tsv")
        vcf_cur = open(p + ".cur.vcfpath").read().strip()
        keys = set(cur) | set(old)
        cv_new, cv_old = clnsig(vcf_cur, keys), clnsig(p + ".oldcv.vcf", keys)
        expected = {k for k in keys if k[0] == "var" and cv_new.get(k[1:], "") != cv_old.get(k[1:], "")}
        expected |= {k for k in keys if cur.get(k, set()) != old.get(k, set())}
        diff = [r for r in rows(p + ".diff.tsv") if r["category"].startswith("change_")]
        found = {key(r) for r in diff}
        wrong = n_cv = 0
        for r in diff:
            m = re.findall(r"ClinVar ([^;()]+?) \(was ([^;()]+?)\)", r["note"])
            if m and key(r)[0] == "var":
                n_cv += 1
                k = key(r)[1:]
                now, was = (x.replace("not listed", "").replace(" ", "_") for x in m[-1])
                if now != cv_new.get(k, "") or was != cv_old.get(k, ""):
                    wrong += 1
                    print("  wrong reason:", s, k, r["note"][:200], cv_new.get(k), cv_old.get(k), file=sys.stderr)
        same = [r for r in rows(p + ".same.tsv") if r["category"].startswith("change_")]
        tiers = collections.Counter((int(re.match(r"tier (\d)", r["note"]).group(1)), r["category"][7:]) for r in diff)
        missed, extra = expected - found, found - expected
        for k in sorted(missed)[:5]:
            print("  missed:", s, k, cur.get(k), old.get(k), file=sys.stderr)
        for k in sorted(extra)[:5]:
            print("  extra:", s, k, cur.get(k), old.get(k), file=sys.stderr)
        vals = [len(expected), len(found & expected), len(missed), len(extra), n_cv, wrong, len(same),
                tiers[(3, "new")], tiers[(3, "changed")], tiers[(3, "removed")],
                sum(v for (t, _), v in tiers.items() if t == 2), sum(v for (t, _), v in tiers.items() if t == 1)]
        print(s, *vals, sep="\t")
        tot.update(dict(zip(range(len(vals)), vals)))
    print("total", *(tot[i] for i in range(12)), sep="\t")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
