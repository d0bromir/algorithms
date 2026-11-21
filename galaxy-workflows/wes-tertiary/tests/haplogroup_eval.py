#!/usr/bin/env python3
"""Evaluate the mtDNA (wes_mito --haplogroup) and Y (wes_yhaplo) haplogroups against truth haplogroups.

The call is compared with the truth on the tree itself, not by name:
  same          the same node
  more derived  a descendant of the truth node (a newer tree, or more markers than the truth data had)
  less resolved an ancestor of the truth node (correct, but fewer markers covered)
  other branch  neither (an error)
"Consistent" = same + more derived + less resolved. Major lineage: mtDNA letters (L with its digit), Y first letter.

1000 Genomes panel (--panel): truth Y = Poznik et al. 2016 (most derived SNP; low-coverage genomes + exomes), mapped
onto the YFull tree of Yleaf through its marker file; truth mtDNA = the curated phase 3 mtDNA sequences
(low-coverage genomes) classified with the same HaploGrep3 tree (--mt-truth, HaploGrep3 output).
GIAB (--giab, --giab-truth): truth = the same tools on the deep (100-300x) GIAB genomes; plus maternal (mtDNA) and
paternal (Y) consistency in the trios.

    haplogroup_eval.py --mt-tree tree.xml --y-tree tree.json --y-positions new_positions.txt
        [--panel panel.tsv --mt-truth kg_mt_truth.txt --panel-res DIR] [--giab DIR --giab-truth DIR]
"""
import argparse
import collections
import csv
import json
import os
import re
import statistics
import xml.etree.ElementTree as ET

MATERNAL = [("HG001", "NA12892"), ("HG002", "HG004"), ("HG005", "HG007")]
PATERNAL = [("HG002", "HG003"), ("HG005", "HG006"), ("HG001", "NA12891")]  # HG001 is female: no Y


def mt_tree(path):
    parent = {}

    def walk(node, up):
        for h in node.findall("haplogroup"):
            parent[h.get("name")] = up
            walk(h, h.get("name"))
    walk(ET.parse(path).getroot(), None)
    # the rCRS trees hang from the reference (H2a2a1): re-root the same edges at mt-MRCA
    adj = collections.defaultdict(set)
    for k, v in parent.items():
        if v:
            adj[k].add(v)
            adj[v].add(k)
    if "mt-MRCA" not in adj:
        return parent
    rooted, todo = {"mt-MRCA": None}, ["mt-MRCA"]
    while todo:
        n = todo.pop()
        for m in adj[n]:
            if m not in rooted:
                rooted[m] = n
                todo.append(m)
    return rooted


def y_tree(path):
    parent = {}
    for node, kids in json.load(open(path)).items():
        for k in kids:
            parent[k] = node
    return parent


def ancestors(parent, n):
    out = []
    while n in parent and parent[n] is not None:
        n = parent[n]
        out.append(n)
    return out


def node(parent, name):
    """Tree node of a call: 'H1e1a+16278' -> H1e1a, 'I1*' -> I1."""
    if not name:
        return None
    if name in parent:
        return name
    for cand in (name.rstrip("*'"), re.split(r"[+*]", name)[0].rstrip("'")):
        if cand in parent:
            return cand
    return None


def relation(parent, call, truth):
    c, t = node(parent, call), node(parent, truth)
    if not c or not t:
        return "not in tree"
    if c == t:
        return "same"
    if t in ancestors(parent, c):
        return "more derived"
    if c in ancestors(parent, t):
        return "less resolved"
    return "other branch"


def mt_major(h):
    m = re.match(r"(L\d|[A-Z]+)", h or "")
    return m.group(1) if m else ""


def load(d, s):
    def j(k):
        p = os.path.join(d, f"{s}.{k}.json")
        return json.load(open(p)) if os.path.exists(p) else {}

    def t(k):
        p = os.path.join(d, f"{s}.{k}.time")
        if not os.path.exists(p):
            return None
        sec, kb = open(p).read().split()[-2:]
        return float(sec), int(kb)
    return j("mito").get("haplogroup", {}), j("mito"), j("y"), t("mito"), t("y")


def tally(name, rels):
    c = collections.Counter(rels)
    n = len(rels)
    ok = c["same"] + c["more derived"] + c["less resolved"]
    pct = lambda k: f"{k}/{n} ({100 * k / n:.1f} %)" if n else "0/0"
    print(f"{name}: consistent {pct(ok)}; same {c['same']}, more derived {c['more derived']}, "
          f"less resolved {c['less resolved']}, other branch {c['other branch']}, not in tree {c['not in tree']}")
    return ok, n


def times(label, ts):
    ts = [t for t in ts if t]
    if ts:
        print(f"{label}: median {statistics.median(t[0] for t in ts):.1f} s (max {max(t[0] for t in ts):.1f} s), "
              f"peak memory median {statistics.median(t[1] for t in ts) / 1024:.0f} MB "
              f"(max {max(t[1] for t in ts) / 1024:.0f} MB), n={len(ts)}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--mt-tree", required=True, help="HaploGrep3 tree.xml (phylotree-fu-rcrs 1.2)")
    ap.add_argument("--y-tree", required=True, help="Yleaf hg_prediction_tables/tree.json (YFull)")
    ap.add_argument("--y-positions", required=True, help="Yleaf hg38/new_positions.txt (marker -> YFull node)")
    ap.add_argument("--panel", help="panel.tsv: sample, population, truth Y (hg-SNP), major Y lineage, CRAM URL")
    ap.add_argument("--mt-truth", help="HaploGrep3 output for the 1000 Genomes mtDNA sequences")
    ap.add_argument("--panel-res", help="directory of {sample}.mito.json / .y.json / .time for the panel")
    ap.add_argument("--giab", help="directory of the GIAB exome results")
    ap.add_argument("--giab-truth", help="directory of the GIAB genome (truth) results")
    ap.add_argument("--table", help="write the per-sample table (TSV)")
    a = ap.parse_args(argv)
    mtp, yp = mt_tree(a.mt_tree), y_tree(a.y_tree)
    snp = collections.defaultdict(set)  # marker (upper case) -> YFull nodes (names repeat across lineages)
    for line in open(a.y_positions):
        t = line.rstrip("\n").split("\t")
        if len(t) > 2:
            snp[t[1].upper()].add(t[2])
    ynodes = {n.upper(): n for n in yp} | {n.upper(): n for n in yp.values() if n}

    def y_truth_node(hg_snp):
        """Poznik 'E1a-L133' -> YFull node: 'E-L133' if it exists, else a node of lineage E carrying marker L133."""
        letter, markers = hg_snp[0], hg_snp.split("-")[-1].upper().split("/")
        for m in markers:
            if f"{letter}-{m}" in ynodes:
                return ynodes[f"{letter}-{m}"]
        for m in markers:
            hits = sorted(n for n in snp.get(m, ()) if n.startswith(letter))
            if hits:
                return hits[0]
        return None
    table = []

    if a.panel:
        mt_truth = {r["SampleID"]: r["Haplogroup"] for r in csv.DictReader(open(a.mt_truth), delimiter="\t")}
        mt_rel, y_rel, mt_maj, y_maj, tm, ty, conf = [], [], [], [], [], [], collections.Counter()
        for i, r in enumerate(csv.reader(open(a.panel), delimiter="\t")):
            s, pop, ytruth, ymajor = r[:4]
            half = "tuning" if i % 2 == 0 else "held-out"  # odd panel lines set MIN_COVERED in wes_mito
            hg, mito, y, t1, t2 = load(a.panel_res, s)
            if not mito:
                continue
            tm.append(t1)
            ty.append(t2)
            mrel = relation(mtp, hg.get("haplogroup"), mt_truth.get(s)) if hg.get("haplogroup") else "no call"
            mt_rel.append(mrel)
            mt_maj.append(mt_major(hg.get("haplogroup")) == mt_major(mt_truth.get(s)))
            conf[(half, "mtDNA " + str(hg.get("confidence")), mrel)] += 1
            tnode = y_truth_node(ytruth)
            ycall = y.get("haplogroup", "")
            yrel = relation(yp, ycall, tnode) if ycall else "no call"
            y_rel.append(yrel)
            letter = ytruth[0]  # major haplogroup (A-T): the call is that node or below it
            cn = node(yp, ycall)
            y_maj.append(bool(cn) and (cn == letter or letter in ancestors(yp, cn) if letter in yp
                                       else cn.startswith(letter)))  # A has no single node (A00, A0, A1...)
            conf[(half, "Y " + str(y.get("confidence")), yrel)] += 1
            table.append([s, pop, mt_truth.get(s, ""), hg.get("haplogroup", ""), str(hg.get("quality", "")),
                          str(hg.get("covered_fraction", "")), mrel, ytruth, tnode or "", ycall,
                          str(y.get("valid_markers", "")), str(y.get("qc_score", "")), yrel])
        print(f"== 1000 Genomes phase 3 exomes ({len(mt_rel)} males, {len(set(r[1] for r in table))} populations)")
        tally("mtDNA (calls)", [r for r in mt_rel if r != "no call"])
        print(f"mtDNA no call: {mt_rel.count('no call')}; same major lineage: {sum(mt_maj)}/{len(mt_maj)}")
        tally("Y (calls)", [r for r in y_rel if r != "no call"])
        print(f"Y no call: {y_rel.count('no call')}; resolved to the major haplogroup (A-T) or below: "
              f"{sum(y_maj)}/{len(y_maj)}")
        for half in ("tuning", "held-out"):
            print(f"{half} half, by confidence:")
            for kind in sorted({k[1] for k in conf if k[0] == half}):
                c = {k[2]: v for k, v in conf.items() if k[:2] == (half, kind)}
                print(f"  {kind}: {sum(c.values())} ({', '.join(f'{r} {n}' for r, n in sorted(c.items()))})")
        times("wes_mito --haplogroup", tm)
        times("wes_yhaplo", ty)

    if a.giab:
        print("== GIAB IDT exomes vs the same tools on deep genomes")
        samples = ["HG001", "HG002", "HG003", "HG004", "HG005", "HG006", "HG007", "NA12891", "NA12892"]
        res, mt_rel, y_rel = {}, [], []
        for s in samples:
            hg, mito, y, t1, t2 = load(a.giab, s)
            thg, _, ty_, _, _ = load(a.giab_truth, s) if a.giab_truth else ({}, {}, {}, None, None)
            res[s] = (hg.get("haplogroup", ""), y.get("haplogroup", ""))
            mrel = relation(mtp, hg.get("haplogroup"), thg.get("haplogroup")) if thg.get("haplogroup") else "no truth"
            yrel = (relation(yp, y.get("haplogroup"), ty_.get("haplogroup")) if ty_.get("haplogroup") and
                    y.get("haplogroup") else "no truth" if not ty_.get("haplogroup") else "no call")
            mt_rel.append(mrel)
            y_rel.append(yrel)
            print(f"{s}\tmt {hg.get('haplogroup', '')} (q {hg.get('quality', '')}, {hg.get('confidence', '')}; "
                  f"genome {thg.get('haplogroup', '-')}: {mrel})\tY {y.get('haplogroup', '') or y.get('status', '')}"
                  f" ({y.get('valid_markers', '')} markers; genome {ty_.get('haplogroup', '-')}: {yrel})\t"
                  f"{t1[0] if t1 else ''} s / {t2[0] if t2 else ''} s, {max(t1[1], t2[1]) // 1024 if t1 and t2 else ''} MB")
            table.append([s, "GIAB", thg.get("haplogroup", ""), hg.get("haplogroup", ""), str(hg.get("quality", "")),
                          str(hg.get("covered_fraction", "")), mrel, ty_.get("haplogroup", ""), "",
                          y.get("haplogroup", ""), str(y.get("valid_markers", "")), str(y.get("qc_score", "")), yrel])
        tally("mtDNA vs genome", [r for r in mt_rel if r != "no truth"])
        tally("Y vs genome", [r for r in y_rel if r not in ("no truth", "no call")])
        for kind, pairs, i in (("maternal mtDNA", MATERNAL, 0), ("paternal Y", PATERNAL, 1)):
            for c, p in pairs:
                if res.get(c, ("", ""))[i] or res.get(p, ("", ""))[i]:
                    same = res[c][i] == res[p][i]
                    print(f"{kind}: {c} {res[c][i] or '-'} / {p} {res[p][i] or '-'}: "
                          f"{'identical' if same else relation(mtp if i == 0 else yp, res[c][i], res[p][i])}")

    if a.table:
        with open(a.table, "w") as f:
            f.write("\t".join(["sample", "population", "mt_truth", "mt_call", "mt_quality", "mt_covered", "mt_result",
                               "y_truth", "y_truth_node", "y_call", "y_markers", "y_qc", "y_result"]) + "\n")
            for r in table:
                f.write("\t".join(r) + "\n")


if __name__ == "__main__":
    main()
