#!/usr/bin/env python3
"""Benchmark of the report's panel suggestions (wes_clinical_report suggest_panels) on simulated patients.

Patients: the diseases of phenotype.hpoa (OMIM, ORPHA) that phenotype_to_genes.txt links to exactly one gene; each
patient gets --terms of the disease's own phenotype terms (aspect P, drawn at random) and --noise terms of another
disease (imprecise or unrelated findings). A suggestion is correct when the panel has the disease gene rated green.
Diseases whose gene is green on no panel are counted apart (no panel could be right). As in the report, panels with
more than PANEL_SUGGEST_MAX_GENES green genes (superpanels) are not suggested.

The ranking method (PANEL_SUGGEST) is chosen on the tuning half (diseases at even positions of the sorted list) by
the share of patients with a correct panel among the first three, and every method is reported on the held-out half.
The terms' information content and the gene matching are those of the report (propagated phenotype_to_genes.txt).
The scores are computed from an inverted index (gene -> panels) for speed and checked against the report's
suggest_panels on the first 30 patients.

    panel_eval.py --report tools/wes_clinical_report/wes_clinical_report.py --panels panelapp_panels.tsv
        --hpo-genes phenotype_to_genes.txt --hpoa phenotype.hpoa [--terms 5] [--noise 1] [--max 1500] [--seed 1]
"""
import argparse
import collections
import importlib.util
import math
import random
import statistics
import time

METHODS = ("max", "sum", "sum_sqrt", "sum_log", "share")


def load_report(path):
    spec = importlib.util.spec_from_file_location("wes_clinical_report", path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--report", required=True)
    ap.add_argument("--panels", required=True)
    ap.add_argument("--hpo-genes", required=True)
    ap.add_argument("--hpoa", required=True)
    ap.add_argument("--terms", type=int, default=5)
    ap.add_argument("--noise", type=int, default=1)
    ap.add_argument("--max", type=int, default=1500, help="patients per half at most")
    ap.add_argument("--seed", type=int, default=2)
    a = ap.parse_args(argv)
    rep = load_report(a.report)
    t0 = time.time()
    panels, _src = rep.load_panels(a.panels)
    green_of = collections.defaultdict(set)  # gene -> panels with the gene green
    size = {}
    for k, p in panels.items():
        g = [x for x, (r, _m) in p["genes"].items() if r == "green"]
        size[k] = len(g)
        if len(g) > rep.PANEL_SUGGEST_MAX_GENES:  # superpanels are not suggested (as in the report)
            continue
        for x in g:
            green_of[x].add(k)

    term_genes, term_dis, dis_genes, all_dis = collections.defaultdict(set), collections.defaultdict(set), \
        collections.defaultdict(set), set()
    with open(a.hpo_genes) as f:
        f.readline()
        for line in f:
            t = line.rstrip("\n").split("\t")
            if len(t) < 5:
                continue
            term_genes[t[0]].add(t[3])
            if t[4]:
                term_dis[t[0]].add(t[4])
                dis_genes[t[4]].add(t[3])
                all_dis.add(t[4])
    n_dis = len(all_dis)
    ic = {t: math.log(n_dis / len(d)) for t, d in term_dis.items() if d}
    direct = collections.defaultdict(set)
    with open(a.hpoa) as f:
        for line in f:
            if line.startswith("#") or line.startswith("database_id"):
                continue
            t = line.rstrip("\n").split("\t")
            if len(t) > 10 and t[10] == "P" and t[2] != "NOT":
                direct[t[0]].add(t[3])
    diseases = sorted(d for d in direct if len(dis_genes.get(d, ())) == 1 and len(direct[d]) >= a.terms)
    print(f"loaded in {time.time() - t0:.0f} s: {len(panels)} panels, {len(diseases)} single-gene diseases with "
          f">= {a.terms} phenotype terms")
    rng = random.Random(a.seed)
    halves = {"tuning": diseases[0::2], "held-out": diseases[1::2]}

    def patient(d):
        terms = set(rng.sample(sorted(direct[d]), a.terms))
        for _ in range(a.noise):
            other = rng.choice(diseases)
            if other != d:
                terms.add(rng.choice(sorted(direct[other])))
        return terms

    def orders(terms):
        """method -> panels in suggestion order (same scores and tie-break as the report)."""
        tot = sum(ic.get(t, 0.0) for t in terms)
        share = collections.Counter()
        for t in terms:
            for g in term_genes.get(t, ()):
                share[g] += ic.get(t, 0.0)
        agg = {}
        for g, s in share.items():
            s /= tot
            for k in green_of.get(g, ()):
                x = agg.setdefault(k, [0.0, 0.0, 0])
                x[0] += s
                x[1] = max(x[1], s)
                x[2] += s >= 0.5
        out = {}
        for m in METHODS:
            sc = {k: {"max": x[1], "sum": x[0], "sum_sqrt": x[0] / math.sqrt(size[k]),
                      "sum_log": x[0] / math.log(size[k] + math.e), "share": x[2] / size[k]}[m]
                  for k, x in agg.items() if x[1] > 0}
            out[m] = sorted(sc, key=lambda k: (-sc[k], len(panels[k]["genes"]), k))
        return out

    def ranks(terms, gene):
        """method -> (rank (1-based) of the first panel with the gene green or None, green genes of the top panel)."""
        return {m: (next((i + 1 for i, k in enumerate(o) if k in green_of.get(gene, ())), None),
                    size[o[0]] if o else 0) for m, o in orders(terms).items()}

    # the fast ranking must equal the report's suggest_panels
    for d in halves["tuning"][:30]:
        terms = patient(d)
        ph = rep.Phenotype()
        ph.terms = terms
        ph.ic = {t: ic.get(t, 0.0) for t in terms}
        for t in terms:
            for g in term_genes.get(t, ()):
                ph.gene_terms[g].add(t)
        fast = orders(terms)
        for m in METHODS:
            want = [k for _s, k, _g in rep.suggest_panels(panels, ph, top=10, method=m)]
            assert fast[m][:10] == want, (d, m, fast[m][:10], want)
    print("fast ranking = report's suggest_panels on 30 patients, all methods")
    rng.seed(a.seed)

    res = {}
    for name, ds in halves.items():
        ds = ds[:a.max]
        per = {m: [] for m in METHODS}
        no_panel = 0
        t1 = time.time()
        for d in ds:
            gene = next(iter(dis_genes[d]))
            terms = patient(d)
            if not green_of.get(gene):
                no_panel += 1
                continue
            r = ranks(terms, gene)
            for m in METHODS:
                per[m].append(r[m])
        res[name] = (per, no_panel, len(ds), (time.time() - t1) / max(len(ds), 1))
    tune = res["tuning"][0]
    best = max(METHODS, key=lambda m: sum(1 for r, _n in tune[m] if r and r <= 3))
    print(f"chosen on the tuning half (correct panel in the first three): {best} (report default {rep.PANEL_SUGGEST})")
    for name, (per, no_panel, n, sec) in res.items():
        print(f"== {name}: {n} patients, {no_panel} with the gene green on no panel ({100 * no_panel / n:.1f} %); "
              f"{1000 * sec:.1f} ms per patient")
        for m in METHODS:
            rr = [x[0] for x in per[m]]
            sizes = [x[1] for x in per[m]]
            k = len(rr)
            hit = lambda c: sum(1 for r in rr if r and r <= c)
            med = statistics.median(r if r else 10 ** 6 for r in rr) if rr else 0
            print(f"  {m:9s} top-1 {hit(1)}/{k} ({100 * hit(1) / k:.1f} %), top-3 {100 * hit(3) / k:.1f} %, "
                  f"top-5 {100 * hit(5) / k:.1f} %, median rank {med:g}, green genes of the first panel: median "
                  f"{statistics.median(sizes):g}" + ("   <- chosen" if m == best else ""))


if __name__ == "__main__":
    main()
