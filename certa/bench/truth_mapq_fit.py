#!/usr/bin/env python3
"""Fits the MAPQ table of `certa map --intervals --mapq cal`
(bench/VALIDATION_RESULTS.md, "A MAPQ fitted at the locus") and writes it as
include/certa/mapq_table.h.

  samtools view RUN.bam | truth_mapq_fit.py feat truth.npz REF.fai OUT.npz
  truth_mapq_fit.py fit OUT.npz ... > include/certa/mapq_table.h

RUN is paired interval mode with --iv-mapq cert, --mapq-limit 10 and
--feature-tags (XO: the best other candidate's distance, XN: clusters
verified). `feat` keeps per mate: placed, at the locus (truth_eval.py's rule
at 150 bp), MAPQ, tier, XE, XD, XO, XN. `fit` pools the runs given and
writes, as Phred of the misplaced share (+0.5 / +1), one value per cell:
- a record with MAPQ >= 40 and tier S0, S1 or PR: tier x edits (XE, 0..5) x
  gap (second-best distance minus XE, 0..10, or none within the limit) x
  clusters verified (1, 2-4, 5-16, 17 or more);
- any other record: tier x MAPQ.
A cell with fewer than MIN_N mates takes its parent's value: clusters
pooled, then edits, then gaps, then the tier.
"""
import re
import sys

import numpy as np

sys.path.insert(0, __import__("os").path.dirname(__file__))
from truth_eval import read_sam  # noqa: E402
from truth_lift import CLS_TRUTH, Ref  # noqa: E402

TIERS = ("S0", "S1", "SR", "S2", "SL", "PR", "R", "IV", "MR")
CONF = (0, 1, 5)  # S0, S1, PR
NE, NG, NX = 6, 12, 4  # edits 0..5; gaps 0..10 and none (11); cluster buckets
MIN_N = 1000
TG = re.compile(r"\t(XE|XD|XO|XN):i:(-?\d+)")


def bucket_x(xn):
    return np.digitize(xn, [2, 5, 17])


def cmd_feat(a):
    t = np.load(a[0])
    n = t["cls"].shape[0]
    ref = Ref(a[1])
    f = {k: np.full((n, 2), -1, np.int32) for k in ("XE", "XD", "XO", "XN")}
    tier = np.full((n, 2), -1, np.int8)

    def tee():
        for line in sys.stdin:
            if line[0] != "@":
                c = line.split("\t", 11)
                flag = int(c[1])
                if not flag & (256 | 2048) and len(c) > 11:
                    i, m = int(c[0][1:]), 1 if flag & 128 else 0
                    rest = "\t" + c[11]
                    for k, v in TG.findall(rest):
                        f[k][i, m] = int(v)
                    x = re.search(r"\tXT:Z:(\w+)", rest)
                    if x and x.group(1) in TIERS:
                        tier[i, m] = TIERS.index(x.group(1))
            yield line
    r = read_sam(tee(), n, ref.index)
    ref.par_to_x(r["rc"], r["left"])
    truth = t["cls"] == CLS_TRUTH
    ok = r["placed"] & truth & (r["rc"] == t["tc"]) & (r["strand"] == t["ts"]) & (np.abs(r["left"] + 1 - t["tp"]) <= 150)
    np.savez(a[2], placed=r["placed"] & truth, ok=ok, mapq=r["mapq"], tier=tier, **f)


def cells(z):
    """Per placed mate: (confident cell index or -1, low cell index or -1)."""
    mq = z["mapq"].astype(np.int64)
    tier = z["tier"].astype(np.int64)
    second = np.where(z["XD"] >= 0, z["XD"], z["XO"])
    gap = np.where(second >= 0, np.minimum(np.maximum(second - z["XE"], 0), NG - 2), NG - 1)
    edits = np.minimum(np.maximum(z["XE"], 0), NE - 1)
    ti = np.full(tier.shape, -1)
    for k, tv in enumerate(CONF):
        ti[tier == tv] = k
    conf = (mq >= 40) & (ti >= 0)
    ci = np.where(conf, ((ti * NE + edits) * NG + gap) * NX + bucket_x(z["XN"]), -1)
    li = np.where(~conf & (tier >= 0), tier * 40 + np.minimum(mq, 39), -1)
    return ci, li


def phred(e, n):
    return min(60, int(round(-10 * np.log10((e + 0.5) / (n + 1)))))


def cmd_fit(paths):
    nc, ec = np.zeros(len(CONF) * NE * NG * NX), np.zeros(len(CONF) * NE * NG * NX)
    nl, el = np.zeros(len(TIERS) * 40), np.zeros(len(TIERS) * 40)
    for p in paths:
        z = np.load(p)
        ci, li = cells(z)
        pl, bad = z["placed"], z["placed"] & ~z["ok"]
        for idx, n_, e_ in ((ci, nc, ec), (li, nl, el)):
            s = pl & (idx >= 0)
            n_ += np.bincount(idx[s], minlength=len(n_))
            e_ += np.bincount(idx[s & bad[...]], minlength=len(e_)) if (s & bad).any() else 0
    n4 = nc.reshape(len(CONF), NE, NG, NX)
    e4 = ec.reshape(len(CONF), NE, NG, NX)
    q = np.zeros(n4.shape, int)
    for t in range(len(CONF)):
        for d in range(NE):
            for g in range(NG):
                for x in range(NX):
                    for n_, e_ in ((n4[t, d, g, x], e4[t, d, g, x]), (n4[t, d, g].sum(), e4[t, d, g].sum()),
                                   (n4[t, :, g].sum(), e4[t, :, g].sum()), (n4[t].sum(), e4[t].sum())):
                        if n_ >= MIN_N:
                            q[t, d, g, x] = phred(e_, n_)
                            break
                    else:
                        q[t, d, g, x] = phred(e4[t].sum(), n4[t].sum())
    nl2, el2 = nl.reshape(len(TIERS), 40), el.reshape(len(TIERS), 40)
    ql = np.zeros(nl2.shape, int)
    for t in range(len(TIERS)):
        for m in range(40):
            n_, e_ = nl2[t, m], el2[t, m]
            if n_ < MIN_N:  # pool the tier's MAPQ band of 10
                b = slice(m // 10 * 10, m // 10 * 10 + 10)
                n_, e_ = nl2[t, b].sum(), el2[t, b].sum()
            if n_ < MIN_N:
                n_, e_ = nl2[t].sum(), el2[t].sum()
            ql[t, m] = phred(e_, n_) if n_ > 0 else 0
    out = sys.stdout
    out.write("// Generated by bench/truth_mapq_fit.py from " + ", ".join(p.rsplit("/", 1)[-1] for p in paths) +
              f": {int(nc.sum() + nl.sum()):,} mates with a truth.\n")
    out.write("// The MAPQ of `certa map --intervals --mapq cal` (paired): Phred of the share placed\n"
              "// more than 150 bp from the truth (bench/VALIDATION_RESULTS.md, \"A MAPQ fitted at the locus\").\n"
              "#pragma once\n#include <cstdint>\n\nnamespace certa {\n\n")
    out.write(f"// Records with MAPQ >= 40 in tiers S0, S1, PR: [tier][edits 0..{NE - 1}][gap 0..{NG - 2}, none]"
              "[clusters 1, 2-4, 5-16, 17+].\n")
    out.write(f"constexpr uint8_t kCalConf[{len(CONF)}][{NE}][{NG}][{NX}] = {{\n")
    for t in range(len(CONF)):
        out.write("  {  // " + TIERS[CONF[t]] + "\n")
        for d in range(NE):
            out.write("    {" + ", ".join("{" + ", ".join(str(v) for v in q[t, d, g]) + "}" for g in range(NG)) + "},\n")
        out.write("  },\n")
    out.write("};\n\n// Every other record: [tier S0..MR][its MAPQ 0..39].\n")
    out.write(f"constexpr uint8_t kCalLow[{len(TIERS)}][40] = {{\n")
    for t in range(len(TIERS)):
        out.write("  {" + ", ".join(str(v) for v in ql[t]) + "},  // " + TIERS[t] + "\n")
    out.write("};\n\n}  // namespace certa\n")


if __name__ == "__main__":
    {"feat": cmd_feat, "fit": cmd_fit}[sys.argv[1]](sys.argv[2:])
