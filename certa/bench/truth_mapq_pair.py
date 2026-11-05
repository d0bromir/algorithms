#!/usr/bin/env python3
"""The pair-aware MAPQ of `certa map --intervals --mapq pair`
(bench/VALIDATION_RESULTS.md, "A MAPQ from both mates"): gradient-boosted
trees over both mates' records, fitted at the locus on HG002, written as
include/certa/mapq_model.h.

  samtools view RUN.bam | truth_mapq_pair.py feat truth.npz REF.fai OUT.npz
  truth_mapq_pair.py fit [--cv] OUT.npz ... > include/certa/mapq_model.h
  truth_mapq_pair.py fit --t2t OUT.npz ... > include/certa/mapq_model_t2t.h
  truth_mapq_pair.py eval [--half] --model MODEL.h ... truth.npz OUT.npz NAME=RUN.npz ...
  samtools view PAIR.bam | truth_mapq_pair.py replay MODEL.h FEAT.npz

RUN is paired interval mode with --iv-mapq cert, --iv-pair all,
--mapq-limit 10 and --feature-tags. `feat` keeps per mate the SAM fields the
model reads (flag, MAPQ, |TLEN|, soft clips, XT, XE, XD, XO, XN, XB, XQ, AS,
NM, XS; -1 when absent and for an unmapped record; with --t2t also XC, the
mate's placement on T2T-CHM13), whether it is placed with a truth, and
whether it is at the locus (truth_eval.py's rule at 150 bp).

The model's inputs (`features`, mirrored by pair_model_features in
certa/interval.h) are, per mate, its own record's
  tier, MAPQ, gap to the second-best locus (XD, else XO, minus XE; 11: none),
  XE, log2 of the clusters verified, XB, XQ, NM, soft clips / 5, the pair
  margin of a PR mate (AS + the other mate's AS - XS), proper / mate-unmapped
  flags, and |TLEN| in classes;
and the other mate's tier, MAPQ, gap, XE, XQ, NM, soft clips and margin.
With --t2t (t2t_model_features) also the mate's T2T placement: whether it
lifts to the record's locus, elsewhere, or not at all; its gap to the
second-best locus there, XE minus its XE there, its MAPQ there / 5 and its
tier; and the other mate's class and gap.

`fit` trains trees of depth --depth (4) by Newton boosting on the log loss of
"placed more than 150 bp from the truth" (--rounds 300, rate 0.1, L2 1, every
misplaced mate and 1 in 8 of the others, weighted 8; seed 7), so the score is
the log odds of a misplacement. The committed T2T model is --depth 6
--rounds 400 (truth_pair.sh fit). MAPQ is -10 log10 of that probability, rounded, at most 60;
the header holds the log-odds thresholds of every MAPQ value, so that CERTA
and this script compute it with the same additions and comparisons. --cv fits
on the even pairs and reports the odd pairs at matched N.
"""
import argparse
import re
import sys

import numpy as np

sys.path.insert(0, __import__("os").path.dirname(__file__))
from truth_eval import read_sam  # noqa: E402
from truth_lift import CLS_TRUTH, Ref  # noqa: E402

TIERS = ("S0", "S1", "SR", "S2", "SL", "PR", "R", "IV", "MR")
PR = TIERS.index("PR")
KEYS = ("XE", "XD", "XO", "XN", "XB", "XQ", "AS", "NM", "XS")
XC = ("t_cls", "t_tier", "t_mq", "t_xe", "t_xd", "t_xo")  # --t2t: XC:Z:cls,tier,mapq,xe,xd,xo
TG = re.compile(r"\t(" + "|".join(KEYS) + r"):i:(-?\d+)")
CLIP = re.compile(r"(\d+)S")
TLEN_BINS = (1, 200, 300, 400, 500, 700, 1000)
OWN = ("tier", "mq", "gap", "xe", "xn", "xb", "xq", "nm", "clip", "margin", "pp", "tlen")
PART = ("tier", "mq", "gap", "xe", "xq", "nm", "clip", "margin")
NAMES = OWN + tuple("p_" + k for k in PART)
T2T = ("t_cls", "t_gap", "t_dx", "t_mq", "t_tier", "p_t_cls", "p_t_gap")
DEPTH, ROUNDS, RATE, LAM, MIN_H, KEEP, SEED = 4, 300, 0.1, 1.0, 20.0, 8, 7  # --depth, --rounds


def cmd_feat(a):
    t = np.load(a.truth)
    n = t["cls"].shape[0]
    ref = Ref(a.fai)
    f = {k: np.full((n, 2), -1, np.int32) for k in KEYS + XC}
    f["t_cls"][:] = 0
    tier = np.full((n, 2), -1, np.int8)
    flag2, clip, tlen, mq = (np.zeros((n, 2), d) for d in (np.int16, np.int16, np.int32, np.int16))

    def tee():
        for line in sys.stdin:
            if line[0] != "@":
                c = line.split("\t", 11)
                flag = int(c[1])
                if not flag & (256 | 2048):
                    i, m = int(c[0][1:]), 1 if flag & 128 else 0
                    flag2[i, m] = flag
                    if not flag & 4:
                        mq[i, m] = int(c[4])
                        tlen[i, m] = abs(int(c[8]))
                        clip[i, m] = sum(int(x) for x in CLIP.findall(c[5]))
                        if len(c) > 11:
                            rest = "\t" + c[11]
                            for k, v in TG.findall(rest):
                                f[k][i, m] = int(v)
                            x = re.search(r"\tXT:Z:(\w+)", rest)
                            if x and x.group(1) in TIERS:
                                tier[i, m] = TIERS.index(x.group(1))
                            x = re.search(r"\tXC:Z:([-\d,]+)", rest)
                            if x:
                                for k, v in zip(XC, x.group(1).split(",")):
                                    f[k][i, m] = int(v)
            yield line
    r = read_sam(tee(), n, ref.index)
    ref.par_to_x(r["rc"], r["left"])
    truth = t["cls"] == CLS_TRUTH
    ok = r["placed"] & truth & (r["rc"] == t["tc"]) & (r["strand"] == t["ts"]) & (np.abs(r["left"] + 1 - t["tp"]) <= 150)
    np.savez(a.out, truth=truth, placed=r["placed"], ok=ok, mapq=mq, tier=tier, flag=flag2, clip=clip, tlen=tlen,
             **f)


def log2i(v):
    """floor(log2(v)) for integers v >= 1 (exact)."""
    return np.frexp(v.astype(np.float64))[1].astype(np.int32) - 1


def features(z, t2t=False):
    """(n, 2, F) int32: each mate's own inputs, then the other mate's (pair_model_features); with t2t
    also the T2T-CHM13 inputs (t2t_model_features)."""
    mapped = (z["flag"] & 4) == 0
    g = lambda k: np.where(mapped, z[k], -1)  # noqa: E731
    tier = np.where(mapped, z["tier"].astype(np.int32), -1)
    xe, xd, xo = g("XE"), g("XD"), g("XO")
    second = np.where(xd >= 0, xd, xo)
    own = {
        "tier": tier + 1,
        "mq": np.where(mapped, np.clip(z["mapq"], 0, 60), 0),
        "gap": np.where(mapped & (second >= 0), np.clip(second - xe, 0, 10), 11),
        "xe": np.clip(xe, -1, 10) + 1,
        "xn": np.minimum(log2i(np.maximum(g("XN"), 0) + 1), 15),
        "xb": np.clip(g("XB"), 0, 3),
        "xq": np.clip(g("XQ"), 0, 60),
        "nm": np.clip(g("NM"), 0, 15),
        "clip": np.where(mapped, np.minimum(z["clip"] // 5, 15), 0),
        "pp": ((z["flag"] & 2) > 0) * 2 + ((z["flag"] & 8) > 0),
        "tlen": np.where(mapped, np.digitize(z["tlen"], TLEN_BINS), 0),
    }
    AS, xs = g("AS"), g("XS")
    pm = mapped[:, ::-1]
    own["margin"] = np.where((tier == PR) & (xs > -1000) & pm, np.clip(AS + AS[:, ::-1] - xs, 0, 40) + 1, 0)
    cols = [own[k] for k in OWN] + [own[k][:, ::-1] for k in PART]
    if t2t:
        cls = np.where(mapped, z["t_cls"], 0)
        on = cls > 0
        tsec = np.where(z["t_xd"] >= 0, z["t_xd"], z["t_xo"])
        tgap = np.where(~on, 12, np.where(tsec >= 0, np.clip(tsec - z["t_xe"], 0, 10), 11))
        tdx = np.where(on & (z["t_xe"] >= 0) & (xe >= 0), np.clip(xe - z["t_xe"], -3, 5) + 3, 9)
        cols += [cls, tgap, tdx, np.where(on, np.minimum(z["t_mq"], 60) // 5, 0), np.where(on, z["t_tier"] + 1, 0),
                 cls[:, ::-1], tgap[:, ::-1]]
    return np.stack(cols, axis=-1).astype(np.int32)


class Trees:
    """Trees of depth `depth` in level order: node k of level d at (1 << d) - 1 + k
    goes right when x[feat] > thr (feat -1: always left)."""

    def __init__(self, base, feat, thr, leaf, depth):
        self.base, self.feat, self.thr, self.leaf, self.depth = base, feat, thr, leaf, depth

    def score(self, X):
        n = X.shape[0]
        f = np.full(n, self.base)
        rows = np.arange(n)
        for t in range(len(self.leaf)):
            node = np.zeros(n, np.int64)
            for d in range(self.depth):
                k = (1 << d) - 1 + node
                fj = self.feat[t][k]
                go = (fj >= 0) & (X[rows, np.maximum(fj, 0)] > self.thr[t][k])
                node = node * 2 + go
            f += self.leaf[t][node]
        return f


def fit(X, y, w, base, depth=DEPTH, rounds=ROUNDS, log=sys.stderr):
    n, F = X.shape
    nb = X.max(axis=0) + 1
    f = np.full(n, base)
    nodes = (1 << depth) - 1
    feat = np.full((rounds, nodes), -1, np.int64)
    thr = np.zeros((rounds, nodes), np.int64)
    leaf = np.zeros((rounds, 1 << depth))
    rows = np.arange(n)
    for r in range(rounds):
        p = 1 / (1 + np.exp(-f))
        gr = (p - y) * w
        h = np.maximum(p * (1 - p), 1e-6) * w
        node = np.zeros(n, np.int64)
        for d in range(depth):
            nn = 1 << d
            best_gain = np.zeros(nn)
            for j in range(F):
                b = int(nb[j])
                key = node * b + X[:, j]
                G = np.bincount(key, gr, nn * b).reshape(nn, b)
                H = np.bincount(key, h, nn * b).reshape(nn, b)
                GL, HL = np.cumsum(G, 1), np.cumsum(H, 1)
                GT, HT = GL[:, -1:], HL[:, -1:]
                GR, HR = GT - GL, HT - HL
                gain = GL ** 2 / (HL + LAM) + GR ** 2 / (HR + LAM) - GT ** 2 / (HT + LAM)
                gain[(HL < MIN_H) | (HR < MIN_H)] = 0
                gain[:, -1] = 0
                bi = np.argmax(gain, axis=1)
                bg = gain[np.arange(nn), bi]
                better = bg > best_gain
                k = (1 << d) - 1 + np.arange(nn)
                feat[r, k[better]] = j
                thr[r, k[better]] = bi[better]
                best_gain = np.where(better, bg, best_gain)
            k = (1 << d) - 1 + node
            fj = feat[r, k]
            go = (fj >= 0) & (X[rows, np.maximum(fj, 0)] > thr[r, k])
            node = node * 2 + go
        G = np.bincount(node, gr, 1 << depth)
        H = np.bincount(node, h, 1 << depth)
        leaf[r] = -G / (H + LAM) * RATE
        f += leaf[r][node]
        if r % 50 == 0 or r == rounds - 1:
            ll = -np.sum(w * (y * np.log(p + 1e-12) + (1 - y) * np.log(1 - p + 1e-12))) / w.sum()
            print(f"[fit] round {r}: log loss {ll:.6f}", file=log, flush=True)
    return Trees(base, feat, thr, leaf, depth)


def q_thresholds():
    """MAPQ q >= 1 when the log odds are <= T[q]: -10 log10 p >= q - 0.5."""
    p = 10.0 ** (-(np.arange(1, 61) - 0.5) / 10)
    return np.log(p / (1 - p))


def to_mapq(score, T):
    return (score[:, None] <= T[None, :]).sum(axis=1)


def training_set(paths, half, t2t):
    Xs, ys = [], []
    for p in paths:
        z = np.load(p)
        F = features(z, t2t)
        n = F.shape[0]
        sel = z["placed"] & z["truth"]
        if half is not None:
            sel &= (np.arange(n) % 2 == half)[:, None]
        Xs.append(F[sel])
        ys.append((~z["ok"][sel]).astype(np.float64))
    X, y = np.concatenate(Xs), np.concatenate(ys)
    rng = np.random.default_rng(SEED)
    keep = (y == 1) | (rng.random(len(y)) < 1 / KEEP)
    X, y = X[keep], y[keep]
    w = np.where(y == 1, 1.0, float(KEEP))
    return X, y, w, float(np.log(np.sum(w * y) / np.sum(w * (1 - y))))


def write_header(m, paths, t2t, out=sys.stdout):
    T = q_thresholds()
    nodes = (1 << m.depth) - 1
    k = "kT2tModel" if t2t else "kModel"
    nf = len(NAMES) + (len(T2T) if t2t else 0)
    what = ("both mates' records and their placements on T2T-CHM13 (--t2t)" if t2t else "both mates' records")
    out.write("// Generated by bench/truth_mapq_pair.py from " + ", ".join(p.rsplit("/", 1)[-1] for p in paths) +
              f".\n// The MAPQ of `certa map --intervals --mapq pair`: boosted trees over {what}\n"
              "// (bench/VALIDATION_RESULTS.md, \"A MAPQ from both mates\"). Inputs: certa/pair_model.h.\n"
              "#pragma once\n#include <cstdint>\n\nnamespace certa {\n\n")
    out.write(f"constexpr int {k}Features = {nf}, {k}Trees = {len(m.leaf)}, {k}Depth = {m.depth};\n")
    out.write(f"constexpr double {k}Base = {m.base!r};\n")
    out.write(f"// MAPQ q >= 1 when the score is <= {k}Q[q - 1].\n")
    out.write(f"constexpr double {k}Q[60] = {{" + ", ".join(repr(float(v)) for v in T) + "};\n")
    out.write("// Per tree, node k of level d at (1 << d) - 1 + k: right when x[feat] > thr (feat -1: left).\n")
    for name, arr, typ, width in (("Feat", m.feat, "int8_t", nodes), ("Thr", m.thr, "int8_t", nodes),
                                  ("Leaf", m.leaf, "double", 1 << m.depth)):
        out.write(f"constexpr {typ} {k}{name}[{len(m.leaf)}][{width}] = {{\n")
        for t in range(len(m.leaf)):
            vals = (repr(float(v)) for v in arr[t]) if typ == "double" else (str(int(v)) for v in arr[t])
            out.write("  {" + ", ".join(vals) + "},\n")
        out.write("};\n")
    out.write("\n}  // namespace certa\n")


def read_header(path):
    """(trees, MAPQ thresholds, whether the model reads the T2T inputs)."""
    s = open(path).read()
    t2t = "kT2tModelBase" in s
    k = "kT2tModel" if t2t else "kModel"

    def arr(name):
        body = s[s.index(name):]
        body = body[body.index("{") + 1:body.index("};")]
        return [float(v) for v in re.findall(r"-?\d+(?:\.\d*)?(?:e-?\d+)?", body)]
    nt = int(re.search(k + r"Trees = (\d+)", s).group(1))
    depth = int(re.search(k + r"Depth = (\d+)", s).group(1))
    nodes = (1 << depth) - 1
    base = float(re.search(k + r"Base = (\S+);", s).group(1))
    feat = np.array(arr(k + "Feat["), np.int64).reshape(nt, nodes)
    thr = np.array(arr(k + "Thr["), np.int64).reshape(nt, nodes)
    leaf = np.array(arr(k + "Leaf["), np.float64).reshape(nt, 1 << depth)
    T = np.array(arr(k + "Q["), np.float64)
    return Trees(base, feat, thr, leaf, depth), T, t2t


def at_n(mq, bad, n):
    """Misplaced among the n most confident (the boundary MAPQ pro rata)."""
    order = np.argsort(-mq, kind="stable")
    m, b = mq[order], bad[order]
    if n >= len(m):
        return float(b.sum())
    q = m[n - 1]
    above = m > q
    at = m == q
    return float(b[above].sum()) + (n - int(above.sum())) * float(b[at].sum()) / float(at.sum())


def cmd_fit(a):
    X, y, w, base = training_set(a.npz, 0 if a.cv else None, a.t2t)
    print(f"[fit] {len(y):,} rows, {int(y.sum()):,} misplaced", file=sys.stderr)
    m = fit(X, y, w, base, a.depth, a.rounds)
    write_header(m, a.npz, a.t2t)


def mapq_of(z, m, T, t2t, mask=None):
    """The model's MAPQ of every mapped mate (of the mates in mask, if given; 0 elsewhere)."""
    F = features(z, t2t)
    sel = (z["flag"] & 4) == 0
    if mask is not None:
        sel = sel & mask
    q = np.zeros(z["flag"].shape, np.int64)
    q[sel] = to_mapq(m.score(F[sel]), T)
    return q


def cmd_eval(a):
    """Misplaced per million at matched N: the models on FEAT.npz (--half: odd pairs only), the run's own
    MAPQ, and other runs (truth_t2t.py place npz files) on the same mates."""
    z = np.load(a.npz)
    n = z["truth"].shape[0]
    half = (np.arange(n) % 2 == 1)[:, None] if a.half else np.ones((n, 1), bool)
    t = np.load(a.truth)
    truth = (t["cls"] == CLS_TRUTH) & half
    nt = int(truth.sum())
    runs = {"the run's MAPQ": (z["mapq"].astype(np.int64), z["placed"], z["ok"])}
    for path in a.model:
        m, T, t2t = read_header(path)
        runs[f"--mapq pair ({path.rsplit('/', 1)[-1]})"] = (mapq_of(z, m, T, t2t, z["placed"] & truth), z["placed"],
                                                           z["ok"])
    for spec in a.runs:
        name, p = spec.split("=", 1)
        r = np.load(p)
        ok = r["placed"] & (r["rc"] == t["tc"]) & (r["strand"] == t["ts"]) & \
            (np.abs(r["left"].astype(np.int64) - (t["tp"].astype(np.int64) - 1)) <= 150)
        runs[name] = (r["mapq"].astype(np.int64), r["placed"], ok)
    shares = [float(s) for s in a.shares.split(",")]
    print("| run | " + " | ".join(f"{100 * s:.2f} %" for s in shares) + " | MAPQ >= 40 |")
    print("|---|" + "---|" * (len(shares) + 1))
    for name, (mq, pl, ok) in runs.items():
        sel = pl & truth
        row = [at_n(mq[sel], ~ok[sel], int(round(s * nt))) / int(round(s * nt)) * 1e6 for s in shares]
        print(f"| {name} | " + " | ".join(f"{v:,.0f}" for v in row) + f" | {100 * np.mean(mq[sel] >= 40):.2f} % |")
    print(f"\n{nt:,} mates with a truth; misplaced per million among each run's N most confident.")


def cmd_replay(a):
    """CERTA's --mapq pair MAPQ (stdin) against this script's on the same records (FEAT.npz: the same run's
    features, from --feature-tags)."""
    m, T, t2t = read_header(a.model)
    z = np.load(a.npz)
    want = mapq_of(z, m, T, t2t)
    n = want.shape[0]
    got = np.full((n, 2), -1, np.int64)
    for line in sys.stdin:
        if line[0] == "@":
            continue
        c = line.split("\t", 6)
        flag = int(c[1])
        if flag & (256 | 2048 | 4):
            continue
        got[int(c[0][1:]), 1 if flag & 128 else 0] = int(c[4])
    mapped = (z["flag"] & 4) == 0
    bad = mapped & (got != want)
    print(f"[replay] {int(mapped.sum()):,} mapped mates, {int(bad.sum()):,} differ")
    sys.exit(1 if bad.any() else 0)


def main():
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    p = sp.add_parser("feat")
    p.add_argument("truth"), p.add_argument("fai"), p.add_argument("out")
    p = sp.add_parser("fit")
    p.add_argument("--cv", action="store_true")
    p.add_argument("--t2t", action="store_true")
    p.add_argument("--depth", type=int, default=DEPTH)
    p.add_argument("--rounds", type=int, default=ROUNDS)
    p.add_argument("npz", nargs="+")
    p = sp.add_parser("eval")
    p.add_argument("--half", action="store_true")
    p.add_argument("--shares", default="0.9975,0.9929,0.98,0.95")
    p.add_argument("--model", action="append", default=[])
    p.add_argument("truth"), p.add_argument("npz"), p.add_argument("runs", nargs="*")
    p = sp.add_parser("replay")
    p.add_argument("model"), p.add_argument("npz")
    a = ap.parse_args()
    {"feat": cmd_feat, "fit": cmd_fit, "eval": cmd_eval, "replay": cmd_replay}[a.cmd](a)


if __name__ == "__main__":
    main()
