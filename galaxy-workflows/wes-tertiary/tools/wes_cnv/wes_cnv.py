#!/usr/bin/env python3
"""Exome copy-number calls from per-target read depth against a panel of samples of the same kit (read-depth HMM).

Input are mosdepth "regions" BEDs (chrom, start, end, [name,] mean depth) of the capture targets: one for the sample,
three or more for the reference panel, all from the same capture kit, sequencer and pipeline. The steps:

1. Targets: kept when the panel's median depth is >= --min-depth and the panel varies little there (relative
   spread <= --max-panel-cv); low-mappability and failing targets drop out on their own.
2. Normalisation: every sample is divided by its median autosomal target depth, then corrected for GC content
   (median ratio per GC bin, from the reference FASTA) and target length.
3. Sex chromosomes: each sample's X (outside the pseudoautosomal regions) is brought to a two-copy scale by its X
   copy number, so female and male panel samples serve both sexes. Y is not called.
4. Ratio and noise: ratio = sample / panel median; each target's error combines the panel's spread there and the
   read-sampling error, scaled per sample so that the median |z| of the autosomes is that of a standard normal.
5. Calls: Viterbi path of a hidden Markov model with copy number 0..4 (0..2 on a male X), entry probability
   --p-cnv per target, expected length --mean-targets targets, reset across gaps > 1 Mb. Every segment off the normal
   copy number is reported with its log10 Bayes factor against the normal copy number.

Also: chromosome-level dosage (median ratio per chromosome: trisomy, monosomy, sex-chromosome complement).

Outputs: --calls TSV (one row per call), --targets TSV (per-target ratios, optional), --summary JSON.
Needs Python >= 3.10 and numpy.
"""
import argparse
import gzip
import json
import math
import os
import sys

import numpy as np

VERSION = "0.1.0"
PAR = {"X": [(10000, 2781479), (155701382, 156030895)], "Y": [(10000, 2781479), (56887902, 57217415)]}
AUTOSOMES = [str(i) for i in range(1, 23)]


def norm_chrom(c):
    c = c[3:] if c.startswith("chr") else c
    return "MT" if c == "M" else c


def opener(path):
    with open(path, "rb") as f:
        gz = f.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if gz else open(path)


def read_depth(path):
    """mosdepth regions BED -> {(chrom, start, end): mean depth}, plus the region name if present."""
    d, names = {}, {}
    with opener(path) as f:
        for line in f:
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue
            t = line.rstrip("\n").split("\t")
            k = (norm_chrom(t[0]), int(t[1]), int(t[2]))
            d[k] = float(t[-1])
            if len(t) >= 5:
                names[k] = t[3]
    return d, names


def in_par(chrom, pos):
    return any(a <= pos <= b for a, b in PAR.get(chrom, ()))


class Fasta:
    """Random access to an uncompressed FASTA through its .fai (built here if missing)."""

    def __init__(self, path):
        fai = path + ".fai"
        lines = open(fai).read().splitlines() if os.path.exists(fai) else self._build_fai(path)
        self.idx = {}
        for line in lines:
            n, ln, off, lb, lw = line.split("\t")[:5]
            self.idx[norm_chrom(n)] = (int(ln), int(off), int(lb), int(lw))
        self.f = open(path, "rb")

    @staticmethod
    def _build_fai(path):
        out, name, ln, off, lb, lw, pos = [], None, 0, 0, 0, 0, 0
        with open(path, "rb") as f:
            for line in f:
                if line.startswith(b">"):
                    if name is not None:
                        out.append(f"{name}\t{ln}\t{off}\t{lb}\t{lw}")
                    name, ln, off, lb, lw = line[1:].split()[0].decode(), 0, pos + len(line), 0, 0
                else:
                    if not lb:
                        lb, lw = len(line.rstrip(b"\r\n")), len(line)
                    ln += len(line.rstrip(b"\r\n"))
                pos += len(line)
        if name is not None:
            out.append(f"{name}\t{ln}\t{off}\t{lb}\t{lw}")
        return out

    def gc(self, chrom, start, end):
        if chrom not in self.idx:
            return float("nan")
        ln, off, lb, lw = self.idx[chrom]
        start, end = max(0, start), min(end, ln)
        if end <= start:
            return float("nan")
        a = off + (start // lb) * lw + start % lb
        b = off + (end // lb) * lw + end % lb
        self.f.seek(a)
        s = self.f.read(b - a).replace(b"\n", b"").replace(b"\r", b"").upper()
        gc, at = s.count(b"G") + s.count(b"C"), s.count(b"A") + s.count(b"T")
        return gc / (gc + at) if gc + at else float("nan")


def binned_correction(ratio, covariate, nbins=40, min_n=200):
    """Median of ratio in quantile bins of the covariate; returns the per-target correction factor."""
    ok = np.isfinite(ratio) & np.isfinite(covariate) & (ratio > 0)
    corr = np.ones_like(ratio)
    if ok.sum() < nbins * min_n // 4:
        return corr
    edges = np.unique(np.quantile(covariate[ok], np.linspace(0, 1, nbins + 1)))
    if len(edges) < 3:  # (nearly) constant covariate, e.g. targets of one length: nothing to correct
        return corr
    which = np.clip(np.searchsorted(edges, covariate, side="right") - 1, 0, len(edges) - 2)
    med = np.ones(len(edges) - 1)
    for b in range(len(edges) - 1):
        sel = ok & (which == b)
        if sel.sum() >= 20:
            med[b] = np.median(ratio[sel])
    good = np.isfinite(covariate)
    corr[good] = med[which[good]]
    corr[~np.isfinite(corr) | (corr <= 0)] = 1.0
    return corr


def fixed_width_correction(ratio, covariate, width=0.01, min_n=60):
    """Median ratio in fixed-width covariate bins (sparse bins merged with their neighbours): GC extremes keep their
    own bins instead of being pooled into one quantile bin."""
    ok = np.isfinite(ratio) & np.isfinite(covariate) & (ratio > 0)
    corr = np.ones_like(ratio)
    if ok.sum() < 1000:
        return corr
    b = np.floor(covariate / width).astype(int)
    lo, hi = int(b[ok].min()), int(b[ok].max())
    counts = np.bincount(b[ok] - lo, minlength=hi - lo + 1)
    groups, cur, n = [], [], 0  # consecutive bins merged until each group has min_n targets
    for i, c in enumerate(counts):
        cur.append(i)
        n += c
        if n >= min_n:
            groups.append(cur)
            cur, n = [], 0
    if cur:
        if groups:
            groups[-1] += cur
        else:
            groups.append(cur)
    group_of = np.empty(len(counts), dtype=int)
    for gi, g in enumerate(groups):
        group_of[g] = gi
    med = np.array([np.median(ratio[ok & np.isin(b - lo, g)]) for g in groups])
    fin = np.isfinite(covariate)
    corr[fin] = med[group_of[np.clip(b[fin] - lo, 0, len(counts) - 1)]]
    corr[~np.isfinite(corr) | (corr <= 0)] = 1.0
    return corr


def mad(x, axis=None):
    m = np.median(x, axis=axis, keepdims=True)
    return 1.4826 * np.median(np.abs(x - m), axis=axis)


def x_copies(norm_x):
    """X copy number of a sample from its median non-PAR X depth relative to its autosomes (1, 2, 3 ...)."""
    r = float(np.median(norm_x)) if len(norm_x) else float("nan")
    return max(1, int(round(2 * r))) if np.isfinite(r) else 2, r


def viterbi(ratio, sd, states_mean, chrom_break, p_cnv, mean_len, normal, df=4.0):
    """Log-space Viterbi over targets. states_mean: (T, S) expected ratio per target and state."""
    T, S = states_mean.shape
    stay = 1.0 - 1.0 / mean_len
    logA = np.full((S, S), -np.inf)
    for i in range(S):
        for j in range(S):
            if i == normal:
                logA[i, j] = math.log(1 - p_cnv) if j == normal else math.log(p_cnv / (S - 1))
            else:
                logA[i, j] = math.log(stay) if j == i else (math.log(1 - stay) if j == normal else math.log(1e-9))
    logpi = np.full(S, math.log(p_cnv / (S - 1)))
    logpi[normal] = math.log(1 - p_cnv)
    # emission: normal in ratio space, sd grows with the expected ratio (multiplicative error), floor for CN0
    sds = sd[:, None] * np.sqrt(np.maximum(states_mean, 0.25))
    # Student-t (df) emission: single outlying targets carry little weight, runs of targets a lot
    logB = -0.5 * (df + 1) * np.log1p(((ratio[:, None] - states_mean) / sds) ** 2 / df) - np.log(sds)
    back = np.zeros((T, S), dtype=np.int8)
    v = logpi + logB[0]
    for t in range(1, T):
        if chrom_break[t]:
            prev = np.full(S, v.max()) + logpi
            back[t] = int(v.argmax())
            v = prev + logB[t]
            continue
        cand = v[:, None] + logA
        back[t] = cand.argmax(axis=0)
        v = cand.max(axis=0) + logB[t]
    path = np.empty(T, dtype=np.int8)
    path[-1] = int(v.argmax())
    for t in range(T - 1, 0, -1):
        path[t - 1] = back[t, path[t]]
    return path, logB


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--sample", required=True, help="mosdepth regions BED of the sample (capture targets)")
    ap.add_argument("--sample-name", default=None)
    ap.add_argument("--panel", nargs="+", required=True, help="mosdepth regions BEDs of reference samples, same kit")
    ap.add_argument("--panel-names", nargs="*", default=None, help="names of the panel files (to drop the sample)")
    ap.add_argument("--reference", help="reference FASTA (GC correction)")
    ap.add_argument("--genes", help="BED of coding exons with the gene symbol in column 4 (annotation)")
    ap.add_argument("--segdups", help="BED of segmental duplications (e.g. GIAB GRCh38_segdups.bed.gz): calls with most "
                                      "of their targets inside are marked low confidence")
    ap.add_argument("--min-depth", type=float, default=20.0)
    ap.add_argument("--max-panel-cv", type=float, default=0.35)
    ap.add_argument("--p-cnv", type=float, default=1e-4)
    ap.add_argument("--mean-targets", type=float, default=8.0)
    ap.add_argument("--min-bf", type=float, default=3.0, help="minimum log10 Bayes factor of a reported call")
    ap.add_argument("--read-length", type=float, default=150.0)
    ap.add_argument("--shrink", type=float, default=10.0, help="prior weight (in samples) of the typical spread")
    ap.add_argument("--pcs", type=int, default=1, help="principal components of the panel removed (batch effects)")
    ap.add_argument("--df", type=float, default=10.0, help="degrees of freedom of the Student-t emission")
    ap.add_argument("--calls", required=True)
    ap.add_argument("--targets", default=None)
    ap.add_argument("--summary", default=None)
    ap.add_argument("--version", action="version", version=f"wes_cnv {VERSION}")
    a = ap.parse_args(argv)

    name = a.sample_name or os.path.basename(a.sample).split(".")[0]
    panel_files = list(a.panel)
    pnames = a.panel_names if a.panel_names and len(a.panel_names) == len(panel_files) else \
        [os.path.basename(p).split(".")[0] for p in panel_files]
    keep = [(p, n) for p, n in zip(panel_files, pnames) if n != name and os.path.realpath(p) != os.path.realpath(a.sample)]
    cols = ["sample", "chrom", "start", "end", "length_bp", "type", "copy_number", "normal_copies", "targets",
            "mean_ratio", "log10_bf", "segdup_fraction", "panel_carriers", "confidence", "genes"]
    if len(keep) < 3:  # not an error: a run with too few samples gets an empty, explained result
        msg = f"not analysed: {len(keep)} panel samples other than {name} (at least 3 of the same kit are needed)"
        with open(a.calls, "w") as f:
            f.write(f"##wes_cnv {VERSION}; {msg}\n#" + "\t".join(cols) + "\n")
        if a.summary:
            json.dump({"version": VERSION, "sample": name, "status": msg}, open(a.summary, "w"), indent=1)
        if a.targets:
            open(a.targets, "w").write("chrom\tstart\tend\tdepth\tpanel_depth\tratio\tsd\tused\tcall_copy\tgc\n")
        print(msg)
        return
    sd_ = read_depth(a.sample)[0]
    panel = [read_depth(p)[0] for p, _ in keep]
    keys = [k for k in sd_ if all(k in p for p in panel)]
    keys.sort(key=lambda k: ((AUTOSOMES + ["X", "Y"]).index(k[0]) if k[0] in AUTOSOMES + ["X", "Y"] else 99, k[1]))
    # chromosome Y (not called): its depth relative to the autosomes, for the sex-chromosome complement
    ykeys = [k for k in keys if k[0] == "Y" and not in_par("Y", (k[1] + k[2]) // 2)]
    keys = [k for k in keys if k[0] in AUTOSOMES or k[0] == "X"]
    if len(keys) < 1000:
        sys.exit(f"wes_cnv: only {len(keys)} targets shared by the sample and the panel: different capture BEDs?")
    T, P = len(keys), len(panel)
    chrom = np.array([k[0] for k in keys])
    start = np.array([k[1] for k in keys])
    end = np.array([k[2] for k in keys])
    length = (end - start).astype(float)
    M = np.empty((P + 1, T))
    M[0] = [sd_[k] for k in keys]
    for i, p in enumerate(panel):
        M[i + 1] = [p[k] for k in keys]
    auto = np.isin(chrom, AUTOSOMES)
    xnp = (chrom == "X") & ~np.array([in_par("X", (s + e) // 2) for s, e in zip(start, end)])
    raw_scale = np.median(M[:, auto], axis=1)  # median autosomal target depth per sample
    N = M / raw_scale[:, None]

    # X copies per sample; X brought to a two-copy scale
    xc = []
    for i in range(P + 1):
        sel = xnp & (np.median(N[1:], axis=0) > 0.2)
        c, r = x_copies(N[i, sel])
        xc.append((c, r))
        N[i, xnp] *= 2.0 / c

    # 75th percentile: about half of the Y targets have no uniquely mapped reads, so the median is unstable
    y_frac = [float(np.percentile([dd[k] for k in ykeys], 75)) / raw_scale[i] if len(ykeys) >= 20 else float("nan")
              for i, dd in enumerate([sd_] + panel)]

    # GC and length correction (two rounds against the panel median)
    gc = np.full(T, np.nan)
    if a.reference:
        fa = Fasta(a.reference)
        gc = np.array([fa.gc(c, s, e) for c, s, e in zip(chrom, start, end)])
    for _ in range(2):
        ref = np.median(N[1:], axis=0)
        okr = ref > 0.1
        for i in range(P + 1):
            r = np.where(okr, N[i] / np.where(okr, ref, 1), np.nan)
            if np.isfinite(gc).any():
                N[i] /= fixed_width_correction(r, gc)
                r = np.where(okr, N[i] / np.where(okr, ref, 1), np.nan)
            N[i] /= binned_correction(r, np.log(length))
        scale = np.median(N[:, auto], axis=1)
        N /= scale[:, None]

    ref = np.median(N[1:], axis=0)
    ref_depth = ref * raw_scale[0]  # in the sample's depth units
    pre = (ref_depth >= a.min_depth) & (ref > 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        L = np.log2(np.maximum(N[:, pre] / ref[pre], 0.05))
    # batch effects: the panel's main axes of variation in log ratio (library preparation, sequencing run) are
    # projected out of every sample (as in XHMM/CoNIFER); a panel from mixed batches otherwise calls the batch
    n_pc = max(0, min(a.pcs, P - 3))
    if n_pc:
        Lp = L[1:] - np.median(L[1:], axis=1, keepdims=True)
        vt = np.linalg.svd(Lp, full_matrices=False)[2]
        V = vt[:n_pc]
        L = L - (L @ V.T) @ V
        L -= np.median(L[:, auto[pre]], axis=1, keepdims=True)
    R = np.full((P + 1, T), np.nan)
    R[:, pre] = 2.0 ** L
    pcv = np.full(T, np.nan)
    pcv[pre] = mad(R[1:, pre], axis=0) * math.sqrt(P / max(P - n_pc, 1))  # panel relative spread per target
    use = pre & np.isfinite(pcv) & (pcv <= a.max_panel_cv)
    ratio = np.where(use, R[0], np.nan)
    reads = np.maximum(ref_depth * length / a.read_length, 1.0)
    # per-target error: read sampling plus the panel's extra spread there, shrunk towards the typical spread
    # (a panel of a few samples estimates one target's spread poorly)
    extra = np.maximum(pcv ** 2 - 1.0 / reads, 0.0)
    extra_typ = float(np.nanmedian(extra[use])) if use.any() else 0.0
    extra = np.where(np.isfinite(extra), (P * extra + a.shrink * extra_typ) / (P + a.shrink), extra_typ)
    sigma = np.sqrt(1.0 / reads + extra * (1 + 1.0 / P))
    sigma = np.maximum(sigma, 0.02)
    z = (ratio - 1) / sigma
    k = float(np.nanmedian(np.abs(z[use & auto])) / 0.6745) if (use & auto).any() else 1.0
    sd = sigma * k

    # chromosome dosage
    chrom_dosage = {}
    for c in AUTOSOMES + ["X"]:
        sel = use & (chrom == c) & (xnp if c == "X" else True)
        if sel.sum() >= 30:
            chrom_dosage[c] = round(float(np.median(ratio[sel])), 3)
    sex_c, sex_r = xc[0]
    male_x = sex_c == 1
    male_y = [y for (c_, _), y in zip(xc[1:], y_frac[1:]) if c_ == 1 and np.isfinite(y)]
    y_ref = float(np.median(male_y)) if male_y else 0.45
    y_rel = y_frac[0] / y_ref if np.isfinite(y_frac[0]) and y_ref > 0.05 else None
    y_copies = None if y_rel is None else (0 if y_rel < 0.2 else (1 if y_rel < 1.75 else 2))
    sex_complement = ("X" * sex_c + "Y" * y_copies) if y_copies is not None else None
    aneuploidy = []
    for c, r_ in chrom_dosage.items():
        if c == "X":
            continue
        if r_ >= 1.35:
            aneuploidy.append(f"chromosome {c}: gain of a whole copy (trisomy), median ratio {r_}")
        elif r_ >= 1.10:
            aneuploidy.append(f"chromosome {c}: partial gain (mosaic trisomy or a large duplication), median ratio {r_}")
        elif r_ <= 0.65:
            aneuploidy.append(f"chromosome {c}: loss of a whole copy (monosomy), median ratio {r_}")
        elif r_ <= 0.90:
            aneuploidy.append(f"chromosome {c}: partial loss (mosaic monosomy or a large deletion), median ratio {r_}")
    if sex_complement and sex_complement not in ("XX", "XY"):
        aneuploidy.append(f"sex chromosomes {sex_complement} (X ratio {sex_r:.2f}, Y ratio {y_rel:.2f} of a male)")

    # HMM over the targets in use: autosomes, PAR and a female X with two copies as normal; a male X with one
    idx_all = np.where(use)[0]
    male_block = (chrom[idx_all] == "X") & xnp[idx_all] & male_x
    blocks = [(idx_all[~male_block], 2), (idx_all[male_block], 1)]
    segs = []  # (indices, path, logB, normal)
    for idx, normal in blocks:
        if not len(idx):
            continue
        cn = np.arange(5, dtype=float)
        means = np.tile(np.where(cn == 0, 0.03, cn / normal), (len(idx), 1))
        brk = np.ones(len(idx), dtype=bool)
        brk[1:] = (chrom[idx][1:] != chrom[idx][:-1]) | (start[idx][1:] - end[idx][:-1] > 1_000_000)
        path, logB = viterbi(ratio[idx], sd[idx], means, brk, a.p_cnv, a.mean_targets, normal, a.df)
        segs.append((idx, path, logB, normal, brk))

    # genes
    gene_iv = {}
    if a.genes:
        with opener(a.genes) as f:
            for line in f:
                t = line.rstrip("\n").split("\t")
                if len(t) >= 4 and not line.startswith(("#", "track")):
                    gene_iv.setdefault(norm_chrom(t[0]), []).append((int(t[1]), int(t[2]), t[3]))
        for c in gene_iv:
            gene_iv[c].sort()

    exons_total = {}
    for c, v in gene_iv.items():
        for _, _, g in v:
            exons_total[(c, g)] = exons_total.get((c, g), 0) + 1

    def genes_in(c, s, e):
        """GENE:k/n per gene with coding exons in [s, e): k exons overlapped of the gene's n."""
        hit = {}
        for gs, ge, g in gene_iv.get(c, ()):
            if gs >= e:
                break
            if ge > s and g:
                hit[g] = hit.get(g, 0) + 1
        return [f"{g}:{k}/{exons_total[(c, g)]}" for g, k in hit.items()]

    in_sd = np.zeros(T, dtype=bool)
    if a.segdups:
        sdi = {}
        with opener(a.segdups) as f:
            for line in f:
                t = line.split("\t")
                if len(t) >= 3 and not line.startswith(("#", "track")):
                    sdi.setdefault(norm_chrom(t[0]), []).append((int(t[1]), int(t[2])))
        for c, v in sdi.items():
            v.sort()
            st = np.array([x[0] for x in v])
            en = np.maximum.accumulate(np.array([x[1] for x in v]))
            sel = np.where(chrom == c)[0]
            mid = (start[sel] + end[sel]) // 2
            j = np.searchsorted(st, mid, side="right") - 1
            in_sd[sel] = (j >= 0) & (en[np.maximum(j, 0)] > mid)

    def panel_carriers(tidx, is_del):
        """Panel samples whose depth over the call's targets moves the same way (common CNV or systematic artefact)."""
        n = 0
        for i_ in range(1, P + 1):
            others = np.delete(R[1:, tidx], i_ - 1, axis=0)
            rr = float(np.mean(R[i_, tidx] / np.maximum(np.median(others, axis=0), 1e-6)))
            n += (rr < 0.75) if is_del else (rr > 1.25)
        return n

    def confidence(c):
        """high: >= 3 targets, mostly outside segmental duplications, not seen in the panel, log10 BF >= 5."""
        why = []
        if c["targets"] < 3:
            why.append(f"{c['targets']} target{'s' if c['targets'] > 1 else ''}")
        if c["segdup_fraction"] >= 0.5:
            why.append("segmental duplication")
        if c["panel_carriers"] >= 2:
            why.append(f"also in {c['panel_carriers']} panel samples")
        if c["log10_bf"] < 5:
            why.append("weak evidence")
        return "high" if not why else "low: " + ", ".join(why)

    calls, copy_at = [], {}
    for idx, path, logB, normal, brk in segs:
        for u, t in enumerate(idx):
            copy_at[t] = int(path[u])
        i, n = 0, len(idx)
        while i < n:
            j = i
            while j + 1 < n and not brk[j + 1] and path[j + 1] == path[i]:
                j += 1
            if path[i] != normal:
                lbf = float(np.sum(logB[i:j + 1, path[i]] - logB[i:j + 1, normal])) / math.log(10)
                if lbf >= a.min_bf:
                    c, s_, e_ = chrom[idx[i]], int(start[idx[i]]), int(end[idx[j]])
                    calls.append({
                        "sample": name, "chrom": c, "start": s_, "end": e_, "length_bp": e_ - s_,
                        "type": "DEL" if path[i] < normal else "DUP", "copy_number": int(path[i]),
                        "normal_copies": normal, "targets": j - i + 1,
                        "mean_ratio": round(float(np.mean(ratio[idx[i:j + 1]])), 3), "log10_bf": round(lbf, 1),
                        "segdup_fraction": round(float(np.mean(in_sd[idx[i:j + 1]])), 2),
                        "panel_carriers": panel_carriers(idx[i:j + 1], path[i] < normal),
                        "genes": ",".join(genes_in(c, s_, e_)),
                    })
                    calls[-1]["confidence"] = confidence(calls[-1])
            i = j + 1
    calls.sort(key=lambda c: ((AUTOSOMES + ["X"]).index(c["chrom"]), c["start"]))

    with open(a.calls, "w") as f:
        f.write(f"##wes_cnv {VERSION}; panel {P} samples: {','.join(n for _, n in keep)}; noise {k:.2f}\n")
        f.write("#" + "\t".join(cols) + "\n")
        for c in calls:
            f.write("\t".join(str(c[x]) for x in cols) + "\n")
    if a.targets:
        with open(a.targets, "w") as f:
            f.write("chrom\tstart\tend\tdepth\tpanel_depth\tratio\tsd\tused\tcall_copy\tgc\n")
            for t in range(T):
                u = copy_at.get(t)
                f.write(f"{chrom[t]}\t{start[t]}\t{end[t]}\t{M[0, t]:.1f}\t{ref_depth[t]:.1f}\t"
                        + (f"{ratio[t]:.3f}" if np.isfinite(ratio[t]) else "NA") + f"\t{sd[t]:.3f}\t{int(use[t])}\t"
                        + (str(u) if u is not None else "NA") + f"\t{gc[t]:.3f}\n")
    summ = {
        "version": VERSION, "status": "analysed", "sample": name, "panel": [n for _, n in keep], "targets_total": T,
        "targets_used": int(use.sum()), "noise_inflation": round(k, 3),
        "ratio_sd_autosomes": round(float(mad(ratio[use & auto])), 4),
        "x_copies": sex_c, "x_ratio": round(sex_r, 3), "chromosome_dosage": chrom_dosage, "pcs_removed": n_pc,
        "calls": len(calls), "deletions": sum(c["type"] == "DEL" for c in calls),
        "duplications": sum(c["type"] == "DUP" for c in calls),
        "panel_x_copies": [c for c, _ in xc[1:]], "y_ratio": round(y_rel, 3) if y_rel is not None else None,
        "sex_complement": sex_complement, "aneuploidy": aneuploidy,
        "high_confidence": sum(c["confidence"] == "high" for c in calls),
    }
    if a.summary:
        json.dump(summ, open(a.summary, "w"), indent=1)
    print(json.dumps({k2: summ[k2] for k2 in ("sample", "targets_used", "ratio_sd_autosomes", "sex_complement",
                                               "aneuploidy", "calls", "high_confidence")}))


if __name__ == "__main__":
    main()
