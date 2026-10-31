#!/usr/bin/env python3
"""Truth placements on GRCh38 for reads that come from the HG002 Q100
haplotypes (bench/VALIDATION_PLAN.md, tests A and B).

A read's origin on a haplotype contig is projected to GRCh38 through GIAB's
dipcall alignment of that haplotype (the alignment the v5.0q benchmark was
built from). The read has a GRCh38 truth when both of its ends fall in one
alignment record, outside query stretches two records cover, and the
projected span is within 50 bases of the read's; otherwise it has no linear
truth (an SV, sequence missing from GRCh38, or a break in the alignment).
An end inside an insertion (haplotype bases absent from GRCh38) moves to the
nearest aligned base inward, at most 10 bases.

  truth_lift.py blocks HAP1.bam HAP2.bam REF.fai blocks.npz [--samtools samtools]
  truth_lift.py sim SIM_1.fq REF.fai blocks.npz STRAT_DIR BENCH_BED truth.npz [--length 150]
  truth_lift.py check SIM_1.fq SIM_2.fq truth.npz REF.fa OUT_PREFIX [--n 100000]

`blocks` reads the dipcall BAMs into aligned blocks. `sim` writes the truth
of every pair tools/simulate.py wrote (names r<i>_<contig>_<pos1>_<s1>_<pos2>_<s2>_<edits>;
the file's k-th pair is renamed a<k> before mapping). Truth arrays have one
row per pair and one column per mate: contig (index in REF.fai, -1 none),
1-based leftmost position, strand (1 reverse), strata bits (STRATA order),
and class (0 truth, 1 no linear truth, 2 excluded: chrM/chrEBV, 3 not
labelled).
"""
import argparse
import gzip
import re
import subprocess
import sys

import numpy as np

SNAP = 10        # an end inside an insertion moves at most this far inward
SPAN_SLACK = 50  # projected span against read span (the --bio gap cap)
PAD = 1 << 20    # gap between contigs in global coordinates
CLS_TRUTH, CLS_NOLIN, CLS_EXCL, CLS_UNLAB = 0, 1, 2, 3

# GRCh38 pseudo-autosomal regions, 0-based half-open. Ensembl's chrY carries
# the X sequence there and GIAB's analysis set masks it, so placements in a
# Y PAR are moved to X before scoring.
PAR_Y = ((10000, 2781479, 10000), (56887902, 57217415, 155701382))  # (y0, y1, x0)
PAR_X = ((10000, 2781479), (155701382, 156030895))

# Strata: (name, file under STRAT_DIR, mode). "in": the read's span lies in
# one interval; "ov": it overlaps one. "bench" is the v5.0q benchmark BED.
STRATA = (
    ("bench", None, "in"),
    ("easy", "Union/GRCh38_notinalldifficultregions.bed.gz", "in"),
    ("difficult", "Union/GRCh38_alldifficultregions.bed.gz", "ov"),
    ("lowmap_segdup", "Union/GRCh38_alllowmapandsegdupregions.bed.gz", "ov"),
    ("segdup", "SegmentalDuplications/GRCh38_segdups.bed.gz", "ov"),
    ("lowmap", "Mappability/GRCh38_lowmappabilityall.bed.gz", "ov"),
    ("tandem", "LowComplexity/GRCh38_AllTandemRepeatsandHomopolymers_slop5.bed.gz", "ov"),
    ("mhc", "OtherDifficult/GRCh38_MHC.bed.gz", "ov"),
    ("cmrg_falsedup", "FunctionalTechnicallyDifficult/GRCh38_CMRGv1.00_falselyduplicatedgenes.bed.gz", "ov"),
    ("par", "PAR", "ov"),
)

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def ensembl(name):
    """GIAB/UCSC contig name to Ensembl's (galaxy's GRCh38)."""
    if name == "chrM":
        return "MT"
    if name.startswith("chrUn_"):
        return name[6:].replace("v", ".")
    if name.endswith("_random"):
        return name.split("_")[1].replace("v", ".")
    return name[3:] if name.startswith("chr") else name


class Ref:
    """GRCh38 contigs from a FASTA index, and global coordinates."""

    def __init__(self, fai):
        self.names, lens = [], []
        with open(fai) as f:
            for line in f:
                t = line.split("\t")
                self.names.append(t[0])
                lens.append(int(t[1]))
        self.lengths = np.array(lens, dtype=np.int64)
        self.offsets = np.concatenate(([0], np.cumsum(self.lengths + PAD)[:-1])).astype(np.int64)
        self.index = {n: i for i, n in enumerate(self.names)}
        self.x, self.y = self.index.get("X", -1), self.index.get("Y", -1)

    def par_to_x(self, contig, pos):
        """Moves 0-based positions in chrY's PARs to chrX (arrays, in place)."""
        for y0, y1, x0 in PAR_Y:
            m = (contig == self.y) & (pos >= y0) & (pos < y1)
            contig[m] = self.x
            pos[m] += x0 - y0


def cmd_blocks(a):
    ref = Ref(a.fai)
    hname, hlen = [], []
    cols = {k: [] for k in ("h", "q0", "q1", "r0", "dirn", "rec", "rc")}
    ext = []  # record extents on the haplotype contig: (contig, start, end)
    nrec = 0
    for bam in (a.hap1, a.hap2):
        p = subprocess.Popen(f"{a.samtools} view {bam} | cut -f1-6", shell=True, stdout=subprocess.PIPE, text=True)
        for line in p.stdout:
            name, flag, rname, pos, _mapq, cig = line.rstrip("\n").split("\t")
            flag = int(flag)
            if flag & (4 | 256):
                continue
            c = ref.index.get(ensembl(rname))
            if c is None:
                continue
            ops = [(int(n), o) for n, o in CIGAR.findall(cig)]
            qlen = sum(n for n, o in ops if o in "MIS=XH")
            if name not in hname:
                hname.append(name)
                hlen.append(qlen)
            h = hname.index(name)
            rev = bool(flag & 16)
            qa, ra = 0, int(pos) - 1
            first = last = None
            run = None  # the current M/=/X run: [qa, ra, n]
            blocks = []
            for n, o in ops:
                if o in "M=X":
                    if run and run[0] + run[2] == qa and run[1] + run[2] == ra:
                        run[2] += n
                    else:
                        run = [qa, ra, n]
                        blocks.append(run)
                    first = qa if first is None else first
                    last = qa + n
                    qa += n
                    ra += n
                elif o in "ISH":
                    qa += n
                elif o in "DN":
                    ra += n
            for bq, br, n in blocks:
                # Haplotype coordinates; reverse records align the reverse
                # complement, so their query offsets count from the contig's end.
                s = qlen - bq - n if rev else bq
                cols["h"].append(h)
                cols["q0"].append(s)
                cols["q1"].append(s + n)
                cols["r0"].append(br + n - 1 if rev else br)  # GRCh38 position of base s
                cols["dirn"].append(-1 if rev else 1)
                cols["rec"].append(nrec)
                cols["rc"].append(c)
            if first is not None:
                ext.append((h, qlen - last, qlen - first) if rev else (h, first, last))
            nrec += 1
        if p.wait():
            sys.exit(f"samtools failed on {bam}")
    hlen = np.array(hlen, dtype=np.int64)
    hoff = np.concatenate(([0], np.cumsum(hlen + PAD)[:-1])).astype(np.int64)
    h = np.array(cols["h"], dtype=np.int64)
    q0 = hoff[h] + np.array(cols["q0"], dtype=np.int64)
    order = np.argsort(q0, kind="stable")
    out = {k: np.array(v)[order] for k, v in cols.items() if k != "h"}
    out["q0"] = q0[order]
    out["q1"] = (hoff[h] + np.array(cols["q1"], dtype=np.int64))[order]
    # Query stretches two or more records cover: no 1:1 truth there.
    ev = []
    for hh, lo, hi in ext:
        ev.append((hoff[hh] + lo, 1))
        ev.append((hoff[hh] + hi, -1))
    ev.sort()
    ms, me, depth = [], [], 0
    for x, d in ev:
        if depth < 2 <= depth + d:
            ms.append(x)
        elif depth >= 2 > depth + d:
            me.append(x)
        depth += d
    print(f"[blocks] {len(hname)} haplotype contigs, {nrec} records, {len(q0)} blocks, "
          f"{sum(e - s for s, e in zip(ms, me))} bases covered twice", file=sys.stderr)
    np.savez(a.out, hname=np.array(hname), hlen=hlen, hoff=hoff, q0=out["q0"], q1=out["q1"],
             r0=out["r0"].astype(np.int64), dirn=out["dirn"].astype(np.int8), rec=out["rec"].astype(np.int32),
             rc=out["rc"].astype(np.int16), ms=np.array(ms, dtype=np.int64), me=np.array(me, dtype=np.int64))


class Lift:
    """Projects haplotype intervals to GRCh38 through the blocks."""

    def __init__(self, path):
        d = np.load(path)
        self.hname = [str(x) for x in d["hname"]]
        self.hidx = {n: i for i, n in enumerate(self.hname)}
        self.hoff = d["hoff"]
        for k in ("q0", "q1", "r0", "dirn", "rec", "rc", "ms", "me"):
            setattr(self, k, d[k])

    def _point(self, x, start):
        n = len(self.q0)
        b = np.searchsorted(self.q0, x, side="right") - 1
        bc = np.maximum(b, 0)
        inside = (b >= 0) & (x < self.q1[bc])
        # A snapped end may land in another record than the other end; project()
        # then finds no truth.
        if start:  # snap forward to the next block's first base
            nb = np.minimum(b + 1, n - 1)
            snap = ~inside & (b + 1 < n) & (self.q0[nb] - x <= SNAP) & (self.q0[nb] > x)
            xs = np.where(snap, self.q0[nb], x)
            bs = np.where(snap, nb, bc)
        else:  # snap back to the previous block's last base
            snap = ~inside & (b >= 0) & (x - (self.q1[bc] - 1) <= SNAP)
            xs = np.where(snap, self.q1[bc] - 1, x)
            bs = bc
        r = self.r0[bs] + (xs - self.q0[bs]) * self.dirn[bs]
        return inside | snap, bs, r

    def _multi(self, x):
        if not len(self.ms):
            return np.zeros(len(x), bool)
        i = np.searchsorted(self.ms, x, side="right") - 1
        return (i >= 0) & (x < self.me[np.maximum(i, 0)])

    def project(self, h, s, e, strand):
        """Haplotype contig index h, 0-based [s, e), strand (1 reverse), as
        arrays. Returns (valid, GRCh38 contig, 0-based leftmost, strand)."""
        xs, xe = self.hoff[h] + s, self.hoff[h] + e - 1
        g1, b1, r1 = self._point(xs, True)
        g2, b2, r2 = self._point(xe, False)
        span = np.abs(r2 - r1) + 1
        valid = g1 & g2 & (self.rec[b1] == self.rec[b2]) & ~self._multi(xs) & ~self._multi(xe)
        valid &= np.abs(span - (e - s)) <= SPAN_SLACK
        left = np.minimum(r1, r2)
        rstrand = strand ^ (self.dirn[b1] < 0)
        return valid, self.rc[b1].astype(np.int16), left, rstrand.astype(np.int8)


class Strata:
    """GIAB stratifications on GRCh38 (Ensembl names, global coordinates)."""

    def __init__(self, ref, strat_dir, bench_bed):
        self.ref = ref
        self.sets = []
        for name, path, mode in STRATA:
            if path is None:
                iv = self._bed(bench_bed)
            elif path == "PAR":
                iv = [(ref.x, s, e) for s, e in PAR_X if ref.x >= 0]
            else:
                iv = self._bed(f"{strat_dir}/{path}")
            g = np.array(sorted((ref.offsets[c] + s, ref.offsets[c] + e) for c, s, e in iv), dtype=np.int64).reshape(-1, 2)
            # Merge, so the intervals are disjoint and sorted by both ends.
            merged = []
            for s, e in g:
                if merged and s <= merged[-1][1]:
                    merged[-1][1] = max(merged[-1][1], e)
                else:
                    merged.append([s, e])
            m = np.array(merged, dtype=np.int64).reshape(-1, 2)
            self.sets.append((name, mode, m[:, 0], m[:, 1]))

    def _bed(self, path):
        out = []
        with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as f:
            for line in f:
                if line.startswith(("#", "track", "browser")):
                    continue
                t = line.split("\t", 3)
                c = self.ref.index.get(ensembl(t[0]))
                if c is not None:
                    out.append((c, int(t[1]), int(t[2])))
        return out

    def bits(self, contig, left, span):
        """Strata bits of 0-based [left, left + span) on contig (arrays)."""
        out = np.zeros(len(contig), dtype=np.uint16)
        ok = contig >= 0
        qs = np.where(ok, self.ref.offsets[np.maximum(contig, 0)] + left, -1)
        qe = qs + span
        for bit, (_name, mode, S, E) in enumerate(self.sets):
            if not len(S):
                continue
            if mode == "ov":
                i = np.searchsorted(E, qs, side="right")
                hit = (i < len(S)) & (S[np.minimum(i, len(S) - 1)] < qe)
            else:
                j = np.searchsorted(S, qs, side="right") - 1
                hit = (j >= 0) & (E[np.maximum(j, 0)] >= qe)
            out |= (hit & ok).astype(np.uint16) << bit
        return out


def truth_arrays(n):
    return dict(tc=np.full((n, 2), -1, np.int16), tp=np.zeros((n, 2), np.int32), ts=np.zeros((n, 2), np.int8),
                st=np.zeros((n, 2), np.uint16), cls=np.full((n, 2), CLS_UNLAB, np.int8))


def fill_truth(t, mate, rows, lift, strata, ref, hname, s, e, strand):
    """Projects haplotype placements (hname: contig names) into truth rows."""
    h = np.array([lift.hidx.get(x, -1) for x in hname], dtype=np.int64)
    known = h >= 0
    t["cls"][rows[~known], mate] = CLS_EXCL
    rows, s, e, strand, h = rows[known], s[known], e[known], strand[known], h[known]
    valid, rc, left, rs = lift.project(h, s, e, strand)
    rc, left = rc.astype(np.int64), left.astype(np.int64)
    ref.par_to_x(rc, left)
    t["cls"][rows[~valid], mate] = CLS_NOLIN
    v = rows[valid]
    t["cls"][v, mate] = CLS_TRUTH
    t["tc"][v, mate] = rc[valid]
    t["tp"][v, mate] = left[valid] + 1
    t["ts"][v, mate] = rs[valid]
    t["st"][v, mate] = strata.bits(rc[valid], left[valid], (e - s)[valid])


def cmd_sim(a):
    ref, lift = Ref(a.fai), Lift(a.blocks)
    strata = Strata(ref, a.strat, a.bench)
    names, pos, strand = [[], []], [[], []], [[], []]
    with open(a.fq) as f:
        for k, line in enumerate(f):
            if k % 4:
                continue
            fl = line[1:].rstrip("\n").split("_")  # r<i>_<contig>_<pos1>_<s1>_<pos2>_<s2>_<edits>
            contig = "_".join(fl[1:-5])
            for m, (p, st) in enumerate(((fl[-5], fl[-4]), (fl[-3], fl[-2]))):
                names[m].append(contig)
                pos[m].append(int(p) - 1)
                strand[m].append(st == "-")
    n = len(names[0])
    t = truth_arrays(n)
    rows = np.arange(n)
    for m in range(2):
        s = np.array(pos[m], dtype=np.int64)
        fill_truth(t, m, rows, lift, strata, ref, names[m], s, s + a.length,
                   np.array(strand[m], dtype=np.int8))
    np.savez(a.out, names=np.array([x for x, _, _ in STRATA]), ref=np.array(ref.names), **t)
    summary(t, file=sys.stderr)


COMP = str.maketrans("ACGTN", "TGCAN")


def edit_banded(a, b, band=12):
    """Edit distance of equal-length strings within a diagonal band."""
    n = len(a)
    inf = n + band
    prev = [j if j <= band else inf for j in range(n + 1)]
    for i in range(1, n + 1):
        cur = [inf] * (n + 1)
        if i <= band:
            cur[0] = i
        for j in range(max(1, i - band), min(n, i + band) + 1):
            cur[j] = min(prev[j - 1] + (a[i - 1] != b[j - 1]), prev[j] + 1, cur[j - 1] + 1)
        prev = cur
    return prev[n]


def cmd_check(a):
    """Projection self-check: error-free simulated mates against GRCh38 at
    their projected truth (sample variants and indels allow a few edits)."""
    t = dict(np.load(a.truth))  # an NpzFile reads the array again on every access
    refnames = [str(x) for x in t["ref"]]
    reads, regions = [], []
    for m, path in enumerate((a.fq1, a.fq2)):
        with open(path) as f:
            for k, line in enumerate(f):
                pair = k // 4
                if pair >= a.n:
                    break
                if k % 4 == 0:
                    edits = int(line.rstrip("\n").rsplit("_", 1)[1])
                elif k % 4 == 1 and edits == 0 and t["cls"][pair, m] == CLS_TRUTH:
                    seq = line.strip()
                    p = int(t["tp"][pair, m])
                    if t["ts"][pair, m]:
                        seq = seq.translate(COMP)[::-1]
                    reads.append((pair, m, seq))
                    regions.append(f"{refnames[t['tc'][pair, m]]}:{p}-{p + len(seq) - 1}")
    with open(a.out + ".regions", "w") as f:
        f.write("\n".join(regions) + "\n")
    out = subprocess.run(f"{a.samtools} faidx {a.ref} -r {a.out}.regions", shell=True, check=True,
                         capture_output=True, text=True).stdout
    seqs, cur = [], []
    for line in out.splitlines():
        if line.startswith(">"):
            if cur or seqs:
                seqs.append("".join(cur))
            cur = []
        else:
            cur.append(line.strip().upper())
    seqs.append("".join(cur))
    hist = {}
    worst = []
    for (pair, m, seq), r in zip(reads, seqs):
        d = sum(x != y for x, y in zip(seq, r))
        if d > 3:
            d = min(d, edit_banded(seq, r))
        hist[min(d, 11)] = hist.get(min(d, 11), 0) + 1
        if d > 10:
            worst.append((pair, m, d))
    n = len(reads)
    print(f"projection check: {n} error-free mates with a truth")
    for d in sorted(hist):
        label = f"{d:2d}" if d < 11 else ">10"
        print(f"  edits {label}: {hist[d]:8d} ({100.0 * hist[d] / n:.3f}%)")
    for pair, m, d in worst[:10]:
        print(f"  far: pair {pair} mate {m + 1}: {d} edits at {refnames[t['tc'][pair, m]]}:{t['tp'][pair, m]}")


def summary(t, file=sys.stdout):
    n = t["cls"].size
    for c, label in ((CLS_TRUTH, "truth"), (CLS_NOLIN, "no linear truth"), (CLS_EXCL, "excluded"),
                     (CLS_UNLAB, "not labelled")):
        k = int((t["cls"] == c).sum())
        print(f"  {label:16s} {k:10d} mates ({100.0 * k / max(n, 1):.3f}%)", file=file)
    tr = t["cls"] == CLS_TRUTH
    for bit, (name, _, _) in enumerate(STRATA):
        k = int(((t["st"] >> bit) & 1)[tr].sum())
        print(f"  stratum {name:14s} {k:10d} mates ({100.0 * k / max(int(tr.sum()), 1):.3f}% of truth)", file=file)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("blocks")
    b.add_argument("hap1")
    b.add_argument("hap2")
    b.add_argument("fai")
    b.add_argument("out")
    b.add_argument("--samtools", default="samtools")
    s = sub.add_parser("sim")
    s.add_argument("fq")
    s.add_argument("fai")
    s.add_argument("blocks")
    s.add_argument("strat")
    s.add_argument("bench")
    s.add_argument("out")
    s.add_argument("--length", type=int, default=150)
    c = sub.add_parser("check")
    c.add_argument("fq1")
    c.add_argument("fq2")
    c.add_argument("truth")
    c.add_argument("ref")
    c.add_argument("out")
    c.add_argument("--n", type=int, default=100000)
    c.add_argument("--samtools", default="samtools")
    a = ap.parse_args()
    {"blocks": cmd_blocks, "sim": cmd_sim, "check": cmd_check}[a.cmd](a)


if __name__ == "__main__":
    main()
