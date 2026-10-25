"""Fault injection for certa check (Aim 4).

  aim4_tamper.py in.sam out.sam faults.tsv [rate] [seed]

Copies a SAM file of `certa map --intervals --witness`, changing one claim
in a random fraction `rate` (default 0.01) of eligible records (or, for the
pair faults, of eligible pairs: both records), and lists each changed
record in faults.tsv (record number, fault). Faults that make a claim false
by construction (certa check must reject every one):
  as+1        AS raised by 1 (the CIGAR no longer scores AS)
  clip        the last aligned base clipped, AS kept (scores AS - 6 at most)
  xu<xl       XU set to XL - 1
  xq+10       XQ raised by 10 (above MAPQ or the formula at AS, XH)
  nowitness   XW removed from a record with XU < read length
  pr-xq+10    a PR mate's XQ raised by 10 (above MAPQ or the pairing rule)
  nopairwitness  XY removed from both mates of a PR pair; sure for a mate
              whose XQ rests on XS: above q_se, BWA-MEM's formula at (AS, XH)
              (0 without XH), when it has no XV (unsure otherwise:
              nopairwitness?)
  dist-1      D in both mates' XY set one below the pair's distance (the
              pair is no longer proper); sure as above (dist-1?)
Faults whose effect on truth is unknown (the checker may accept a claim that
stays true and provable, or reject one it can no longer prove):
  pos+1       POS moved by one base, AS kept
  xh-10       XH lowered by 10, XQ kept
  probe+1     one witness probe moved by one base
  xs-10       XS lowered by 10 in both mates of a PR pair, XQ kept
  mate:F      an untouched PR mate of a record with fault F: its XQ rests on
              the pair, so it may no longer be provable either
"""
import math
import random
import re
import sys

SURE = ["as+1", "clip", "xu<xl", "xq+10", "nowitness", "pr-xq+10"]
UNSURE = ["pos+1", "xh-10", "probe+1"]
PAIR_SURE = ["nopairwitness", "dist-1"]
PAIR_UNSURE = ["xs-10"]


def tag(fields, name):
    for i in range(11, len(fields)):
        if fields[i].startswith(name + ":"):
            return i
    return -1


def eligible(f, fault):
    mapped = not int(f[1]) & 4
    xu, xl, xq, xh, xw = (tag(f, t) for t in ("XU", "XL", "XQ", "XH", "XW"))
    L = len(f[9])
    if fault in ("as+1", "pos+1"):
        return mapped
    if fault == "clip":
        return mapped and re.fullmatch(r"(\d+S)?\d+M", f[5]) is not None and int(re.findall(r"(\d+)M", f[5])[0]) > 1
    if fault == "xu<xl":
        return xu >= 0 and xl >= 0
    if fault in ("xq+10", "pr-xq+10"):
        pr = tag(f, "XT") >= 0 and f[tag(f, "XT")] == "XT:Z:PR"
        return xq >= 0 and 0 < int(f[xq][5:]) <= 50 and pr == (fault == "pr-xq+10")
    if fault == "nowitness":
        return xw >= 0 and xu >= 0 and int(f[xu][5:]) < L
    if fault == "xh-10":
        return xh >= 0
    if fault == "probe+1":
        return xw >= 0 and re.search(r"\d+-\d+", f[xw]) is not None
    return False


def apply(f, fault):
    if fault == "as+1":
        i = tag(f, "AS")
        f[i] = "AS:i:%d" % (int(f[i][5:]) + 1)
    elif fault == "pos+1":
        f[3] = str(int(f[3]) + 1)
    elif fault == "clip":
        m = re.fullmatch(r"(\d+S)?(\d+)M", f[5])
        f[5] = (m.group(1) or "") + "%dM1S" % (int(m.group(2)) - 1)
    elif fault == "xu<xl":
        i, j = tag(f, "XU"), tag(f, "XL")
        f[i] = "XU:i:%d" % (int(f[j][5:]) - 1)
    elif fault in ("xq+10", "pr-xq+10"):
        i = tag(f, "XQ")
        f[i] = "XQ:i:%d" % (int(f[i][5:]) + 10)
    elif fault == "nowitness":
        del f[tag(f, "XW")]
    elif fault == "xh-10":
        i = tag(f, "XH")
        f[i] = "XH:i:%d" % (int(f[i][5:]) - 10)
    elif fault == "probe+1":
        i = tag(f, "XW")
        probes = list(re.finditer(r"(\d+)-(\d+)", f[i]))
        m = random.choice(probes)
        a, b = int(m.group(1)) + 1, int(m.group(2)) + 1
        f[i] = f[i][: m.start()] + "%d-%d" % (a, b) + f[i][m.end():]


def mapq_se(score, sub, L):
    """BWA-MEM's single-end MAPQ formula, as certa check computes it."""
    sub = max(sub, 19)
    if score <= 0 or sub >= score:
        return 0
    identity = 1.0 - (L - score) / 5.0 / L
    t = 1.0 if L < 50 else math.log(50.0) / math.log(L)
    t *= identity * identity
    return max(0, min(60, int(6.02 * (score - sub) * t * t + 0.499)))


def rests_on_xs(f):
    """True if this PR record's XQ exceeds its single-end term q_se (so it
    rests on XS); None if unknown (XV: its score needs the reference)."""
    xq = int(f[tag(f, "XQ")][5:])
    if xq <= 0:
        return False
    if tag(f, "XV") >= 0:
        return None
    h = tag(f, "XH")
    qse = mapq_se(int(f[tag(f, "AS")][5:]), int(f[h][5:]), len(f[9])) if h >= 0 else 0
    return xq > qse


def pair_eligible(a, b, fault):
    """A PR pair: both records tier PR with XS and XY, and some mate whose
    XQ rests on XS (for xs-10: some mate with XQ > 0)."""
    for f in (a, b):
        if tag(f, "XT") < 0 or f[tag(f, "XT")] != "XT:Z:PR" or tag(f, "XS") < 0 or tag(f, "XY") < 0:
            return False
    if fault == "xs-10":
        return any(int(f[tag(f, "XQ")][5:]) > 0 for f in (a, b))
    return any(rests_on_xs(f) for f in (a, b))


def pair_apply(a, b, fault):
    if fault == "nopairwitness":
        for f in (a, b):
            del f[tag(f, "XY")]
    elif fault == "dist-1":
        dist = abs(int(a[3]) - int(b[3]))  # the reverse mate starts after the forward one
        for f in (a, b):
            i = tag(f, "XY")
            rest = f[i][5:].split(";", 1)[1]
            f[i] = "XY:Z:%d;%s" % (dist - 1, rest)
    elif fault == "xs-10":
        for f in (a, b):
            i = tag(f, "XS")
            f[i] = "XS:i:%d" % (int(f[i][5:]) - 10)


def main():
    src, dst, log = sys.argv[1:4]
    rate = float(sys.argv[4]) if len(sys.argv) > 4 else 0.01
    random.seed(int(sys.argv[5]) if len(sys.argv) > 5 else 1)
    n = 0
    with open(src) as fi, open(dst, "w") as fo, open(log, "w") as fl:
        fl.write("record\tfault\tsure\n")
        pending = None  # (record number, fields) of a record whose mate may follow

        def single(num, f, write=True):
            fault = None
            if random.random() < rate:
                faults = [x for x in SURE + UNSURE if eligible(f, x)]
                if faults:
                    fault = random.choice(faults)
                    apply(f, fault)
                    fl.write("%d\t%s\t%d\n" % (num, fault, fault in SURE))
            if write:
                fo.write("\t".join(f) + "\n")
            return fault

        def is_pr(f):
            return tag(f, "XT") >= 0 and f[tag(f, "XT")] == "XT:Z:PR"

        for line in fi:
            if line[0] == "@":
                fo.write(line)
                continue
            n += 1
            f = line.rstrip("\n").split("\t")
            if pending and pending[1][0] == f[0]:  # a pair: both records
                num0, f0 = pending
                pending = None
                faults = [x for x in PAIR_SURE + PAIR_UNSURE if pair_eligible(f0, f, x)]
                if faults and random.random() < rate:
                    fault = random.choice(faults)
                    rest = [rests_on_xs(f0), rests_on_xs(f)]
                    pair_apply(f0, f, fault)
                    for num, r in zip((num0, n), rest):
                        sure = fault in PAIR_SURE and r is True
                        name = fault if (sure or fault not in PAIR_SURE) else fault + "?"
                        fl.write("%d\t%s\t%d\n" % (num, name, sure))
                    fo.write("\t".join(f0) + "\n")
                    fo.write("\t".join(f) + "\n")
                else:
                    g0 = single(num0, f0, False)
                    g1 = single(n, f, False)
                    if g0 and not g1 and is_pr(f):
                        fl.write("%d\tmate:%s\t0\n" % (n, g0))
                    if g1 and not g0 and is_pr(f0):
                        fl.write("%d\tmate:%s\t0\n" % (num0, g1))
                    fo.write("\t".join(f0) + "\n")
                    fo.write("\t".join(f) + "\n")
                continue
            if pending:
                single(*pending)
            pending = (n, f)
        if pending:
            single(*pending)


if __name__ == "__main__":
    main()
