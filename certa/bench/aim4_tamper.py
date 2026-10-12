"""Fault injection for certa check (Aim 4).

  aim4_tamper.py in.sam out.sam faults.tsv [rate] [seed]

Copies a SAM file of `certa map --intervals --witness`, changing one claim
in a random fraction `rate` (default 0.01) of eligible records, and lists
each change in faults.tsv (record number, fault). Faults that make a claim
false by construction (certa check must reject every one):
  as+1        AS raised by 1 (the CIGAR no longer scores AS)
  clip        the last aligned base clipped, AS kept (scores AS - 6 at most)
  xu<xl       XU set to XL - 1
  xq+10       XQ raised by 10 (above MAPQ or the formula at AS, XH)
  nowitness   XW removed from a record with XU < read length
Faults whose effect on truth is unknown (the checker may accept a claim that
stays true and provable, or reject one it can no longer prove):
  pos+1       POS moved by one base, AS kept
  xh-10       XH lowered by 10, XQ kept
  probe+1     one witness probe moved by one base
"""
import random
import re
import sys

SURE = ["as+1", "clip", "xu<xl", "xq+10", "nowitness"]
UNSURE = ["pos+1", "xh-10", "probe+1"]


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
    if fault == "xq+10":
        return xq >= 0 and 0 < int(f[xq][5:]) <= 50 and tag(f, "XT") >= 0 and f[tag(f, "XT")] != "XT:Z:PR"
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
    elif fault == "xq+10":
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


def main():
    src, dst, log = sys.argv[1:4]
    rate = float(sys.argv[4]) if len(sys.argv) > 4 else 0.01
    random.seed(int(sys.argv[5]) if len(sys.argv) > 5 else 1)
    n = 0
    with open(src) as fi, open(dst, "w") as fo, open(log, "w") as fl:
        fl.write("record\tfault\tsure\n")
        for line in fi:
            if line[0] == "@":
                fo.write(line)
                continue
            n += 1
            if random.random() < rate:
                f = line.rstrip("\n").split("\t")
                faults = [x for x in SURE + UNSURE if eligible(f, x)]
                if faults:
                    fault = random.choice(faults)
                    apply(f, fault)
                    line = "\t".join(f) + "\n"
                    fl.write("%d\t%s\t%d\n" % (n, fault, fault in SURE))
            fo.write(line)


if __name__ == "__main__":
    main()
