"""Join tier R output (TSV) with a reference aligner's SAM (BWA-MEM2).

  aim2_join.py tierR.tsv aligner.sam [--fasta ref.fa]

- Validity: the aligner's score (BWA-MEM defaults, -5 per clipped end) must
  be <= hi for every read, and <= lo for reads claimed exact. Violations are
  explained if the aligner's window contains reference N (BWA replaces N with
  random bases); others are reported.
- Breakdown by the aligner's MAPQ (< 20: repeat-like; >= 20: unique-like).
"""
import collections
import re
import sys

CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def score(cig, nm):
    ops = [(int(n), o) for n, o in CIG.findall(cig)]
    m = sum(n for n, o in ops if o in "M=X")
    gaps = [n for n, o in ops if o in "ID"]
    mism = max(0, nm - sum(gaps))
    clips = sum(1 for n, o in ops if o in "SH")
    return (m - mism) - 4 * mism - sum(6 + g for g in gaps) - 5 * clips


def ref_span(cig):
    return sum(int(n) for n, o in CIG.findall(cig) if o in "MD=XN")


class Fasta:
    def __init__(self, path):
        self.f = open(path, "rb")
        self.idx = {}
        for line in open(path + ".fai"):
            name, length, offset, bases, width = line.split("\t")[:5]
            self.idx[name] = (int(length), int(offset), int(bases), int(width))

    def has_n(self, name, start1, n):
        length, offset, bases, width = self.idx[name]
        for p in range(max(start1 - 1, 0), min(start1 - 1 + n, length)):
            self.f.seek(offset + (p // bases) * width + p % bases)
            if self.f.read(1).upper() not in b"ACGT":
                return True
        return False


def main():
    tsv, sam = sys.argv[1], sys.argv[2]
    fasta = Fasta(sys.argv[sys.argv.index("--fasta") + 1]) if "--fasta" in sys.argv else None
    rows = {}
    for line in open(tsv):
        if line.startswith("read\t"):
            continue
        t = line.rstrip("\n").split("\t")
        rows[t[0]] = dict(exact=int(t[1]), lo=int(t[2]), hi=int(t[3]), copies=int(t[4]), parts=int(t[8]),
                          absent=int(t[9]))
    al = {}
    for line in open(sam):
        if line[0] == "@":
            continue
        t = line.split("\t")
        if t[0] not in rows or int(t[1]) & 0x900:
            continue
        f = int(t[1])
        if f & 4:
            al[t[0]] = None
            continue
        nm = int(re.search(r"NM:i:(\d+)", line).group(1))
        al[t[0]] = (int(t[4]), score(t[5], nm), t[2], int(t[3]), ref_span(t[5]))
    cats = collections.defaultdict(collections.Counter)
    bad = []
    for name, r in rows.items():
        a = al.get(name)
        cat = "unmapped by aligner" if a is None else ("aligner MAPQ < 20" if a[0] < 20 else "aligner MAPQ >= 20")
        c = cats[cat]
        c["n"] += 1
        if r["lo"] == -1:
            c["no lower bound"] += 1
        elif r["exact"]:
            c["exact"] += 1
            c["exact + proven tie"] += r["copies"] >= 2
        else:
            w = r["hi"] - r["lo"]
            c["interval <= 5"] += w <= 5
            c["proven ambiguous (2 loci, interval <= 5)"] += r["copies"] >= 2 and w <= 5
        if a is not None:
            bound = r["lo"] if r["exact"] else r["hi"]
            if a[1] > bound:
                why = "UNEXPLAINED"
                if fasta and fasta.has_n(a[2], a[3], a[4]):
                    why = "explained (reference N)"
                bad.append((name, r, a, why))
    order = ["n", "exact", "exact + proven tie", "interval <= 5", "proven ambiguous (2 loci, interval <= 5)",
             "no lower bound"]
    for cat in ["aligner MAPQ < 20", "aligner MAPQ >= 20", "unmapped by aligner"]:
        c = cats[cat]
        if not c["n"]:
            continue
        print(f"{cat}: {c['n']} reads")
        for k in order[1:]:
            print(f"   {k:45s} {c[k]:7d}  {100 * c[k] / c['n']:6.2f}%")
        good = c["exact"] + c["proven ambiguous (2 loci, interval <= 5)"]
        print(f"   {'exact or proven ambiguous':45s} {good:7d}  {100 * good / c['n']:6.2f}%")
    unexpl = [b for b in bad if b[3] == "UNEXPLAINED"]
    print(f"validity: aligner score above the bound for {len(bad)} reads; unexplained {len(unexpl)}")
    for b in unexpl[:10]:
        print("  ", b)


if __name__ == "__main__":
    main()
