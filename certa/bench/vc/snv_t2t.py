#!/usr/bin/env python3
"""snv_t2t.py SAMTOOLS REF SITES.tsv CERTA.cram FEAT.tsv.zst...: at CERTA-only false SNV calls, where T2T-CHM13
puts the reads CERTA counts there with MAPQ >= 5: those with the false allele (the called alt base; the ref base
where the truth is homalt) against the others. FEAT: per primary record of a --t2t --feature-tags run of the same
reads (t2t_feat.sh: qname, mate, contig, pos, MAPQ, XT, XC). XC's class: T2T's placement lifts to within 150 bp of
the record (here), lifts elsewhere in GRCh38, is CHM13 sequence that does not lift (CHM13 only), or there is none.
A class from a run that placed the read more than 150 bp from CERTA.cram's record refers to another locus and is
counted apart. Reads in CHM13-only sequence are those a sink contig could take; reads T2T puts elsewhere in
GRCh38 need the placement itself to change."""
import re
import subprocess
import sys
from collections import Counter, defaultdict

st, ref, sites_path, certa = sys.argv[1:5]
feat_paths = sys.argv[5:]
OP = re.compile(r"(\d+)([MIDNSHP=X])")
CLS = {0: "none on T2T", 1: "here", 2: "elsewhere in GRCh38", 3: "CHM13 only"}


def view(cram, region):
    cmd = [st, "view", "-T", ref, "-F", "0x904", cram, region]
    return subprocess.run(cmd, capture_output=True, text=True, check=True).stdout.splitlines()


def base_at(c, pos):
    r, q = int(c[3]), 0
    for ln, op in OP.findall(c[5]):
        ln = int(ln)
        if op in "M=X":
            if r <= pos < r + ln:
                return c[9][q + pos - r]
            r += ln
            q += ln
        elif op in "IS":
            q += ln
        elif op in "DN":
            if r <= pos < r + ln:
                return None
            r += ln
    return None


def tier(rec):
    m = re.search(r"\tXT:Z:(\w+)", rec)
    return m.group(1) if m else "-"


# (qname, mate 1|2) -> [(site, allele, contig, pos, tier)]: one entry per site the read covers.
want = defaultdict(list)
sites = [l.rstrip("\n").split("\t") for l in open(sites_path) if l.strip()]
for chrom, pos, kind, rb, ab, tgt, qgt in sites:
    pos = int(pos)
    excess = rb if (kind == "geno" and tgt == "homalt") else ab.split(",")[0]
    for rec in view(certa, f"{chrom}:{pos}-{pos}"):
        c = rec.split("\t", 11)
        if int(c[4]) < 5:
            continue
        b = base_at(c, pos)
        if b is None:
            continue
        mate = 1 if int(c[1]) & 64 else 2
        want[(c[0], mate)].append((f"{chrom}:{pos}", "false" if b == excess else "other", c[2], int(c[3]), tier(rec)))

found = {}
for path in feat_paths:
    z = subprocess.Popen(["zstd", "-dc", path], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    for line in z.stdout:
        q, mate, contig, p, mq, xt, xc = line.rstrip("\n").split("\t")
        k = (q, int(mate))
        if k in want:
            if xc == "-":
                cls, tmq = -1, -1
            else:
                f = xc.split(",")
                cls, tmq = int(f[0]), int(f[2])
            found[k] = (contig, int(p), cls, tmq)
    z.wait()


def state(k, contig, pos):
    f = found.get(k)
    if f is None:
        return "not in the T2T run"
    if f[0] != contig or abs(f[1] - pos) > 150:
        return "placed elsewhere in that run"
    if f[2] < 0:
        return "no XC"
    return CLS[f[2]] + (", T2T MAPQ >= 20" if f[2] and f[3] >= 20 else (", T2T MAPQ < 20" if f[2] else ""))


table, by_tier = Counter(), Counter()
for k, lst in want.items():
    for site, allele, contig, pos, t in lst:
        s = state(k, contig, pos)
        table[(allele, s)] += 1
        if allele == "false":
            by_tier[(t, s)] += 1
n = Counter(a for (a, _), v in table.items() for _ in range(v))
print(f"{len(sites)} sites; reads with MAPQ >= 5: {n['false']} with the false allele, {n['other']} others")
for allele in ("false", "other"):
    print(f"\n{allele} allele ({n[allele]}):")
    for (a, s), v in sorted(table.items(), key=lambda x: -x[1]):
        if a == allele:
            print(f"   {s:42} {v:6}  {100 * v / max(n[allele], 1):5.1f} %")
print("\nfalse allele by CERTA's tier:")
for (t, s), v in sorted(by_tier.items(), key=lambda x: (x[0][0], -x[1])):
    print(f"   {t:3} {s:42} {v:6}")
