#!/usr/bin/env python3
"""Compare primary alignments of several mappers on the same single-end reads.

  compare_sam.py REFERENCE_TOOL=ref.sam TOOL=x.sam [TOOL=y.sam ...] --certa certa.sam

For every tool: fraction of reads mapped, and agreement with the reference
tool (same contig and unclipped start within --window bp), overall and for
reads the reference tool maps with MAPQ >= 30.

For CERTA (certified reads only; uncertified reads went to the fallback):
agreement with the reference tool by CERTA tier and MAPQ, and a real-data
check of the certificate. A certified read's NM is the minimum end-to-end
edit distance over the whole reference, so the reference tool must never
report an unclipped alignment with fewer edits. Violations are counted and
the first few are printed.

With --fasta (indexed by samtools faidx), a violation whose reference window
contains a non-ACGT base is reported separately. BWA-family indexes replace
reference N/IUPAC bases with random nucleotides, so they can report a match
there; CERTA always scores them as mismatches. Only violations without such
an explanation fail the check.
"""
import argparse
import re
import sys

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def primary(path):
    """name -> (contig, unclipped_start, strand, mapq, nm, clipped)"""
    out = {}
    with open(path) as f:
        for line in f:
            if line.startswith("@"):
                continue
            t = line.rstrip("\n").split("\t", 11)  # t[11] = all optional tags
            flag = int(t[1])
            if flag & 0x900:
                continue
            name = t[0]
            if flag & 4:
                out[name] = None
                continue
            ops = CIGAR.findall(t[5])
            lead = 0
            for n, op in ops:
                if op in "SH":
                    lead += int(n)
                else:
                    break
            clipped = any(op in "SH" for _, op in ops)
            nm = -1
            for tag in t[11].split("\t") if len(t) > 11 else []:
                if tag.startswith("NM:i:"):
                    nm = int(tag[5:])
            tags = t[11] if len(t) > 11 else ""
            out[name] = (t[2], int(t[3]) - lead, "-" if flag & 16 else "+", int(t[4]), nm, clipped, tags)
    return out


class Fasta:
    """Random access to a FASTA indexed with samtools faidx."""

    def __init__(self, path):
        self.f = open(path, "rb")
        self.idx = {}
        with open(path + ".fai") as fai:
            for line in fai:
                name, length, offset, bases, width = line.split("\t")[:5]
                self.idx[name] = (int(length), int(offset), int(bases), int(width))

    def fetch(self, name, start1, n):
        """n bases starting at 1-based position start1."""
        length, offset, bases, width = self.idx[name]
        out = []
        for p in range(max(start1 - 1, 0), min(start1 - 1 + n, length)):
            self.f.seek(offset + (p // bases) * width + p % bases)
            out.append(self.f.read(1).decode())
        return "".join(out)


def agree(a, b, window):
    return a is not None and b is not None and a[0] == b[0] and abs(a[1] - b[1]) <= window


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ref", help="REFERENCE_TOOL=file.sam")
    ap.add_argument("tools", nargs="*", help="TOOL=file.sam")
    ap.add_argument("--certa", help="CERTA certified SAM")
    ap.add_argument("--reads", type=int, required=True, help="number of input reads")
    ap.add_argument("--window", type=int, default=10)
    ap.add_argument("--fasta", help="reference FASTA (with .fai) to explain violations")
    a = ap.parse_args()

    ref_name, ref_path = a.ref.split("=", 1)
    ref = primary(ref_path)
    confident = {k for k, v in ref.items() if v is not None and v[3] >= 30}
    print(f"reference: {ref_name}; {a.reads} reads; {sum(v is not None for v in ref.values())} mapped; "
          f"{len(confident)} with MAPQ >= 30")
    print(f"\n{'tool':14s} {'mapped %':>9s} {'agree %':>9s} {'agree % (ref MAPQ>=30)':>24s}")
    for spec in [a.ref] + a.tools:
        name, path = spec.split("=", 1)
        aln = ref if name == ref_name else primary(path)
        mapped = sum(v is not None for v in aln.values())
        ag = sum(agree(v, ref.get(k), a.window) for k, v in aln.items())
        agc = sum(agree(aln.get(k), ref[k], a.window) for k in confident)
        print(f"{name:14s} {100 * mapped / a.reads:9.3f} {100 * ag / a.reads:9.3f} "
              f"{100 * agc / max(len(confident), 1):24.3f}")

    if not a.certa:
        return
    cer = primary(a.certa)
    print(f"\nCERTA certified reads: {len(cer)} ({100 * len(cer) / a.reads:.2f}% of input)")
    groups = {}
    violations, checked = [], 0
    for k, v in cer.items():
        tier = re.search(r"XT:Z:(S\d)", v[6]).group(1)
        g = groups.setdefault((tier, v[3]), [0, 0, 0])
        g[0] += 1
        r = ref.get(k)
        g[1] += agree(v, r, a.window)
        g[2] += r is not None and r[3] >= 30
        if r is not None and not r[5] and r[4] >= 0:  # unclipped reference alignment
            checked += 1
            if r[4] < v[4]:
                violations.append((k, v, r))
    print(f"{'tier':5s} {'MAPQ':>4s} {'reads':>10s} {'agree with ref %':>17s} {'ref MAPQ>=30 %':>15s}")
    for (tier, mq) in sorted(groups, key=lambda x: (x[0], -x[1])):
        n, ag, rc = groups[(tier, mq)]
        print(f"{tier:5s} {mq:4d} {n:10d} {100 * ag / n:17.3f} {100 * rc / n:15.3f}")
    print(f"\ncertificate check: {checked} certified reads with an unclipped {ref_name} alignment; "
          f"{len(violations)} where {ref_name} reports fewer edits than CERTA's certified minimum")
    fasta = Fasta(a.fasta) if a.fasta else None
    unexplained = 0
    for k, v, r in violations:
        why = "UNEXPLAINED"
        if fasta:
            window = fasta.fetch(r[0], r[1], 160).upper()
            bad = sum(c not in "ACGT" for c in window)
            if bad:
                why = f"explained: {bad} non-ACGT reference base(s) in window (random base in {ref_name} index)"
        unexplained += why == "UNEXPLAINED"
        print(f"  {k}: CERTA {v[0]}:{v[1]}{v[2]} NM={v[4]}  {ref_name} {r[0]}:{r[1]}{r[2]} NM={r[4]}  -> {why}")
    print(f"unexplained violations: {unexplained}" + ("" if fasta else " (pass --fasta to classify)"))
    if unexplained:
        sys.exit(1)


if __name__ == "__main__":
    main()
