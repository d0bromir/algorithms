"""MAPQ calibration on simulated reads (Aim 3).

  samtools view -h aln.bam | aim3_calibration.py LABEL [TAG]

Reads named by simulate_pairs.py carry their origin. A primary record is
correct if it lies on the true contig and strand and its unclipped start is
within 30 bases of the true start. Records are binned by MAPQ (or by the
integer tag TAG, e.g. XQ); for each bin: reads, wrong placements, the
observed error rate and the rate MAPQ promises (mean of 10^(-q/10)).
"""
import re
import sys

BINS = [(0, 0), (1, 9), (10, 19), (20, 29), (30, 39), (40, 49), (50, 59), (60, 60)]
LEAD = re.compile(r"^(\d+)S")


def main():
    label = sys.argv[1]
    tag = sys.argv[2] + ":i:" if len(sys.argv) > 2 else None
    n = mapped = wrong = 0
    reads = [0] * len(BINS)
    bad = [0] * len(BINS)
    promise = [0.0] * len(BINS)
    for line in sys.stdin:
        if line[0] == "@":
            continue
        t = line.split("\t", 12)
        flag = int(t[1])
        if flag & 0x900:
            continue
        n += 1
        if flag & 4:
            continue
        mapped += 1
        f = t[0].split(":")  # sim<i>:contig:start1:strand1:start2:strand2
        k = 2 if flag & 64 or not flag & 128 else 4
        contig, start, strand = f[1], int(f[k]), f[k + 1]
        m = LEAD.match(t[5])
        pos = int(t[3]) - (int(m.group(1)) if m else 0)
        ok = t[2] == contig and ("-" if flag & 16 else "+") == strand and abs(pos - start) <= 30
        wrong += not ok
        q = int(t[4])
        if tag:
            i = line.find("\t" + tag)
            q = int(line[i + len(tag) + 1:].split("\t", 1)[0]) if i >= 0 else 0
        b = next(j for j, (lo, hi) in enumerate(BINS) if lo <= q <= hi) if q <= 60 else len(BINS) - 1
        reads[b] += 1
        bad[b] += not ok
        promise[b] += 10 ** (-q / 10)
    what = tag[:-3] if tag else "MAPQ"
    print(f"{label} ({what}): {n} reads, mapped {100 * mapped / n:.2f}%, wrong {wrong} "
          f"({100 * wrong / max(mapped, 1):.3f}% of mapped)")
    print(f"  {what:>7} {'reads':>9} {'wrong':>7} {'observed':>10} {'promised':>10}")
    hi_reads = hi_bad = 0
    for j, (lo, hi) in enumerate(BINS):
        if not reads[j]:
            continue
        name = f"{lo}" if lo == hi else f"{lo}-{hi}"
        print(f"  {name:>7} {reads[j]:9d} {bad[j]:7d} {bad[j] / reads[j]:10.2e} {promise[j] / reads[j]:10.2e}")
        if lo >= 20:
            hi_reads += reads[j]
            hi_bad += bad[j]
    print(f"  {what} >= 20: {hi_reads} reads ({100 * hi_reads / n:.2f}%), wrong {hi_bad} "
          f"({hi_bad / max(hi_reads, 1):.2e})")


if __name__ == "__main__":
    main()
