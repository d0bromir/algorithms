#!/bin/bash
# Exact-duplicate read sequences in a FASTQ (single-end) and pair-level
# (R1 + R2 sequence identical), with the reverse complement counted as the
# same read (a mapper's search is strand-symmetric): distinct sequences,
# reads that repeat an earlier one, and the duplicate-group size histogram.
#   dup_scan.sh R1.fq R2.fq OUT_DIR
set -eu
R1=$1; R2=$2; O=$3
mkdir -p "$O"
S="sort -S 24G --parallel=16 -T $O"
report() {  # label: reads a "count" per distinct sequence on stdin
  awk -v L="$1" '{n+=$1; d++; if ($1>1) {dup+=$1-1; g[$1>=10?10:$1]++}}
    END {printf "%s: %d reads, %d distinct, %d (%.2f %%) repeat an earlier read\n", L, n, d, dup, 100*dup/n;
         for (k=2;k<=10;k++) printf "  groups of %s%d: %d\n", (k==10?">=":""), k, g[k]}'
}
paste <(awk 'NR % 4 == 2' "$R1") <(awk 'NR % 4 == 2' "$R1" | rev | tr ACGT TGCA) | awk '{print ($1 < $2 ? $1 : $2)}' |
  $S | uniq -c | report "single-end (R1, strand-canonical)"
awk 'NR % 4 == 2' "$R1" | $S | uniq -c | report "single-end (R1, as read)"
paste <(awk 'NR % 4 == 2' "$R1") <(awk 'NR % 4 == 2' "$R2") | $S | uniq -c | report "pairs (R1 + R2 as read)"
rm -rf "$O"/sort*
