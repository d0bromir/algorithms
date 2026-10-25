#!/bin/bash
# Read lengths and N patterns of a FASTQ: the length histogram (lengths
# holding >= 0.01 % of reads), reads with any N, a trailing N run, all N,
# and N count classes; also how many reads of length 35 end in N (BCL
# Convert pads reads trimmed below 35 bp with N up to 35).
#   len_n_scan.sh FASTQ...
for f in "$@"; do
  awk -v F="$(basename "$f")" 'NR % 4 == 2 {
      n++; L = length($0); len[L]++
      c = gsub(/N/, "N"); if (c > 0) anyn++
      if (c == L) alln++
      if ($0 ~ /N+$/) { tail++; if (L == 35) tail35++ }
      cls = c == 0 ? "0" : c == 1 ? "1" : c <= 5 ? "2-5" : c <= 15 ? "6-15" : ">15"; nc[cls]++
    }
    END {
      printf "== %s: %d reads\n", F, n
      printf "any N %d (%.3f %%), trailing N run %d (%.3f %%), all N %d (%.4f %%), length 35 ending in N %d\n",
        anyn, 100 * anyn / n, tail, 100 * tail / n, alln, 100 * alln / n, tail35
      printf "N per read:"; split("0 1 2-5 6-15 >15", K, " "); for (i = 1; i <= 5; i++) printf " %s:%d", K[i], nc[K[i]]; printf "\n"
      printf "lengths (>= 0.01 %% of reads):"
      for (L in len) if (len[L] >= n / 10000) printf " %d:%d", L, len[L]
      printf "\n"
    }' "$f"
done
