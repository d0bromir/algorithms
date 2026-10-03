#!/usr/bin/env bash
# Uniform random subsample of a FASTQ(.gz) in one pass (seeded, reproducible).
#   bench/subsample.sh IN.fastq.gz TOTAL_READS WANTED_READS SEED OUT.fastq
set -euo pipefail
in="$1"; total="$2"; want="$3"; seed="$4"; out="$5"
frac=$(awk -v w="$want" -v t="$total" 'BEGIN { printf "%.10f", w / t }')
dec=cat
case "$in" in *.gz) dec="gzip -dc"; command -v pigz >/dev/null && dec="pigz -dc" ;; esac
$dec "$in" | awk -v p="$frac" -v seed="$seed" '
  BEGIN { srand(seed) }
  NR % 4 == 1 { keep = (rand() < p) }
  keep' > "$out"
echo "$(( $(wc -l < "$out") / 4 )) reads written to $out (fraction $frac, seed $seed)"
