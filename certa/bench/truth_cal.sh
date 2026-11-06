#!/bin/bash
# truth_cal.sh: re-paired certified mates (--iv-pair all), the MAPQ fitted at
# the locus (--mapq cal) and CERTA on T2T-CHM13 as the second opinion, on
# HG002 tests A and B (bench/VALIDATION_RESULTS.md, "Re-paired mates and a
# fitted MAPQ"). Aligners on cores 0-63, the rest on 64-127.
#   truth_cal.sh BUILD_DIR T2T_DIR HS1_CIDX CHAIN OUT_DIR
# T2T_DIR is truth_t2t.sh's OUT_DIR (every run's placements, X/*.npz);
# HS1_CIDX is `certa index hs1.fa -q 22 -s 1`; CHAIN is UCSC's
# hs1ToHg38.over.chain.gz. Steps per test:
# 1. --iv-pair all --mapq-limit 10 --feature-tags: features for the fit;
# 2. the fit on both tests, which must reproduce include/certa/mapq_table.h;
# 3. --iv-pair all --mapq cal, and CERTA on hs1, as placements in T2T_DIR;
# 4. the tables, with BWA-MEM2's and CERTA's T2T placements as second opinion.
set -euo pipefail
BUILD=$1; TD=$(cd "$2" && pwd); HS1_CIDX=$3; CHAIN=$4; O=$(mkdir -p "$5" && cd "$5" && pwd)
C=$BUILD/certa; T=$(cd "$(dirname "$0")" && pwd)
W=~/certa-truth; FAI=~/certa-bench/idx/ref.fa.fai; IDX=~/certa-bench/idx/grch38.s1v3.cidx
export PATH=~/tools/env/truth/bin:$PATH
PIN="taskset -c 0-63"; AUX="taskset -c 64-127"
TV="/usr/bin/time -f %e\t%U\t%S\t%M"
IV="--intervals --iv-mapq cert --iv-pair all -t 64 -o /dev/stdout"
for X in A B; do
  R="$W/$X/${X}_1.fq $W/$X/${X}_2.fq"
  echo "== $(date +%T) test $X: features"
  [ -s $O/$X.feat.npz ] || $PIN $C map $IDX $R $IV --mapq-limit 10 --feature-tags 2> $O/$X.feat.log \
    | $AUX python3 $T/truth_mapq_fit.py feat $W/$X/truth.npz $FAI $O/$X.feat.npz
  echo "== $(date +%T) test $X: --mapq cal"
  [ -s $TD/$X/ivcal.npz ] || $PIN $TV -o $O/$X.cal.time $C map $IDX $R $IV --mapq cal 2> $O/$X.cal.log \
    | $AUX python3 $T/truth_t2t.py place $W/$X/truth.npz $FAI $TD/$X/ivcal.npz
  echo "== $(date +%T) test $X: CERTA on hs1"
  [ -s $TD/$X/certa_t2t.npz ] || $PIN $TV -o $O/$X.hs1.time $C map $HS1_CIDX $R $IV 2> $O/$X.hs1.log \
    | $AUX python3 $T/truth_t2t.py place --chain $CHAIN $W/$X/truth.npz $FAI $TD/$X/certa_t2t.npz
done
echo "== $(date +%T) fit"
$AUX python3 $T/truth_mapq_fit.py fit $O/A.feat.npz $O/B.feat.npz | tail -n +2 > $O/mapq_table.fit.h
tail -n +2 $T/../include/certa/mapq_table.h | cmp - $O/mapq_table.fit.h && echo "the fit reproduces include/certa/mapq_table.h"
for X in A B; do
  echo "== $(date +%T) test $X: tables"
  {
    echo "## Test $X: BWA-MEM2 on T2T-CHM13 as the second opinion"; echo
    $AUX python3 $T/truth_t2t.py eval $W/$X/truth.npz $TD/$X
    echo; echo "## Test $X: CERTA on T2T-CHM13 as the second opinion"; echo
    $AUX python3 $T/truth_t2t.py eval --t2t certa_t2t $W/$X/truth.npz $TD/$X
  } > $O/cal_$X.md
done
echo "== $(date +%T) finished"
