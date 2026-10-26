#!/bin/bash
# Columba 2.0 (lossless, best mode, >= 96 % identity: every co-optimal
# alignment within 6 edits of a 150-bp read) against CERTA's certificates on
# the same 2 M HG002 reads / pairs, as slow_bench.sh does for Yara and
# RazerS 3 (bench/weaknesses.py, bench/within_rate.py). galaxy, 64 cores.
#   k_reliability.sh CERTA_BIN OUT_DIR
set -u
C=$1; O=$2
cd ~/certa-bench
mkdir -p "$O"
B=$(dirname "$C")/../bench
COL=$HOME/tools/src/columba/build_Vanilla/columba
IDX=idx/grch38.s1v3.cidx
R1=reads/hg002_2M.fq; R2=reads/hg002_2M_2.fq
P="taskset -c 0-63"
for m in se pe; do
  in="-f $R1"; cin=$R1
  [ $m = pe ] && { in="-f $R1 -F $R2"; cin="$R1 $R2"; }
  $P "$COL" -r idx/columba/grch38 $in -t 64 -a best -I 96 -o "$O/columba_$m.sam" -l "$O/columba_$m.log" > /dev/null
  $P "$C" map $IDX $cin -t 64 --full-cert --intervals -o "$O/certa.sam" 2> "$O/certa_$m.log"
  echo "-- $m: CERTA --full-cert --intervals against Columba (best, >= 96 %), 2 M"
  grep -E "reads: certified|pairs .*both mates" "$O/certa_$m.log"
  python3 "$B/weaknesses.py" "$O/certa.sam" Columba:all="$O/columba_$m.sam"
  python3 "$B/within_rate.py" "$O/certa.sam" "$O/columba_$m.sam"
  rm -f "$O/certa.sam" "$O/columba_$m.sam"
done
echo ALL_DONE
