#!/bin/bash
# Peak memory and time at full scale, the build before the memory changes
# against the build with them, interleaved: 20,001,073 HG002 reads (single-
# end, R1) or pairs, the final protocol's conditions (bench/final_bench.sh:
# the same 64 cores, SAM to /dev/null, interval mode), default and
# --full-cert, three rounds; CPU, and on galaxy the A100 too.
#   memory_bench.sh HOST OLD_CERTA NEW_CERTA OUT_DIR
set -u
host=$1; OLD=$2; NEW=$3; O=$4
cd ~/certa-bench
mkdir -p "$O"
IDX=idx/grch38.s1v3.cidx
R1=reads/hg002_20M.fq
if [ "$host" = a2 ]; then R2=/data/certa-bench/reads/hg002_20M_2.fq; PIN="numactl --interleave=all"; GPU=0
else R2=reads/hg002_20M_2.fq; PIN="taskset -c 0-63"; GPU=1; fi
TSV=$O/timings.tsv
[ -f "$TSV" ] || printf "host\tmode\tpipeline\tconfig\tbuild\trep\twall_s\tuser_s\tsys_s\tmax_rss_kb\n" > "$TSV"
run() {  # mode pipeline config build rep certa-args...
  local mode=$1 pipe=$2 cfg=$3 b=$4 rep=$5; shift 5
  local C=$OLD x=""
  [ "$b" = new ] && C=$NEW
  [ "$cfg" != default ] && x=$cfg
  /usr/bin/time -f "%e\t%U\t%S\t%M" -o "$O/time.txt" $PIN "$C" map $IDX "$@" $x -t 64 --intervals -o /dev/null 2> "$O/log"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$host" "$mode" "$pipe" "$cfg" "$b" "$rep" "$(cat "$O/time.txt")" | tee -a "$TSV"
}
# Warm-up: index and reads into the page cache.
cat $R1 $R2 > /dev/null
$PIN "$NEW" map $IDX reads/hg002_2M.fq -t 64 --intervals -o /dev/null 2>/dev/null
for rep in 1 2 3; do
  for cfg in default --full-cert; do
    for b in old new; do
      run single-end certa-cpu-iv $cfg $b $rep $R1
      run paired certa-cpu-iv $cfg $b $rep $R1 $R2
      if [ $GPU = 1 ]; then
        run single-end certa-gpu-iv $cfg $b $rep $R1 --gpu --device 1
        run paired certa-gpu-iv $cfg $b $rep $R1 $R2 --gpu --device 1
      fi
    done
  done
done
rm -f "$O/time.txt" "$O/log"
echo ALL_DONE
