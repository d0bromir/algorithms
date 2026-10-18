#!/bin/bash
# Final timing protocol, one session per host: every competitor, CERTA on the
# CPU and (with a GPU) CERTA on the GPU, on 20,001,073 HG002 reads (single-end,
# R1) and the same reads' pairs. Every pipeline gets 64 cores: on galaxy
# (128 cores) everything runs under taskset -c 0-63; a2 has 64 cores. CERTA
# and its fallback (minibwa) run concurrently; "64+64" starts 64 threads in
# each on those 64 cores. SAM goes to /dev/null. Three rounds, the tools in
# turn (paired BWA-MEM2, at over four minutes a run, and the paired fallback
# with minibwa's default batch: one run each).
#   final_bench.sh HOST CERTA_BIN OUT_DIR
set -u
host=$1; C=$2; O=$3
cd ~/certa-bench
mkdir -p "$O"
BIN=$HOME/micromamba/envs/certa-bench/bin
AWK=$(dirname "$C")/../bench/drop_certified_mates.awk
IDX=idx/grch38.s1v3.cidx
R1=reads/hg002_20M.fq
if [ "$host" = a2 ]; then R2=/data/certa-bench/reads/hg002_20M_2.fq; PIN=""; GPU=0
else R2=reads/hg002_20M_2.fq; PIN="taskset -c 0-63"; GPU=1; fi
TSV=$O/timings.tsv
[ -f "$TSV" ] || printf "host\tmode\tpipeline\tconfig\trep\twall_s\tuser_s\tsys_s\tmax_rss_kb\n" > "$TSV"
run() {  # mode pipeline config rep command
  local mode=$1 name=$2 cfg=$3 rep=$4; shift 4
  /usr/bin/time -f "%e\t%U\t%S\t%M" -o "$O/time.txt" $PIN bash -c "$*"
  printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$host" "$mode" "$name" "$cfg" "$rep" "$(cat "$O/time.txt")" | tee -a "$TSV"
}
se_pipe() {  # name config certa-args certa-threads minibwa-threads rep
  run single-end "$1" "$2" "$6" "set -o pipefail; '$C' map $IDX $R1 -t $4 $3 -o /dev/null -u /dev/stdout 2>/dev/null | '$BIN/minibwa' map -t $5 -K 10m idx/minibwa /dev/stdin > /dev/null 2>&1"
}
pe_pipe() {  # name config certa-args certa-threads minibwa-threads minibwa-K rep
  rm -f "$O/u1" "$O/u2"; mkfifo "$O/u1" "$O/u2"
  run paired "$1" "$2" "$7" "'$C' map $IDX $R1 $R2 -t $4 $3 -o /dev/null -u '$O/u1' -U '$O/u2' 2>/dev/null & '$BIN/minibwa' map -t $5 $6 idx/minibwa '$O/u1' '$O/u2' 2>/dev/null | awk -f '$AWK' > /dev/null; wait"
}
GPUARG="--gpu --device 1"
# Warm-up: indexes and reads into the page cache.
cat $R1 $R2 > /dev/null
for f in reads/hg002_2M.fq; do
  "$C" map $IDX $f -t 64 -o /dev/null -u /dev/null 2>/dev/null
  "$BIN/minibwa" map -t 64 idx/minibwa $f > /dev/null 2>&1
  "$BIN/bwa-mem2" mem -t 64 idx/bwa-mem2 $f > /dev/null 2>&1
  "$BIN/strobealign" -t 64 --use-index idx/ref.fa $f > /dev/null 2>&1
  "$BIN/minimap2" -ax sr -t 64 idx/minimap2.sr.mmi $f > /dev/null 2>&1
  "$BIN/bowtie2" -p 64 -x idx/bowtie2 -U $f > /dev/null 2>&1
done
for rep in 1 2 3; do
  run single-end bwa-mem2 - $rep "'$BIN/bwa-mem2' mem -t 64 idx/bwa-mem2 $R1 > /dev/null 2>&1"
  run single-end minibwa - $rep "'$BIN/minibwa' map -t 64 idx/minibwa $R1 > /dev/null 2>&1"
  run single-end strobealign - $rep "'$BIN/strobealign' -t 64 --use-index idx/ref.fa $R1 > /dev/null 2>&1"
  run single-end minimap2 - $rep "'$BIN/minimap2' -ax sr -t 64 idx/minimap2.sr.mmi $R1 > /dev/null 2>&1"
  run single-end bowtie2 - $rep "'$BIN/bowtie2' -p 64 -x idx/bowtie2 -U $R1 > /dev/null 2>&1"
  se_pipe certa-cpu 64+64 "" 64 64 $rep
  run single-end certa-cpu-iv 64 $rep "'$C' map $IDX $R1 -t 64 --intervals -o /dev/null 2>/dev/null"
  if [ $GPU = 1 ]; then
    se_pipe certa-gpu 16+48 "$GPUARG" 16 48 $rep
    se_pipe certa-gpu 64+64 "$GPUARG" 64 64 $rep
    run single-end certa-gpu-iv 64 $rep "'$C' map $IDX $R1 -t 64 $GPUARG --intervals -o /dev/null 2>/dev/null"
  fi
done
run single-end certa-cpu-tier 64 1 "'$C' map $IDX $R1 -t 64 -o /dev/null -u /dev/null 2>/dev/null"
[ $GPU = 1 ] && run single-end certa-gpu-tier 16 1 "'$C' map $IDX $R1 -t 16 $GPUARG -o /dev/null -u /dev/null 2>/dev/null"
run paired bwa-mem2 - 1 "'$BIN/bwa-mem2' mem -t 64 idx/bwa-mem2 $R1 $R2 > /dev/null 2>&1"
pe_pipe certa-cpu "64+64 Kdefault" "" 64 64 "" 1
for rep in 1 2 3; do
  run paired minibwa - $rep "'$BIN/minibwa' map -t 64 idx/minibwa $R1 $R2 > /dev/null 2>&1"
  pe_pipe certa-cpu "64+64 K10m" "" 64 64 "-K 10m" $rep
  run paired certa-cpu-iv 64 $rep "'$C' map $IDX $R1 $R2 -t 64 --intervals -o /dev/null 2>/dev/null"
  if [ $GPU = 1 ]; then
    run paired certa-gpu-iv 64 $rep "'$C' map $IDX $R1 $R2 -t 64 $GPUARG --intervals -o /dev/null 2>/dev/null"
    pe_pipe certa-gpu "20+44 K10m" "$GPUARG" 20 44 "-K 10m" $rep
    pe_pipe certa-gpu "64+64 K10m" "$GPUARG" 64 64 "-K 10m" $rep
  fi
done
echo done > "$O/done"
