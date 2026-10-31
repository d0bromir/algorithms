#!/bin/bash
# --full-cert at full scale, in the final protocol's conditions
# (bench/final_bench.sh): 20,001,073 HG002 reads (single-end, R1) or pairs,
# every run on the same 64 cores (galaxy: taskset -c 0-63; a2: numactl
# --interleave=all), SAM to /dev/null, interval mode (every read a record),
# three rounds. Configurations: default; --full-cert --adaptive-parts 0
# (passes 3 and L); --full-cert (with pass V). CPU, and on galaxy the A100
# too. The first round's logs give the certified shares (certified.txt).
#   full_cert_bench.sh HOST CERTA_BIN OUT_DIR
set -u
if [ $# -ne 3 ]; then
  echo "usage: full_cert_bench.sh HOST CERTA_BIN OUT_DIR  (HOST: galaxy or a2)" >&2; exit 2
fi
host=$1; C=$2; O=$3
cd ~/certa-bench
IDX=idx/grch38.s1v3.cidx
R1=reads/hg002_20M.fq
if [ "$host" = a2 ]; then R2=/data/certa-bench/reads/hg002_20M_2.fq; PIN="numactl --interleave=all"; GPU=0
else R2=reads/hg002_20M_2.fq; PIN="taskset -c 0-63"; GPU=1; fi
# On another host, or without an input, the timings would mean nothing.
if [ "$(hostname -s)" != "$host" ]; then
  echo "full_cert_bench.sh: this host is $(hostname -s), not $host" >&2; exit 2
fi
for f in "$C" $IDX $R1 $R2; do
  [ -e "$f" ] || { echo "full_cert_bench.sh: $f not found" >&2; exit 2; }
done
mkdir -p "$O"
TSV=$O/timings.tsv
[ -f "$TSV" ] || printf "host\tmode\tpipeline\tconfig\trep\twall_s\tuser_s\tsys_s\tmax_rss_kb\n" > "$TSV"
CFGS=("default" "--full-cert --adaptive-parts 0" "--full-cert")
run() {  # mode pipeline config rep certa-args...
  local mode=$1 pipe=$2 cfg=$3 rep=$4; shift 4
  local x=""; [ "$cfg" != default ] && x=$cfg
  /usr/bin/time -f "%e\t%U\t%S\t%M" -o "$O/time.txt" $PIN "$C" map $IDX "$@" $x -t 64 --intervals -o /dev/null 2> "$O/log"
  printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$host" "$mode" "$pipe" "$cfg" "$rep" "$(cat "$O/time.txt")" | tee -a "$TSV"
  if [ "$rep" = 1 ]; then
    { echo "== $host $mode $pipe $cfg"; grep -E "reads: certified|pairs .*both mates|pass V \(" "$O/log"; } >> "$O/certified.txt"
  fi
}
# Warm-up: index and reads into the page cache.
cat $R1 $R2 > /dev/null
$PIN "$C" map $IDX reads/hg002_2M.fq -t 64 --full-cert --intervals -o /dev/null 2>/dev/null
for rep in 1 2 3; do
  for cfg in "${CFGS[@]}"; do
    run single-end certa-cpu-iv "$cfg" $rep $R1
    run paired certa-cpu-iv "$cfg" $rep $R1 $R2
    if [ $GPU = 1 ]; then
      run single-end certa-gpu-iv "$cfg" $rep $R1 --gpu --device 1
      run paired certa-gpu-iv "$cfg" $rep $R1 $R2 --gpu --device 1
    fi
  done
done
rm -f "$O/time.txt" "$O/log"
echo ALL_DONE
