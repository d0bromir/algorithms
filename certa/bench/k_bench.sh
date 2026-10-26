#!/bin/bash
# Competitors added from SWOT_ANALYSIS.md section K, in the final protocol's
# conditions (bench/final_bench.sh): 64 cores (galaxy: taskset -c 0-63; a2:
# numactl --interleave=all), warmed page cache, three rounds, output to
# /dev/null where the tool allows.
#   k_bench.sh HOST OUT_DIR TOOL...
# Tools:
#   dragmap   DRAGMAP 1.3.0 (Illumina's open-source DRAGEN mapper), 20 M reads
#             and pairs (galaxy: aarch64 build; a2: x86-64 build)
#   columba   Columba 2.0 (lossless, best mode, >= 96 % identity, as Yara -e 4),
#             2 M reads and pairs (galaxy; built with -fsigned-char: with
#             aarch64's unsigned char its FASTQ reader misses EOF and drops the
#             last batch)
#   mm2fast   mm2-fast (AVX-512 minimap2 2.24) against minimap2 2.24 and 2.31
#             -ax sr on 20 M reads, and its output against 2.24's (a2)
# Not runnable here: Parabricks 4.8 fq2bam (its arm64 build stops with SIGILL
# on galaxy's Ampere Altra, ARMv8.2; a2 has no GPU).
set -u
host=$1; O=$2; shift 2
cd ~/certa-bench
mkdir -p "$O"
TSV=$O/timings.tsv
[ -f "$TSV" ] || printf "host\tmode\tpipeline\tconfig\trep\twall_s\tuser_s\tsys_s\tmax_rss_kb\n" > "$TSV"
R1=reads/hg002_20M.fq
if [ "$host" = a2 ]; then PIN="numactl --interleave=all"; R2=/data/certa-bench/reads/hg002_20M_2.fq
else PIN="taskset -c 0-63"; R2=reads/hg002_20M_2.fq; fi
E=$HOME/micromamba/envs/certa-k
run() {  # mode pipeline config rep command
  local mode=$1 name=$2 cfg=$3 rep=$4; shift 4
  /usr/bin/time -f "%e\t%U\t%S\t%M" -o "$O/time.txt" $PIN bash -c "$*"
  printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$host" "$mode" "$name" "$cfg" "$rep" "$(cat "$O/time.txt")" | tee -a "$TSV"
}
for tool in "$@"; do
  case $tool in
    dragmap)
      echo "dragmap $($E/bin/dragen-os --version 2>&1 | head -1)" >> "$O/versions.txt"
      cat $R1 $R2 idx/dragmap/* > /dev/null
      for rep in 1 2 3; do
        run single-end dragmap - $rep "'$E/bin/dragen-os' -r idx/dragmap -1 $R1 --num-threads 64 > /dev/null 2>&1"
        run paired dragmap - $rep "'$E/bin/dragen-os' -r idx/dragmap -1 $R1 -2 $R2 --num-threads 64 > /dev/null 2>&1"
      done ;;
    columba)
      COL=$HOME/tools/src/columba/build_Vanilla/columba
      echo "columba $(git -C $HOME/tools/src/columba rev-parse --short HEAD) (-fsigned-char)" >> "$O/versions.txt"
      for rep in 1 2 3; do
        # It writes only to a .sam/.rhs file: a scratch SAM, deleted after each run.
        run single-end-2M columba "best, >= 96 %" $rep "'$COL' -r idx/columba/grch38 -f reads/hg002_2M.fq -t 64 -a best -I 96 -o $O/columba.sam -l $O/columba.log; rm -f $O/columba.sam"
        run paired-2M columba "best, >= 96 %" $rep "'$COL' -r idx/columba/grch38 -f reads/hg002_2M.fq -F reads/hg002_2M_2.fq -t 64 -a best -I 96 -o $O/columba.sam -l $O/columba.log; rm -f $O/columba.sam"
      done ;;
    mm2fast)
      MM2=$HOME/micromamba/envs/certa-bench/bin/minimap2   # 2.31, as in final_a2.tsv
      OLD=$HOME/tools/src/minimap2-2.24/minimap2           # upstream 2.24, mm2-fast's base
      FAST=$HOME/tools/src/mm2-fast/minimap2
      echo "minimap2 $($MM2 --version), minimap2 $($OLD --version), mm2-fast $($FAST --version)" >> "$O/versions.txt"
      cat $R1 > /dev/null
      for b in minimap2-2.24 mm2-fast; do  # one index for all; outputs without @PG
        B=$OLD; [ $b = mm2-fast ] && B=$FAST
        $PIN $B -ax sr -t 64 idx/minimap2.sr.mmi reads/hg002_2M.fq 2>/dev/null | grep -v '^@PG' | md5sum | cut -c1-32 > "$O/$b.2M.md5"
      done
      for rep in 1 2 3; do
        run single-end minimap2 "2.31 -ax sr" $rep "'$MM2' -ax sr -t 64 idx/minimap2.sr.mmi $R1 > /dev/null 2>&1"
        run single-end minimap2 "2.24 -ax sr" $rep "'$OLD' -ax sr -t 64 idx/minimap2.sr.mmi $R1 > /dev/null 2>&1"
        run single-end mm2-fast "-ax sr" $rep "'$FAST' -ax sr -t 64 idx/minimap2.sr.mmi $R1 > /dev/null 2>&1"
      done ;;
    *) echo "unknown tool $tool" >&2 ;;
  esac
done
echo ALL_DONE
