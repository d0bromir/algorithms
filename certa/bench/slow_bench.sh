#!/bin/bash
# CERTA's full certification (--full-cert) against the slow but reliable
# mappers, galaxy, 64 threads on cores 0-63, index loading included (the
# indexes are built beforehand by bench/slow_index.sh).
#  - Short reads, 2 M HG002 reads (single-end) and pairs: the full-sensitivity
#    read mappers, each set to CERTA's edit radius (certificate radius <= 5
#    edits; 4 % and 96 % of 150 bases allow 6):
#      Yara 1.0.5      -e 4 -y full: every alignment with <= 4 % errors;
#      RazerS 3 3.5.12 -i 96 -rr 100: 100 % recognition of every match with
#                      >= 96 % identity (pairs: -ll 450 -le 400, i.e. 50-850).
#    CERTA in interval mode (every read gets a record), CPU and GPU, with and
#    without the q = 15 second index.
#    RazerS 3's paired mode aborts (std::bad_alloc) at the end of every run
#    in this build (aarch64 Bioconda), with -le 150 too: its single-end time
#    is the comparison, the paired run is kept to record that.
#  - Long reads, 242,534 HG002 HiFi reads: Winnowmap 2.03 (-x map-pb,
#    repetitive 15-mers from meryl; it has no short-read mode) against
#    minimap2 2.31 + lr_ceiling --certify, CERTA's long-read pipeline.
#    Winnowmap's dumped index (-d) crashes on load in this build, so it
#    indexes the FASTA in every run; its log gives the time that took.
#    lr_ceiling also runs on Winnowmap's alignments: it proves the optimum and
#    counts the reads whose proven optimum scores above the tool's alignment
#    (CERTA's scores, +1/-4, gaps -6-1k).
# Reliability (bench/weaknesses.py): each short-read tool against CERTA's
# certificates (missed better locus, unmapped although <= 5 edits exist,
# MAPQ >= 20 on a proven exact tie).
#   slow_bench.sh CERTA_BIN LR_CEILING_BIN OUT_DIR [RazerS 3 reads = 2000000] [reps = 3]
#                 [sections = "certa yara razers hifi winnowmap check"]
# The sections may run in separate calls (same OUT_DIR); "check" uses the
# SAMs the others left and deletes them.
set -u
C=$1; LR=$2; O=$3; NRZ=${4:-2000000}; REPS=${5:-3}; SECTIONS=${6:-certa yara razers hifi winnowmap check}
cd ~/certa-bench
mkdir -p "$O"
E=$HOME/micromamba/envs/certa-slow/bin
BIN=$HOME/micromamba/envs/certa-bench/bin
S=idx/slow
IDX=idx/grch38.s1v3.cidx
IDX2=idx/grch38.s1q15.cidx
REF=idx/ref.fa
R1=reads/hg002_2M.fq
R2=reads/hg002_2M_2.fq
HIFI=$HOME/hifi_real/hifi_1x.fa
WEAK=$(dirname "$C")/../bench/weaknesses.py
PIN="taskset -c 0-63"
GPU="--gpu --device 1"
TSV=$O/timings.tsv
[ -f "$TSV" ] || printf "data\ttool\tconfig\trep\twall_s\tuser_s\tsys_s\tmax_rss_kb\texit\n" > "$TSV"
run() {  # data tool config rep command; returns the command's exit status
  local data=$1 tool=$2 cfg=$3 rep=$4; shift 4
  /usr/bin/time -f "%e\t%U\t%S\t%M\t%x" -o "$O/time.txt" $PIN bash -c "$*"
  local st=$?
  printf "%s\t%s\t%s\t%s\t%s\n" "$data" "$tool" "$cfg" "$rep" "$(grep -v '^Command' "$O/time.txt")" | tee -a "$TSV"
  return $st
}
has() { [[ " $SECTIONS " == *" $1 "* ]]; }
# RazerS 3 on the first NRZ reads (pairs) when the whole set takes too long;
# its reliability is then measured against CERTA on the same reads.
Z1=$R1; Z2=$R2
if [ "$NRZ" -lt 2000000 ]; then Z1=$O/rz_1.fq; Z2=$O/rz_2.fq; fi

if has certa; then
  echo "== $(date +%T) CERTA, short reads"
  cat $R1 $R2 > /dev/null
  $PIN "$C" map $IDX $R1 -t 64 --full-cert --index2 $IDX2 --intervals -o /dev/null 2>/dev/null  # warm-up
  for rep in $(seq 1 $REPS); do
    for g in cpu gpu; do
      G=""; [ $g = gpu ] && G=$GPU
      run single-end certa-$g "default --intervals" $rep "'$C' map $IDX $R1 -t 64 $G --intervals -o /dev/null 2>/dev/null"
      run single-end certa-$g "--full-cert --intervals" $rep "'$C' map $IDX $R1 -t 64 $G --full-cert --intervals -o /dev/null 2>/dev/null"
      run single-end certa-$g "--full-cert --index2 q15 --intervals" $rep "'$C' map $IDX $R1 -t 64 $G --full-cert --index2 $IDX2 --intervals -o /dev/null 2>/dev/null"
      run paired certa-$g "default --intervals" $rep "'$C' map $IDX $R1 $R2 -t 64 $G --intervals -o /dev/null 2>/dev/null"
      run paired certa-$g "--full-cert --intervals" $rep "'$C' map $IDX $R1 $R2 -t 64 $G --full-cert --intervals -o /dev/null 2>/dev/null"
      run paired certa-$g "--full-cert --index2 q15 --intervals" $rep "'$C' map $IDX $R1 $R2 -t 64 $G --full-cert --index2 $IDX2 --intervals -o /dev/null 2>/dev/null"
    done
  done
fi
# The slow tools once each; their SAM is kept for the reliability check.
if has yara; then
  echo "== $(date +%T) Yara"
  $PIN "$E/yara_mapper" -e 4 -y full -t 64 -o "$O/warm.sam" $S/yara_grch38 reads/hg002_2M.fq > /dev/null 2>&1  # warm-up
  rm -f "$O/warm.sam"
  run single-end yara "-e 4 -y full" 1 "'$E/yara_mapper' -e 4 -y full -t 64 -o '$O/yara_se.sam' $S/yara_grch38 $R1 > '$O/yara_se.log' 2>&1" ||
    rm -f "$O/yara_se.sam"
  run paired yara "-e 4 -y full" 1 "'$E/yara_mapper' -e 4 -y full -t 64 -o '$O/yara_pe.sam' $S/yara_grch38 $R1 $R2 > '$O/yara_pe.log' 2>&1" ||
    rm -f "$O/yara_pe.sam"
fi
if has razers; then
  echo "== $(date +%T) RazerS 3 on $NRZ reads"
  if [ "$NRZ" -lt 2000000 ]; then
    head -n $((4 * NRZ)) $R1 > "$Z1"
    head -n $((4 * NRZ)) $R2 > "$Z2"
  fi
  cat $REF > /dev/null
  run single-end-$NRZ razers3 "-i 96 -rr 100" 1 "'$E/razers3' -i 96 -rr 100 -tc 64 -o '$O/razers_se.sam' $REF $Z1 > '$O/razers_se.log' 2>&1" ||
    rm -f "$O/razers_se.sam"
  run paired-$NRZ razers3 "-i 96 -rr 100 -ll 450 -le 400" 1 "'$E/razers3' -i 96 -rr 100 -ll 450 -le 400 -tc 64 -o '$O/razers_pe.sam' $REF $Z1 $Z2 > '$O/razers_pe.log' 2>&1" ||
    rm -f "$O/razers_pe.sam"
  if [ "$NRZ" -lt 2000000 ]; then
    for rep in $(seq 1 $REPS); do
      run single-end-$NRZ certa-cpu "--full-cert --index2 q15 --intervals" $rep "'$C' map $IDX $Z1 -t 64 --full-cert --index2 $IDX2 --intervals -o /dev/null 2>/dev/null"
      run paired-$NRZ certa-cpu "--full-cert --index2 q15 --intervals" $rep "'$C' map $IDX $Z1 $Z2 -t 64 --full-cert --index2 $IDX2 --intervals -o /dev/null 2>/dev/null"
    done
  fi
fi
if has hifi; then
  echo "== $(date +%T) long reads (HiFi): minimap2 + lr_ceiling"
  cat $HIFI $S/mm2_hifi.mmi > /dev/null
  run hifi minimap2 "-x map-hifi -c --eqx" 1 "'$BIN/minimap2' -x map-hifi -c --eqx -t 64 $S/mm2_hifi.mmi $HIFI > '$O/mm2_hifi.paf' 2> '$O/mm2_hifi.log'" &&
    run hifi lr_ceiling "--certify (minimap2 alignments)" 1 "'$LR' $IDX $HIFI '$O/mm2_hifi.paf' -t 64 --certify > /dev/null 2> '$O/lr_mm2.log'"
  echo "minimap2: $(grep -c tp:A:P "$O/mm2_hifi.paf") primary alignments, $(cut -f1 "$O/mm2_hifi.paf" | sort -u | wc -l) reads with one" | tee -a "$O/hifi_counts.txt"
  rm -f "$O/mm2_hifi.paf"
fi
if has winnowmap; then
  echo "== $(date +%T) long reads (HiFi): Winnowmap"
  cat $HIFI $REF > /dev/null
  run hifi winnowmap "-x map-pb -W repetitive_k15 -c --eqx (indexes the FASTA)" 1 "'$E/winnowmap' -W $S/repetitive_k15.txt -x map-pb -c --eqx -t 64 $REF $HIFI > '$O/wm_hifi.paf' 2> '$O/wm_hifi.log'" &&
    run hifi lr_ceiling "--certify (Winnowmap alignments)" 1 "'$LR' $IDX $HIFI '$O/wm_hifi.paf' -t 64 --certify > /dev/null 2> '$O/lr_wm.log'"
  echo "winnowmap: $(grep -m1 'loaded/built the index' "$O/wm_hifi.log" | cut -d: -f5 | cut -d'*' -f1) s of its run to index the FASTA" | tee -a "$O/hifi_counts.txt"
  echo "winnowmap: $(grep -c tp:A:P "$O/wm_hifi.paf") primary alignments, $(cut -f1 "$O/wm_hifi.paf" | sort -u | wc -l) reads with one" | tee -a "$O/hifi_counts.txt"
  rm -f "$O/wm_hifi.paf"
fi
if has check; then
  echo "== $(date +%T) certified fractions and reliability"
  log() {  # certa arguments: the run's summary lines
    $PIN "$C" map $IDX "$@" -t 64 -o "$O/certa.sam" 2> "$O/certa.log"
    grep -E "reads: certified|pairs .*both mates|second index" "$O/certa.log"
  }
  tools() {  # the SAMs present, as weaknesses.py arguments: name=file ...
    local n f
    for n in "$@"; do f=${n#*=}; [ -s "$f" ] && echo "$n"; done
  }
  # RazerS 3 flags its first match primary, not its best: RazerS3:all takes
  # each read's best record (weaknesses.py).
  echo "-- CERTA single-end --full-cert --index2 --intervals, 2 M reads"
  log $R1 --full-cert --index2 $IDX2 --intervals
  rz=""; [ "$NRZ" -ge 2000000 ] && rz=RazerS3:all=$O/razers_se.sam
  python3 "$WEAK" "$O/certa.sam" $(tools Yara=$O/yara_se.sam $rz)
  echo "-- CERTA paired --full-cert --index2 --intervals, 2 M pairs"
  log $R1 $R2 --full-cert --index2 $IDX2 --intervals
  rz=""; [ "$NRZ" -ge 2000000 ] && rz=RazerS3:all=$O/razers_pe.sam
  python3 "$WEAK" "$O/certa.sam" $(tools Yara=$O/yara_pe.sam $rz)
  if [ "$NRZ" -lt 2000000 ]; then
    for m in se:"$Z1" pe:"$Z1 $Z2"; do
      [ -s "$O/razers_${m%%:*}.sam" ] || continue
      echo "-- CERTA on RazerS 3's $NRZ reads (${m%%:*})"
      log ${m#*:} --full-cert --index2 $IDX2 --intervals
      python3 "$WEAK" "$O/certa.sam" RazerS3:all="$O/razers_${m%%:*}.sam" RazerS3-primary="$O/razers_${m%%:*}.sam"
    done
  fi
  echo "-- HiFi"
  for f in mm2:minimap2 wm:winnowmap; do
    [ -s "$O/lr_${f%%:*}.log" ] && grep -E "^reads|^certify" "$O/lr_${f%%:*}.log" | sed "s/^/${f#*:}: /"
  done
  cat "$O/hifi_counts.txt" 2>/dev/null
  rm -f "$O"/*.sam "$O"/rz_?.fq
fi
rm -f "$O/time.txt"
echo ALL_DONE
