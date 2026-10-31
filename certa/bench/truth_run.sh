#!/bin/bash
# truth_run.sh: phase 1 runs (VALIDATION_PLAN.md): truth for tests A and B, the aligners on
# GRCh38, scoring and reports. Aligners get cores 0-63 and 64 threads; BAM
# compression, scoring and truth work run on 64-127.
#   truth_run.sh [STEP...]   steps: truthA truthB align score report (default: all)
set -eu
W=~/certa-truth; D=~/certa-data/truth; B=~/certa-bench
C=~/certa-bio/certa/build/certa
S=~/certa-bio/certa/build/bio_score
T=~/certa-bio/certa/bench
BIN=~/micromamba/envs/certa-bench/bin
DRAG=~/micromamba/envs/certa-k/bin/dragen-os
IDX=$B/idx/grch38.s1v3.cidx
FAI=$B/idx/ref.fa.fai
STRAT=$D/strat
BENCH=$D/giab/HG002_GRCh38_v5.0q_smvar.benchmark.bed
INSERT=429,99
export PATH=~/tools/env/truth/bin:$PATH
PIN="taskset -c 0-63"; AUX="taskset -c 64-127"
TV="/usr/bin/time -f %e\t%U\t%S\t%M"
cd $W
step() { echo "== $(date +%T) $*"; }
bam() { $AUX samtools view -@ 16 -b --output-fmt-option level=1 -o "$1" -; }
STEPS="${*:-truthA truthB align score report}"
has() { case " $STEPS " in *" $1 "*) return 0;; esac; return 1; }

if has truthA; then
  step "test A truth"
  $AUX python3 $T/truth_lift.py sim A/orig_1.fq $FAI lift/blocks.npz $STRAT $BENCH A/truth.npz
  $AUX python3 $T/truth_lift.py check A/orig_1.fq A/orig_2.fq A/truth.npz $B/idx/ref.fa A/check --n 100000 | tee A/check.txt
fi

if has truthB; then
  step "test B: CERTA default and minibwa on each haplotype"
  NB=$(( $(wc -l < B/B_1.fq) / 4 ))
  for h in mat_Y_EBV_MT pat_X_EBV_MT; do
    $PIN $TV -o B/hap_$h.certa.time $C map hap/$h.cidx B/B_1.fq B/B_2.fq -t 64 -o /dev/stdout -u /dev/null -U /dev/null \
      --stats B/hap_$h.certa.json 2> B/hap_$h.certa.log | bam B/hap_$h.certa.bam
    $PIN $TV -o B/hap_$h.minibwa.time $BIN/minibwa map -t 64 hap/$h B/B_1.fq B/B_2.fq 2> B/hap_$h.minibwa.log | bam B/hap_$h.minibwa.bam
  done
  step "test B truth ($NB pairs)"
  $AUX python3 $T/truth_label.py <(samtools view B/hap_mat_Y_EBV_MT.certa.bam) <(samtools view B/hap_pat_X_EBV_MT.certa.bam) \
    $NB $FAI lift/blocks.npz $STRAT $BENCH B/truth.npz --margin 2 \
    --check <(samtools view B/hap_mat_Y_EBV_MT.minibwa.bam) <(samtools view B/hap_pat_X_EBV_MT.minibwa.bam) 2>&1 | tee B/truth.txt
fi

ALIGNERS="certa-bio certa minibwa bwa-mem2 strobealign minimap2 dragmap"
if has align; then
  for X in A B; do
    R1=$X/${X}_1.fq; R2=$X/${X}_2.fq
    for a in $ALIGNERS; do
      [ -s $X/$a.bam ] && continue
      step "test $X: $a"
      case $a in
        certa-bio) cmd="$C map $IDX $R1 $R2 --bio --bio-insert $INSERT -t 64 -o /dev/stdout -u /dev/null -U /dev/null --stats $X/$a.json" ;;
        certa)     cmd="$C map $IDX $R1 $R2 -t 64 -o /dev/stdout -u /dev/null -U /dev/null --stats $X/$a.json" ;;
        minibwa)   cmd="$BIN/minibwa map -t 64 $B/idx/minibwa $R1 $R2" ;;
        bwa-mem2)  cmd="$BIN/bwa-mem2 mem -t 64 $B/idx/bwa-mem2 $R1 $R2" ;;
        strobealign) cmd="$BIN/strobealign -t 64 --use-index $B/idx/ref.fa $R1 $R2" ;;
        minimap2)  cmd="$BIN/minimap2 -ax sr -t 64 $B/idx/minimap2.sr.mmi $R1 $R2" ;;
        dragmap)   cmd="$DRAG -r $B/idx/dragmap -1 $R1 -2 $R2 --num-threads 64" ;;
      esac
      $PIN $TV -o $X/$a.time $cmd 2> $X/$a.log | bam $X/$a.bam.tmp && mv $X/$a.bam.tmp $X/$a.bam
      printf "%s\t%s\t%s\n" "$X" "$a" "$(cat $X/$a.time)" | tee -a timings.tsv
    done
  done
fi

if has score; then
  for X in A B; do
    step "test $X: scoring"
    for a in $ALIGNERS; do
      [ -s $X/$a.npz ] && continue
      ( samtools view $X/$a.bam | tee >($S $IDX > $X/$a.score.tsv 2> $X/$a.score.log) \
          | python3 $T/truth_eval.py score $X/truth.npz $FAI $X/$a.npz 2> $X/$a.eval.log ) &
      while [ "$(jobs -rp | wc -l)" -ge 7 ]; do sleep 5; done
    done
    wait
  done
fi

if has report; then
  for X in A B; do
    step "test $X: report"
    runs=""; scores=""
    for a in $ALIGNERS; do
      n=$a; [ $a = certa-bio ] && n="certa --bio"
      runs="$runs $(printf '%q' "$n=$X/$a.npz")"
      [ $a = certa-bio ] || scores="$scores $(printf '%q' "$n=$X/$a.score.tsv")"
    done
    eval python3 $T/truth_eval.py report $X/truth.npz $runs --scores $scores > $X/report.md
    cat $X/report.md
  done
fi
step finished
