#!/bin/bash
# truth_hg005.sh: the held-out test (VALIDATION_PLAN.md, phase 2;
# bench/VALIDATION_RESULTS.md, "Held out: HG005"). HG005's truth is built as
# test B's from its HPRC release 2 haplotypes, then the frozen configurations
# and the other aligners are scored at the locus, with the T2T second
# opinions. Aligners on cores 0-63, the rest on 64-127.
#   truth_hg005.sh BUILD_DIR HS1_FA HS1_CIDX CHAIN OUT_DIR [STEP...]   steps: prep label align score
# HS1_FA is indexed for bwa-mem2, HS1_CIDX is `certa index HS1_FA -q 22 -s 1`,
# CHAIN is UCSC's hs1ToHg38.over.chain.gz.
# - Haplotypes: HG005_{pat,mat}_hprc_r2_v1.0.1.fa.gz, aligned to GRCh38 with
#   dipcall's minimap2 settings (-x asm5 -z200000,10000 --cs -r2k; -K 4g puts
#   every chromosome in one batch, one thread each); GIAB has no dipcall
#   output for HG005.
# - Reads: the first 20 M pairs of SRR14724528 (NovaSeq 6000, 2x151),
#   trimmed as test B's, named h<pair index>.
# - Labels: truth_label.py on CERTA's default mode and minibwa on each
#   haplotype; each haplotype index is deleted once its reads are mapped.
set -euo pipefail
BUILD=$1; HS1_FA=$2; HS1_CIDX=$3; CHAIN=$4; O=$(mkdir -p "$5" && cd "$5" && pwd); shift 5
C=$BUILD/certa; T=$(cd "$(dirname "$0")" && pwd)
W=~/certa-truth; X=H; D=~/certa-data; B=~/certa-bench
BIN=~/micromamba/envs/certa-bench/bin; BWA_HS1=~/micromamba/envs/refmappers/bin/bwa-mem2
FAI=$B/idx/ref.fa.fai; IDX=$B/idx/grch38.s1v3.cidx
STRAT=$D/truth/strat; BENCH=$D/giab/HG005_GRCh38_1_22_v4.2.1_benchmark.bed
export PATH=~/tools/env/truth/bin:$PATH
PIN="taskset -c 0-63"; AUX="taskset -c 64-127"
TV="/usr/bin/time -f %e\t%U\t%S\t%M"
bam() { $AUX samtools view -@ 16 -b --output-fmt-option level=1 -o "$1" -; }
step() { echo "== $(date +%T) $*"; }
STEPS="${*:-prep label align score}"
has() { case " $STEPS " in *" $1 "*) return 0;; esac; return 1; }
mkdir -p $W/$X $W/hap5 $W/lift5 $O/place
cd $W
R1=$X/${X}_1.fq; R2=$X/${X}_2.fq

if has prep; then
  step "haplotypes"
  for h in pat mat; do
    [ -s hap5/$h.fa ] || pigz -dc -p 8 $D/truth/hg005/HG005_${h}_hprc_r2_v1.0.1.fa.gz > hap5/$h.fa
    [ -s hap5/$h.fa.fai ] || samtools faidx hap5/$h.fa
  done
  step "haplotypes to GRCh38 (minimap2, dipcall's settings), in the background"
  for h in pat mat; do
    [ -s hap5/$h.grch38.bam ] || ( $AUX minimap2 -a -x asm5 -z200000,10000 --cs -r2k -K 4g -t 32 $B/idx/ref.fa hap5/$h.fa \
      2> hap5/$h.mm2.log | samtools view -b -o hap5/$h.grch38.bam.tmp - && mv hap5/$h.grch38.bam.tmp hap5/$h.grch38.bam ) &
  done
  if [ ! -s $R1 ]; then
    step "reads: 20 M pairs of SRR14724528, trimmed"
    for m in 1 2; do ( pigz -dc -p 4 $D/SRR14724528_$m.fastq.gz | head -n 80000000 > $X/raw_$m.fq ) || true; done
    [ "$(wc -l < $X/raw_1.fq)" = 80000000 ] && [ "$(wc -l < $X/raw_2.fq)" = 80000000 ]
    $PIN $TV -o $X/trim.time ~/tools/bin/fastp_d0bromir_cpu_v5 -i $X/raw_1.fq -I $X/raw_2.fq -o $X/t_1.fq -O $X/t_2.fq \
      -2 -g -Q -l 30 -w 32 -j $X/trim.json -h $X/trim.html > $X/trim.log 2>&1
    for m in 1 2; do awk 'NR % 4 == 1 {print "@h" (NR - 1) / 4; next} 1' $X/t_$m.fq > $X/${X}_$m.fq && rm $X/t_$m.fq $X/raw_$m.fq; done
  fi
  echo "pairs: $(( $(wc -l < $R1) / 4 ))"
  wait
  step "projection blocks"
  [ -s lift5/blocks.npz ] || $AUX python3 $T/truth_lift.py blocks hap5/pat.grch38.bam hap5/mat.grch38.bam $FAI lift5/blocks.npz
fi

if has label; then
  NB=$(( $(wc -l < $R1) / 4 ))
  for h in pat mat; do
    if [ ! -s $X/hap_$h.certa.bam ]; then
      step "CERTA index and default mode on $h"
      $AUX $TV -o hap5/$h.cidx.time $C index hap5/$h.fa -o hap5/$h.cidx -q 22 -s 1 -t 32 > hap5/$h.cidx.log 2>&1
      $PIN $C map hap5/$h.cidx $R1 $R2 -t 64 -o /dev/stdout -u /dev/null -U /dev/null 2> $X/hap_$h.certa.log \
        | bam $X/hap_$h.certa.bam
      rm -f hap5/$h.cidx
    fi
    if [ ! -s $X/hap_$h.minibwa.bam ]; then
      step "minibwa index and map on $h"
      $AUX $TV -o hap5/$h.mbw.time $BIN/minibwa index -t 32 hap5/$h.fa hap5/$h > hap5/$h.mbw.log 2>&1
      $PIN $BIN/minibwa map -t 64 hap5/$h $R1 $R2 2> $X/hap_$h.minibwa.log | bam $X/hap_$h.minibwa.bam
      rm -f hap5/$h.mbw hap5/$h.l2b
    fi
  done
  step "truth ($NB pairs)"
  [ -s $X/truth.npz ] || $AUX python3 $T/truth_label.py <(samtools view $X/hap_mat.certa.bam) <(samtools view $X/hap_pat.certa.bam) \
    $NB $FAI lift5/blocks.npz $STRAT $BENCH $X/truth.npz --margin 2 \
    --check <(samtools view $X/hap_mat.minibwa.bam) <(samtools view $X/hap_pat.minibwa.bam) 2>&1 | tee $X/truth.txt
fi

if has align; then
  for a in ivcert ivcal minibwa bwa-mem2 strobealign minimap2 t2t certa_t2t; do
    [ -s $O/place/$a.npz ] && continue
    step "$a"
    chain=""
    case $a in
      ivcert) cmd="$C map $IDX $R1 $R2 --intervals --iv-mapq cert -t 64 -o /dev/stdout" ;;
      ivcal) cmd="$C map $IDX $R1 $R2 --intervals --iv-mapq cert --iv-pair all --mapq cal -t 64 -o /dev/stdout" ;;
      minibwa) cmd="$BIN/minibwa map -t 64 $B/idx/minibwa $R1 $R2" ;;
      bwa-mem2) cmd="$BIN/bwa-mem2 mem -t 64 $B/idx/bwa-mem2 $R1 $R2" ;;
      strobealign) cmd="$BIN/strobealign -t 64 --use-index $B/idx/ref.fa $R1 $R2" ;;
      minimap2) cmd="$BIN/minimap2 -ax sr -t 64 $B/idx/minimap2.sr.mmi $R1 $R2" ;;
      t2t) cmd="$BWA_HS1 mem -t 64 $HS1_FA $R1 $R2"; chain="--chain $CHAIN" ;;
      certa_t2t) cmd="$C map $HS1_CIDX $R1 $R2 --intervals --iv-mapq cert --iv-pair all -t 64 -o /dev/stdout"; chain="--chain $CHAIN" ;;
    esac
    $PIN $TV -o $O/place/$a.time $cmd 2> $O/place/$a.log \
      | $AUX python3 $T/truth_t2t.py place $chain $W/$X/truth.npz $FAI $O/place/$a.npz.tmp.npz 2> $O/place/$a.place.log
    mv $O/place/$a.npz.tmp.npz $O/place/$a.npz
  done
fi

if has score; then
  step "tables"
  {
    echo "## HG005: BWA-MEM2 on T2T-CHM13 as the second opinion"; echo
    $AUX python3 $T/truth_t2t.py eval $W/$X/truth.npz $O/place
    echo; echo "## HG005: CERTA on T2T-CHM13 as the second opinion, certified placements only"; echo
    $AUX python3 $T/truth_t2t.py eval --t2t certa_t2t --certified $W/$X/truth.npz $O/place
    echo; echo "## HG005, control: BWA-MEM2 on GRCh38 as the second opinion"; echo
    $AUX python3 $T/truth_t2t.py eval --second bwa-mem2 $W/$X/truth.npz $O/place
  } > $O/t2t_H.md
fi
step finished
