#!/bin/bash
# vc_map.sh GENOME OUT_DIR (galaxy): Google's 30x NovaSeq PCR-free reads of
# GENOME (HG002, HG005), downloaded once and streamed to five mappings at the
# same time (fanout.py; nothing stored), each keeping the records that overlap
# ~/vc/targets.bed (chr20 and the CMRG genes, 1 kb padding) as a sorted BAM.
# Aligners on cores 0-63, the rest on 64-127; CERTA is the frozen build of
# PR 19 (~/certa-pm).
set -uo pipefail
G=$1; O=$(mkdir -p "$2" && cd "$2" && pwd)
U=https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/30x/$G.novaseq.pcr-free.30x
B=~/certa-bench; IDX=$B/idx/grch38.s1v3.cidx; BIN=~/micromamba/envs/certa-bench/bin
C=~/certa-pm/build/certa; HS1=~/certa-t2t/hs1.cidx; CHAIN=~/headroom/hs1ToHg38.over.chain.gz
PIN="taskset -c 0-63"; AUX="taskset -c 64-127"
TV="/usr/bin/time -f %e\t%U\t%S\t%M"
IV="--intervals --iv-mapq cert --iv-pair all -t 64 -o /dev/stdout"
T=~/vc/targets.bed
names="ivcal pair pairt2t minibwa bwa-mem2"
(cd ~/certa-pm && find include src -type f | sort | xargs sha256sum | sha256sum | cut -c1-16) > $O/certa_source.sha
F=$O/fifo; rm -rf $F; mkdir $F
for n in $names; do for m in 1 2; do mkfifo $F/$n.$m.fq; done; done
echo "== $(date +%T) $G: stream to $names"
$AUX python3 ~/vc/scripts/fanout.py $U.R1.fastq.gz $U.R2.fastq.gz $F $names 2> $O/fanout.log &
FP=$!
for n in $names; do
  r="$F/$n.1.fq $F/$n.2.fq"
  case $n in
    ivcal) cmd="$C map $IDX $r $IV --mapq cal" ;;
    pair) cmd="$C map $IDX $r $IV --mapq pair" ;;
    pairt2t) cmd="$C map $IDX $r $IV --mapq pair --t2t $HS1 --t2t-chain $CHAIN" ;;
    minibwa) cmd="$BIN/minibwa map -t 64 $B/idx/minibwa $r" ;;
    bwa-mem2) cmd="$BIN/bwa-mem2 mem -t 64 $B/idx/bwa-mem2 $r" ;;
  esac
  ( $PIN $TV -o $O/$n.time $cmd 2> $O/$n.log \
      | $AUX $BIN/samtools view -u -L $T - \
      | $AUX $BIN/samtools sort -@ 4 -m 1500M -T $O/$n.sorttmp -o $O/$n.bam - \
    && $AUX $BIN/samtools index $O/$n.bam && echo ok > $O/$n.ok ) &
done
wait $FP; echo "fanout exit $?"
wait
for n in $names; do echo "$n $(cat $O/$n.ok 2>/dev/null || echo FAILED) $(cat $O/$n.time | tr '\t' ' ')"; done
rm -rf $F
echo "== $(date +%T) $G finished"
