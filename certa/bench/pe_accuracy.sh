#!/bin/bash
# Paired-end accuracy on GIAB HG002 chr20 (all 444M pairs mapped, chr20 kept):
# CERTA PE (A100) + minibwa paired fallback via FIFOs, and BWA-MEM2 PE.
# Then GATK HaplotypeCaller + bcftools + rtg vcfeval (eval_bam.sh).
set -euo pipefail
cd ~/certa-bench/acc2
E=~/micromamba/envs/certa-bench/bin
C=~/certa-bench/certa-sl/certa/build/certa
D=~/certa-bench/certa-sl/certa/bench
R1=~/certa-data/SRR37356339_1.fastq.gz
R2=~/certa-data/SRR37356339_2.fastq.gz
RG='@RG\tID:hg002\tSM:HG002'

run_certa() {
  echo "== $(date +%T) CERTA PE on all pairs (A100 + 32 threads; minibwa paired fallback 32 threads)"
  rm -f pe_u1 pe_u2 pe_r1 pe_r2; mkfifo pe_u1 pe_u2 pe_r1 pe_r2
  local start=$(date +%s)
  pigz -dc -p 4 $R1 > pe_r1 &
  pigz -dc -p 4 $R2 > pe_r2 &
  $C map ../idx/grch38.s1v3.cidx pe_r1 pe_r2 -t 32 --gpu --device 1 -u pe_u1 -U pe_u2 \
      --stats ourspe.certa.json 2> ourspe.certa.log \
    | awk -f chr20.awk | $E/samtools sort -@ 4 -m 2G -o ourspe.certa.chr20.bam - &
  local cp=$!
  $E/minibwa map -t 32 ../idx/minibwa pe_u1 pe_u2 2> ourspe.minibwa.log \
    | awk -f $D/drop_certified_mates.awk | awk -f chr20.awk \
    | $E/samtools sort -@ 4 -m 2G -o ourspe.fallback.chr20.bam -
  wait $cp
  wait
  echo "CERTA PE pipeline wall $(( $(date +%s) - start )) s"
  grep -E "reads: cert|pairs|time:" ourspe.certa.log
  $E/samtools merge -f -@ 8 ourspe.merged.bam ourspe.certa.chr20.bam ourspe.fallback.chr20.bam
  $E/samtools addreplacerg -r "$RG" -o ourspe.chr20.bam ourspe.merged.bam
  $E/samtools index ourspe.chr20.bam
  rm -f ourspe.merged.bam pe_u1 pe_u2 pe_r1 pe_r2
  echo "== $(date +%T) CERTA PE done"
}

run_bwa() {
  echo "== $(date +%T) BWA-MEM2 PE on all pairs (64 threads)"
  local start=$(date +%s)
  $E/bwa-mem2 mem -t 64 -R "$RG" ../idx/bwa-mem2 $R1 $R2 2> bwape.log \
    | awk -f chr20.awk | $E/samtools sort -@ 4 -m 2G -o bwape.chr20.bam -
  echo "BWA-MEM2 PE pipeline wall $(( $(date +%s) - start )) s"
  $E/samtools index bwape.chr20.bam
  echo "== $(date +%T) BWA-MEM2 PE done"
}

case "${1:-all}" in
  certa) run_certa; ./eval_bam.sh ourspe ourspe.chr20.bam > eval_ourspe.tsv 2> eval_ourspe.err ;;
  bwa)   run_bwa;   ./eval_bam.sh bwape bwape.chr20.bam > eval_bwape.tsv 2> eval_bwape.err ;;
esac
echo "== $(date +%T) all done"
