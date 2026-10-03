#!/bin/bash
# Paired-end GIAB accuracy for one sample (chr20 vs GIAB v4.2.1):
#   bench/sample_accuracy.sh SAMPLE R1.fq.gz R2.fq.gz STEP
# STEP: prep (chr20 truth), certa (CERTA PE + minibwa paired fallback),
#       bwa (BWA-MEM2 PE), eval-certa, eval-bwa.
# Works in ~/certa-bench/acc_SAMPLE; truth from ~/certa-data/giab.
set -euo pipefail
S=$1; R1=$2; R2=$3; STEP=$4
B=~/certa-bench
E=~/micromamba/envs/certa-bench/bin
export PATH=$E:$PATH
C=${CERTA:-$B/certa-sl/certa/build/certa}
D=$(cd "$(dirname "$0")" && pwd)
W=$B/acc_$S
mkdir -p $W && cd $W
RG="@RG\tID:$S\tSM:$S"
RTG="java -Xmx8g -jar $(ls $E/../share/rtg-tools-*/RTG.jar)"

case $STEP in
prep)
  G=~/certa-data/giab
  printf "chr20\t20\n" > rename.txt
  bcftools view -r chr20 $G/${S}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz | bcftools annotate --rename-chrs rename.txt -Oz -o truth20.vcf.gz
  bcftools index -t -f truth20.vcf.gz
  awk -v OFS="\t" '$1 == "chr20" { $1 = "20"; print }' $G/${S}_GRCh38_1_22_v4.2.1_benchmark.bed > truth20.bed
  for t in snps indels; do bcftools view -v $t truth20.vcf.gz -Oz -o truth20.$t.vcf.gz && bcftools index -t -f truth20.$t.vcf.gz; done
  ln -sfn $B/acc2/ref20.sdf ref20.sdf
  echo "truth chr20 $S: $(bcftools view -H truth20.snps.vcf.gz | wc -l) SNVs, $(bcftools view -H truth20.indels.vcf.gz | wc -l) indels"
  ;;
certa)
  rm -f u1 u2 r1 r2; mkfifo u1 u2 r1 r2
  start=$(date +%s)
  pigz -dc -p 4 $R1 > r1 &
  pigz -dc -p 4 $R2 > r2 &
  $C map $B/idx/grch38.s1v3.cidx r1 r2 -t 32 --gpu --device 1 -u u1 -U u2 --stats certa.json 2> certa.log \
    | awk -f $B/acc2/chr20.awk | samtools sort -@ 4 -m 2G -o certa.only.chr20.bam - &
  cp=$!
  minibwa map -t 32 $B/idx/minibwa u1 u2 2> minibwa.log \
    | awk -f $D/drop_certified_mates.awk | awk -f $B/acc2/chr20.awk \
    | samtools sort -@ 4 -m 2G -o fallback.chr20.bam -
  wait $cp; wait
  echo "CERTA PE pipeline wall $(( $(date +%s) - start )) s"
  grep -E "reads: cert|pairs|time:" certa.log
  samtools merge -f -@ 8 merged.bam certa.only.chr20.bam fallback.chr20.bam
  samtools addreplacerg -r "$RG" -o certa.chr20.bam merged.bam
  samtools index certa.chr20.bam
  rm -f merged.bam u1 u2 r1 r2
  ;;
bwa)
  start=$(date +%s)
  bwa-mem2 mem -t 64 -R "$RG" $B/idx/bwa-mem2 $R1 $R2 2> bwa.log \
    | awk -f $B/acc2/chr20.awk | samtools sort -@ 4 -m 2G -o bwa.chr20.bam -
  echo "BWA-MEM2 PE pipeline wall $(( $(date +%s) - start )) s"
  samtools index bwa.chr20.bam
  ;;
eval-certa|eval-bwa)
  L=${STEP#eval-}
  gatk --java-options "-Xmx16g" HaplotypeCaller -R $B/idx/ref.fa -I $L.chr20.bam -L 20 -O $L.gatk.vcf.gz \
    --native-pair-hmm-threads 16 > $L.gatk.log 2>&1
  bcftools mpileup -f $B/idx/ref.fa -r 20 -a AD,DP -Ou $L.chr20.bam 2> /dev/null | bcftools call -mv -Oz -o $L.bcftools.vcf.gz 2> /dev/null
  bcftools index -t -f $L.bcftools.vcf.gz
  for caller in gatk bcftools; do
    for t in snps indels; do
      bcftools view -v $t $L.$caller.vcf.gz -Oz -o $L.$caller.$t.vcf.gz && bcftools index -t -f $L.$caller.$t.vcf.gz
      rm -rf ve.$L.$caller.$t
      $RTG vcfeval -b truth20.$t.vcf.gz -c $L.$caller.$t.vcf.gz -e truth20.bed -t ref20.sdf -o ve.$L.$caller.$t --threads 16 > /dev/null 2>&1
      printf "%s\t%s\t%s\t%s\t" $S $L $caller $t; awk 'NR == 3 { print $0 }' ve.$L.$caller.$t/summary.txt
    done
  done | tee eval_$L.tsv
  ;;
*) echo "unknown step $STEP" >&2; exit 1 ;;
esac
echo "== $(date +%T) $S $STEP done"
