#!/bin/bash
# Does minibwa's batch size (-K 10m vs its default) change variant calling?
# HG005 paired, chr20 vs GIAB v4.2.1 (protocol of bench/sample_accuracy.sh).
# CERTA runs once (its fallback pairs saved to files); minibwa then maps the
# fallback twice, with each -K; both results are merged with the same
# certified records and evaluated. Only the fallback's batching differs.
set -euo pipefail
S=HG005
B=~/certa-bench
E=~/micromamba/envs/certa-bench/bin
export PATH=$E:$PATH
C=${CERTA:-$HOME/certa-cpuopt2/build/certa}
AWKD=$(cd "$(dirname "$0")" && pwd)/drop_certified_mates.awk
W=$B/acc_$S
cd $W
RG="@RG\tID:$S\tSM:$S"
RTG="java -Xmx8g -jar $(ls $E/../share/rtg-tools-*/RTG.jar)"
R1=~/certa-data/SRR14724528_1.fastq.gz
R2=~/certa-data/SRR14724528_2.fastq.gz
echo "== $(date +%T) CERTA"
rm -f r1 r2; mkfifo r1 r2
pigz -dc -p 8 $R1 > r1 &
pigz -dc -p 8 $R2 > r2 &
$C map $B/idx/grch38.s1v3.cidx r1 r2 -t 32 --gpu --device 1 -u kfb_1.fq -U kfb_2.fq --stats k.json 2> k.log \
  | awk -f $B/acc2/chr20.awk | samtools sort -@ 4 -m 2G -o k.certa.chr20.bam -
wait
rm -f r1 r2
grep -E "reads: cert|pairs|time:" k.log
for K in def 10m; do
  opt=""; [ $K = 10m ] && opt="-K 10m"
  echo "== $(date +%T) minibwa $opt"
  minibwa map -t 64 $opt $B/idx/minibwa kfb_1.fq kfb_2.fq 2> k$K.minibwa.log \
    | awk -f $AWKD | awk -f $B/acc2/chr20.awk | samtools sort -@ 4 -m 2G -o k$K.fallback.chr20.bam -
  samtools merge -f -@ 8 k$K.merged.bam k.certa.chr20.bam k$K.fallback.chr20.bam
  samtools addreplacerg -r "$RG" -o k$K.chr20.bam k$K.merged.bam
  samtools index k$K.chr20.bam
  rm -f k$K.merged.bam
done
evalb() {
  L=$1
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
  done > eval_$L.tsv
}
echo "== $(date +%T) evaluation"
evalb kdef &
evalb k10m &
wait
cat eval_kdef.tsv eval_k10m.tsv
echo "records differing: $(diff <(samtools view kdef.chr20.bam | cut -f1-9) <(samtools view k10m.chr20.bam | cut -f1-9) | grep -c '^<' || true)"
rm -f kfb_1.fq kfb_2.fq
echo "== $(date +%T) done"
echo done > k.done
