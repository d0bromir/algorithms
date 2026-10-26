#!/bin/bash
# Is the current build's interval-mode output on HG005 the one whose F1 was
# measured (acc_HG005/iv.chr20.bam, mapper 906161c, aim5_2026-10/hg005_rerun_f1.tsv)?
# Same steps as bench/sample_accuracy.sh certa-iv; the chr20 records compared.
set -euo pipefail
S=HG005
B=~/certa-bench
E=~/micromamba/envs/certa-bench/bin
export PATH=$E:$PATH
C=${CERTA:-$HOME/certa-cpuopt2/build/certa}
cd $B/acc_$S
RG="@RG\tID:$S\tSM:$S"
R1=~/certa-data/SRR14724528_1.fastq.gz
R2=~/certa-data/SRR14724528_2.fastq.gz
echo "== $(date +%T) interval mode, current build"
rm -f r1 r2; mkfifo r1 r2
start=$(date +%s)
pigz -dc -p 8 $R1 > r1 &
pigz -dc -p 8 $R2 > r2 &
$C map $B/idx/grch38.s1v3.cidx r1 r2 -t 64 --gpu --device 1 --intervals --stats ivnow.json 2> ivnow.log \
  | awk -f $B/acc2/chr20.awk | samtools sort -@ 4 -m 2G -o ivnow.only.chr20.bam -
wait
rm -f r1 r2
echo "wall $(( $(date +%s) - start )) s"
grep -E "reads: cert|pairs|intervals|time:" ivnow.log
samtools addreplacerg -r "$RG" -o ivnow.chr20.bam ivnow.only.chr20.bam
rm -f ivnow.only.chr20.bam
a=$(samtools view iv.chr20.bam | md5sum | cut -d' ' -f1)
b=$(samtools view ivnow.chr20.bam | md5sum | cut -d' ' -f1)
echo "records: evaluated $a, current $b: $([ "$a" = "$b" ] && echo identical || echo DIFFERENT)"
if [ "$a" != "$b" ]; then
  echo "differing records: $(diff <(samtools view iv.chr20.bam | cut -f1-9) <(samtools view ivnow.chr20.bam | cut -f1-9) | grep -c '^<' || true)"
fi
echo done > ivnow.done
