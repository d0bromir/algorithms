#!/bin/bash
# galaxy: dv_map.sh BUILD_DIR TAG [OPTIONS...]: the vc subsets of HG002 and
# HG005 (vc_subset.sh's extraction, from galaxy's copies of the five target
# CRAMs) mapped by BUILD_DIR/certa with --intervals --iv-pair all -t 64
# OPTIONS, the records overlapping the targets sorted into ~/vc/G/S_TAG.bam
# (+ .bai, .md5), timed into ~/vc/G/S_TAG.time, for call.sh on a2.
set -euo pipefail
C=$1/certa; TAG=$2; shift 2
V=~/vc; ST=~/micromamba/envs/certa-bench/bin/samtools; REF=~/certa-bench/idx/ref.fa
IDX=~/certa-bench/idx/grch38.s1v3.cidx; T=$V/targets.bed
TV="/usr/bin/time -f %e\t%U\t%S\t%M"
step() { echo "== $(date +%T) $*"; }
for G in HG002 HG005; do
  O=$V/$G; S=$V/subset_$G; mkdir -p $S; n=S_$TAG
  [ -s $O/$n.bam.md5 ] && continue
  if [ ! -s $S/R2.fq ]; then
    step "$G reads"
    { $ST view -H $O/minibwa.cram
      for m in ivcal pair pairt2t minibwa bwa-mem2; do taskset -c 64-127 $ST view -T $REF -F 0x900 $O/$m.cram; done \
        | awk 'BEGIN{FS=OFS="\t"} { k = $1 "\t" int($2 / 64) % 4; if (!(k in s)) { s[k] = 1; print } }'
    } | taskset -c 64-127 $ST collate -u -O - $S/col | taskset -c 64-127 $ST fastq -n -1 $S/R1.fq -2 $S/R2.fq \
        -s $S/single.fq -0 /dev/null - 2> $S/fastq.log
  fi
  step "$G map $n: $*"
  taskset -c 0-63 $TV -o $O/$n.time $C map $IDX $S/R1.fq $S/R2.fq --intervals --iv-pair all -t 64 -o /dev/stdout "$@" \
      2> $O/$n.log | taskset -c 64-127 $ST view -@ 8 -u -L $T - \
    | taskset -c 64-127 $ST sort -@ 16 -m 2G -T $S/sort.$n -o $O/$n.bam -
  taskset -c 64-127 $ST index $O/$n.bam
  md5sum < $O/$n.bam > $O/$n.bam.md5
  cat $O/$n.time
done
rm -f $V/subset_HG00?/R1.fq $V/subset_HG00?/R2.fq $V/subset_HG00?/single.fq
step "done"
