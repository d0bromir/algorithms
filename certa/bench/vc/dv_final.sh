#!/bin/bash
# galaxy: dv_final.sh BUILD_DIR TAG=OPTIONS...: the vc subsets of HG002 and
# HG005 (dv_map.sh's extraction, done once for every TAG) mapped by
# BUILD_DIR/certa with --intervals --iv-pair all -t 64 and each TAG's
# OPTIONS; the records overlapping the targets, sorted, as
# ~/vc/G/S_TAG.cram (+ .crai, .md5; timed into .time) for dv_relay_cram.sh.
# CRAM is about a third of the BAM, which the relay through a third machine
# (3 MB/s) needs.
set -euo pipefail
C=$1/certa; shift
V=~/vc; ST=~/micromamba/envs/certa-bench/bin/samtools; REF=~/certa-bench/idx/ref.fa
IDX=~/certa-bench/idx/grch38.s1v3.cidx; T=$V/targets.bed
TV="/usr/bin/time -f %e\t%U\t%S\t%M"
step() { echo "== $(date +%T) $*"; }
for G in HG002 HG005; do
  O=$V/$G; S=$V/subset_$G; mkdir -p $S
  if [ ! -s $S/R2.fq ]; then
    step "$G reads"
    { $ST view -H $O/minibwa.cram
      for m in ivcal pair pairt2t minibwa bwa-mem2; do taskset -c 64-127 $ST view -T $REF -F 0x900 $O/$m.cram; done \
        | awk 'BEGIN{FS=OFS="\t"} { k = $1 "\t" int($2 / 64) % 4; if (!(k in s)) { s[k] = 1; print } }'
    } | taskset -c 64-127 $ST collate -u -O - $S/col | taskset -c 64-127 $ST fastq -n -1 $S/R1.fq -2 $S/R2.fq \
        -s $S/single.fq -0 /dev/null - 2> $S/fastq.log
  fi
  for spec in "$@"; do
    tag=${spec%%=*}; opts=${spec#*=}; n=S_$tag
    [ -s $O/$n.cram.md5 ] && continue
    step "$G map $n: $opts"
    # shellcheck disable=SC2086  # the options are words
    taskset -c 0-63 $TV -o $O/$n.time $C map $IDX $S/R1.fq $S/R2.fq --intervals --iv-pair all -t 64 -o /dev/stdout \
        $opts 2> $O/$n.log | taskset -c 64-127 $ST view -@ 8 -u -L $T - \
      | taskset -c 64-127 $ST sort -@ 16 -m 2G -T $S/sort.$n -O cram --reference $REF -o $O/$n.cram -
    taskset -c 64-127 $ST index $O/$n.cram
    md5sum < $O/$n.cram > $O/$n.cram.md5
    cat $O/$n.time
  done
done
rm -f $V/subset_HG00?/R1.fq $V/subset_HG00?/R2.fq $V/subset_HG00?/single.fq
step "done"
