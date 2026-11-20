#!/bin/bash
# a2: bwscale.sh GENOME IN OUT TABLE: a copy of ~/vc/GENOME/IN.bam (or .cram)
# as OUT.bam with every record's MAPQ m replaced by the m-th value of TABLE
# (comma-separated, from MAPQ 0; the last value for larger m), for call.sh.
# The diagnostic of VALIDATION_RESULTS.md "The fitted MAPQ on BWA's scale":
# the fitted MAPQ (a Phred) rewritten on BWA-MEM's scale. The code's table
# is kBwaScale in include/certa/pair_model.h (--mapq-scale bwa):
#   0,1,4,10,16,22,27,32,38,43,45,46,48,49,50,51,52,53,54,54,55,55,55,55,
#   56,56,56,56,57,57,57,57,58,58,58,58,59,59,59,59,60
set -euo pipefail
G=$1; IN=$2; OUT=$3; TAB=$4
V=~/vc; O=$V/$G; ST=$V/dv/rootfs/opt/conda/envs/bio/bin/samtools; REF=~/certa-bench/idx/ref.fa
SRC=$O/$IN.bam; [ -s $SRC ] || SRC=$O/$IN.cram
$ST view -@8 -h -T $REF $SRC \
  | awk -v tab="$TAB" 'BEGIN{FS=OFS="\t"; n=split(tab,t,",")} /^@/{print; next} {m=$5+1; if (m>n) m=n; $5=t[m]; print}' \
  | $ST view -@8 -b -o $O/$OUT.bam -
$ST index -@8 $O/$OUT.bam
[ -e $O/$IN.time ] && cp $O/$IN.time $O/$OUT.time || true
echo "== $(date +%T) $OUT written"
