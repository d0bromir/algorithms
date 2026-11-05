#!/bin/bash
# a2: vc_opts.sh GENOME TAG OPTIONS...: DeepVariant on CERTA with other map
# options. GENOME's subset of read pairs (vc_subset.sh; `vc_subset.sh GENOME`
# with no mode only extracts them) is mapped with --intervals --iv-pair all
# -t 64 OPTIONS by $CERTA (default: the frozen PR 19 build), timed into
# GENOME/S_TAG.time (vc_table.py reads it); the records overlapping the
# targets go to call.sh as GENOME/S_TAG.bam, and lost_reads.py compares all
# records (kept as BAM while it runs) with minibwa's (lost_TAG.txt).
set -euo pipefail
G=$1; TAG=$2; shift 2
V=~/vc; O=$V/$G; S=$V/subset_$G
ST=$V/dv/rootfs/opt/conda/envs/bio/bin/samtools
IDX=~/certa-bench/idx/grch38.s1v3.cidx; T=$V/targets.bed; C=${CERTA:-~/certa-pm/build/certa}
TV="/usr/bin/time -f %e\t%U\t%S\t%M"
step() { echo "== $(date +%T) $*"; }
[ -s $S/R2.fq ] || { echo "no subset reads in $S"; exit 2; }
n=S_$TAG
step "map $n with $C: $*"
$TV -o $O/$n.time $C map $IDX $S/R1.fq $S/R2.fq --intervals --iv-pair all -t 64 -o /dev/stdout "$@" 2> $S/$n.log \
  | $ST view -@ 8 -b -o $S/$n.all.bam -
$ST view -u -L $T $S/$n.all.bam | $ST sort -@ 8 -m 1500M -T $S/sort.$n -o $O/$n.bam -
$ST index $O/$n.bam
python3 $V/scripts/lost_reads.py <($ST view -F 0x900 $O/minibwa.cram) <($ST view $S/$n.all.bam) > $S/lost_$TAG.txt 2>&1
rm -f $S/$n.all.bam
grep "^| all" $S/lost_$TAG.txt
cd $V && scripts/call.sh $G $n
step "all finished"
