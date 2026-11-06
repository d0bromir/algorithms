#!/bin/bash
# a2: vc_opts.sh GENOME TAG OPTIONS...: DeepVariant on CERTA with other map
# options (the frozen PR 19 build). GENOME's subset of read pairs
# (vc_subset.sh) is mapped with --intervals --iv-pair all -t 64 OPTIONS; the
# records overlapping the targets go to call.sh as GENOME/S_TAG.bam, and
# lost_reads.py compares the full output with minibwa's (lost_TAG.txt).
set -euo pipefail
G=$1; TAG=$2; shift 2
V=~/vc; O=$V/$G; S=$V/subset_$G
ST=$V/dv/rootfs/opt/conda/envs/bio/bin/samtools
IDX=~/certa-bench/idx/grch38.s1v3.cidx; T=$V/targets.bed; C=~/certa-pm/build/certa
step() { echo "== $(date +%T) $*"; }
[ -s $S/R2.fq ] || { echo "no subset reads in $S"; exit 2; }
n=S_$TAG
step "map $n: $*"
$C map $IDX $S/R1.fq $S/R2.fq --intervals --iv-pair all -t 64 -o $S/$n.all.sam "$@" 2> $S/$n.log
$ST view -u -L $T $S/$n.all.sam | $ST sort -@ 8 -m 1500M -T $S/sort.$n -o $O/$n.bam -
$ST index $O/$n.bam
$ST view -F 0x900 $O/minibwa.cram > $S/mb.$TAG.sam
python3 $V/scripts/lost_reads.py $S/mb.$TAG.sam $S/$n.all.sam > $S/lost_$TAG.txt 2>&1
rm -f $S/mb.$TAG.sam $S/$n.all.sam
grep "^| all" $S/lost_$TAG.txt
cd $V && scripts/call.sh $G $n
step "all finished"
