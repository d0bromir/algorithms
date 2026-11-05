#!/bin/bash
# a2: vc_model.sh GENOME TAG MODEL.h: DeepVariant on another --mapq pair
# model. The frozen PR 19 source (~/certa-pm) with include/certa/mapq_model.h
# replaced by MODEL.h is built as ~/certa-TAG; GENOME's subset of read pairs
# (vc_subset.sh) is mapped with --iv-mapq cert --iv-pair all --mapq pair, and
# the records overlapping the targets go to call.sh as GENOME/S_pair_TAG.bam.
set -euo pipefail
G=$1; TAG=$2; MODEL=$3
V=~/vc; O=$V/$G; S=$V/subset_$G
ST=$V/dv/rootfs/opt/conda/envs/bio/bin/samtools
IDX=~/certa-bench/idx/grch38.s1v3.cidx; T=$V/targets.bed; P=~/certa-$TAG; C=$P/build/certa
step() { echo "== $(date +%T) $*"; }
[ -s $S/R2.fq ] || { echo "no subset reads in $S"; exit 2; }
if [ ! -x $C ]; then
  rm -rf $P; mkdir -p $P
  cp -r ~/certa-pm/include ~/certa-pm/src $P/
  cp $MODEL $P/include/certa/mapq_model.h
  mkdir -p $P/build/obj
  for f in index seqio mapper_cpu check chain main cmd_map map_config map_pipeline map_stats sam_output gpu_stub; do
    echo "g++ -O3 -DNDEBUG -march=native -std=c++17 -I$P/include -I$P/src -c $P/src/$f.cpp -o $P/build/obj/$f.o"
  done | xargs -P 12 -I{} sh -c {}
  g++ -O3 -march=native $P/build/obj/*.o -lpthread -o $C
fi
step "build $TAG: $(cd $P && find include src -type f | LC_ALL=C sort | xargs sha256sum | sha256sum | cut -c1-16), model $(md5sum < $MODEL | cut -c1-12)"
n=S_pair_$TAG
if [ ! -s $O/$n.bam.bai ]; then
  step "map $n"
  $C map $IDX $S/R1.fq $S/R2.fq --intervals --iv-pair all -t 64 -o /dev/stdout --iv-mapq cert --mapq pair 2> $S/$n.log \
    | $ST view -u -L $T - | $ST sort -@ 8 -m 1500M -T $S/sort.$n -o $O/$n.bam -
  $ST index $O/$n.bam
fi
cd $V && scripts/call.sh $G $n
step "all finished"
