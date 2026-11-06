#!/bin/bash
# a2: mate_rescue.sh BUILD_DIR: --iv-rescue local in variant calling
# (VALIDATION_RESULTS.md, "A local mate rescue"). BUILD_DIR holds the source
# (include/, src/), built here with vc_subset.sh's compile lines (a2 has no
# cmake). For HG002 and HG005, the subset of read pairs is mapped by
# vc_opts.sh with CERTA=BUILD_DIR/build/certa, without the local rescue
# (S_e2e) and with it (S_local), and called. On HG002 also the swap
# diagnostic (S_swap, swap_lost.py): CERTA's records with the local rescue,
# minibwa's for the reads CERTA still loses. HG005's subset reads are
# extracted first (vc_subset.sh with no mode) and removed after (disk).
# Run from ~/vc/scripts, where the scripts of this folder are installed.
set -uo pipefail
P=$(cd "$1" && pwd); export CERTA=$P/build/certa
V=~/vc; ST=$V/dv/rootfs/opt/conda/envs/bio/bin/samtools; T=$V/targets.bed
IDX=~/certa-bench/idx/grch38.s1v3.cidx; H=$V/scripts
step() { echo "== $(date +%T) $*"; }
step "build $(cd $P && find include src -type f | LC_ALL=C sort | xargs sha256sum | sha256sum | cut -c1-16)"
if [ ! -x $CERTA ]; then
  mkdir -p $P/build/obj
  for f in index seqio mapper_cpu check chain main cmd_map map_config map_pipeline map_stats sam_output gpu_stub; do
    echo "g++ -O3 -DNDEBUG -march=native -std=c++17 -I$P/include -I$P/src -c $P/src/$f.cpp -o $P/build/obj/$f.o"
  done | xargs -P 12 -I{} sh -c {}
  g++ -O3 -march=native $P/build/obj/*.o -lpthread -o $CERTA
fi
for G in HG002 HG005; do
  if [ ! -s $V/subset_$G/R2.fq ]; then
    step "$G subset reads"
    $H/vc_subset.sh $G || { step "$G reads failed"; continue; }
  fi
  for m in e2e local; do
    $H/vc_opts.sh $G $m --iv-rescue $m || step "$G $m failed"
  done
  if [ $G = HG002 ]; then
    step "HG002 swap"
    S=$V/subset_HG002; O=$V/HG002
    python3 $H/swap_lost.py <($ST view -F 0x900 $O/minibwa.cram) \
        <($CERTA map $IDX $S/R1.fq $S/R2.fq --intervals --iv-pair all -t 64 --iv-rescue local -o /dev/stdout \
          2> $S/S_swap.log) 2> $S/S_swap.swap.log | $ST view -u -L $T - | $ST sort -@ 8 -m 1500M -T $S/sort.swap \
        -o $O/S_swap.bam -
    $ST index $O/S_swap.bam
    (cd $V && scripts/call.sh HG002 S_swap) || step "HG002 swap failed"
  fi
done
rm -f $V/subset_HG005/R1.fq $V/subset_HG005/R2.fq $V/subset_HG005/single.fq
step "done"
