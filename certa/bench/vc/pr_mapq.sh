#!/bin/bash
# a2: pr_mapq.sh BUILD_DIR: PR mates' point MAPQ from the best other pair
# found (VALIDATION_RESULTS.md, "The MAPQ of mates placed by the pair").
# BUILD_DIR holds the source (include/, src/), built here with
# mate_rescue.sh's compile lines. For HG002 and HG005, the subset of read
# pairs is mapped by vc_opts.sh with CERTA=BUILD_DIR/build/certa and
# --iv-rescue local (S_prq), as S_mr3_local was by the build before, and
# called; same_records.py then checks that only MAPQ and XQ differ from
# S_mr3_local. HG005's subset reads are extracted first and removed after.
# Run from ~/vc/scripts, where the scripts of this folder are installed.
set -uo pipefail
P=$(cd "$1" && pwd); export CERTA=$P/build/certa
V=~/vc; ST=$V/dv/rootfs/opt/conda/envs/bio/bin/samtools; H=$V/scripts
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
  $H/vc_opts.sh $G prq --iv-rescue local || step "$G prq failed"
  python3 $H/same_records.py <($ST view $V/$G/S_mr3_local.bam) <($ST view $V/$G/S_prq.bam) || step "$G records differ"
done
rm -f $V/subset_HG005/R1.fq $V/subset_HG005/R2.fq $V/subset_HG005/single.fq
step "done"
