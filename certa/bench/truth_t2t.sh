#!/bin/bash
# truth_t2t.sh: a MAPQ informed by T2T-CHM13 on tests A and B
# (bench/VALIDATION_RESULTS.md, "A MAPQ from T2T-CHM13"). BWA-MEM2 maps each
# test to hs1 unless OUT_DIR/X/t2t.hs1.bam exists; truth_t2t.py lifts its
# records to GRCh38, keeps every run's placements and scores the rule.
#   truth_t2t.sh COMPLETE_DIR OUT_DIR HS1_FA CHAIN
# COMPLETE_DIR is truth_complete.sh's OUT_DIR; HS1_FA is indexed for
# bwa-mem2; CHAIN is UCSC's hs1ToHg38.over.chain.gz.
set -euo pipefail
CD=$(cd "$1" && pwd); O=$(mkdir -p "$2" && cd "$2" && pwd); HS1=$3; CHAIN=$(cd "$(dirname "$4")" && pwd)/$(basename "$4")
W=~/certa-truth; FAI=~/certa-bench/idx/ref.fa.fai
T=$(cd "$(dirname "$0")" && pwd)
BWA=~/micromamba/envs/refmappers/bin/bwa-mem2
export PATH=~/tools/env/truth/bin:$PATH
PIN="taskset -c 0-63"; AUX="taskset -c 64-127"
TV="/usr/bin/time -f %e\t%U\t%S\t%M"
place() {  # out.npz bam...: the union of the BAMs' primary placements
  local out=$1; shift
  [ -s "$out" ] && return
  for f in "$@"; do $AUX samtools view "$f"; done | $AUX python3 $T/truth_t2t.py place $W/$X/truth.npz $FAI "$out"
}
for X in A B; do
  mkdir -p $O/$X; cd $O/$X
  if [ ! -s t2t.hs1.bam ]; then
    echo "== $(date +%T) test $X: bwa-mem2 on hs1"
    $PIN $TV -o t2t.time $BWA mem -t 64 $HS1 $W/$X/${X}_1.fq $W/$X/${X}_2.fq 2> t2t.bwa.log \
      | $AUX samtools view -@ 8 -b --output-fmt-option level=1 -o t2t.hs1.bam.tmp -
    mv t2t.hs1.bam.tmp t2t.hs1.bam
  fi
  echo "== $(date +%T) test $X: placements"
  pids=()
  [ -s t2t.npz ] || ( $AUX samtools view t2t.hs1.bam \
    | $AUX python3 $T/truth_t2t.py place --chain $CHAIN $W/$X/truth.npz $FAI t2t.npz 2> t2t.place.log ) & pids+=($!)
  for m in ivcert def bio; do place $m.npz $CD/$X/$m.bam & pids+=($!); done
  place def.complete.npz $CD/$X/def.bam $CD/$X/def.fb.bam & pids+=($!)
  place bio.complete.npz $CD/$X/bio.bam $CD/$X/bio.fb.bam & pids+=($!)
  for a in minibwa bwa-mem2 strobealign minimap2; do place $a.npz $W/$X/$a.bam & pids+=($!); done
  for p in "${pids[@]}"; do wait $p; done
  echo "== $(date +%T) test $X: tables"
  {
    echo "## Test $X: BWA-MEM2 on T2T-CHM13 as the second opinion"; echo
    $AUX python3 $T/truth_t2t.py eval $W/$X/truth.npz .
    echo; echo "## Test $X, control: BWA-MEM2 on GRCh38 as the second opinion"; echo
    $AUX python3 $T/truth_t2t.py eval --second bwa-mem2 $W/$X/truth.npz .
  } > $O/t2t_$X.md
done
echo "== $(date +%T) finished"
