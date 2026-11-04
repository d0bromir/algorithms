#!/bin/bash
# truth_headroom.sh: the runs of truth_complete.sh and phase 1's other
# aligners at matched sensitivity, at a 20 bp window and at the locus (150
# bp), with runs against other references on test B
# (bench/VALIDATION_RESULTS.md, "Matched sensitivity").
#   truth_headroom.sh COMPLETE_DIR OUT_DIR [REF_DIR]
# COMPLETE_DIR is truth_complete.sh's OUT_DIR. The 20 bp tables reuse its
# npz files and phase 1's; for 150 bp every run is rescored into OUT_DIR.
# REF_DIR holds test B against other references:
# - t2t.hs1.bam: bwa-mem2 mem -t 64 on T2T-CHM13 (hs1), lifted here through
#   hs1ToHg38.over.chain.gz (UCSC) by lift_sam.py;
# - giraffe.bam: vg giraffe 1.78.0 on HPRC v2.1 hprc-v2.1-mc-grch38-eval.d46
#   (HG002 not in the graph; its published distance index is format v3, so
#   rebuilt by vg autoindex -w sr-giraffe), -o SAM --ref-name GRCh38 -N HG002.
set -euo pipefail
CD=$(cd "$1" && pwd); O=$(mkdir -p "$2" && cd "$2" && pwd); RD=${3:+$(cd "$3" && pwd)}
W=~/certa-truth; FAI=~/certa-bench/idx/ref.fa.fai
T=$(cd "$(dirname "$0")" && pwd)
export PATH=~/tools/env/truth/bin:$PATH
AUX="taskset -c 64-127"
score() {  # window out.npz filter bam...: the union of the BAMs' records, through filter, scored
  local w=$1 out=$2 filter=$3; shift 3
  [ -s "$out" ] && return
  for f in "$@"; do $AUX samtools view "$f"; done | $filter \
    | $AUX python3 $T/truth_eval.py score --window $w $W/$X/truth.npz $FAI "$out" 2> "${out%.npz}.eval.log"
}
for X in A B; do
  echo "== $(date +%T) test $X: scoring"
  L=$O/$X/w150
  mkdir -p $L/truth $L/runs $O/$X/w20
  ln -sf $W/$X/truth.npz $L/truth/truth.npz
  pids=()
  for m in ivcert iv def bio; do score 150 $L/runs/$m.npz cat $CD/$X/$m.bam & pids+=($!); done
  score 150 $L/runs/def.complete.npz cat $CD/$X/def.bam $CD/$X/def.fb.bam & pids+=($!)
  score 150 $L/runs/bio.complete.npz cat $CD/$X/bio.bam $CD/$X/bio.fb.bam & pids+=($!)
  for a in minibwa bwa-mem2 strobealign minimap2; do score 150 $L/truth/$a.npz cat $W/$X/$a.bam & pids+=($!); done
  ext20=(); ext150=()
  if [ $X = B ] && [ -n "$RD" ]; then
    T2T="BWA-MEM2 on T2T-CHM13, lifted"; GIR="vg giraffe, HPRC v2.1 (HG002 held out)"
    for d in $O/B/w20:20 $L/runs:150; do
      score ${d#*:} ${d%:*}/t2t.npz "$AUX python3 $T/lift_sam.py $RD/hs1ToHg38.over.chain.gz" $RD/t2t.hs1.bam & pids+=($!)
      score ${d#*:} ${d%:*}/giraffe.npz "$AUX python3 $T/lift_sam.py --names" $RD/giraffe.bam & pids+=($!)
    done
    ext20=("$T2T=$O/B/w20/t2t.npz" "$GIR=$O/B/w20/giraffe.npz")
    ext150=("$T2T=$L/runs/t2t.npz" "$GIR=$L/runs/giraffe.npz")
  fi
  for p in "${pids[@]}"; do wait $p; done
  echo "== $(date +%T) test $X: tables"
  pids=()
  for r in ivcert:$CD/$X:$L/runs bio:$CD/$X:$L/runs minibwa:$W/$X:$L/truth bwa-mem2:$W/$X:$L/truth \
           strobealign:$W/$X:$L/truth minimap2:$W/$X:$L/truth; do
    IFS=: read -r m d20 d150 <<< "$r"  # the BAM and its 20 bp scores are in d20
    $AUX samtools view $d20/$m.bam | $AUX python3 $T/truth_nearmiss.py $W/$X/truth.npz $d20/$m.npz $d150/$m.npz > $O/$X/nearmiss.$m.txt &
    pids+=($!)
  done
  for p in "${pids[@]}"; do wait $p; done
  {
    echo "## Test $X, window 20 bp"; echo
    $AUX python3 $T/truth_headroom.py $W/$X $CD/$X "${ext20[@]}"
    echo; echo "## Test $X, window 150 bp (the locus)"; echo
    $AUX python3 $T/truth_headroom.py $L/truth $L/runs "${ext150[@]}"
    echo; echo "## Test $X, near misses (truth_nearmiss.py)"
    for m in ivcert bio minibwa bwa-mem2 strobealign minimap2; do echo; echo "**$m**: $(cat $O/$X/nearmiss.$m.txt)"; done
  } > $O/headroom_$X.md
done
echo "== $(date +%T) finished"
