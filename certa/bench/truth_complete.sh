#!/bin/bash
# truth_complete.sh: CERTA's complete pipelines on tests A and B, scored
# against phase 1's truth next to phase 1's runs of the other aligners
# (bench/VALIDATION_RESULTS.md, "Complete pipelines"). The default mode and
# --bio write their uncertified pairs to FASTQ, which minibwa maps as in
# final_bench.sh (-K 10m; drop_certified_mates.awk keeps CERTA's record of a
# certified mate). Interval mode maps every read itself, with the default
# MAPQ, --iv-mapq cons and --iv-mapq cert.
#   truth_complete.sh BUILD_DIR OUT_DIR [TEST...]   (tests: A B, the default)
set -euo pipefail
BUILD=$1; O=$(mkdir -p "$2" && cd "$2" && pwd); shift 2
W=~/certa-truth; B=~/certa-bench
C=$BUILD/certa; T=$(cd "$(dirname "$0")" && pwd)
MB=~/micromamba/envs/certa-bench/bin/minibwa
IDX=$B/idx/grch38.s1v3.cidx; FAI=$B/idx/ref.fa.fai
export PATH=~/tools/env/truth/bin:$PATH
PIN="taskset -c 0-63"; AUX="taskset -c 64-127"
TV="/usr/bin/time -f %e\t%U\t%S\t%M"
bam() { $AUX samtools view -@ 16 -b --output-fmt-option level=1 -o "$1" -; }
score() {  # out.npz bam...: the union of the BAMs' records, scored
  local out=$1; shift
  [ -s "$out" ] && return
  for f in "$@"; do $AUX samtools view "$f"; done \
    | $AUX python3 $T/truth_eval.py score $W/$X/truth.npz $FAI "$out" 2> "${out%.npz}.eval.log"
}
for X in ${*:-A B}; do
  mkdir -p $O/$X; cd $O/$X
  R1=$W/$X/${X}_1.fq; R2=$W/$X/${X}_2.fq
  for m in def iv ivc ivcert bio; do
    [ -s $m.bam ] && continue
    case $m in
      def) opt="-u $m.u1.fq -U $m.u2.fq" ;;
      iv)  opt="--intervals" ;;
      ivc) opt="--intervals --iv-mapq cons" ;;
      ivcert) opt="--intervals --iv-mapq cert" ;;
      bio) opt="--bio --bio-insert 429,99 -u $m.u1.fq -U $m.u2.fq" ;;
    esac
    echo "== $(date +%T) test $X: certa $m"
    $PIN $TV -o $m.time $C map $IDX $R1 $R2 $opt -t 64 -o /dev/stdout --stats $m.json 2> $m.log \
      | bam $m.bam.tmp
    mv $m.bam.tmp $m.bam
  done
  for m in def bio; do
    [ -s $m.fb.bam ] && continue
    echo "== $(date +%T) test $X: minibwa on certa $m's fallback, $(( $(wc -l < $m.u1.fq) / 4 )) pairs"
    $PIN $TV -o $m.fb.time $MB map -t 64 -K 10m $B/idx/minibwa $m.u1.fq $m.u2.fq 2> $m.fb.log \
      | awk -f $T/drop_certified_mates.awk | bam $m.fb.bam.tmp
    mv $m.fb.bam.tmp $m.fb.bam
  done
  echo "== $(date +%T) test $X: identity"
  for m in def bio; do md5sum < $m.u1.fq > $m.u1.md5; md5sum < $m.u2.fq > $m.u2.md5; done
  cmp -s <($AUX samtools view def.bam | md5sum) <($AUX samtools view $W/$X/certa.bam | md5sum) \
    && echo "test $X: default records identical to phase 1's" || echo "test $X: default records DIFFER from phase 1's"
  mates() { $AUX samtools view -F 2304 "$1" | awk '{print $1 "\t" and($2, 192)}' | LC_ALL=C sort -u; }
  echo "test $X: mates with a primary record in both CERTA's SAM and the fallback's:" \
    "$(LC_ALL=C comm -12 <(mates def.bam) <(mates def.fb.bam) | wc -l) (default)," \
    "$(LC_ALL=C comm -12 <(mates bio.bam) <(mates bio.fb.bam) | wc -l) (--bio)"
  echo "== $(date +%T) test $X: scoring"
  pids=()
  score def.npz def.bam & pids+=($!)
  score def.complete.npz def.bam def.fb.bam & pids+=($!)
  score iv.npz iv.bam & pids+=($!)
  score ivc.npz ivc.bam & pids+=($!)
  score ivcert.npz ivcert.bam & pids+=($!)
  score bio.npz bio.bam & pids+=($!)
  score bio.complete.npz bio.bam bio.fb.bam & pids+=($!)
  for p in "${pids[@]}"; do wait $p; done
  echo "== $(date +%T) test $X: report"
  runs=("certa --bio + minibwa=bio.complete.npz" "certa + minibwa=def.complete.npz" "certa interval mode=iv.npz"
        "certa interval mode, --iv-mapq cons=ivc.npz"
        "certa interval mode, --iv-mapq cert=ivcert.npz" "certa --bio=bio.npz" "certa=def.npz")
  scores=()
  for a in minibwa bwa-mem2 strobealign minimap2; do runs+=("$a=$W/$X/$a.npz"); done
  for a in certa minibwa bwa-mem2 strobealign minimap2; do scores+=("$a=$W/$X/$a.score.tsv"); done
  $AUX python3 $T/truth_eval.py report $W/$X/truth.npz "${runs[@]}" --scores "${scores[@]}" > report.md
  $AUX python3 $T/truth_caps.py $W/$X . > caps.md
  for f in def iv ivc ivcert bio def.fb bio.fb; do printf "%s\t%s\t%s\n" $X $f "$(cat $f.time)"; done | tee -a $O/timings.tsv
done
echo "== $(date +%T) finished"
