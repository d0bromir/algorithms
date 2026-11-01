#!/bin/bash
# truth_rerun_bio.sh: tests A and B again for a new --bio build only, scored
# against phase 1's truth and reported next to phase 1's runs of every
# aligner (bench/VALIDATION_RESULTS.md, "After phase 1"). Also checks that
# bio_score reproduces YP on every --bio record.
#   truth_rerun_bio.sh BUILD_DIR OUT_DIR [TEST...]   (tests: A B, the default)
set -eu
BUILD=$1; O=$2; shift 2
W=~/certa-truth; B=~/certa-bench
C=$BUILD/certa; S=$BUILD/bio_score; T=$(cd "$(dirname "$0")" && pwd)
IDX=$B/idx/grch38.s1v3.cidx; FAI=$B/idx/ref.fa.fai
export PATH=~/tools/env/truth/bin:$PATH
PIN="taskset -c 0-63"; AUX="taskset -c 64-127"
for X in ${*:-A B}; do
  mkdir -p $O/$X
  if [ ! -s $O/$X/bio.bam ]; then
    echo "== $(date +%T) test $X: certa --bio"
    $PIN /usr/bin/time -f "%e\t%U\t%S\t%M" -o $O/$X/bio.time $C map $IDX $W/$X/${X}_1.fq $W/$X/${X}_2.fq --bio \
      --bio-insert 429,99 -t 64 -o /dev/stdout -u /dev/null -U /dev/null --stats $O/$X/bio.json 2> $O/$X/bio.log \
      | $AUX samtools view -@ 16 -b --output-fmt-option level=1 -o $O/$X/bio.bam -
  fi
  echo "== $(date +%T) test $X: scoring"
  $AUX samtools view $O/$X/bio.bam | $AUX $S $IDX > $O/$X/bio.score.tsv 2> $O/$X/bio.score.log
  $AUX samtools view -F 2308 $O/$X/bio.bam | $AUX awk -F '\t' '{m = and($2, 64) ? 1 : (and($2, 128) ? 2 : 0);
    yp = -1; yr = 0; for (i = 12; i <= NF; i++) { if ($i ~ /^YP:i:/) yp = substr($i, 6); if ($i == "YR:i:1") yr = 1 }
    print $1 "\t" m "\t" yp "\t" yr}' > $O/$X/bio.yp
  $AUX samtools view $O/$X/bio.bam | $AUX python3 $T/truth_eval.py score $W/$X/truth.npz $FAI $O/$X/bio.npz \
    2> $O/$X/bio.eval.log
  echo "== $(date +%T) test $X: report"
  runs=("certa --bio (new)=$O/$X/bio.npz" "certa --bio=$W/$X/certa-bio.npz"); scores=()
  for a in certa minibwa bwa-mem2 strobealign minimap2; do
    runs+=("$a=$W/$X/$a.npz"); scores+=("$a=$W/$X/$a.score.tsv")
  done
  $AUX python3 $T/truth_eval.py report $W/$X/truth.npz "${runs[@]}" --bio "certa --bio (new)" --scores "${scores[@]}" \
    > $O/$X/report.md
  $AUX python3 - $O/$X/bio.yp $O/$X/bio.score.tsv <<'PY'
import sys
yp = {}
for line in open(sys.argv[1]):
    n, m, p, r = line.rstrip("\n").split("\t")
    if p != "-1":
        yp[(n, m)] = (int(p), r == "1")
same = diff = 0
for line in open(sys.argv[2]):
    n, m, v = line.rstrip("\n").split("\t")
    if (n, m) in yp:
        same += int(v) == yp[(n, m)][0]
        diff += int(v) != yp[(n, m)][0]
print(f"YP check: {len(yp)} --bio records ({sum(r for _, r in yp.values())} placed only by the pair): "
      f"bio_score equal on {same}, different on {diff}")
PY
done
echo "== $(date +%T) finished"
