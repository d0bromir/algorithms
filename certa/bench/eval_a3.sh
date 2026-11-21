#!/bin/bash
# A3 accuracy: single-end --bio with and without --bio-phantom on a 4M-read
# subsample of HG002 test B (mate 0), scored against truth. The phantom term
# is in bio_place's single-end posterior. Disk-safe: subsample + SAM stream +
# npz/reports all in /dev/shm (RAM); no BAM.
set -uo pipefail
mkdir -p /dev/shm/a2/tmp
export TMPDIR=/dev/shm/a2/tmp
cd /dev/shm/a2/src/certa
cmake --build build -j 32 --target certa > /dev/shm/a3_cli.log 2>&1 \
  || { echo CLI_BUILD_FAIL; tail -30 /dev/shm/a3_cli.log; exit 1; }
echo CLI_OK
C=/dev/shm/a2/src/certa/build/certa
B=~/certa-bench; W=~/certa-truth/B; T=/dev/shm/a2/src/certa/bench
IDX=$B/idx/grch38.s1v3.cidx; FAI=$B/idx/ref.fa.fai
cd /dev/shm/a2
head -n 16000000 "$W/B_1.fq" > /dev/shm/a2/sub_1.fq  # 4M reads (mate 0 of the first 4M pairs)
echo "== $(date +%T) bio baseline (single-end, 4M reads)"
taskset -c 0-63 "$C" map "$IDX" sub_1.fq --bio -t 64 -o /dev/stdout -u /dev/null 2>b0.log \
  | python3 "$T/truth_eval.py" score "$W/truth.npz" "$FAI" b0.npz 2>b0.eval.log
echo "== $(date +%T) bio --bio-phantom 0.01 (single-end, 4M reads)"
taskset -c 0-63 "$C" map "$IDX" sub_1.fq --bio --bio-phantom 0.01 -t 64 -o /dev/stdout -u /dev/null 2>b1.log \
  | python3 "$T/truth_eval.py" score "$W/truth.npz" "$FAI" b1.npz 2>b1.eval.log
echo "== $(date +%T) report"
python3 "$T/truth_eval.py" report "$W/truth.npz" bio=b0.npz bio-phantom=b1.npz > a3_report.md
echo DONE
sed -n '1,32p' a3_report.md
