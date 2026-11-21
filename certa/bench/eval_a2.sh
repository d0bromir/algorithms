#!/bin/bash
# A2 accuracy: default-mode CPU runs with and without --mapq-local on HG002
# test B, scored against truth. The certificate is unchanged, so the placed
# set is identical; only MAPQ differs. Disk-safe: SAM streamed into the scorer
# (no BAM), npz/reports in /dev/shm.
set -uo pipefail
mkdir -p /dev/shm/a2/tmp
export TMPDIR=/dev/shm/a2/tmp  # compiler temporaries in RAM, not on the full disk
cd /dev/shm/a2/src/certa
cmake --build build -j 32 --target certa > /dev/shm/a2/build_cli.log 2>&1 \
  || { echo CLI_BUILD_FAIL; tail -30 /dev/shm/a2/build_cli.log; exit 1; }
echo CLI_OK
C=/dev/shm/a2/src/certa/build/certa
B=~/certa-bench; W=~/certa-truth/B; T=/dev/shm/a2/src/certa/bench
IDX=$B/idx/grch38.s1v3.cidx; FAI=$B/idx/ref.fa.fai
cd /dev/shm/a2
echo "== $(date +%T) baseline (no flag)"
taskset -c 0-63 "$C" map "$IDX" "$W/B_1.fq" "$W/B_2.fq" -t 64 -o /dev/stdout -u /dev/null -U /dev/null 2>r0.log \
  | python3 "$T/truth_eval.py" score "$W/truth.npz" "$FAI" r0.npz 2>r0.eval.log
echo "== $(date +%T) --mapq-local"
taskset -c 0-63 "$C" map "$IDX" "$W/B_1.fq" "$W/B_2.fq" --mapq-local -t 64 -o /dev/stdout -u /dev/null -U /dev/null 2>r1.log \
  | python3 "$T/truth_eval.py" score "$W/truth.npz" "$FAI" r1.npz 2>r1.eval.log
echo "== $(date +%T) report"
python3 "$T/truth_eval.py" report "$W/truth.npz" baseline=r0.npz mapq-local=r1.npz > a2_report.md
echo DONE
sed -n '1,40p' a2_report.md
