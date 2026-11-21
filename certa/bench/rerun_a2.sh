#!/bin/bash
# A2 re-run after the cli.h flag fix: pull, rebuild certa, re-run only the
# --mapq-local run (baseline r0.npz is unchanged by the flag-parsing fix).
set -uo pipefail
export TMPDIR=/dev/shm/a2/tmp
cd /dev/shm/a2/src
git fetch --depth 1 origin claude/a2-mapq-local > /dev/shm/a2/pull.log 2>&1 \
  && git reset --hard FETCH_HEAD >> /dev/shm/a2/pull.log 2>&1 \
  || { echo PULL_FAIL; cat /dev/shm/a2/pull.log; exit 1; }
cd certa
cmake --build build -j 32 --target certa > /dev/shm/a2/build_cli2.log 2>&1 \
  || { echo BUILD_FAIL; tail -30 /dev/shm/a2/build_cli2.log; exit 1; }
echo CLI_OK
C=/dev/shm/a2/src/certa/build/certa
B=~/certa-bench; W=~/certa-truth/B; T=/dev/shm/a2/src/certa/bench
IDX=$B/idx/grch38.s1v3.cidx; FAI=$B/idx/ref.fa.fai
cd /dev/shm/a2
echo "== $(date +%T) --mapq-local (rerun)"
taskset -c 0-63 "$C" map "$IDX" "$W/B_1.fq" "$W/B_2.fq" --mapq-local -t 64 -o /dev/stdout -u /dev/null -U /dev/null 2>r1.log \
  | python3 "$T/truth_eval.py" score "$W/truth.npz" "$FAI" r1.npz 2>r1.eval.log
echo "== $(date +%T) report"
python3 "$T/truth_eval.py" report "$W/truth.npz" baseline=r0.npz mapq-local=r1.npz > a2_report.md
echo DONE
sed -n '1,30p' a2_report.md
