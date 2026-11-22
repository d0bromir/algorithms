#!/bin/bash
# Paired production: build + certificate tests, then paired --bio with/without
# the calibrated phantom prior (A5 frozen delta=0.006, pi=0.003, segdup track)
# on a 4M-pair HG002 subsample. Disk-safe: RAM, streamed SAM, no BAM.
set -uo pipefail
mkdir -p /dev/shm/a2/tmp
export TMPDIR=/dev/shm/a2/tmp
cd /dev/shm/a2/src
git fetch --depth 1 origin claude/paired-phantom > /dev/shm/pp_pull.log 2>&1 \
  && git reset --hard FETCH_HEAD >> /dev/shm/pp_pull.log 2>&1 \
  || { echo PULL_FAIL; cat /dev/shm/pp_pull.log; exit 1; }
cd certa
cmake --build build -j 32 --target test_certa verify_lemmas > /dev/shm/pp_build.log 2>&1 \
  || { echo BUILD_FAIL; tail -40 /dev/shm/pp_build.log; exit 1; }
echo BUILD_OK
./build/verify_lemmas > /dev/shm/pp_lemmas.log 2>&1; LEMMAS_EXIT=$?; echo LEMMAS_EXIT=$LEMMAS_EXIT
./build/test_certa > /dev/shm/pp_test.log 2>&1; CERTA_EXIT=$?; echo CERTA_EXIT=$CERTA_EXIT
echo ---PAIR-TEST---; grep bio-pair /dev/shm/pp_test.log
echo ---FAILCOUNT---; grep -c FAIL /dev/shm/pp_test.log
[ "$LEMMAS_EXIT" = 0 ] && [ "$CERTA_EXIT" = 0 ] || { echo TESTS_FAILED; tail -6 /dev/shm/pp_test.log; exit 1; }
cmake --build build -j 32 --target certa > /dev/shm/pp_cli.log 2>&1 \
  || { echo CLI_BUILD_FAIL; tail -30 /dev/shm/pp_cli.log; exit 1; }
echo CLI_OK
C=/dev/shm/a2/src/certa/build/certa
B=~/certa-bench; W=~/certa-truth/B; T=/dev/shm/a2/src/certa/bench
IDX=$B/idx/grch38.s1v3.cidx; FAI=$B/idx/ref.fa.fai
cd /dev/shm/a2
[ -f segdups.bed ] || zcat ~/certa-data/truth/strat/SegmentalDuplications/GRCh38_segdups.bed.gz > segdups.bed
[ -f sub_2.fq ] || head -n 16000000 "$W/B_2.fq" > sub_2.fq
echo "== $(date +%T) paired bio baseline (4M pairs)"
taskset -c 0-63 "$C" map "$IDX" sub_1.fq sub_2.fq --bio -t 64 -o /dev/stdout -u /dev/null -U /dev/null 2>p0.log \
  | python3 "$T/truth_eval.py" score "$W/truth.npz" "$FAI" p0.npz 2>p0.eval.log
echo "== $(date +%T) paired bio phantom (delta 0.006 pi 0.003, 4M pairs)"
taskset -c 0-63 "$C" map "$IDX" sub_1.fq sub_2.fq --bio --bio-phantom 0.006 --bio-phantom-bed segdups.bed \
  --bio-phantom-pi 0.003 --bio-phantom-pi-unique 0 -t 64 -o /dev/stdout -u /dev/null -U /dev/null 2>p1.log \
  | python3 "$T/truth_eval.py" score "$W/truth.npz" "$FAI" p1.npz 2>p1.eval.log
python3 "$T/truth_eval.py" report "$W/truth.npz" paired-bio=p0.npz paired-phantom=p1.npz > paired_report.md
echo DONE
sed -n '1,18p' paired_report.md
