#!/bin/bash
# B1-impl: build + certificate tests, then paired interval --t2t with/without
# --t2t-commit (selective lift-and-commit) on a 4M-pair HG002 subsample. The
# commit should rescue missing-copy placements (raise at-truth %) and cut
# confident errors. Disk-safe: RAM, streamed SAM, no BAM.
set -uo pipefail
mkdir -p /dev/shm/a2/tmp
export TMPDIR=/dev/shm/a2/tmp
cd /dev/shm/a2/src
git fetch --depth 1 origin claude/b1-lift-commit > /dev/shm/b1i_pull.log 2>&1 \
  && git reset --hard FETCH_HEAD >> /dev/shm/b1i_pull.log 2>&1 \
  || { echo PULL_FAIL; cat /dev/shm/b1i_pull.log; exit 1; }
cd certa
cmake --build build -j 32 --target test_certa verify_lemmas certa > /dev/shm/b1i_build.log 2>&1 \
  || { echo BUILD_FAIL; tail -40 /dev/shm/b1i_build.log; exit 1; }
echo BUILD_OK
./build/verify_lemmas > /dev/shm/b1i_lemmas.log 2>&1; LEMMAS_EXIT=$?; echo LEMMAS_EXIT=$LEMMAS_EXIT
./build/test_certa > /dev/shm/b1i_test.log 2>&1; CERTA_EXIT=$?; echo CERTA_EXIT=$CERTA_EXIT
echo ---FAILCOUNT---; grep -c FAIL /dev/shm/b1i_test.log
[ "$LEMMAS_EXIT" = 0 ] && [ "$CERTA_EXIT" = 0 ] || { echo TESTS_FAILED; tail -6 /dev/shm/b1i_test.log; exit 1; }
C=/dev/shm/a2/src/certa/build/certa
B=~/certa-bench; W=~/certa-truth/B; T=/dev/shm/a2/src/certa/bench
IDX=$B/idx/grch38.s1v3.cidx; FAI=$B/idx/ref.fa.fai
HS1=~/certa-t2t/hs1.cidx; CHAIN=~/headroom/hs1ToHg38.over.chain.gz
cd /dev/shm/a2
[ -f sub_2.fq ] || head -n 16000000 "$W/B_2.fq" > sub_2.fq
echo "== $(date +%T) interval --t2t baseline (4M pairs)"
taskset -c 0-63 "$C" map "$IDX" sub_1.fq sub_2.fq --intervals --iv-pair all --iv-rescue local --mapq pair \
  --t2t "$HS1" --t2t-chain "$CHAIN" -t 64 -o /dev/stdout 2>c0.log \
  | python3 "$T/truth_eval.py" score "$W/truth.npz" "$FAI" c0.npz 2>c0.eval.log
echo "== $(date +%T) interval --t2t --t2t-commit (4M pairs)"
taskset -c 0-63 "$C" map "$IDX" sub_1.fq sub_2.fq --intervals --iv-pair all --iv-rescue local --mapq pair \
  --t2t "$HS1" --t2t-chain "$CHAIN" --t2t-commit -t 64 -o /dev/stdout 2>c1.log \
  | python3 "$T/truth_eval.py" score "$W/truth.npz" "$FAI" c1.npz 2>c1.eval.log
python3 "$T/truth_eval.py" report "$W/truth.npz" t2t=c0.npz t2t-commit=c1.npz > b1impl_report.md
echo DONE
sed -n '1,18p' b1impl_report.md
