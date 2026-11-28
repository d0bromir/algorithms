#!/bin/bash
# Verify the CI speedup: rebuild test_certa with OpenMP, run seed 11 with 4
# threads (CI-runner-like), check correctness + wall time. Disk-safe (RAM build).
set -uo pipefail
export TMPDIR=/dev/shm/a2/tmp
cd /dev/shm/a2/src
git fetch --depth 1 origin claude/ci-speedup > /dev/shm/a2/ci_pull.log 2>&1 \
  && git reset --hard FETCH_HEAD >> /dev/shm/a2/ci_pull.log 2>&1 \
  || { echo PULL_FAIL; cat /dev/shm/a2/ci_pull.log; exit 1; }
cd certa
cmake -S . -B build -DCERTA_CUDA=OFF > /dev/shm/a2/ci_cmake.log 2>&1
cmake --build build -j 32 --target test_certa > /dev/shm/a2/ci_build.log 2>&1 \
  || { echo BUILD_FAIL; tail -30 /dev/shm/a2/ci_build.log; exit 1; }
echo BUILD_OK
grep -i "openmp" /dev/shm/a2/ci_cmake.log || echo "(OpenMP not found by cmake)"
echo "== seed 11, OMP_NUM_THREADS=4 (CI-like)"
t0=$(date +%s)
OMP_NUM_THREADS=4 ./build/test_certa 11 > /dev/shm/a2/ci_test11.log 2>&1
t1=$(date +%s)
echo "seed 11 wall: $((t1 - t0)) s"
tail -2 /dev/shm/a2/ci_test11.log
echo "FAIL count: $(grep -c FAIL /dev/shm/a2/ci_test11.log)"
echo "== seed 11, OMP_NUM_THREADS=1 (serial baseline)"
t0=$(date +%s)
OMP_NUM_THREADS=1 ./build/test_certa 11 > /dev/shm/a2/ci_test11s.log 2>&1
t1=$(date +%s)
echo "seed 11 serial wall: $((t1 - t0)) s"
echo "FAIL count: $(grep -c FAIL /dev/shm/a2/ci_test11s.log)"
