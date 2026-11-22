#!/bin/bash
# Build A4 (--bio-phantom-bed) from the pushed branch in RAM and run the
# certificate tests (incl. RegionTrack + locus-prior unit tests). Disk-safe.
set -uo pipefail
mkdir -p /dev/shm/a2/tmp
export TMPDIR=/dev/shm/a2/tmp
cd /dev/shm/a2/src
git fetch --depth 1 origin claude/a4-locus-prior > /dev/shm/a4_pull.log 2>&1 \
  && git reset --hard FETCH_HEAD >> /dev/shm/a4_pull.log 2>&1 \
  || { echo PULL_FAIL; cat /dev/shm/a4_pull.log; exit 1; }
cd certa
cmake --build build -j 32 --target test_certa verify_lemmas > /dev/shm/a4_build.log 2>&1 \
  || { echo BUILD_FAIL; tail -40 /dev/shm/a4_build.log; exit 1; }
echo BUILD_OK
./build/verify_lemmas > /dev/shm/a4_lemmas.log 2>&1; echo LEMMAS_EXIT=$?
./build/test_certa > /dev/shm/a4_test.log 2>&1; echo CERTA_EXIT=$?
echo ---TAIL---; tail -4 /dev/shm/a4_test.log
echo ---A4---; grep bio-phantom-bed /dev/shm/a4_test.log
echo ---FAILCOUNT---; grep -c FAIL /dev/shm/a4_test.log
