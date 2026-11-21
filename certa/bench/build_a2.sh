#!/bin/bash
# Build A2 (--mapq-local) from the pushed branch in RAM and run the certificate
# tests. Disk-safe: clone + build live in /dev/shm; the production checkout
# ~/certa-bio/certa is untouched.
set -uo pipefail
mkdir -p /dev/shm/a2
cd /dev/shm/a2
if [ ! -d src ]; then
  git clone --depth 1 --branch claude/a2-mapq-local https://github.com/d0bromir/algorithms.git src \
    > /dev/shm/a2/clone.log 2>&1 || { echo CLONE_FAIL; tail -20 /dev/shm/a2/clone.log; exit 1; }
fi
cd src/certa
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release > /dev/shm/a2/cmake.log 2>&1 \
  || { echo CMAKE_FAIL; tail -25 /dev/shm/a2/cmake.log; exit 1; }
cmake --build build -j 32 --target test_certa verify_lemmas > /dev/shm/a2/build.log 2>&1 \
  || { echo BUILD_FAIL; tail -40 /dev/shm/a2/build.log; exit 1; }
echo BUILD_OK
./build/verify_lemmas > /dev/shm/a2/lemmas.log 2>&1; echo LEMMAS_EXIT=$?
./build/test_certa > /dev/shm/a2/test_certa.log 2>&1; echo CERTA_EXIT=$?
echo ---TEST-TAIL---
tail -8 /dev/shm/a2/test_certa.log
echo ---MAPQ-LOCAL-LINE---
grep mapq-local /dev/shm/a2/test_certa.log
