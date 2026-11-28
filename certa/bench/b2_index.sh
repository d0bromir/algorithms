#!/bin/bash
# B2-build step 2a: build the augmented (GRCh38 + CHM13 sinks) CERTA index in
# /dev/shm (RAM) so the ~20.5 Gb index never touches the nearly-full root disk.
# q=22 s=1 matches grch38.s1v3. Disk/RAM guarded.
set -uo pipefail
export TMPDIR=/dev/shm/a2/tmp
cd /dev/shm/a2/src
git fetch --depth 1 origin main > /dev/shm/a2/b2_pull.log 2>&1 \
  && git reset --hard FETCH_HEAD >> /dev/shm/a2/b2_pull.log 2>&1 \
  || { echo PULL_FAIL; cat /dev/shm/a2/b2_pull.log; exit 1; }
cd certa
cmake --build build -j 32 --target certa > /dev/shm/a2/b2_build.log 2>&1 \
  || { echo BUILD_FAIL; tail -30 /dev/shm/a2/b2_build.log; exit 1; }
echo BUILD_OK
C=/dev/shm/a2/src/certa/build/certa
mkdir -p /dev/shm/a2/b2
avail=$(df --output=avail -BG /dev/shm | tail -1 | tr -dc '0-9')
echo "/dev/shm available: ${avail}G"
[ "$avail" -ge 30 ] || { echo LOW_RAM_ABORT; exit 1; }
echo "== $(date +%T) augmented fasta (GRCh38 + sinks)"
cat ~/certa-bench/idx/ref.fa /dev/shm/a2/sinks.fa > /dev/shm/a2/b2/augmented.fa
ls -la /dev/shm/a2/b2/augmented.fa
echo "== $(date +%T) index build (q=22 s=1)"
taskset -c 0-63 "$C" index /dev/shm/a2/b2/augmented.fa -o /dev/shm/a2/b2/augmented.cidx -q 22 -s 1 -t 64 2>/dev/shm/a2/b2/index.log
echo "== $(date +%T) index done"
tail -4 /dev/shm/a2/b2/index.log
ls -la /dev/shm/a2/b2/augmented.cidx
df -h /dev/shm /
echo DONE_INDEX
