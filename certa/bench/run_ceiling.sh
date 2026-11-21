#!/bin/bash
# A1 (research guideline Pillar A): the offline gap-graded MAPQ ceiling on
# HG002 test B. Disk-safe: BWA-MEM2's SAM is streamed straight into the
# extractor (no BAM written) and every intermediate lives in /dev/shm (RAM);
# only small .npz/.md are produced. Fit on even pairs, report odd pairs held
# out (calibrate-then-freeze), plus the full-data ceiling.
set -euo pipefail
cd /dev/shm/ceiling
BWA=~/micromamba/envs/certa-bench/bin/bwa-mem2
IDX=~/certa-bench/idx/bwa-mem2
W=~/certa-truth/B
echo "== $(date +%T) feat: stream bwa-mem2 -> gap, join CERTA default + truth"
taskset -c 0-63 "$BWA" mem -t 64 "$IDX" "$W/B_1.fq" "$W/B_2.fq" 2>/dev/shm/ceiling/bwa.log \
  | python3 truth_ceiling.py feat "$W/truth.npz" "$W/certa.npz" /dev/shm/ceiling/joined.npz
echo "== $(date +%T) fit (even pairs) + fit (all)"
python3 truth_ceiling.py fit --cv /dev/shm/ceiling/joined.npz /dev/shm/ceiling/model_cv.npz
python3 truth_ceiling.py fit /dev/shm/ceiling/joined.npz /dev/shm/ceiling/model.npz
echo "== $(date +%T) eval: odd pairs held out, then full"
python3 truth_ceiling.py eval --half /dev/shm/ceiling/model_cv.npz /dev/shm/ceiling/joined.npz > /dev/shm/ceiling/ceiling_holdout.md
python3 truth_ceiling.py eval /dev/shm/ceiling/model.npz /dev/shm/ceiling/joined.npz > /dev/shm/ceiling/ceiling_full.md
echo "== $(date +%T) DONE"
cat /dev/shm/ceiling/ceiling_holdout.md
