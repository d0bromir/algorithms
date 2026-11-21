#!/bin/bash
# Re-run A1 eval only (joined.npz + models already in /dev/shm/ceiling).
set -euo pipefail
cd /dev/shm/ceiling
python3 truth_ceiling.py eval --half model_cv.npz joined.npz > ceiling_holdout.md
python3 truth_ceiling.py eval model.npz joined.npz > ceiling_full.md
echo DONE
cat ceiling_holdout.md
