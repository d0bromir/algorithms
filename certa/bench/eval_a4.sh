#!/bin/bash
# A4 accuracy: single-end --bio with the locus-specific phantom prior
# (--bio-phantom-bed segdups, pi inside 0.01, pi_unique 0) on the same 4M-read
# subsample as A3; compared to the A3 baseline (b0) and global phantom (b1).
# Disk-safe: BED decompressed to RAM, SAM streamed, npz/reports in RAM.
set -uo pipefail
mkdir -p /dev/shm/a2/tmp
export TMPDIR=/dev/shm/a2/tmp
cd /dev/shm/a2/src/certa
cmake --build build -j 32 --target certa > /dev/shm/a4_cli.log 2>&1 \
  || { echo CLI_BUILD_FAIL; tail -30 /dev/shm/a4_cli.log; exit 1; }
echo CLI_OK
C=/dev/shm/a2/src/certa/build/certa
B=~/certa-bench; W=~/certa-truth/B; T=/dev/shm/a2/src/certa/bench
IDX=$B/idx/grch38.s1v3.cidx; FAI=$B/idx/ref.fa.fai
cd /dev/shm/a2
zcat ~/certa-data/truth/strat/SegmentalDuplications/GRCh38_segdups.bed.gz > /dev/shm/a2/segdups.bed
echo "== $(date +%T) bio locus prior (single-end, 4M reads)"
taskset -c 0-63 "$C" map "$IDX" sub_1.fq --bio --bio-phantom 0.01 --bio-phantom-bed segdups.bed \
  --bio-phantom-pi 0.01 --bio-phantom-pi-unique 0 -t 64 -o /dev/stdout -u /dev/null 2>b2.log \
  | python3 "$T/truth_eval.py" score "$W/truth.npz" "$FAI" b2.npz 2>b2.eval.log
echo "== $(date +%T) report"
python3 "$T/truth_eval.py" report "$W/truth.npz" bio=b0.npz global=b1.npz locus=b2.npz > a4_report.md
echo DONE
sed -n '1,34p' a4_report.md
