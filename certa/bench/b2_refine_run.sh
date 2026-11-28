#!/bin/bash
# B2 refinement test: build the augmented index with the Cen/Sat-refined sinks
# (174 contigs, 1.5 Mb) in RAM, map the 4M-pair subsample, score vs the GRCh38
# baseline (base.npz from step 2), and analyse sink stealing. Disk-safe (RAM).
set -uo pipefail
export TMPDIR=/dev/shm/a2/tmp
cd /dev/shm/a2/src
git fetch --depth 1 origin main > /dev/shm/a2/b2r_pull.log 2>&1 \
  && git reset --hard FETCH_HEAD >> /dev/shm/a2/b2r_pull.log 2>&1 \
  || { echo PULL_FAIL; cat /dev/shm/a2/b2r_pull.log; exit 1; }
cd certa
cmake --build build -j 32 --target certa > /dev/shm/a2/b2r_build.log 2>&1 \
  || { echo BUILD_FAIL; tail -30 /dev/shm/a2/b2r_build.log; exit 1; }
echo BUILD_OK
C=/dev/shm/a2/src/certa/build/certa
B=~/certa-bench; W=~/certa-truth/B; T=/dev/shm/a2/src/certa/bench
FAI=$B/idx/ref.fa.fai
cd /dev/shm/a2
echo "== $(date +%T) refined augmented fasta + index"
cat ~/certa-bench/idx/ref.fa /dev/shm/a2/sinks_refined.fa > /dev/shm/a2/b2/aug_ref.fa
taskset -c 0-63 "$C" index /dev/shm/a2/b2/aug_ref.fa -o /dev/shm/a2/b2/aug_ref.cidx -q 22 -s 1 -t 64 2>/dev/shm/a2/b2/aug_ref.index.log
tail -2 /dev/shm/a2/b2/aug_ref.index.log
echo "== $(date +%T) map refined (4M pairs)"
taskset -c 0-63 "$C" map /dev/shm/a2/b2/aug_ref.cidx sub_1.fq sub_2.fq --intervals --iv-pair all --iv-rescue local --mapq pair \
  -t 64 -o /dev/shm/a2/b2/aug_ref.sam 2>/dev/shm/a2/b2/aug_ref.log
python3 "$T/truth_eval.py" score "$W/truth.npz" "$FAI" /dev/shm/a2/b2/aug_ref.npz < /dev/shm/a2/b2/aug_ref.sam 2>/dev/shm/a2/b2/aug_ref.eval.log
echo "== $(date +%T) report + sink analysis"
python3 "$T/truth_eval.py" report "$W/truth.npz" grch38=/dev/shm/a2/b2/base.npz refined=/dev/shm/a2/b2/aug_ref.npz > /dev/shm/a2/b2/b2_refine_report.md
python3 "$T/b2_sink_analysis.py" /dev/shm/a2/b2/aug_ref.sam "$W/truth.npz" | tee /dev/shm/a2/b2/b2_refine_sinks.md
echo DONE
sed -n '1,16p' /dev/shm/a2/b2/b2_refine_report.md
