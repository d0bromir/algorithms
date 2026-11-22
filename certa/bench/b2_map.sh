#!/bin/bash
# B2-build step 2b: map a 4M-pair HG002 subsample with the GRCh38 index
# (baseline) and the augmented GRCh38+sinks index, score against truth, and
# analyse sink placements (stolen real reads vs correct missing-copy abstentions).
# SAMs + npz in /dev/shm (RAM). Sink contigs (sink_*) are absent from the GRCh38
# truth reference, so truth_eval counts them as not-placed (abstain) -- the B2
# semantics. Disk-safe.
set -uo pipefail
export TMPDIR=/dev/shm/a2/tmp
C=/dev/shm/a2/src/certa/build/certa
B=~/certa-bench; W=~/certa-truth/B; T=/dev/shm/a2/src/certa/bench
GR=$B/idx/grch38.s1v3.cidx; AUG=/dev/shm/a2/b2/augmented.cidx; FAI=$B/idx/ref.fa.fai
cd /dev/shm/a2
[ -f sub_2.fq ] || head -n 16000000 "$W/B_2.fq" > sub_2.fq
echo "== $(date +%T) baseline GRCh38 (4M pairs)"
taskset -c 0-63 "$C" map "$GR" sub_1.fq sub_2.fq --intervals --iv-pair all --iv-rescue local --mapq pair \
  -t 64 -o /dev/shm/a2/b2/base.sam 2>/dev/shm/a2/b2/base.log
python3 "$T/truth_eval.py" score "$W/truth.npz" "$FAI" /dev/shm/a2/b2/base.npz < /dev/shm/a2/b2/base.sam 2>/dev/shm/a2/b2/base.eval.log
echo "== $(date +%T) augmented GRCh38+sinks (4M pairs)"
taskset -c 0-63 "$C" map "$AUG" sub_1.fq sub_2.fq --intervals --iv-pair all --iv-rescue local --mapq pair \
  -t 64 -o /dev/shm/a2/b2/aug.sam 2>/dev/shm/a2/b2/aug.log
python3 "$T/truth_eval.py" score "$W/truth.npz" "$FAI" /dev/shm/a2/b2/aug.npz < /dev/shm/a2/b2/aug.sam 2>/dev/shm/a2/b2/aug.eval.log
echo "== $(date +%T) report"
python3 "$T/truth_eval.py" report "$W/truth.npz" grch38=/dev/shm/a2/b2/base.npz augmented=/dev/shm/a2/b2/aug.npz > /dev/shm/a2/b2/b2_report.md
echo "== $(date +%T) sink analysis"
python3 "$T/b2_sink_analysis.py" /dev/shm/a2/b2/aug.sam "$W/truth.npz" | tee /dev/shm/a2/b2/b2_sinks.md
echo DONE
sed -n '1,18p' /dev/shm/a2/b2/b2_report.md
