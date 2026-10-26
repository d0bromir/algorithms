#!/usr/bin/env bash
# End-to-end run on synthetic data: simulate -> index -> map (CPU, and GPU if
# available) -> check CPU/GPU outputs are identical -> placement accuracy.
#   scripts/demo.sh [genome_bp] [n_reads]          (defaults 50 Mbp, 2 M reads)
# Environment: BUILD_DIR (default build-<host>-<arch>), WORK (default ./demo-<host>),
#              K (edit budget, default 2), THREADS (default all cores),
#              DEVICE (CUDA device index, default 0).
set -euo pipefail
here="$(cd "$(dirname "$0")/.." && pwd)"
bin="${BUILD_DIR:-$here/build-$(hostname -s)-$(uname -m)}/certa"
work="${WORK:-$PWD/demo-$(hostname -s)}"
gbp="${1:-50000000}"
nreads="${2:-2000000}"
k="${K:-2}"
threads="${THREADS:-$(nproc)}"
mkdir -p "$work"
[[ -x "$bin" ]] || { echo "build first: scripts/build.sh"; exit 1; }

[[ -s "$work/genome.fa" ]] || python3 "$here/tools/simulate.py" genome "$work/genome.fa" --length "$gbp"
[[ -s "$work/reads.fq" ]] || python3 "$here/tools/simulate.py" reads "$work/genome.fa" "$work/reads.fq" --n "$nreads"
[[ -s "$work/genome.cidx" ]] || "$bin" index "$work/genome.fa" -o "$work/genome.cidx" -t "$threads"

echo "== CPU ($threads threads)"
"$bin" map "$work/genome.cidx" "$work/reads.fq" -k "$k" -t "$threads" \
  -o "$work/cpu.sam" -u "$work/cpu.uncertified.fq" --stats "$work/cpu.stats.json"
python3 "$here/tools/simulate.py" eval "$work/cpu.sam"

usage="$("$bin" 2>&1 || true)"  # the usage screen reports GPU support
if [[ "$usage" == *"GPU support compiled in: yes"* ]] && command -v nvidia-smi >/dev/null; then
  echo "== GPU"
  "$bin" map "$work/genome.cidx" "$work/reads.fq" -k "$k" -t "$threads" --gpu --device "${DEVICE:-0}" \
    -o "$work/gpu.sam" -u "$work/gpu.uncertified.fq" --stats "$work/gpu.stats.json"
  if cmp -s <(grep -v '^@PG' "$work/cpu.sam") <(grep -v '^@PG' "$work/gpu.sam") &&
     cmp -s "$work/cpu.uncertified.fq" "$work/gpu.uncertified.fq"; then
    echo "== CPU and GPU outputs are identical"
  else
    echo "== ERROR: CPU and GPU outputs differ"; exit 1
  fi
fi
echo "== outputs in $work"
