#!/usr/bin/env bash
# Build CERTA on this host (x86-64 or aarch64; CUDA auto-detected) and run the
# brute-force certificate test.
#   scripts/build.sh            # auto: GPU back-end if nvcc is available
#   CERTA_CUDA=OFF scripts/build.sh
#   CUDA_ARCHS=80 scripts/build.sh   # A100 only (default "80;86")
set -euo pipefail
here="$(cd "$(dirname "$0")/.." && pwd)"
build="${BUILD_DIR:-$here/build-$(hostname -s)-$(uname -m)}"
cuda="${CERTA_CUDA:-AUTO}"

# Pick up a CUDA toolkit that is installed but not on PATH.
if [[ "$cuda" != OFF ]] && ! command -v nvcc >/dev/null; then
  for d in /usr/local/cuda/bin /usr/local/cuda-*/bin; do
    [[ -x "$d/nvcc" ]] && export PATH="$d:$PATH" && break
  done
fi

echo "== host $(hostname -s), $(uname -m), $(nproc) cores"
command -v nvcc >/dev/null && nvcc --version | tail -1 || echo "== nvcc not found (CPU-only build)"
command -v nvidia-smi >/dev/null && nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

cmake -S "$here" -B "$build" -DCERTA_CUDA="$cuda" -DCERTA_CUDA_ARCHS="${CUDA_ARCHS:-80;86}"
cmake --build "$build" -j"$(nproc)"
"$build/test_certa"
echo "== built: $build/certa"
