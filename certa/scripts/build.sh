#!/usr/bin/env bash
# Build CERTA on this host (x86-64 or aarch64; CUDA auto-detected) and run the
# brute-force certificate test. Uses CMake when available, otherwise calls
# the compilers directly (force that with CERTA_NO_CMAKE=1).
#   scripts/build.sh            # auto: GPU back-end if nvcc is available
#   CERTA_CUDA=OFF scripts/build.sh
#   CUDA_ARCHS=80 scripts/build.sh   # A100 only (default "80;86")
set -euo pipefail
here="$(cd "$(dirname "$0")/.." && pwd)"
build="${BUILD_DIR:-$here/build-$(hostname -s)-$(uname -m)}"
cuda="${CERTA_CUDA:-AUTO}"
archs="${CUDA_ARCHS:-80;86}"

# Pick up a CUDA toolkit that is installed but not on PATH.
if [[ "$cuda" != OFF ]] && ! command -v nvcc >/dev/null; then
  for d in /usr/local/cuda/bin /usr/local/cuda-*/bin; do
    [[ -x "$d/nvcc" ]] && export PATH="$d:$PATH" && break
  done
fi

echo "== host $(hostname -s), $(uname -m), $(nproc) cores"
command -v nvcc >/dev/null && nvcc --version | tail -1 || echo "== nvcc not found (CPU-only build)"
command -v nvidia-smi >/dev/null && nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

if [[ -z "${CERTA_NO_CMAKE:-}" ]] && command -v cmake >/dev/null; then
  cmake -S "$here" -B "$build" -DCERTA_CUDA="$cuda" -DCERTA_CUDA_ARCHS="$archs"
  cmake --build "$build" -j"$(nproc)"
else
  echo "== cmake not found: compiling directly"
  mkdir -p "$build"
  cxx="${CXX:-g++}"
  flags=(-O3 -std=c++17 -pthread -I"$here/include" -I"$here/src")
  libs=()
  if "$cxx" -march=native -x c++ -E /dev/null >/dev/null 2>&1; then flags+=(-march=native)
  elif "$cxx" -mcpu=native -x c++ -E /dev/null >/dev/null 2>&1; then flags+=(-mcpu=native); fi
  if echo '#include <zlib.h>' | "$cxx" -x c++ -E - >/dev/null 2>&1; then
    flags+=(-DCERTA_WITH_ZLIB); libs+=(-lz)
  else
    echo "== zlib headers not found: gzip input disabled (plain FASTA/FASTQ only)"
  fi
  objs=()
  for f in index seqio mapper_cpu main; do
    "$cxx" "${flags[@]}" -c "$here/src/$f.cpp" -o "$build/$f.o" &
    objs+=("$build/$f.o")
  done
  "$cxx" "${flags[@]}" -c "$here/tests/test_certa.cpp" -o "$build/test_certa.o" &
  "$cxx" "${flags[@]}" -c "$here/src/gpu_stub.cpp" -o "$build/gpu_stub.o" &
  wait
  core=("$build/index.o" "$build/seqio.o" "$build/mapper_cpu.o")
  if [[ "$cuda" != OFF ]] && command -v nvcc >/dev/null; then
    gencode=()
    for a in ${archs//;/ }; do gencode+=(-gencode "arch=compute_$a,code=sm_$a"); done
    nvcc -O3 -std=c++17 --expt-relaxed-constexpr "${gencode[@]}" -I"$here/include" -I"$here/src" \
      -c "$here/src/mapper_gpu.cu" -o "$build/mapper_gpu.o"
    nvcc -o "$build/certa" "$build/main.o" "$build/mapper_gpu.o" "${core[@]}" -Xcompiler -pthread "${libs[@]}"
  elif [[ "$cuda" == ON ]]; then
    echo "CERTA_CUDA=ON but nvcc was not found" >&2; exit 1
  else
    "$cxx" -pthread -o "$build/certa" "$build/main.o" "$build/gpu_stub.o" "${core[@]}" "${libs[@]}"
  fi
  "$cxx" -pthread -o "$build/test_certa" "$build/test_certa.o" "$build/gpu_stub.o" "${core[@]}" "${libs[@]}"
fi
"$build/test_certa"
echo "== built: $build/certa"
