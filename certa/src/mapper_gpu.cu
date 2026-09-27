// CUDA back-end: one thread per read, running the same certa::process_read
// as the CPU back-end. Reference and index stay resident on the device;
// reads are streamed in batches.
#include <cuda_runtime.h>

#include <stdexcept>
#include <string>

#include "mapper.h"

namespace certa {

namespace {

void check(cudaError_t e, const char* what) {
  if (e != cudaSuccess)
    throw std::runtime_error(std::string("CUDA error in ") + what + ": " +
                             cudaGetErrorString(e));
}

template <class T>
T* upload(const T* host, size_t n, const char* what) {
  T* dev = nullptr;
  check(cudaMalloc(&dev, (n ? n : 1) * sizeof(T)), what);
  if (n) check(cudaMemcpy(dev, host, n * sizeof(T), cudaMemcpyHostToDevice), what);
  return dev;
}

__global__ void map_kernel(IndexView ix, Params p, const uint8_t* codes,
                           const uint64_t* offs, const uint16_t* lens,
                           const uint64_t* hashes, size_t n, Result* out) {
  size_t i = blockIdx.x * static_cast<size_t>(blockDim.x) + threadIdx.x;
  if (i >= n) return;
  Workspace ws;  // per-thread local memory
  Result r;
  process_read(ix, p, codes + offs[i], lens[i], hashes[i], ws, r);
  out[i] = r;
}

}  // namespace

struct GpuMapper::Impl {
  int device = 0;
  IndexView view{};
  uint8_t* d_ref = nullptr;
  uint64_t* d_keys = nullptr;
  uint32_t* d_pos = nullptr;
  uint64_t* d_dir = nullptr;
  // Batch buffers, grown on demand.
  uint8_t* d_codes = nullptr;
  uint64_t* d_offs = nullptr;
  uint16_t* d_lens = nullptr;
  uint64_t* d_hashes = nullptr;
  Result* d_out = nullptr;
  size_t cap_codes = 0, cap_reads = 0;
  cudaDeviceProp prop{};

  void reserve(size_t codes, size_t reads) {
    if (codes > cap_codes) {
      cudaFree(d_codes);
      cap_codes = codes + codes / 4;
      check(cudaMalloc(&d_codes, cap_codes), "alloc codes");
    }
    if (reads > cap_reads) {
      cudaFree(d_offs); cudaFree(d_lens); cudaFree(d_hashes); cudaFree(d_out);
      cap_reads = reads + reads / 4;
      check(cudaMalloc(&d_offs, cap_reads * sizeof(uint64_t)), "alloc offs");
      check(cudaMalloc(&d_lens, cap_reads * sizeof(uint16_t)), "alloc lens");
      check(cudaMalloc(&d_hashes, cap_reads * sizeof(uint64_t)), "alloc hashes");
      check(cudaMalloc(&d_out, cap_reads * sizeof(Result)), "alloc results");
    }
  }
};

GpuMapper::GpuMapper(const Reference& ref, const Index& ix, int device)
    : impl_(new Impl) {
  Impl& m = *impl_;
  m.device = device;
  check(cudaSetDevice(device), "cudaSetDevice");
  check(cudaGetDeviceProperties(&m.prop, device), "cudaGetDeviceProperties");
  // Workspace lives in per-thread local memory (~26 KB); raise the limit.
  check(cudaDeviceSetLimit(cudaLimitStackSize, 48 * 1024), "stack limit");
  m.d_ref = upload(ref.seq.data(), ref.seq.size(), "upload reference");
  m.d_keys = upload(ix.keys.data(), ix.keys.size(), "upload keys");
  m.d_pos = upload(ix.pos.data(), ix.pos.size(), "upload positions");
  m.d_dir = upload(ix.dir.data(), ix.dir.size(), "upload directory");
  m.view = ix.view(ref);
  m.view.ref = m.d_ref;
  m.view.keys = m.d_keys;
  m.view.pos = m.d_pos;
  m.view.dir = m.d_dir;
}

GpuMapper::~GpuMapper() {
  Impl& m = *impl_;
  cudaFree(m.d_ref); cudaFree(m.d_keys); cudaFree(m.d_pos); cudaFree(m.d_dir);
  cudaFree(m.d_codes); cudaFree(m.d_offs); cudaFree(m.d_lens);
  cudaFree(m.d_hashes); cudaFree(m.d_out);
}

void GpuMapper::map(const Params& p, const ReadBatch& b, std::vector<Result>& out) {
  Impl& m = *impl_;
  const size_t n = b.size();
  out.resize(n);
  if (n == 0) return;
  m.reserve(b.codes.size(), n);
  check(cudaMemcpy(m.d_codes, b.codes.data(), b.codes.size(), cudaMemcpyHostToDevice), "copy codes");
  check(cudaMemcpy(m.d_offs, b.offs.data(), n * sizeof(uint64_t), cudaMemcpyHostToDevice), "copy offs");
  check(cudaMemcpy(m.d_lens, b.lens.data(), n * sizeof(uint16_t), cudaMemcpyHostToDevice), "copy lens");
  check(cudaMemcpy(m.d_hashes, b.hashes.data(), n * sizeof(uint64_t), cudaMemcpyHostToDevice), "copy hashes");
  const int threads = 128;
  const unsigned blocks = static_cast<unsigned>((n + threads - 1) / threads);
  map_kernel<<<blocks, threads>>>(m.view, p, m.d_codes, m.d_offs, m.d_lens,
                                  m.d_hashes, n, m.d_out);
  check(cudaGetLastError(), "kernel launch");
  check(cudaDeviceSynchronize(), "kernel");
  check(cudaMemcpy(out.data(), m.d_out, n * sizeof(Result), cudaMemcpyDeviceToHost), "copy results");
}

std::string GpuMapper::device_name() const {
  const cudaDeviceProp& pr = impl_->prop;
  return std::string(pr.name) + " (sm_" + std::to_string(pr.major) +
         std::to_string(pr.minor) + ", " +
         std::to_string(pr.totalGlobalMem >> 20) + " MiB)";
}

bool GpuMapper::compiled_in() { return true; }

}  // namespace certa
