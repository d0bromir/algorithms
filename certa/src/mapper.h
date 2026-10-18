// CPU and GPU drivers around certa::process_read.
#pragma once
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "certa/core.h"
#include "index.h"
#include "seqio.h"

namespace certa {

// Reads packed for either back-end: base codes back to back.
struct ReadBatch {
  std::vector<uint8_t> codes;
  std::vector<uint64_t> offs;
  std::vector<uint16_t> lens;
  std::vector<uint64_t> hashes;  // FNV-1a of the read name (tie-breaking)
  size_t size() const { return lens.size(); }
};

void encode_batch(const std::vector<FastqRecord>& recs, ReadBatch& out, int threads = 1);

void map_cpu(const IndexView& ix, const Params& p, const ReadBatch& batch,
             std::vector<Result>& out, int threads);

// Maps reads idx[0..n) of the batch on the calling thread into out[0..n),
// with map_cpu's prefetching; each result is the one map_cpu gives.
void map_reads(const IndexView& ix, const Params& p, const ReadBatch& batch, const size_t* idx, size_t n,
               Result* out, Workspace& ws);

// Keeps the reference and index resident on one CUDA device.
class GpuMapper {
 public:
  GpuMapper(const Reference& ref, const Index& ix, int device);
  ~GpuMapper();
  void map(const Params& p, const ReadBatch& batch, std::vector<Result>& out);
  std::string device_name() const;
  static bool compiled_in();

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace certa
