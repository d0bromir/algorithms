#include <algorithm>
#include <atomic>
#include <thread>

#include "mapper.h"

namespace certa {

void encode_batch(const std::vector<FastqRecord>& recs, ReadBatch& out) {
  out.codes.clear();
  out.offs.resize(recs.size());
  out.lens.resize(recs.size());
  out.hashes.resize(recs.size());
  for (size_t i = 0; i < recs.size(); ++i) {
    const FastqRecord& r = recs[i];
    out.offs[i] = out.codes.size();
    // Reads longer than 65535 are clamped; process_read rejects > LMAX anyway.
    out.lens[i] = static_cast<uint16_t>(std::min<size_t>(r.seq.size(), 65535));
    for (size_t j = 0; j < out.lens[i]; ++j) out.codes.push_back(encode_base(r.seq[j]));
    uint64_t h = 1469598103934665603ULL;
    for (char c : r.name) { h ^= static_cast<uint8_t>(c); h *= 1099511628211ULL; }
    out.hashes[i] = h;
  }
}

void map_cpu(const IndexView& ix, const Params& p, const ReadBatch& batch,
             std::vector<Result>& out, int threads) {
  out.resize(batch.size());
  std::atomic<size_t> next{0};
  const size_t chunk = 256;
  auto worker = [&] {
    std::unique_ptr<Workspace> ws(new Workspace);
    for (;;) {
      size_t a = next.fetch_add(chunk);
      if (a >= batch.size()) break;
      size_t b = std::min(batch.size(), a + chunk);
      for (size_t i = a; i < b; ++i)
        process_read(ix, p, batch.codes.data() + batch.offs[i], batch.lens[i],
                     batch.hashes[i], *ws, out[i]);
    }
  };
  std::vector<std::thread> pool;
  for (int t = 1; t < threads; ++t) pool.emplace_back(worker);
  worker();
  for (auto& th : pool) th.join();
}

}  // namespace certa
