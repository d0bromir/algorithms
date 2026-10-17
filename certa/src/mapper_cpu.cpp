#include <algorithm>
#include <atomic>
#include <thread>

#include "mapper.h"

namespace certa {

void encode_batch(const std::vector<FastqRecord>& recs, ReadBatch& out, int threads) {
  const size_t n = recs.size();
  out.offs.resize(n);
  out.lens.resize(n);
  out.hashes.resize(n);
  uint64_t total = 0;
  for (size_t i = 0; i < n; ++i) {
    out.offs[i] = total;
    // Reads longer than 65535 are clamped; process_read rejects > LMAX anyway.
    out.lens[i] = static_cast<uint16_t>(std::min<size_t>(recs[i].seq.size(), 65535));
    total += out.lens[i];
  }
  out.codes.resize(total);
  threads = std::max(1, std::min<int>(threads, static_cast<int>(n / 4096) + 1));
  std::vector<std::thread> pool;
  for (int t = 0; t < threads; ++t) {
    pool.emplace_back([&, t] {
      for (size_t i = n * t / threads; i < n * (t + 1) / threads; ++i) {
        const FastqRecord& r = recs[i];
        uint8_t* dst = out.codes.data() + out.offs[i];
        for (size_t j = 0; j < out.lens[i]; ++j) dst[j] = encode_base(r.seq[j]);
        uint64_t h = 1469598103934665603ULL;
        for (char c : r.name) { h ^= static_cast<uint8_t>(c); h *= 1099511628211ULL; }
        out.hashes[i] = h;
      }
    });
  }
  for (auto& th : pool) th.join();
}

namespace {

// The q-mer keys count_parts looks up for a read: every part, both strands,
// every sampling offset (parts with an N give fewer keys; that only means
// fewer prefetches).
constexpr int kMaxKeys = 2 * PMAX * SMAX;
int part_keys(const IndexView& ix, const uint8_t* codes, int L, uint64_t* keys) {
  const int P = L <= LMAX ? part_count(L, ix.q, ix.s) : 0;
  if (P < 1) return 0;
  uint8_t rc[LMAX];
  for (int i = 0; i < L; ++i) {
    const uint8_t c = codes[L - 1 - i];
    rc[i] = c < 4 ? (uint8_t)(3 - c) : (uint8_t)4;
  }
  int n = 0;
  for (int st = 0; st < 2; ++st) {
    const uint8_t* s = st ? rc : codes;
    for (int j = 0; j < P; ++j) {
      int off, len;
      part_geometry(L, P, j, &off, &len);
      for (int t = 0; t < ix.s; ++t) {
        uint64_t k;
        if (kmer_key(s + off + t, ix.q, &k)) keys[n++] = k;
      }
    }
  }
  return n;
}

}  // namespace

void map_cpu(const IndexView& ix, const Params& p, const ReadBatch& batch,
             std::vector<Result>& out, int threads) {
  out.resize(batch.size());
  std::atomic<size_t> next{0};
  const size_t chunk = 256;
  const int shift = 2 * ix.q - ix.dir_bits;
  // Every lookup costs two dependent random reads of the index (the bucket
  // directory, then the bucket), and enumeration a third (the hits'
  // positions). They are pipelined across reads: while read i is mapped, the
  // buckets and positions of read i + 1 and the directory entries of read
  // i + 2 are prefetched. Results do not change.
  auto worker = [&] {
    std::unique_ptr<Workspace> ws(new Workspace);
    uint64_t keys[3][kMaxKeys];
    int nk[3] = {0, 0, 0};
    auto read_keys = [&](size_t i, int slot) {
      nk[slot] = part_keys(ix, batch.codes.data() + batch.offs[i], batch.lens[i], keys[slot]);
      for (int k = 0; k < nk[slot]; ++k) __builtin_prefetch(ix.dir + (keys[slot][k] >> shift));
    };
    auto prefetch_buckets = [&](int slot) {
      for (int k = 0; k < nk[slot]; ++k) {
        const uint64_t l = ix.dir[keys[slot][k] >> shift];
        if (ix.keys16) __builtin_prefetch(ix.keys16 + l);
        else __builtin_prefetch(ix.keys + l);
        __builtin_prefetch(ix.pos + l);
      }
    };
    for (;;) {
      size_t a = next.fetch_add(chunk);
      if (a >= batch.size()) break;
      size_t b = std::min(batch.size(), a + chunk);
      if (a < b) read_keys(a, a % 3);
      if (a + 1 < b) read_keys(a + 1, (a + 1) % 3);
      if (a < b) prefetch_buckets(a % 3);
      for (size_t i = a; i < b; ++i) {
        if (i + 2 < b) read_keys(i + 2, (i + 2) % 3);
        if (i + 1 < b) prefetch_buckets((i + 1) % 3);
        process_read(ix, p, batch.codes.data() + batch.offs[i], batch.lens[i],
                     batch.hashes[i], *ws, out[i]);
      }
    }
  };
  std::vector<std::thread> pool;
  for (int t = 1; t < threads; ++t) pool.emplace_back(worker);
  worker();
  for (auto& th : pool) th.join();
}

}  // namespace certa
