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

// The q-mer keys count_parts looks up for a read: every part without an N,
// both strands, every sampling offset (rolled as count_parts rolls them).
constexpr int kMaxKeys = 2 * PMAX * SMAX;
int part_keys(const IndexView& ix, const uint8_t* codes, int L, uint64_t* keys) {
  const int P = L <= LMAX ? part_count(L, ix.q, ix.s) : 0;
  if (P < 1) return 0;
  uint8_t rc[LMAX];
  for (int i = 0; i < L; ++i) {
    const uint8_t c = codes[L - 1 - i];
    rc[i] = c < 4 ? (uint8_t)(3 - c) : (uint8_t)4;
  }
  const uint64_t mask = ix.q >= 32 ? ~0ull : (1ull << (2 * ix.q)) - 1;
  int n = 0;
  for (int st = 0; st < 2; ++st) {
    const uint8_t* s = st ? rc : codes;
    for (int j = 0; j < P; ++j) {
      int off, len;
      part_geometry(L, P, j, &off, &len);
      // One pass over the part: the rolled key ends a q-mer at i = q - 1 + t
      // for shift t < s; the keys are kept only if the part has no N.
      uint64_t k = 0;
      bool has_n = false;
      for (int i = 0; i < len; ++i) {
        const uint8_t c = s[off + i];
        has_n |= c > 3;
        k = ((k << 2) | (c & 3)) & mask;
        const int t = i - (ix.q - 1);
        if (t >= 0 && t < ix.s) keys[n + t] = k;
      }
      if (!has_n) n += ix.s;
    }
  }
  return n;
}

}  // namespace

void map_reads(const IndexView& ix, const Params& p, const ReadBatch& batch, const size_t* idx, size_t n,
               Result* out, Workspace& ws) {
  // Every lookup costs two dependent random reads of the index (the bucket
  // directory, then the bucket), and enumeration a third (the hits'
  // positions). They are pipelined across reads: while read k is mapped, the
  // buckets and positions of read k + 1 and the directory entries of read
  // k + 2 are prefetched. Results do not change.
  const int shift = 2 * ix.q - ix.dir_bits;
  uint64_t keys[3][kMaxKeys];
  int nk[3] = {0, 0, 0};
  auto read_keys = [&](size_t k) {
    const size_t i = idx[k];
    const int slot = static_cast<int>(k % 3);
    nk[slot] = part_keys(ix, batch.codes.data() + batch.offs[i], batch.lens[i], keys[slot]);
    for (int j = 0; j < nk[slot]; ++j) __builtin_prefetch(ix.dir + (keys[slot][j] >> shift));
  };
  auto prefetch_buckets = [&](size_t k) {
    const int slot = static_cast<int>(k % 3);
    for (int j = 0; j < nk[slot]; ++j) {
      const uint64_t l = ix.dir[keys[slot][j] >> shift];
      if (ix.keys16) __builtin_prefetch(ix.keys16 + l);
      else __builtin_prefetch(ix.keys + l);
      __builtin_prefetch(ix.pos + l);
    }
  };
  if (n > 0) read_keys(0);
  if (n > 1) read_keys(1);
  if (n > 0) prefetch_buckets(0);
  for (size_t k = 0; k < n; ++k) {
    if (k + 2 < n) read_keys(k + 2);
    if (k + 1 < n) prefetch_buckets(k + 1);
    const size_t i = idx[k];
    process_read(ix, p, batch.codes.data() + batch.offs[i], batch.lens[i], batch.hashes[i], ws, out[k]);
  }
}

void map_cpu(const IndexView& ix, const Params& p, const ReadBatch& batch,
             std::vector<Result>& out, int threads) {
  out.resize(batch.size());
  std::atomic<size_t> next{0};
  constexpr size_t chunk = 256;
  auto worker = [&] {
    std::unique_ptr<Workspace> ws(new Workspace);
    size_t idx[chunk];
    for (;;) {
      const size_t a = next.fetch_add(chunk);
      if (a >= batch.size()) break;
      const size_t b = std::min(batch.size(), a + chunk);
      for (size_t i = a; i < b; ++i) idx[i - a] = i;
      map_reads(ix, p, batch, idx, b - a, out.data() + a, *ws);
    }
  };
  std::vector<std::thread> pool;
  for (int t = 1; t < threads; ++t) pool.emplace_back(worker);
  worker();
  for (auto& th : pool) th.join();
}

}  // namespace certa
