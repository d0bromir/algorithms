#include "index.h"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <thread>

#include "seqio.h"

namespace certa {

uint8_t encode_base(char c) {
  switch (c) {
    case 'A': case 'a': return 0;
    case 'C': case 'c': return 1;
    case 'G': case 'g': return 2;
    case 'T': case 't': return 3;
    default: return 4;
  }
}

Reference Reference::load_fasta(const std::string& path) {
  Reference ref;
  LineReader in(path);
  std::string line;
  ref.seq.assign(kContigPad, 4);
  bool open = false;
  while (in.getline(line)) {
    if (line.empty()) continue;
    if (line[0] == '>') {
      if (open) {
        ref.lengths.back() = ref.seq.size() - ref.offsets.back();
        ref.seq.insert(ref.seq.end(), kContigPad, 4);
      }
      size_t sp = line.find_first_of(" \t");
      ref.names.push_back(line.substr(1, sp == std::string::npos ? std::string::npos : sp - 1));
      ref.offsets.push_back(ref.seq.size());
      ref.lengths.push_back(0);
      open = true;
      continue;
    }
    if (!open) throw std::runtime_error("FASTA does not start with '>': " + path);
    for (char c : line) ref.seq.push_back(encode_base(c));
  }
  if (!open) throw std::runtime_error("no sequences in " + path);
  ref.lengths.back() = ref.seq.size() - ref.offsets.back();
  ref.seq.insert(ref.seq.end(), kContigPad, 4);
  if (ref.seq.size() >= (uint64_t(1) << 32))
    throw std::runtime_error("reference longer than 2^32 bases is not supported");
  return ref;
}

int Reference::contig_of(int64_t pos) const {
  auto it = std::upper_bound(offsets.begin(), offsets.end(), static_cast<uint64_t>(pos < 0 ? 0 : pos));
  if (it == offsets.begin()) return -1;
  int c = static_cast<int>(it - offsets.begin()) - 1;
  if (pos < static_cast<int64_t>(offsets[c]) ||
      pos >= static_cast<int64_t>(offsets[c] + lengths[c]))
    return -1;
  return c;
}

namespace {

struct KeyPos {
  uint64_t key;
  uint32_t pos;
  bool operator<(const KeyPos& o) const {
    return key < o.key || (key == o.key && pos < o.pos);
  }
};

template <class F>
void parallel_for(int threads, F&& f) {
  std::vector<std::thread> pool;
  for (int t = 0; t < threads; ++t) pool.emplace_back(f, t);
  for (auto& th : pool) th.join();
}

}  // namespace

Index Index::build(const Reference& ref, int q, int s, int threads) {
  if (q < 8 || q > 32) throw std::runtime_error("q must be in [8, 32]");
  if (s < 1 || s > SMAX) throw std::runtime_error("s must be in [1, 16]");
  threads = std::max(1, threads);
  Index ix;
  ix.q = q;
  ix.s = s;
  const uint64_t len = ref.seq.size();
  const uint64_t nsamp = len >= static_cast<uint64_t>(q) ? (len - q) / s + 1 : 0;
  const uint8_t* seq = ref.seq.data();

  // Pass 1: count valid (N-free) sampled q-mers per thread range.
  std::vector<uint64_t> counts(threads + 1, 0);
  auto range = [&](int t, uint64_t* a, uint64_t* b) {
    *a = nsamp * t / threads;
    *b = nsamp * (t + 1) / threads;
  };
  parallel_for(threads, [&](int t) {
    uint64_t a, b, n = 0, key;
    range(t, &a, &b);
    for (uint64_t i = a; i < b; ++i) n += kmer_key(seq + i * s, q, &key);
    counts[t + 1] = n;
  });
  for (int t = 0; t < threads; ++t) counts[t + 1] += counts[t];
  std::vector<KeyPos> kp(counts[threads]);

  // Pass 2: fill, then sort each thread's chunk.
  parallel_for(threads, [&](int t) {
    uint64_t a, b, key, w = counts[t];
    range(t, &a, &b);
    for (uint64_t i = a; i < b; ++i)
      if (kmer_key(seq + i * s, q, &key))
        kp[w++] = {key, static_cast<uint32_t>(i * s)};
    std::sort(kp.begin() + counts[t], kp.begin() + counts[t + 1]);
  });
  // Merge sorted chunks pairwise in parallel rounds.
  std::vector<uint64_t> bounds(counts.begin(), counts.end());
  while (bounds.size() > 2) {
    std::vector<uint64_t> next{bounds[0]};
    std::vector<std::thread> pool;
    size_t i = 0;
    for (; i + 2 < bounds.size(); i += 2) {
      uint64_t a = bounds[i], m = bounds[i + 1], b = bounds[i + 2];
      pool.emplace_back([&kp, a, m, b] {
        std::inplace_merge(kp.begin() + a, kp.begin() + m, kp.begin() + b);
      });
      next.push_back(b);
    }
    if (i + 1 < bounds.size()) next.push_back(bounds[i + 1]);  // odd chunk
    for (auto& th : pool) th.join();
    bounds.swap(next);
  }

  ix.keys.resize(kp.size());
  ix.pos.resize(kp.size());
  for (size_t i = 0; i < kp.size(); ++i) {
    ix.keys[i] = kp[i].key;
    ix.pos[i] = kp[i].pos;
  }
  std::vector<KeyPos>().swap(kp);

  int bits = 8;
  while (bits < 26 && (uint64_t(1) << (bits + 1)) <= ix.keys.size()) ++bits;
  ix.dir_bits = std::min(bits, 2 * q);
  const int shift = 2 * q - ix.dir_bits;
  const uint64_t nb = uint64_t(1) << ix.dir_bits;
  ix.dir.assign(nb + 1, 0);
  for (uint64_t k : ix.keys) ++ix.dir[(k >> shift) + 1];
  for (uint64_t b = 0; b < nb; ++b) ix.dir[b + 1] += ix.dir[b];
  return ix;
}

IndexView Index::view(const Reference& ref) const {
  IndexView v;
  v.ref = ref.seq.data();
  v.ref_len = ref.seq.size();
  v.keys = keys.data();
  v.pos = pos.data();
  v.dir = dir.data();
  v.n = keys.size();
  v.q = q;
  v.s = s;
  v.dir_bits = dir_bits;
  return v;
}

namespace {

const char kMagic[8] = {'C', 'E', 'R', 'T', 'A', 'I', 'X', '1'};

struct File {
  std::FILE* f;
  File(const std::string& p, const char* mode) : f(std::fopen(p.c_str(), mode)) {
    if (!f) throw std::runtime_error("cannot open " + p);
  }
  ~File() { std::fclose(f); }
  void write(const void* p, size_t n) {
    if (n && std::fwrite(p, 1, n, f) != n) throw std::runtime_error("write failed");
  }
  void read(void* p, size_t n) {
    if (n && std::fread(p, 1, n, f) != n) throw std::runtime_error("index file truncated");
  }
  template <class T> void put(T v) { write(&v, sizeof v); }
  template <class T> T get() { T v; read(&v, sizeof v); return v; }
  template <class T> void put_vec(const std::vector<T>& v) {
    put<uint64_t>(v.size());
    write(v.data(), v.size() * sizeof(T));
  }
  template <class T> void get_vec(std::vector<T>& v) {
    v.resize(get<uint64_t>());
    read(v.data(), v.size() * sizeof(T));
  }
};

}  // namespace

void save_index(const std::string& path, const Reference& ref, const Index& ix) {
  File f(path, "wb");
  f.write(kMagic, sizeof kMagic);
  f.put<int32_t>(ix.q);
  f.put<int32_t>(ix.s);
  f.put<int32_t>(ix.dir_bits);
  f.put<uint64_t>(ref.names.size());
  for (size_t i = 0; i < ref.names.size(); ++i) {
    f.put<uint64_t>(ref.names[i].size());
    f.write(ref.names[i].data(), ref.names[i].size());
    f.put<uint64_t>(ref.offsets[i]);
    f.put<uint64_t>(ref.lengths[i]);
  }
  f.put_vec(ref.seq);
  f.put_vec(ix.keys);
  f.put_vec(ix.pos);
  f.put_vec(ix.dir);
}

void load_index(const std::string& path, Reference& ref, Index& ix) {
  File f(path, "rb");
  char magic[8];
  f.read(magic, sizeof magic);
  if (std::memcmp(magic, kMagic, sizeof magic) != 0)
    throw std::runtime_error(path + " is not a CERTA index (or wrong version)");
  ix.q = f.get<int32_t>();
  ix.s = f.get<int32_t>();
  ix.dir_bits = f.get<int32_t>();
  uint64_t n = f.get<uint64_t>();
  ref.names.resize(n);
  ref.offsets.resize(n);
  ref.lengths.resize(n);
  for (uint64_t i = 0; i < n; ++i) {
    ref.names[i].resize(f.get<uint64_t>());
    f.read(&ref.names[i][0], ref.names[i].size());
    ref.offsets[i] = f.get<uint64_t>();
    ref.lengths[i] = f.get<uint64_t>();
  }
  f.get_vec(ref.seq);
  f.get_vec(ix.keys);
  f.get_vec(ix.pos);
  f.get_vec(ix.dir);
}

}  // namespace certa
