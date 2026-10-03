#include "index.h"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <thread>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

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
  std::vector<uint8_t>& seq = ref.seq.owned;
  seq.assign(kContigPad, 4);
  bool open = false;
  while (in.getline(line)) {
    if (line.empty()) continue;
    if (line[0] == '>') {
      if (open) {
        ref.lengths.back() = seq.size() - ref.offsets.back();
        seq.insert(seq.end(), kContigPad, 4);
      }
      size_t sp = line.find_first_of(" \t");
      ref.names.push_back(line.substr(1, sp == std::string::npos ? std::string::npos : sp - 1));
      ref.offsets.push_back(seq.size());
      ref.lengths.push_back(0);
      open = true;
      continue;
    }
    if (!open) throw std::runtime_error("FASTA does not start with '>': " + path);
    for (char c : line) seq.push_back(encode_base(c));
  }
  if (!open) throw std::runtime_error("no sequences in " + path);
  ref.lengths.back() = seq.size() - ref.offsets.back();
  seq.insert(seq.end(), kContigPad, 4);
  if (seq.size() >= (uint64_t(1) << 32))
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

  std::vector<uint64_t>& keys = ix.keys.owned;
  std::vector<uint32_t>& pos = ix.pos.owned;
  keys.resize(kp.size());
  pos.resize(kp.size());
  for (size_t i = 0; i < kp.size(); ++i) {
    keys[i] = kp[i].key;
    pos[i] = kp[i].pos;
  }
  std::vector<KeyPos>().swap(kp);

  int bits = 8;
  while (bits < 28 && (uint64_t(1) << (bits + 1)) <= keys.size()) ++bits;
  ix.dir_bits = std::min(bits, 2 * q);
  const int shift = 2 * q - ix.dir_bits;
  const uint64_t nb = uint64_t(1) << ix.dir_bits;
  std::vector<uint64_t>& dir = ix.dir.owned;
  dir.assign(nb + 1, 0);
  for (uint64_t k : keys) ++dir[(k >> shift) + 1];
  for (uint64_t b = 0; b < nb; ++b) dir[b + 1] += dir[b];
  return ix;
}

IndexView Index::view(const Reference& ref) const {
  IndexView v;
  v.ref = ref.seq.data();
  v.ref_len = ref.seq.size();
  v.keys = keys16.size() ? nullptr : keys.data();
  v.keys16 = keys16.size() ? keys16.data() : nullptr;
  v.pos = pos.data();
  v.dir = dir.data();
  v.n = pos.size();
  v.q = q;
  v.s = s;
  v.dir_bits = dir_bits;
  return v;
}

namespace {

// v1: arrays packed after their sizes (read with fread).
// v2: every array starts on a 64-byte boundary, so it can be used in place
// from a memory mapping.
const char kMagicV1[8] = {'C', 'E', 'R', 'T', 'A', 'I', 'X', '1'};
const char kMagicV2[8] = {'C', 'E', 'R', 'T', 'A', 'I', 'X', '2'};
// v3: as v2, plus a key-width field; keys are stored as their low 16 bits
// whenever 2q - dir_bits <= 16 (the bucket directory holds the rest).
const char kMagicV3[8] = {'C', 'E', 'R', 'T', 'A', 'I', 'X', '3'};
constexpr uint64_t kAlign = 64;

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
  // v2 layout: size, zero padding to a 64-byte boundary, then the data.
  template <class T> void put_array(const Array<T>& a) {
    put<uint64_t>(a.size());
    static const char zeros[kAlign] = {};
    long at = std::ftell(f);
    write(zeros, (kAlign - static_cast<uint64_t>(at) % kAlign) % kAlign);
    write(a.data(), a.size() * sizeof(T));
  }
  template <class T> void get_vec(std::vector<T>& v) {  // v1 layout
    v.resize(get<uint64_t>());
    read(v.data(), v.size() * sizeof(T));
  }
};

// A read-only file mapping, released when the last Index referring to it goes.
struct Mapping {
  void* base = MAP_FAILED;
  size_t size = 0;
  ~Mapping() {
    if (base != MAP_FAILED) munmap(base, size);
  }
};

// Sequential reader over the mapped v2 file.
struct Cursor {
  const char* base;
  uint64_t at, size;
  void need(uint64_t n) const {
    if (at + n > size) throw std::runtime_error("index file truncated");
  }
  template <class T> T get() {
    need(sizeof(T));
    T v;
    std::memcpy(&v, base + at, sizeof v);
    at += sizeof v;
    return v;
  }
  template <class T> void array(Array<T>& a) {
    a.len = get<uint64_t>();
    at += (kAlign - at % kAlign) % kAlign;
    need(a.len * sizeof(T));
    a.ptr = reinterpret_cast<const T*>(base + at);
    a.owned.clear();
    at += a.len * sizeof(T);
  }
};

}  // namespace

void save_index(const std::string& path, const Reference& ref, const Index& ix) {
  File f(path, "wb");
  f.write(kMagicV3, sizeof kMagicV3);
  f.put<int32_t>(ix.q);
  f.put<int32_t>(ix.s);
  f.put<int32_t>(ix.dir_bits);
  const int shift = 2 * ix.q - ix.dir_bits;
  const int32_t key_bits = shift <= 16 ? 16 : 64;
  f.put<int32_t>(key_bits);
  f.put<uint64_t>(ref.names.size());
  for (size_t i = 0; i < ref.names.size(); ++i) {
    f.put<uint64_t>(ref.names[i].size());
    f.write(ref.names[i].data(), ref.names[i].size());
    f.put<uint64_t>(ref.offsets[i]);
    f.put<uint64_t>(ref.lengths[i]);
  }
  f.put_array(ref.seq);
  if (key_bits == 16) {
    Array<uint16_t> k16;
    if (ix.keys16.size()) {
      k16.ptr = ix.keys16.data();
      k16.len = ix.keys16.size();
    } else {
      const uint64_t mask = (uint64_t(1) << shift) - 1;
      k16.owned.resize(ix.keys.size());
      for (uint64_t i = 0; i < ix.keys.size(); ++i) k16.owned[i] = static_cast<uint16_t>(ix.keys[i] & mask);
    }
    f.put_array(k16);
  } else {
    if (!ix.keys.size() && ix.keys16.size())
      throw std::runtime_error("cannot widen 16-bit keys");
    f.put_array(ix.keys);
  }
  f.put_array(ix.pos);
  f.put_array(ix.dir);
}

void load_index(const std::string& path, Reference& ref, Index& ix) {
  char magic[8];
  {
    File f(path, "rb");
    f.read(magic, sizeof magic);
    if (std::memcmp(magic, kMagicV1, sizeof magic) == 0) {  // legacy: fread
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
      f.get_vec(ref.seq.owned);
      f.get_vec(ix.keys.owned);
      f.get_vec(ix.pos.owned);
      f.get_vec(ix.dir.owned);
      return;
    }
    if (std::memcmp(magic, kMagicV2, sizeof magic) != 0 &&
        std::memcmp(magic, kMagicV3, sizeof magic) != 0)
      throw std::runtime_error(path + " is not a CERTA index");
  }
  int fd = ::open(path.c_str(), O_RDONLY);
  if (fd < 0) throw std::runtime_error("cannot open " + path);
  struct stat st;
  if (fstat(fd, &st) != 0) {
    ::close(fd);
    throw std::runtime_error("cannot stat " + path);
  }
  auto m = std::make_shared<Mapping>();
  m->size = static_cast<size_t>(st.st_size);
  m->base = mmap(nullptr, m->size, PROT_READ, MAP_SHARED | MAP_POPULATE, fd, 0);
  ::close(fd);
  if (m->base == MAP_FAILED) throw std::runtime_error("cannot mmap " + path);
  Cursor c{static_cast<const char*>(m->base), sizeof magic, m->size};
  ix.q = c.get<int32_t>();
  ix.s = c.get<int32_t>();
  ix.dir_bits = c.get<int32_t>();
  const bool v3 = std::memcmp(magic, kMagicV3, sizeof magic) == 0;
  const int key_bits = v3 ? c.get<int32_t>() : 64;
  if (key_bits != 16 && key_bits != 64) throw std::runtime_error(path + ": bad key width");
  uint64_t n = c.get<uint64_t>();
  ref.names.resize(n);
  ref.offsets.resize(n);
  ref.lengths.resize(n);
  for (uint64_t i = 0; i < n; ++i) {
    uint64_t len = c.get<uint64_t>();
    c.need(len);
    ref.names[i].assign(c.base + c.at, len);
    c.at += len;
    ref.offsets[i] = c.get<uint64_t>();
    ref.lengths[i] = c.get<uint64_t>();
  }
  c.array(ref.seq);
  if (key_bits == 16) c.array(ix.keys16);
  else c.array(ix.keys);
  c.array(ix.pos);
  c.array(ix.dir);
  ix.mapping = m;
}

}  // namespace certa
