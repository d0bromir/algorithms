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

  if (kp.size() > UINT32_MAX) throw std::runtime_error("more than 2^32 - 1 index entries");
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
  std::vector<uint32_t>& dir = ix.dir.owned;
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
constexpr uint64_t kChunk = uint64_t(1) << 21;  // directory entries per step (16 MB at 64 bits)

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
    pad();
    write(a.data(), a.size() * sizeof(T));
  }
  // A 32-bit array stored as 64-bit values (the directory in v1 to v3 files).
  void put_array64(const Array<uint32_t>& a) {
    put<uint64_t>(a.size());
    pad();
    std::vector<uint64_t> buf;
    for (uint64_t i = 0; i < a.size(); i += kChunk) {
      buf.assign(a.data() + i, a.data() + std::min(a.size(), i + kChunk));
      write(buf.data(), buf.size() * sizeof(uint64_t));
    }
  }
  void pad() {
    static const char zeros[kAlign] = {};
    long at = std::ftell(f);
    write(zeros, (kAlign - static_cast<uint64_t>(at) % kAlign) % kAlign);
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

// The file's 64-bit directory as 32-bit offsets (each is at most the number
// of entries, < 2^32), converted on several threads. With `drop`, src lies in
// a file mapping whose pages are released as they are read, so the 64-bit
// copy is never resident as a whole.
void narrow_dir(const uint64_t* src, uint64_t len, std::vector<uint32_t>& dst, bool drop,
                const std::string& path) {
  static const uint64_t page = static_cast<uint64_t>(sysconf(_SC_PAGESIZE));
  auto release = [&](uint64_t a, uint64_t b) {  // whole pages of src[a, b)
    const uint64_t lo = (reinterpret_cast<uint64_t>(src + a) + page - 1) / page * page;
    const uint64_t hi = reinterpret_cast<uint64_t>(src + b) / page * page;
    if (drop && hi > lo) madvise(reinterpret_cast<void*>(lo), hi - lo, MADV_DONTNEED);
  };
  dst.resize(len);
  const uint64_t steps = (len + kChunk - 1) / kChunk;
  const int threads = static_cast<int>(std::max<uint64_t>(1, std::min<uint64_t>({16, steps, std::thread::hardware_concurrency()})));
  std::vector<char> wide(threads, 0);
  parallel_for(threads, [&](int t) {
    uint64_t high = 0;  // the high halves, or-ed
    for (uint64_t k = steps * t / threads; k < steps * (t + 1) / threads; ++k) {
      const uint64_t a = k * kChunk, b = std::min(len, a + kChunk);
      for (uint64_t i = a; i < b; ++i) {
        high |= src[i] >> 32;
        dst[i] = static_cast<uint32_t>(src[i]);
      }
      release(a, b);
    }
    wide[t] = high != 0;
  });
  release(0, len);  // the pages split between steps
  for (char w : wide)
    if (w) throw std::runtime_error(path + ": directory offset beyond 2^32 - 1");
}

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
  f.put_array64(ix.dir);
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
      std::vector<uint64_t> dir64;
      f.get_vec(dir64);
      ix.dir = Array<uint32_t>();
      narrow_dir(dir64.data(), dir64.size(), ix.dir.owned, false, path);
      return;
    }
    if (std::memcmp(magic, kMagicV2, sizeof magic) != 0 &&
        std::memcmp(magic, kMagicV3, sizeof magic) != 0)
      throw std::runtime_error(path + " is not a CERTA index");
  }
  struct Fd {
    int fd;
    ~Fd() { ::close(fd); }
  } fd{::open(path.c_str(), O_RDONLY)};
  if (fd.fd < 0) throw std::runtime_error("cannot open " + path);
  struct stat st;
  if (fstat(fd.fd, &st) != 0) throw std::runtime_error("cannot stat " + path);
  auto m = std::make_shared<Mapping>();
  m->size = static_cast<size_t>(st.st_size);
  m->base = mmap(nullptr, m->size, PROT_READ, MAP_SHARED, fd.fd, 0);
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
  Array<uint64_t> dir64;
  c.array(dir64);
  // The arrays before the directory are used in place: one thread faults them
  // in (MAP_POPULATE on the same pages again) while the others read the
  // directory, once, into 32 bits, dropping its pages as they go.
  const size_t keep = static_cast<size_t>(reinterpret_cast<const char*>(dir64.data()) - c.base);
  bool populated = true;
  std::thread populate([&] {
    populated = keep == 0 ||
                mmap(m->base, keep, PROT_READ, MAP_SHARED | MAP_FIXED | MAP_POPULATE, fd.fd, 0) != MAP_FAILED;
  });
  ix.dir = Array<uint32_t>();
  try {
    narrow_dir(dir64.data(), dir64.size(), ix.dir.owned, true, path);
  } catch (...) {
    populate.join();
    throw;
  }
  populate.join();
  if (!populated) throw std::runtime_error("cannot mmap " + path);
  ix.mapping = m;
}

}  // namespace certa
