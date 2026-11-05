// Reference loading and the q-mer index (built in memory, or memory-mapped).
#pragma once
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "certa/core.h"

namespace certa {

// Ns inserted before, between and after contigs so no alignment within the
// certified radius can span two contigs.
constexpr uint64_t kContigPad = 512;

// Read-only array that either owns its storage or points into a mapping (or
// into `held`, storage of its own that is not zero-filled first).
template <class T>
struct Array {
  std::vector<T> owned;
  const T* ptr = nullptr;
  uint64_t len = 0;
  std::shared_ptr<T> held;
  const T* data() const { return ptr ? ptr : owned.data(); }
  uint64_t size() const { return ptr ? len : owned.size(); }
  const T& operator[](uint64_t i) const { return data()[i]; }
};

struct Reference {
  std::vector<std::string> names;
  std::vector<uint64_t> offsets;  // start of each contig in `seq`
  std::vector<uint64_t> lengths;
  Array<uint8_t> seq;             // codes 0..3, 4 = N

  static Reference load_fasta(const std::string& path);
  // Index of the contig containing concatenated position `pos`, or -1.
  int contig_of(int64_t pos) const;
};

struct Index {
  int q = 22, s = 8, dir_bits = 0;
  Array<uint64_t> keys;    // full keys (built in memory, or v2 files)
  Array<uint16_t> keys16;  // low 16 bits (v3 files, when 2q - dir_bits <= 16)
  // v4 files: the high 8 of those bits (keys8), or no key bits at all; the
  // other ref_key_bits are read from the reference (IndexView). 0: not v4.
  Array<uint8_t> keys8;
  int ref_key_bits = 0;
  Array<uint32_t> pos;
  Array<uint32_t> dir;     // in memory as 32 bits (files store 64)
  std::shared_ptr<void> mapping;  // keeps a memory-mapped index file alive

  static Index build(const Reference& ref, int q, int s, int threads);
  IndexView view(const Reference& ref) const;
  // Key bits stored per entry: 64, 16, 8 or 0.
  int key_bits() const;
};

// key_bits: what the file stores per entry. -1: as before (16 when the
// bucket leaves at most 16 key bits, else 64); 8 or 0: a v4 file (16 or fewer
// bits left; 8 needs more than 8), whose other key bits are read from the
// reference. Any index (built, or loaded in any version) can be saved so.
void save_index(const std::string& path, const Reference& ref, const Index& ix, int key_bits = -1);
// Memory-maps the index file; `ref` and `ix` point into the mapping, which
// stays alive while `ix` does. Repeated runs reuse the OS page cache.
void load_index(const std::string& path, Reference& ref, Index& ix);

// A base's code: A, C, G, T (either case) are 0 to 3, anything else 4. A
// table, so encoding a read is one load per base.
struct BaseCodes {
  uint8_t code[256];
  constexpr BaseCodes() : code() {
    for (int i = 0; i < 256; ++i) code[i] = 4;
    code['A'] = code['a'] = 0;
    code['C'] = code['c'] = 1;
    code['G'] = code['g'] = 2;
    code['T'] = code['t'] = 3;
  }
};
inline constexpr BaseCodes kBaseCodes{};
inline uint8_t encode_base(char c) { return kBaseCodes.code[static_cast<unsigned char>(c)]; }

}  // namespace certa
