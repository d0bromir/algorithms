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

// Read-only array that either owns its storage or points into a mapping.
template <class T>
struct Array {
  std::vector<T> owned;
  const T* ptr = nullptr;
  uint64_t len = 0;
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
  Array<uint32_t> pos;
  Array<uint64_t> dir;
  std::shared_ptr<void> mapping;  // keeps a memory-mapped index file alive

  static Index build(const Reference& ref, int q, int s, int threads);
  IndexView view(const Reference& ref) const;
};

void save_index(const std::string& path, const Reference& ref, const Index& ix);
// Memory-maps the index file; `ref` and `ix` point into the mapping, which
// stays alive while `ix` does. Repeated runs reuse the OS page cache.
void load_index(const std::string& path, Reference& ref, Index& ix);

uint8_t encode_base(char c);

}  // namespace certa
