// Reference loading and the sampled q-mer index.
#pragma once
#include <cstdint>
#include <string>
#include <vector>

#include "certa/core.h"

namespace certa {

// Ns inserted before, between and after contigs so no alignment within the
// certified radius can span two contigs.
constexpr uint64_t kContigPad = 512;

struct Reference {
  std::vector<std::string> names;
  std::vector<uint64_t> offsets;  // start of each contig in `seq`
  std::vector<uint64_t> lengths;
  std::vector<uint8_t> seq;       // codes 0..3, 4 = N

  static Reference load_fasta(const std::string& path);
  // Index of the contig containing concatenated position `pos`, or -1.
  int contig_of(int64_t pos) const;
};

struct Index {
  int q = 22, s = 8, dir_bits = 0;
  std::vector<uint64_t> keys;
  std::vector<uint32_t> pos;
  std::vector<uint64_t> dir;

  static Index build(const Reference& ref, int q, int s, int threads);
  IndexView view(const Reference& ref) const;
};

void save_index(const std::string& path, const Reference& ref, const Index& ix);
void load_index(const std::string& path, Reference& ref, Index& ix);

uint8_t encode_base(char c);

}  // namespace certa
