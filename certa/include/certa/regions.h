// A set of reference regions (segmental duplications / copy-deficit loci) as
// merged intervals of the concatenated (global) position, for a per-locus
// prior. Built from a BED (chrom start end); a leading "chr" on chrom is
// ignored so UCSC and Ensembl names both map to the reference's contigs.
#pragma once
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace certa {

class RegionTrack {
 public:
  void add(int64_t lo, int64_t hi) { iv_.push_back({lo, hi}); }

  // Load a plain BED (chrom start end ...), mapping chrom to a contig via
  // `names` and `offsets` to global positions. Contigs absent from the
  // reference are skipped. Call finalize() after add() if not using load_bed.
  void load_bed(const std::string& path, const std::vector<std::string>& names,
                const std::vector<uint64_t>& offsets) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open region BED " + path);
    std::string line;
    while (std::getline(f, line)) {
      if (line.empty() || line[0] == '#' || line.compare(0, 5, "track") == 0 ||
          line.compare(0, 6, "browser") == 0)
        continue;
      std::istringstream ss(line);
      std::string chrom;
      int64_t lo, hi;
      if (!(ss >> chrom >> lo >> hi)) continue;
      if (chrom.compare(0, 3, "chr") == 0) chrom = chrom.substr(3);
      int c = -1;
      for (size_t i = 0; i < names.size(); ++i)
        if (names[i] == chrom) { c = static_cast<int>(i); break; }
      if (c < 0) continue;
      add(static_cast<int64_t>(offsets[c]) + lo, static_cast<int64_t>(offsets[c]) + hi);
    }
    finalize();
  }

  // Sort and merge the intervals, and index their starts for contains().
  void finalize() {
    std::sort(iv_.begin(), iv_.end());
    std::vector<std::pair<int64_t, int64_t>> m;
    for (const auto& iv : iv_) {
      if (!m.empty() && iv.first <= m.back().second)
        m.back().second = std::max(m.back().second, iv.second);
      else
        m.push_back(iv);
    }
    iv_.swap(m);
    starts_.resize(iv_.size());
    for (size_t i = 0; i < iv_.size(); ++i) starts_[i] = iv_[i].first;
  }

  // Is global position pos in any region ([start, end))?
  bool contains(int64_t pos) const {
    auto it = std::upper_bound(starts_.begin(), starts_.end(), pos);
    if (it == starts_.begin()) return false;
    const size_t i = static_cast<size_t>(it - starts_.begin()) - 1;
    return pos < iv_[i].second;
  }

  bool empty() const { return iv_.empty(); }
  size_t size() const { return iv_.size(); }

 private:
  std::vector<std::pair<int64_t, int64_t>> iv_;
  std::vector<int64_t> starts_;
};

}  // namespace certa
