#include "chain.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>

#ifdef CERTA_WITH_ZLIB
#include <zlib.h>
#endif

namespace certa {

std::string ensembl_name(const std::string& n) {
  if (n == "chrM") return "MT";
  auto dotted = [](std::string s) {
    std::replace(s.begin(), s.end(), 'v', '.');
    return s;
  };
  if (n.compare(0, 6, "chrUn_") == 0) return dotted(n.substr(6));
  if (n.size() > 7 && n.compare(n.size() - 7, 7, "_random") == 0) {
    const size_t a = n.find('_'), b = n.find('_', a + 1);
    return dotted(n.substr(a + 1, b - a - 1));
  }
  return n.compare(0, 3, "chr") == 0 ? n.substr(3) : n;
}

namespace {

// Lines of a gzip or plain text file.
class Lines {
 public:
  explicit Lines(const std::string& path) {
#ifdef CERTA_WITH_ZLIB
    gz_ = gzopen(path.c_str(), "rb");
    if (!gz_) throw std::runtime_error("cannot open " + path);
#else
    if (path.size() > 3 && path.compare(path.size() - 3, 3, ".gz") == 0)
      throw std::runtime_error("built without zlib; cannot read " + path);
    f_ = std::fopen(path.c_str(), "r");
    if (!f_) throw std::runtime_error("cannot open " + path);
#endif
  }
  ~Lines() {
#ifdef CERTA_WITH_ZLIB
    if (gz_) gzclose(gz_);
#else
    if (f_) std::fclose(f_);
#endif
  }
  bool next(std::string& line) {
    line.clear();
    char buf[4096];
    for (;;) {
#ifdef CERTA_WITH_ZLIB
      if (!gzgets(gz_, buf, sizeof buf)) return !line.empty();
#else
      if (!std::fgets(buf, sizeof buf, f_)) return !line.empty();
#endif
      line += buf;
      if (!line.empty() && line.back() == '\n') {
        line.pop_back();
        return true;
      }
    }
  }

 private:
#ifdef CERTA_WITH_ZLIB
  gzFile gz_ = nullptr;
#else
  std::FILE* f_ = nullptr;
#endif
};

std::vector<std::string> fields(const std::string& s) {
  std::vector<std::string> out;
  size_t i = 0;
  while (i < s.size()) {
    while (i < s.size() && (s[i] == ' ' || s[i] == '\t' || s[i] == '\r')) ++i;
    size_t j = i;
    while (j < s.size() && s[j] != ' ' && s[j] != '\t' && s[j] != '\r') ++j;
    if (j > i) out.push_back(s.substr(i, j - i));
    i = j;
  }
  return out;
}

}  // namespace

Chain Chain::load(const std::string& path, const Reference& ref) {
  std::unordered_map<std::string, int> idx;
  for (size_t i = 0; i < ref.names.size(); ++i) idx[ref.names[i]] = static_cast<int>(i);
  auto resolve = [&](const std::string& q) {
    auto it = idx.find(q);
    if (it == idx.end()) it = idx.find(ensembl_name(q));
    return it == idx.end() ? -1 : it->second;
  };
  Chain c;
  Lines in(path);
  std::string line;
  Contig* cur = nullptr;
  int64_t tpos = 0, qpos = 0, qsize = 0;
  int qc = -1;
  bool rev = false;
  while (in.next(line)) {
    const std::vector<std::string> f = fields(line);
    if (f.empty()) continue;
    if (f[0] == "chain") {
      // chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
      if (f.size() < 12) throw std::runtime_error(path + ": bad chain header: " + line);
      cur = &c.by_src_[f[2]];
      tpos = std::atoll(f[5].c_str());
      qc = resolve(f[7]);
      qsize = std::atoll(f[8].c_str());
      rev = f[9] == "-";
      qpos = std::atoll(f[10].c_str());
      continue;
    }
    if (!cur) throw std::runtime_error(path + ": block before any chain header");
    const int64_t n = std::atoll(f[0].c_str());
    cur->b.push_back(Block{tpos, tpos + n, qc, rev, qsize, qpos});
    ++c.n_blocks_;
    tpos += n;
    qpos += n;
    if (f.size() == 3) {
      tpos += std::atoll(f[1].c_str());
      qpos += std::atoll(f[2].c_str());
    }
  }
  for (auto& kv : c.by_src_) {
    Contig& ct = kv.second;
    std::stable_sort(ct.b.begin(), ct.b.end(), [](const Block& x, const Block& y) { return x.start < y.start; });
    ct.starts.resize(ct.b.size());
    for (size_t i = 0; i < ct.b.size(); ++i) ct.starts[i] = ct.b[i].start;
  }
  return c;
}

bool Chain::lift(const std::string& src, int64_t pos, int* contig, int64_t* out, bool* reversed) const {
  const auto it = by_src_.find(src);
  if (it == by_src_.end()) return false;
  const Contig& ct = it->second;
  const size_t i = static_cast<size_t>(std::upper_bound(ct.starts.begin(), ct.starts.end(), pos) - ct.starts.begin());
  if (i == 0) return false;
  const Block& b = ct.b[i - 1];
  if (pos >= b.end) return false;
  const int64_t p = b.qstart + (pos - b.start);
  *contig = b.contig;
  *reversed = b.reversed;
  *out = b.reversed ? b.qsize - 1 - p : p;
  return true;
}

T2tView t2t_view(const T2tLift& l, const Placement& p, const Reference& ref) {
  T2tView v;
  if (!p.mapped || l.state == 0) return v;
  v.tier = l.tier;
  v.mapq = l.mapq;
  v.xe = l.xe;
  v.xd = l.xd;
  v.xo = l.xo;
  if (l.state != 1) {
    v.cls = 3;
    return v;
  }
  const int lead = p.r.n_cigar && (p.r.cigar[0] & 0xF) == kOpS ? static_cast<int>(p.r.cigar[0] >> 4) : 0;
  const int64_t left = p.r.ref_pos - lead;
  const bool here = l.contig == ref.contig_of(p.r.ref_pos) && l.strand == p.r.strand && std::llabs(l.left - left) <= 150;
  v.cls = here ? 1 : 2;
  return v;
}

bool Chain::lift_record(const std::string& src, int64_t left, int64_t span, int* contig, int64_t* out,
                        bool* flip) const {
  int ca = -1, cb = -1;
  int64_t pa = 0, pb = 0;
  bool ra = false, rb = false;
  const bool a = lift(src, left, &ca, &pa, &ra), b = lift(src, left + span - 1, &cb, &pb, &rb);
  if (a && !ra) {
    *contig = ca, *out = pa, *flip = false;
  } else if (b && rb) {
    *contig = cb, *out = pb, *flip = true;
  } else if (b && !rb) {
    *contig = cb, *out = pb - (span - 1), *flip = false;
  } else if (a && ra) {
    *contig = ca, *out = pa - (span - 1), *flip = true;
  } else {
    return false;
  }
  return *contig >= 0;
}

}  // namespace certa
