// certa check: see check.h. Written apart from the mapper (core.h, pair.h,
// repeat.h, interval.h): its own lookup, Lemma L, DP and SAM parsing.
#include "check.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <fstream>
#include <functional>
#include <mutex>
#include <thread>

namespace certa {
namespace check {

const char* const kVerdictNames[kNumVerdicts] = {"accepted",    "malformed",   "wrong_score", "wrong_xl",
                                                 "xu_unproven", "xh_unproven", "wrong_xq"};

namespace {

// BWA-MEM's default scores: the objective every claim is about.
constexpr int kM = 1, kX = 4, kO = 6, kE = 1, kC = 5;
constexpr int kNeg = -(1 << 28);
constexpr uint64_t kMaxOccurrences = 1u << 20;  // per record; beyond it the claim counts as unproven

uint8_t code_of(char c) {
  switch (c) {
    case 'A': case 'a': return 0;
    case 'C': case 'c': return 1;
    case 'G': case 'g': return 2;
    case 'T': case 't': return 3;
    default: return 4;
  }
}

std::vector<uint8_t> revcomp(const std::vector<uint8_t>& v) {
  std::vector<uint8_t> o(v.size());
  for (size_t i = 0; i < v.size(); ++i) {
    const uint8_t c = v[v.size() - 1 - i];
    o[i] = c < 4 ? static_cast<uint8_t>(3 - c) : 4;
  }
  return o;
}

// Lemma L: the most an alignment of a read of length L scores when, of the
// probes, only those in `free` and at most `extra` others occur exactly in
// it. Every other probe lying inside its read span [x, y) costs >= 5, and
// the costs are disjoint; x need only be 0 or just inside a probe, y L or
// just inside one.
int lemma_l(int L, const std::vector<Probe>& p, uint32_t free, int extra) {
  const int n = static_cast<int>(p.size());
  int best = kNeg;
  for (int a = -1; a < n; ++a) {
    const int x = a < 0 ? 0 : p[a].off + 1;
    for (int b = -1; b < n; ++b) {
      const int y = b < 0 ? L : p[b].end - 1;
      if (y <= x) continue;
      int k = 0;
      for (int i = 0; i < n; ++i) k += !((free >> i) & 1) && p[i].off >= x && p[i].end <= y;
      k = std::max(0, k - extra);
      best = std::max(best, (y - x) * kM - (kM + kX) * k - (x > 0 ? kC : 0) - (y < L ? kC : 0));
    }
  }
  return best;
}

// BWA-MEM's single-end MAPQ formula (mem_approx_mapq_se, mapQ_coef_len 50).
int mapq_formula(int score, int sub, int L) {
  if (sub < 19) sub = 19;
  if (score <= 0 || sub >= score) return 0;
  const double l = L, identity = 1.0 - (l * kM - score) / (kM + kX) / l;
  double t = l < 50 ? 1.0 : std::log(50.0) / std::log(l);
  t *= identity * identity;
  const int q = static_cast<int>(6.02 * (score - sub) / kM * t * t + 0.499);
  return q < 0 ? 0 : (q > 60 ? 60 : q);
}

bool parse_cigar(const std::string& s, std::vector<uint32_t>& out) {
  out.clear();
  uint32_t len = 0;
  bool digits = false;
  for (char ch : s) {
    if (ch >= '0' && ch <= '9') {
      len = len * 10 + static_cast<uint32_t>(ch - '0');
      digits = true;
      continue;
    }
    if (!digits) return false;
    uint32_t op;
    switch (ch) {
      case 'M': op = 0; break;
      case 'I': op = 1; break;
      case 'D': op = 2; break;
      case 'S': op = 4; break;
      default: return false;
    }
    out.push_back((len << 4) | op);
    len = 0;
    digits = false;
  }
  return !digits && !out.empty();
}

int64_t span_of(const std::vector<uint32_t>& cigar) {
  int64_t s = 0;
  for (uint32_t c : cigar)
    if ((c & 15) == 0 || (c & 15) == 2) s += c >> 4;
  return s;
}

void split(const std::string& s, char sep, std::vector<std::string>& out) {
  out.clear();
  size_t a = 0;
  for (;;) {
    const size_t b = s.find(sep, a);
    out.push_back(s.substr(a, b == std::string::npos ? std::string::npos : b - a));
    if (b == std::string::npos) return;
    a = b + 1;
  }
}

void parallel(int threads, const std::function<void(int)>& f) {
  std::vector<std::thread> pool;
  for (int t = 1; t < threads; ++t) pool.emplace_back(f, t);
  f(0);
  for (auto& th : pool) th.join();
}

}  // namespace

Checker::Checker(const Reference& ref, const Index& ix)
    : ref_(ref.seq.data()),
      ref_len_(static_cast<int64_t>(ref.seq.size())),
      names_(ref),
      q_(ix.q),
      s_(ix.s),
      dir_bits_(ix.dir_bits),
      keys_(ix.keys16.size() ? nullptr : ix.keys.data()),
      keys16_(ix.keys16.size() ? ix.keys16.data() : nullptr),
      pos_(ix.pos.data()),
      dir_(ix.dir.data()) {
  for (size_t c = 0; c < ref.names.size(); ++c) contig_[ref.names[c]] = static_cast<int>(c);
}

int Checker::contig_of(int64_t pos) const {
  const auto& off = names_.offsets;
  auto it = std::upper_bound(off.begin(), off.end(), static_cast<uint64_t>(pos < 0 ? 0 : pos));
  if (pos < 0 || it == off.begin()) return -1;
  const size_t c = static_cast<size_t>(it - off.begin()) - 1;
  return pos < static_cast<int64_t>(off[c] + names_.lengths[c]) ? static_cast<int>(c) : -1;
}

void Checker::range(uint64_t key, uint64_t* lo, uint64_t* hi) const {
  const int shift = 2 * q_ - dir_bits_;
  const uint64_t b = key >> shift, a = dir_[b], e = dir_[b + 1];
  if (keys16_) {
    const uint16_t k = static_cast<uint16_t>(key & ((uint64_t(1) << shift) - 1));
    *lo = static_cast<uint64_t>(std::lower_bound(keys16_ + a, keys16_ + e, k) - keys16_);
    *hi = static_cast<uint64_t>(std::upper_bound(keys16_ + a, keys16_ + e, k) - keys16_);
  } else {
    *lo = static_cast<uint64_t>(std::lower_bound(keys_ + a, keys_ + e, key) - keys_);
    *hi = static_cast<uint64_t>(std::upper_bound(keys_ + a, keys_ + e, key) - keys_);
  }
}

// Every exact occurrence of probe p of read rd, as a diagonal (reference
// start of the read). An occurrence at p0 contains exactly one sampled
// position p0 + i, i < s, so each is found once, through the q-mer at i.
void Checker::occurrences(const uint8_t* rd, const Probe& p, std::vector<int64_t>& out) const {
  const int len = p.end - p.off;
  for (int k = 0; k < len; ++k)
    if (rd[p.off + k] > 3) return;  // N never matches: no occurrence
  for (int i = 0; i < s_ && i + q_ <= len; ++i) {
    uint64_t key = 0;
    for (int k = 0; k < q_; ++k) key = (key << 2) | rd[p.off + i + k];
    uint64_t lo, hi;
    range(key, &lo, &hi);
    for (uint64_t e = lo; e < hi; ++e) {
      const int64_t start = static_cast<int64_t>(pos_[e]) - i;
      if (start < 0 || start + len > ref_len_) continue;
      bool ok = true;
      for (int k = 0; ok && k < len; ++k) ok = ref_[start + k] == rd[p.off + k];
      if (ok) out.push_back(start - p.off);
    }
  }
}

// The best local score of the read (one strand) over alignments whose cells
// all lie on diagonals [dlo, dhi] (reference position minus read position)
// and in reference [j0, j1). Gotoh, alignments start and end with an
// aligned base, -5 per clipped read end. Cell (i, b) aligns read base i - 1
// with reference i - 1 + dlo + b; an insertion keeps the reference position
// (from b + 1 in the row above), a deletion the read position (from b - 1
// in the same row). An indel's cells lie between the diagonals of the
// aligned bases around it, so restricting every cell to the band loses no
// alignment whose aligned bases are in it.
int Checker::best_band(const uint8_t* rd, int L, int64_t dlo, int64_t dhi, int64_t j0, int64_t j1,
                       Work* w) const {
  if (dhi < dlo) return kNeg;
  const int B = static_cast<int>(dhi - dlo + 1);
  std::vector<int> H(B + 1, kNeg), I(B + 1, kNeg), Hn(B + 1, kNeg), In(B + 1, kNeg);
  int best = kNeg;
  for (int i = 1; i <= L; ++i) {
    const uint8_t c = rd[i - 1];
    const int start = i == 1 ? 0 : -kC;
    int d = kNeg, left = kNeg;  // deletion state and H at (i, b - 1)
    for (int b = 0; b < B; ++b) {
      const int64_t j = i - 1 + dlo + b;
      if (j < j0 || j >= j1) {
        Hn[b] = In[b] = left = d = kNeg;
        continue;
      }
      const uint8_t r = (j >= 0 && j < ref_len_) ? ref_[j] : 4;
      const int m = std::max(H[b], start) + ((c == r && r < 4) ? kM : -kX);
      const int ins = i > 1 ? std::max(H[b + 1] - kO - kE, I[b + 1] - kE) : kNeg;
      d = std::max(left - kO - kE, d - kE);
      const int h = std::max(m, std::max(ins, d));
      Hn[b] = h;
      In[b] = ins;
      left = h;
      best = std::max(best, m - (i < L ? kC : 0));
    }
    Hn[B] = In[B] = kNeg;
    std::swap(H, Hn);
    std::swap(I, In);
    w->cells += static_cast<uint64_t>(B);
  }
  return best;
}

// Whether the witness proves that every alignment (sub: every alignment not
// overlapping the record's) scores <= bound. See check.h, (i)-(iii).
bool Checker::proven(const Claim& c, const std::vector<uint8_t>* rd, int bound, bool sub, Work* w) const {
  const int L = static_cast<int>(c.read.size());
  if (bound >= L * kM) return true;  // no alignment scores more than L
  if (c.t <= 0) return false;
  for (int st = 0; st < 2; ++st) {
    if (c.probes[st].size() > 32) return false;
    int prev = 0;
    for (const Probe& p : c.probes[st]) {  // disjoint, increasing, long enough to look up
      if (p.off < prev || p.end > L || p.end - p.off < q_ + s_ - 1) return false;
      prev = p.end;
    }
  }
  int ut = kNeg;
  for (int st = 0; st < 2; ++st) ut = std::max(ut, lemma_l(L, c.probes[st], 0, c.t - 1));
  if (ut > bound) return false;  // (i)
  const int pad = std::max(0, L * kM - bound - kO - kE);
  int64_t p0 = 0, p1 = 0;
  if (sub) {
    p0 = c.aln.pos;
    p1 = p0 + span_of(c.aln.cigar);
  }
  std::vector<std::pair<int64_t, int>> occ;
  std::vector<int64_t> tmp;
  for (int st = 0; st < 2; ++st) {
    occ.clear();
    const std::vector<Probe>& P = c.probes[st];
    for (size_t k = 0; k < P.size(); ++k) {
      tmp.clear();
      occurrences(rd[st].data(), P[k], tmp);
      for (int64_t d : tmp) occ.push_back({d, static_cast<int>(k)});
      if (occ.size() > kMaxOccurrences) return false;
    }
    w->occurrences += occ.size();
    std::sort(occ.begin(), occ.end());
    size_t i = 0;
    while (i < occ.size()) {  // (ii) chains
      size_t j = i + 1;
      uint32_t mask = 1u << occ[i].second;
      while (j < occ.size() && occ[j].first - occ[j - 1].first <= pad) mask |= 1u << occ[j++].second;
      ++w->chains;
      int distinct = 0;
      for (uint32_t m = mask; m; m &= m - 1) ++distinct;
      if (distinct >= c.t && lemma_l(L, P, mask, 0) > bound) {  // (iii)
        const int64_t dlo = occ[i].first - pad, dhi = occ[j - 1].first + pad;
        std::pair<int64_t, int64_t> win[2] = {{0, ref_len_}, {0, 0}};
        if (sub && st == c.aln.strand) {  // wholly left or wholly right of the record's interval
          win[0] = {0, p0};
          win[1] = {p1, ref_len_};
        }
        for (const auto& wd : win) {
          if (wd.second <= wd.first) continue;
          ++w->windows;
          if (best_band(rd[st].data(), L, dlo, dhi, wd.first, wd.second, w) > bound) return false;
        }
      }
      i = j;
    }
  }
  return true;
}

int Checker::rescore(const std::vector<uint8_t>& read, const Aln& a) const {
  const std::vector<uint8_t> rd = a.strand ? revcomp(read) : read;
  const int L = static_cast<int>(rd.size());
  const int n = static_cast<int>(a.cigar.size());
  if (n == 0) return kNone;
  int64_t j = a.pos;
  int i = 0, score = 0, aligned = 0;
  for (int k = 0; k < n; ++k) {
    const int len = static_cast<int>(a.cigar[k] >> 4), op = static_cast<int>(a.cigar[k] & 15);
    if (len <= 0) return kNone;
    switch (op) {
      case 4:  // soft clip: read ends only
        if (k != 0 && k != n - 1) return kNone;
        score -= kC;
        i += len;
        break;
      case 0:
        for (int x = 0; x < len; ++x, ++i, ++j) {
          if (i >= L || j < 0 || j >= ref_len_) return kNone;
          score += (rd[i] == ref_[j] && rd[i] < 4) ? kM : -kX;
        }
        aligned += len;
        break;
      case 1:
        score -= kO + kE * len;
        i += len;
        break;
      case 2:
        score -= kO + kE * len;
        j += len;
        break;
      default:
        return kNone;
    }
  }
  if (i != L || aligned == 0) return kNone;
  const int ctg = contig_of(a.pos);
  if (ctg < 0 || contig_of(j - 1) != ctg) return kNone;
  return score;
}

Verdict Checker::check(const Claim& c, Work* w) const {
  const int L = static_cast<int>(c.read.size());
  if (L == 0) return kMalformed;
  const std::vector<uint8_t> rd[2] = {c.read, revcomp(c.read)};
  if (c.mapped) {
    const int s = rescore(c.read, c.aln);
    if (s == kNone) return kMalformed;
    if (s != c.as) return kWrongScore;
  }
  if (c.xl != kNone) {
    if (c.xu != kNone && c.xl > c.xu) return kWrongXL;
    bool shown = c.mapped && c.xl <= c.as;  // the record's own alignment proves OPT >= AS
    if (!shown && c.has_alt) {
      const int s = rescore(c.read, c.alt);
      if (s == kNone) return kMalformed;
      shown = c.xl <= s;
    }
    if (!shown) return kWrongXL;
  }
  if (c.xu != kNone) {
    if (c.mapped && c.as > c.xu) return kXUUnproven;  // its own alignment refutes it
    w->xu_checked = c.xu < L * kM;
    if (!proven(c, rd, c.xu, false, w)) return kXUUnproven;
  }
  if (c.xq > 0) {
    if (c.xq > c.mapq) return kWrongXQ;
    if (c.pair_xq) {
      w->pair_xq = true;
    } else {
      if (!c.mapped || c.xh == kNone) return kXHUnproven;
      w->xh_checked = true;
      if (!proven(c, rd, c.xh, true, w)) return kXHUnproven;
      if (c.xq > mapq_formula(c.as, c.xh, L)) return kWrongXQ;
    }
  }
  return kAccepted;
}

bool Checker::parse_sam(const std::string& line, Claim* c) const {
  if (line.empty() || line[0] == '@') return false;
  std::vector<std::string> f;
  split(line, '\t', f);
  if (f.size() < 11) return false;
  *c = Claim();
  const int flag = std::atoi(f[1].c_str());
  c->mapped = !(flag & 4);
  c->mapq = std::atoi(f[4].c_str());
  for (char ch : f[9]) c->read.push_back(code_of(ch));
  if (c->mapped) {
    auto it = contig_.find(f[2]);
    if (it == contig_.end()) return false;
    c->aln.pos = static_cast<int64_t>(names_.offsets[it->second]) + std::atoll(f[3].c_str()) - 1;
    c->aln.strand = (flag & 16) ? 1 : 0;
    if (!parse_cigar(f[5], c->aln.cigar)) return false;
    if (c->aln.strand) c->read = revcomp(c->read);  // SEQ is reverse-complemented
  }
  std::vector<std::string> parts, iv;
  for (size_t k = 11; k < f.size(); ++k) {
    const std::string& t = f[k];
    if (t.size() < 5) continue;
    const std::string tag = t.substr(0, 5), v = t.substr(5);
    if (tag == "AS:i:") c->as = std::atoi(v.c_str());
    else if (tag == "XL:i:") c->xl = std::atoi(v.c_str());
    else if (tag == "XU:i:") c->xu = std::atoi(v.c_str());
    else if (tag == "XQ:i:") c->xq = std::atoi(v.c_str());
    else if (tag == "XH:i:") c->xh = std::atoi(v.c_str());
    else if (tag == "XT:Z:") c->pair_xq = v == "PR";
    else if (tag == "XW:Z:") {
      split(v, ';', parts);
      if (parts.size() != 3) return false;
      c->t = std::atoi(parts[0].c_str());
      for (int st = 0; st < 2; ++st) {
        if (parts[st + 1].empty()) continue;
        split(parts[st + 1], ',', iv);
        for (const std::string& x : iv) {
          const size_t d = x.find('-');
          if (d == std::string::npos) return false;
          c->probes[st].push_back({std::atoi(x.substr(0, d).c_str()), std::atoi(x.substr(d + 1).c_str())});
        }
      }
    } else if (tag == "XV:Z:") {
      split(v, ',', parts);
      if (parts.size() != 4) return false;
      auto it = contig_.find(parts[0]);
      if (it == contig_.end() || !parse_cigar(parts[3], c->alt.cigar)) return false;
      c->alt.pos = static_cast<int64_t>(names_.offsets[it->second]) + std::atoll(parts[1].c_str()) - 1;
      c->alt.strand = parts[2] == "-" ? 1 : 0;
      c->has_alt = true;
    }
  }
  if (c->mapped && c->as == kNone) return false;
  return true;
}

bool validate_index(const Reference& ref, const Index& ix, const std::string& fasta, int threads,
                    std::string* why) {
  auto fail = [&](const std::string& m) {
    *why = m;
    return false;
  };
  threads = std::max(1, threads);
  const uint8_t* seq = ref.seq.data();
  const uint64_t len = ref.seq.size();
  const int q = ix.q, s = ix.s;
  if (q < 1 || q > 32 || s < 1) return fail("bad q or s");
  const uint64_t n = ix.pos.size();
  const bool k16 = ix.keys16.size() > 0;
  if ((k16 ? ix.keys16.size() : ix.keys.size()) != n) return fail("key and position arrays differ in length");
  const int shift = 2 * q - ix.dir_bits;
  if (ix.dir_bits < 0 || shift < 0 || (k16 && shift > 16)) return fail("bad directory width");
  const uint64_t nb = uint64_t(1) << ix.dir_bits;
  if (ix.dir.size() != nb + 1 || ix.dir[0] != 0 || ix.dir[nb] != n) return fail("directory does not span the entries");
  for (uint64_t b = 0; b < nb; ++b)
    if (ix.dir[b] > ix.dir[b + 1]) return fail("directory not monotone");
  auto key_at = [&](uint64_t p, uint64_t* k) {
    if (p + static_cast<uint64_t>(q) > len) return false;
    uint64_t v = 0;
    for (int i = 0; i < q; ++i) {
      const uint8_t ch = seq[p + i];
      if (ch > 3) return false;
      v = (v << 2) | ch;
    }
    *k = v;
    return true;
  };
  const uint64_t nsamp = len >= static_cast<uint64_t>(q) ? (len - q) / s + 1 : 0;
  std::vector<uint64_t> seen((nsamp + 63) / 64, 0);
  std::mutex mu;
  std::string reason;
  std::atomic<bool> bad{false};
  auto flag = [&](const std::string& m) {
    std::lock_guard<std::mutex> g(mu);
    if (!bad.exchange(true)) reason = m;
  };
  // Every entry: key in its bucket, keys sorted, position sampled, the
  // reference q-mer there equals the key, no position twice.
  parallel(threads, [&](int t) {
    for (uint64_t b = nb * t / threads; b < nb * (t + 1) / threads && !bad; ++b) {
      uint64_t prev = 0;
      for (uint64_t e = ix.dir[b]; e < ix.dir[b + 1]; ++e) {
        const uint64_t key = k16 ? ((b << shift) | ix.keys16[e]) : ix.keys[e];
        if ((key >> shift) != b) return flag("key outside its bucket");
        if (e > ix.dir[b] && key < prev) return flag("keys not sorted within a bucket");
        prev = key;
        const uint64_t p = ix.pos[e];
        uint64_t rk;
        if (p % s) return flag("position not a multiple of s");
        if (!key_at(p, &rk) || rk != key) return flag("entry does not match the reference");
        const uint64_t bit = uint64_t(1) << ((p / s) & 63);
        if (__atomic_fetch_or(&seen[p / s / 64], bit, __ATOMIC_RELAXED) & bit)
          return flag("a position is listed twice");
      }
    }
  });
  if (bad) return fail(reason);
  // Their number equals the number of N-free sampled q-mers: with distinct
  // matching positions, the entries are exactly those q-mers.
  std::vector<uint64_t> cnt(threads, 0);
  parallel(threads, [&](int t) {
    uint64_t k;
    for (uint64_t i = nsamp * t / threads; i < nsamp * (t + 1) / threads; ++i) cnt[t] += key_at(i * s, &k);
  });
  uint64_t total = 0;
  for (uint64_t x : cnt) total += x;
  if (total != n) return fail("the index misses N-free q-mers of the reference");
  // Contigs and N padding.
  const size_t nc = ref.names.size();
  if (ref.offsets.size() != nc || ref.lengths.size() != nc) return fail("contig tables differ in length");
  uint64_t cur = 0;
  for (size_t c = 0; c < nc; ++c) {
    if (ref.offsets[c] < cur || ref.offsets[c] + ref.lengths[c] > len) return fail("contigs overlap or overrun");
    for (uint64_t p = cur; p < ref.offsets[c]; ++p)
      if (seq[p] != 4) return fail("padding between contigs is not N");
    cur = ref.offsets[c] + ref.lengths[c];
  }
  for (uint64_t p = cur; p < len; ++p)
    if (seq[p] != 4) return fail("padding after the last contig is not N");
  if (fasta.empty()) return true;
  std::ifstream in(fasta);
  if (!in) return fail("cannot open " + fasta);
  std::string line;
  int64_t c = -1;
  uint64_t at = 0;
  auto close = [&]() { return c < 0 || at == ref.lengths[c]; };
  while (std::getline(in, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty()) continue;
    if (line[0] == '>') {
      if (!close()) return fail("contig " + ref.names[c] + " differs in length from the FASTA");
      ++c;
      at = 0;
      const std::string name = line.substr(1, line.find_first_of(" \t") - 1);
      if (static_cast<size_t>(c) >= nc || ref.names[c] != name) return fail("contig names or order differ: " + name);
      continue;
    }
    if (c < 0) return fail("FASTA does not start with a header");
    if (at + line.size() > ref.lengths[c]) return fail("contig " + ref.names[c] + " is longer in the FASTA");
    const uint8_t* r = seq + ref.offsets[c] + at;
    for (size_t i = 0; i < line.size(); ++i)
      if (r[i] != code_of(line[i])) return fail("contig " + ref.names[c] + " differs from the FASTA");
    at += line.size();
  }
  if (!close() || static_cast<size_t>(c + 1) != nc) return fail("the FASTA has a different number of contigs");
  return true;
}

}  // namespace check
}  // namespace certa
