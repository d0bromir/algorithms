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
                                                 "xu_unproven", "xh_unproven", "wrong_xq",    "xs_unproven"};

namespace {

// BWA-MEM's default scores: the objective every claim is about.
constexpr int kM = 1, kX = 4, kO = 6, kE = 1, kC = 5;
constexpr int kNeg = -(1 << 28);
constexpr uint64_t kMaxOccurrences = 1u << 20;  // per record; beyond it the claim counts as unproven
constexpr uint64_t kMaxPairCells = uint64_t(1) << 27;  // DP cells per pair proof, likewise

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

// BWA-MEM's raw MAPQ of a pair score margin d (mem_sam_pe), clamped to [0, 60].
int margin_mapq(int d) {
  const int q = static_cast<int>(6.02 * d / kM + 0.499);
  return q < 0 ? 0 : (q > 60 ? 60 : q);
}

// Probes per strand: at most 32, disjoint, increasing, inside the read, each
// long enough to be looked up (>= min_len).
bool probes_ok(const std::vector<Probe>* P, int L, int min_len) {
  for (int st = 0; st < 2; ++st) {
    if (P[st].size() > 32) return false;
    int prev = 0;
    for (const Probe& p : P[st]) {
      if (p.off < prev || p.end > L || p.end - p.off < min_len) return false;
      prev = p.end;
    }
  }
  return true;
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

// Whether the witness proves that every alignment (excl given: every
// alignment not overlapping excl) scores <= bound. See check.h, (i)-(iii).
bool Checker::proven(const Claim& c, const std::vector<uint8_t>* rd, int bound, const Aln* excl, Work* w) const {
  const int L = static_cast<int>(c.read.size());
  if (bound >= L * kM) return true;  // no alignment scores more than L
  if (c.t <= 0 || !probes_ok(c.probes, L, q_ + s_ - 1)) return false;
  int ut = kNeg;
  for (int st = 0; st < 2; ++st) ut = std::max(ut, lemma_l(L, c.probes[st], 0, c.t - 1));
  if (ut > bound) return false;  // (i)
  const int pad = std::max(0, L * kM - bound - kO - kE);
  const bool sub = excl != nullptr;
  int64_t p0 = 0, p1 = 0;
  if (sub) {
    p0 = excl->pos;
    p1 = p0 + span_of(excl->cigar);
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
        if (sub && st == excl->strand) {  // wholly left or wholly right of excl's interval
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
    if (!proven(c, rd, c.xu, nullptr, w)) return kXUUnproven;
  }
  if (c.xq > 0) {
    if (c.xq > c.mapq) return kWrongXQ;
    if (c.pair_xq) {  // max(q_se, min(q_pe, q_se + 40)): q_se here, q_pe in check_pair
      w->pair_xq = true;
      w->q_se = 0;
      if (c.xh != kNone) {
        const Aln& a = c.has_alt ? c.alt : c.aln;  // the single-end alignment
        const int s = c.has_alt ? rescore(c.read, c.alt) : (c.mapped ? c.as : kNone);
        if (s == kNone) return kXHUnproven;
        w->xh_checked = true;
        if (!proven(c, rd, c.xh, &a, w)) return kXHUnproven;
        w->q_se = mapq_formula(s, c.xh, L);
      }
      if (c.xq > w->q_se) {
        w->needs_pair = true;
        return kXSUnproven;
      }
    } else {
      if (!c.mapped || c.xh == kNone) return kXHUnproven;
      w->xh_checked = true;
      if (!proven(c, rd, c.xh, &c.aln, w)) return kXHUnproven;
      if (c.xq > mapq_formula(c.as, c.xh, L)) return kWrongXQ;
    }
  }
  return kAccepted;
}

// For each reference position j in [j0, j1): the best local score of the
// read (one strand) over alignments whose first aligned base is j and whose
// cells all lie in [j0, j1); kNeg where there is none. The Gotoh recurrences
// of best_band, run on the read and the reference backwards (an alignment's
// first aligned base is then where it ends).
void Checker::starts(const uint8_t* rd, int L, int64_t j0, int64_t j1, std::vector<int>& best, Work* w) const {
  j0 = std::max<int64_t>(j0, 0);
  j1 = std::min(j1, ref_len_);
  best.assign(static_cast<size_t>(std::max<int64_t>(0, j1 - j0)), kNeg);
  std::vector<int> H(L + 1, kNeg), D(L + 1, kNeg), Hn(L + 1), Dn(L + 1);
  H[0] = 0;
  for (int64_t j = j1 - 1; j >= j0; --j) {
    const uint8_t r = ref_[j];
    Hn[0] = 0;
    Dn[0] = kNeg;
    int ins = kNeg;
    for (int i = 1; i <= L; ++i) {
      const uint8_t c = rd[L - i];
      const int start = i == 1 ? 0 : -kC;  // read bases after it are clipped
      const int m = std::max(H[i - 1], start) + ((c == r && r < 4) ? kM : -kX);
      ins = i > 1 ? std::max(Hn[i - 1] - kO - kE, ins - kE) : kNeg;
      Dn[i] = std::max(H[i] - kO - kE, D[i] - kE);
      Hn[i] = std::max(m, std::max(ins, Dn[i]));
      best[j - j0] = std::max(best[j - j0], m - (i < L ? kC : 0));  // read bases before it are clipped
    }
    std::swap(H, Hn);
    std::swap(D, Dn);
    w->cells += static_cast<uint64_t>(L);
  }
}

// Every proper pair that differs from the records' pair in some mate scores
// <= bound, from the pair witness (check.h, P1-P3).
bool Checker::pair_proven(const Claim& a, const Claim& b, int bound, Work* w) const {
  const Claim* m[2] = {&a, &b};
  const int64_t D = a.pd;
  int L[2], U[2][2];
  std::vector<uint8_t> rd[2][2];
  for (int x = 0; x < 2; ++x) {
    L[x] = static_cast<int>(m[x]->read.size());
    rd[x][0] = m[x]->read;
    rd[x][1] = revcomp(m[x]->read);
    if (m[x]->pt <= 0 || !probes_ok(m[x]->pprobes, L[x], q_ + s_ - 1)) return false;
  }
  if (bound >= (L[0] + L[1]) * kM) return true;
  for (int x = 0; x < 2; ++x)
    for (int st = 0; st < 2; ++st) U[x][st] = lemma_l(L[x], m[x]->pprobes[st], 0, m[x]->pt - 1);
  const int u0 = std::max(U[0][0], U[0][1]), u1 = std::max(U[1][0], U[1][1]);
  if (u0 + u1 > bound) return false;  // (P1)
  const int pad = std::max(0, (L[0] + L[1]) * kM - bound - kO - kE);
  // (P2) Chains of each mate's part occurrences, per strand: diagonals
  // [first, last], the Lemma L bound with their parts free.
  struct Chain {
    int64_t first, last;
    int distinct, ub;
  };
  std::vector<Chain> ch[2][2];
  std::vector<std::pair<int64_t, int>> occ;
  std::vector<int64_t> tmp;
  for (int x = 0; x < 2; ++x)
    for (int st = 0; st < 2; ++st) {
      const std::vector<Probe>& P = m[x]->pprobes[st];
      occ.clear();
      for (size_t k = 0; k < P.size(); ++k) {
        tmp.clear();
        occurrences(rd[x][st].data(), P[k], tmp);
        for (int64_t d : tmp) occ.push_back({d, static_cast<int>(k)});
        if (occ.size() > kMaxOccurrences) return false;
      }
      w->occurrences += occ.size();
      std::sort(occ.begin(), occ.end());
      for (size_t i = 0, j; i < occ.size(); i = j) {
        j = i + 1;
        uint32_t mask = 1u << occ[i].second;
        while (j < occ.size() && occ[j].first - occ[j - 1].first <= pad) mask |= 1u << occ[j++].second;
        int distinct = 0;
        for (uint32_t u = mask; u; u &= u - 1) ++distinct;
        if (distinct >= m[x]->pt) ch[x][st].push_back({occ[i].first, occ[j - 1].first, distinct, lemma_l(L[x], P, mask, 0)});
      }
      w->chains += ch[x][st].size();
    }
  // (P3) A pair above bound in which mate x has t exact parts: x lies in a
  // chain C of x (its cells on diagonals [first - pad, last + pad], its first
  // aligned base in [a_lo, a_hi]), y on the other strand starts within D of
  // it. Every such pair is a real proper pair, so the exact maximum over
  // them, with the records' pair excluded, must not exceed bound.
  std::vector<int> xa, xo, ya, yo;
  for (int x = 0; x < 2; ++x) {
    const int y = 1 - x;
    for (int st = 0; st < 2; ++st) {
      const int sy = 1 - st;
      for (const Chain& c : ch[x][st]) {
        const int64_t a_lo = c.first - pad, a_hi = c.last + pad + L[x] - 1;
        const int64_t s_lo = st == 0 ? a_lo : a_lo - D, s_hi = st == 0 ? a_hi + D : a_hi;
        // y in [s_lo, s_hi]: no exact part (<= U), or in a chain of y whose
        // start range meets it (<= its bound).
        int yub = U[y][sy];
        for (const Chain& e : ch[y][sy])
          if (e.first - pad <= s_hi && e.last + pad + L[y] - 1 >= s_lo) yub = std::max(yub, e.ub);
        if (c.ub + yub <= bound) continue;
        if (w->cells > kMaxPairCells) return false;
        ++w->windows;
        const int64_t x0 = a_lo, x1 = a_hi + L[x] + pad;  // x's cells
        starts(rd[x][st].data(), L[x], x0, x1, xa, w);
        int bx = kNeg;
        for (int v : xa) bx = std::max(bx, v);
        if (bx + yub <= bound) continue;
        const int64_t y0 = s_lo, y1 = s_hi + L[y] + pad;  // y's cells
        starts(rd[y][sy].data(), L[y], y0, y1, ya, w);
        // The records' pair excluded: x wholly left or right of its record,
        // or y of its record (when on their strands).
        const bool ox = st == m[x]->aln.strand, oy = sy == m[y]->aln.strand;
        auto outside = [&](const uint8_t* r, int len, int64_t j0, int64_t j1, const Aln& ex, std::vector<int>& out) {
          const int64_t p0 = ex.pos, p1 = p0 + span_of(ex.cigar);
          std::vector<int> left, right;
          starts(r, len, j0, std::min(j1, p0), left, w);
          starts(r, len, std::max(j0, p1), j1, right, w);
          const int64_t b0 = std::max<int64_t>(j0, 0);
          out.assign(static_cast<size_t>(std::max<int64_t>(0, std::min(j1, ref_len_) - b0)), kNeg);
          for (size_t k = 0; k < left.size(); ++k) out[k] = left[k];
          const int64_t rb = std::max<int64_t>(std::max(j0, p1), 0) - b0;
          for (size_t k = 0; k < right.size(); ++k) out[static_cast<size_t>(rb) + k] = right[k];
        };
        if (ox) outside(rd[x][st].data(), L[x], x0, x1, m[x]->aln, xo);
        if (oy) outside(rd[y][sy].data(), L[y], y0, y1, m[y]->aln, yo);
        // For x starting at f, y starts in [f, f + D] (x forward) or
        // [f - D, f] (x reverse): sliding maxima of ya and yo.
        const int64_t xb = std::max<int64_t>(x0, 0), yb = std::max<int64_t>(y0, 0);
        const int64_t ny = static_cast<int64_t>(ya.size());
        std::vector<int> ma, mo;  // per f
        auto slide = [&](const std::vector<int>& v, std::vector<int>& out) {
          out.assign(xa.size(), kNeg);
          std::vector<int64_t> dq;  // indices into v, values decreasing
          size_t head = 0;
          int64_t next = 0;
          for (size_t f = 0; f < xa.size(); ++f) {
            const int64_t pf = xb + static_cast<int64_t>(f);
            const int64_t lo = (st == 0 ? pf : pf - D) - yb, hi = (st == 0 ? pf + D : pf) - yb;
            while (next < ny && next <= hi) {
              while (dq.size() > head && v[dq.back()] <= v[next]) dq.pop_back();
              dq.push_back(next++);
            }
            while (head < dq.size() && dq[head] < lo) ++head;
            if (head < dq.size()) out[f] = v[dq[head]];
          }
        };
        slide(ya, ma);
        if (ox || oy) slide(oy ? yo : ya, mo);
        int pairs = kNeg;
        for (size_t f = 0; f < xa.size(); ++f) {
          if (!ox) pairs = std::max(pairs, xa[f] + ma[f]);
          else pairs = std::max(pairs, std::max(xo[f] + ma[f], xa[f] + mo[f]));
        }
        if (pairs > bound) return false;
      }
    }
  }
  return true;
}

void Checker::check_pair(const Claim& a, const Claim& b, Work* wa, Work* wb, Verdict* va, Verdict* vb) const {
  *va = check(a, wa);
  *vb = check(b, wb);
  const bool need[2] = {*va == kXSUnproven && wa->needs_pair, *vb == kXSUnproven && wb->needs_pair};
  if (!need[0] && !need[1]) return;
  // The shown pair: both alignments verified (check got past AS), proper.
  bool pair = a.mapped && b.mapped && a.pair_xq && b.pair_xq && a.pd >= 0 && a.pd == b.pd &&
              (*va == kAccepted || need[0]) && (*vb == kAccepted || need[1]) &&
              a.aln.strand != b.aln.strand && contig_of(a.aln.pos) >= 0 &&
              contig_of(a.aln.pos) == contig_of(b.aln.pos);
  if (pair) {
    const Aln& f = a.aln.strand ? b.aln : a.aln;
    const Aln& r = a.aln.strand ? a.aln : b.aln;
    pair = r.pos >= f.pos && r.pos - f.pos <= a.pd;
  }
  const Claim* c[2] = {&a, &b};
  Work* w[2] = {wa, wb};
  Verdict* v[2] = {va, vb};
  int done_bound = kNone;
  bool done_ok = false;
  for (int x = 0; x < 2; ++x) {
    if (!need[x]) continue;
    const int bound = c[x]->xs;
    bool ok = false;
    if (pair && bound != kNone) {
      if (bound != done_bound) {  // both mates normally carry the same XS
        done_ok = pair_proven(a, b, bound, w[x]);
        done_bound = bound;
      }
      ok = done_ok;
    }
    w[x]->xs_checked = ok;
    const int q_pe = ok ? margin_mapq(a.as + b.as - bound) : 0;
    const int q = std::max(w[x]->q_se, std::min(q_pe, w[x]->q_se + 40));
    *v[x] = c[x]->xq <= q ? kAccepted : (ok ? kWrongXQ : kXSUnproven);
  }
}

bool Checker::parse_sam(const std::string& line, Claim* c) const {
  if (line.empty() || line[0] == '@') return false;
  std::vector<std::string> f;
  split(line, '\t', f);
  if (f.size() < 11) return false;
  *c = Claim();
  c->qname = f[0];
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
  // "o-e,o-e,..." -> probes
  auto probe_list = [&](const std::string& s, std::vector<Probe>& out) {
    if (s.empty()) return true;
    split(s, ',', iv);
    for (const std::string& x : iv) {
      const size_t d = x.find('-');
      if (d == std::string::npos) return false;
      out.push_back({std::atoi(x.substr(0, d).c_str()), std::atoi(x.substr(d + 1).c_str())});
    }
    return true;
  };
  for (size_t k = 11; k < f.size(); ++k) {
    const std::string& t = f[k];
    if (t.size() < 5) continue;
    const std::string tag = t.substr(0, 5), v = t.substr(5);
    if (tag == "AS:i:") c->as = std::atoi(v.c_str());
    else if (tag == "XL:i:") c->xl = std::atoi(v.c_str());
    else if (tag == "XU:i:") c->xu = std::atoi(v.c_str());
    else if (tag == "XQ:i:") c->xq = std::atoi(v.c_str());
    else if (tag == "XH:i:") c->xh = std::atoi(v.c_str());
    else if (tag == "XS:i:") c->xs = std::atoi(v.c_str());
    else if (tag == "XT:Z:") c->pair_xq = v == "PR";
    else if (tag == "XW:Z:") {  // t;strand 0;strand 1
      split(v, ';', parts);
      if (parts.size() != 3) return false;
      c->t = std::atoi(parts[0].c_str());
      for (int st = 0; st < 2; ++st)
        if (!probe_list(parts[st + 1], c->probes[st])) return false;
    } else if (tag == "XY:Z:") {  // D;t;strand 0;strand 1
      split(v, ';', parts);
      if (parts.size() != 4) return false;
      c->pd = std::atoi(parts[0].c_str());
      c->pt = std::atoi(parts[1].c_str());
      for (int st = 0; st < 2; ++st)
        if (!probe_list(parts[st + 2], c->pprobes[st])) return false;
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
