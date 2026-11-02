// CERTA, biologically weighted placement (certa map --bio; host only).
//
// The default certificates are about an optimization problem: the fewest
// edits (S0/S1) or BWA-MEM's score (SL, PR). Neither models how a read
// arises, so their optimum need not be the read's origin: two mismatches at
// Q2 bases are likelier than one at a Q37 base, a transition likelier than a
// transversion, an indel in a homopolymer far likelier than one in unique
// sequence.
//
// Here the objective is the read's likelihood. The read was sequenced from
// the sample's genome, which differs from the reference by SNVs (rate snv,
// transitions titv times as frequent as transversions) and small indels
// (rate indel, plus sequencing indels at seq_indel, geometric lengths,
// opening cheaper in tandem repeats), and each base call is wrong with the
// probability its quality states. An alignment's penalty is -10 log10 of its
// likelihood relative to a perfect match, in integer Phred units:
//   - a match costs 0; a read N costs 0 against anything (it says nothing);
//   - read base b at quality Q on reference base a != b costs
//     10 log10(P(a | a) / P(b | a)), from the SNV and error rates (at least 1);
//     a reference N costs the cheapest mismatch of that read base (as in the
//     edit certificates, N never matches);
//   - a gap run of k bases costs go(ctx) + (k - 1) ge, where ctx is the
//     length of the longest tandem run (period 1-4) the gap touches; each
//     inserted base also costs 10 log10(4 P(match)), its emission.
// The best alignment over the whole reference, both strands, end to end, is
// the maximum-likelihood placement. A locus's posterior uses its best
// alignment (Viterbi) and a uniform prior over loci.
//
// Theorem B (certified maximum-likelihood placement; THEORY.md §12). Let S
// be parts of the read with every exact occurrence enumerated (Lemma I), on
// each strand. An alignment that leaves a part of S exactly occurring has a
// candidate diagonal there (Theorem 1). One that touches every part of S
// costs at least cover(S), the cheapest set of events that touches them all
// (bio_cover: a mismatch or a gap inside a part touches that part; only an
// insertion run touches several, and it inserts every base between them).
// One whose gaps total more than G costs at least go_min + G ge. So with
// H = min(cover(S) per strand, go_min + G ge), every alignment below H lies
// in the band of a cluster of candidates (diagonals +- G); an alignment in
// a cluster costs at least cover(S \ parts with a hit there). Bands are
// evaluated exactly by DP, cheapest bound first; those not evaluated (bound
// at or above p1 + delta, or wider than max_band) lower the horizon to their
// bound. Then every end position whose best alignment costs less than the
// horizon YH is known exactly, and if the best of them, p1, is below YH,
// p1 is the global optimum and every locus below YH is listed.
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <utility>
#include <vector>
#if defined(__aarch64__)
#include <arm_neon.h>
#elif defined(__SSSE3__)
#include <tmmintrin.h>
#endif

#include "certa/core.h"

namespace certa {

constexpr int kBioInf = 1 << 28;
constexpr int kBioCtxCap = 50;  // tandem runs are measured up to this many bases
constexpr int kBioQMin = 2, kBioQMax = 60;

struct BioParams {
  double snv = 1e-3;        // SNVs per base between sample and reference
  double titv = 2.1;        // transitions per transversion
  double indel = 1e-4;      // small indels per base (half insertions, half deletions)
  double seq_indel = 1e-5;  // sequencing indel errors per base
  double gap_ext = 0.5;     // an indel continues one more base with this probability
  int str_min = 4;          // tandem runs of this many bases or more make opening cheaper,
  int str_slope = 3;        // by this many Phred per base,
  int str_cap = 24;         // at most this much
  int gap_cap = 50;         // largest band reach G (gap total) a read may need
  int delta = 60;           // clusters bounded at p1 + delta or more are not aligned
  int max_band = 1024;      // widest band aligned (wider ones are bounded only)
  int max_loci = 32;        // loci kept per read for pairing
  double insert_mean = 400, insert_sd = 100;  // fragment length (pairs)
  int improper = 40;        // a pair that is not proper, in Phred
  int rescue_loci = 8;      // a certified mate's loci whose proper windows a pair rescue aligns
  int lanes = 16;           // bands aligned at once, one SIMD lane each (bio_bands; 1: one at a time)
};

// The model's integer costs, from BioParams.
class BioModel {
 public:
  explicit BioModel(const BioParams& p = BioParams()) : p_(p) {
    const double th = p.snv;
    const double wti = p.titv / (1 + p.titv), wtv = (1 - wti) / 2;
    for (int q = 0; q < 64; ++q) {
      const int qq = std::min(kBioQMax, std::max(kBioQMin, q));
      const double e = std::pow(10.0, -qq / 10.0);
      const double match = (1 - th) * (1 - e) + th * e / 3;
      int lo = kBioInf;
      for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b) {
          if (a == b) { mis_[q][a][b] = 0; continue; }
          const double w = (a ^ b) == 2 ? wti : wtv;  // A<->G and C<->T are transitions
          const double pm = (1 - th) * e / 3 + th * (w * (1 - e) + (1 - w) * e / 3);
          const int c = std::max(1, static_cast<int>(std::lround(10 * std::log10(match / pm))));
          mis_[q][a][b] = static_cast<int16_t>(c);
          lo = std::min(lo, c);
        }
      mmin_[q] = static_cast<int16_t>(lo);
      emit_[q] = static_cast<int16_t>(std::max(0L, std::lround(10 * std::log10(4 * match))));
    }
    const double open = (p.indel + p.seq_indel) / 2 * (1 - p.gap_ext);
    const int go = static_cast<int>(std::lround(-10 * std::log10(open)));
    ge_ = std::max(1, static_cast<int>(std::lround(-10 * std::log10(p.gap_ext))));
    // Opening never costs less than extending, so one run is the cheapest way
    // to a given gap total (Theorem B's band reach).
    go_min_ = kBioInf;
    for (int c = 0; c <= kBioCtxCap; ++c) {
      const int disc = c >= p.str_min ? std::min(p.str_cap, p.str_slope * (c - p.str_min + 1)) : 0;
      go_[c] = std::max(ge_, go - disc);
      go_min_ = std::min(go_min_, go_[c]);
    }
    // go never rises with c, so it stays at go(kBioCtxCap) from ctx_sat_ on.
    ctx_sat_ = kBioCtxCap;
    while (ctx_sat_ > 1 && go_[ctx_sat_ - 1] == go_[kBioCtxCap]) --ctx_sat_;
    lut_ok_ = true;
    for (int c = 0; c < 16; ++c) {
      lut_ok_ = lut_ok_ && go_[c] <= 255;
      go_lut_[c] = static_cast<uint8_t>(std::min(go_[c], 255));
    }
    max_add_ = std::max(ge_, go_[0]);
    for (int q = 0; q < 64; ++q) {
      max_add_ = std::max<int>(max_add_, emit_[q]);
      for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b) max_add_ = std::max<int>(max_add_, mis_[q][a][b]);
    }
  }
  const BioParams& params() const { return p_; }
  // Read base b (0-3) at quality q on reference base a (0-3), a != b.
  int mis(int q, int a, int b) const { return mis_[q][a][b]; }
  int mmin(int q) const { return mmin_[q]; }  // the cheapest mismatch at quality q
  int emit(int q) const { return emit_[q]; }  // an inserted base at quality q
  int go(int ctx) const { return go_[ctx]; }  // opening a gap that touches a tandem run of ctx bases
  int ge() const { return ge_; }
  int go_min() const { return go_min_; }
  // The tandem-run length from which a gap opens at go_min: a run length
  // clamped there costs every gap the same.
  int ctx_sat() const { return ctx_sat_; }
  // Optional: the reference's tandem runs, clamped at ctx_sat() (bio_pack_contexts);
  // bio_band reads them there instead of scanning its window. Null: no table.
  void set_ctx_table(const uint8_t* t) { tab_ = t; }
  const uint8_t* ctx_table() const { return tab_; }
  // go(c) for c < 16 as bytes (a run the table holds is below 16), or null
  // if one does not fit a byte: bio_opens8's lookup table.
  const uint8_t* go_lut() const { return lut_ok_ ? go_lut_ : nullptr; }
  // The largest single cost the DP adds (a base, an emission, a gap open or
  // extension).
  int max_add() const { return max_add_; }

 private:
  BioParams p_;
  int16_t mis_[64][4][4];
  int16_t mmin_[64], emit_[64];
  int go_[kBioCtxCap + 1];
  uint8_t go_lut_[16];
  bool lut_ok_;
  int ge_, go_min_, ctx_sat_, max_add_;
  const uint8_t* tab_ = nullptr;
};

// The read's costs on both strands (strand 1: the reverse complement, its
// qualities reversed). cost[st][i][a]: read base i on reference code a
// (4 = N); emit: as an inserted base; mmin: its cheapest mismatch (0 for N).
// lane[i]: read base i's costs for bio_bands, both strands in one 16-byte
// table: [5 st + a] = cost[st][i][a], [10 + st] = emit[st][i] (each at most
// 255).
struct BioRead {
  int L = 0;
  uint8_t seq[2][LMAX];
  int16_t cost[2][LMAX][5];
  int16_t emit[2][LMAX];
  int16_t mmin[2][LMAX];
  uint8_t lane[LMAX][16];
};

// `qual`: FASTQ qualities (Phred + 33), or null for Q30 everywhere.
inline void bio_read(const BioModel& m, const uint8_t* codes, const char* qual, int L, BioRead& br) {
  br.L = L;
  auto fill = [&](int st, int i, uint8_t b, int q) {
    br.seq[st][i] = b;
    int16_t* c = br.cost[st][i];
    if (b > 3) {
      for (int a = 0; a < 5; ++a) c[a] = 0;
      br.emit[st][i] = br.mmin[st][i] = 0;
      return;
    }
    for (int a = 0; a < 4; ++a) c[a] = static_cast<int16_t>(a == b ? 0 : m.mis(q, a, b));
    c[4] = static_cast<int16_t>(m.mmin(q));
    br.emit[st][i] = static_cast<int16_t>(m.emit(q));
    br.mmin[st][i] = static_cast<int16_t>(m.mmin(q));
  };
  for (int i = 0; i < L; ++i) {
    const uint8_t b = codes[i];
    const int q = qual ? std::min(63, std::max(0, qual[i] - 33)) : 30;
    fill(0, i, b, q);
    fill(1, L - 1 - i, b < 4 ? static_cast<uint8_t>(3 - b) : uint8_t(4), q);
  }
  for (int i = 0; i < L; ++i) {
    uint8_t* t = br.lane[i];
    for (int st = 0; st < 2; ++st) {
      for (int a = 0; a < 5; ++a) t[5 * st + a] = static_cast<uint8_t>(std::min<int>(255, br.cost[st][i][a]));
      t[10 + st] = static_cast<uint8_t>(std::min<int>(255, br.emit[st][i]));
    }
    for (int x = 12; x < 16; ++x) t[x] = 0;
  }
}

// A position's tandem run depends on the reference this far on each side.
constexpr int kBioCtxReach = kBioCtxCap + 4;

// buf[0, n) = the reference from lo (N outside it).
inline void bio_fetch(const IndexView& ix, int64_t lo, int n, uint8_t* buf) {
  if (lo >= 0 && static_cast<uint64_t>(lo) + static_cast<uint64_t>(n) <= ix.ref_len) {
    std::memcpy(buf, ix.ref + lo, static_cast<size_t>(n));
    return;
  }
  for (int x = 0; x < n; ++x) buf[x] = ref_at(ix, lo + x);
}

// bio_str_context below on bases already fetched: buf holds the reference
// from a - kBioCtxReach to b + kBioCtxReach (n = b - a).
inline void bio_str_context_buf(const uint8_t* buf, int n, uint8_t* ctx) {
  const int E = kBioCtxReach, hi = n + 2 * E;
  for (int x = 0; x < n; ++x) ctx[x] = 1;
  for (int p = 1; p <= 4; ++p) {
    int s = -1;  // start of the current stretch of matches
    for (int x = p; x <= hi; ++x) {
      if (x < hi && buf[x] < 4 && buf[x] == buf[x - p]) {
        if (s < 0) s = x;
        continue;
      }
      if (s >= 0) {  // stretch [s, x) of matches: run [s - p, x)
        const int r0 = s - p, len = std::min(kBioCtxCap, x - r0);
        for (int y = std::max(r0, E), e = std::min(x, E + n); y < e; ++y)
          if (len > ctx[y - E]) ctx[y - E] = static_cast<uint8_t>(len);
        s = -1;
      }
    }
  }
}

// ctx[x - a] for reference positions x in [a, b): the length (capped at
// kBioCtxCap) of the longest tandem run of period 1-4 containing x, at least
// 1. A run of period p is a stretch where every base equals the one p before
// it, plus those p; N ends it. Depends only on the reference within
// kBioCtxReach of x.
inline void bio_str_context(const IndexView& ix, int64_t a, int64_t b, uint8_t* ctx) {
  const int n = static_cast<int>(b - a), len = n + 2 * kBioCtxReach;
  uint8_t local[512];
  std::vector<uint8_t> big;
  uint8_t* buf = local;
  if (len > static_cast<int>(sizeof local)) {
    big.resize(static_cast<size_t>(len));
    buf = big.data();
  }
  bio_fetch(ix, a - kBioCtxReach, len, buf);
  bio_str_context_buf(buf, n, ctx);
}

// The table bio_band can read tandem runs from: position x's run length
// (bio_str_context), clamped at `clamp` (<= kBioCtxPacked), four bits a base,
// two to a byte (x even in the low half), (ref_len + 1) / 2 bytes in all.
// Since a run length depends only on the reference within kBioCtxReach,
// stretches [a, b) can be packed independently (a even): together they give
// the same table.
constexpr int kBioCtxPacked = 15;
static_assert(kBioCtxPacked < 16 && kBioCtxPacked <= kBioCtxCap, "a packed run must index BioModel::go_lut");

inline void bio_pack_contexts(const IndexView& ix, int64_t a, int64_t b, int clamp, uint8_t* tab) {
  constexpr int64_t kStep = 1 << 16;
  std::vector<uint8_t> ctx(static_cast<size_t>(kStep));
  for (int64_t x0 = a; x0 < b; x0 += kStep) {
    const int64_t x1 = std::min(b, x0 + kStep);
    bio_str_context(ix, x0, x1, ctx.data());
    for (int64_t x = x0; x < x1; ++x) {
      const uint8_t c = static_cast<uint8_t>(std::min<int>(ctx[x - x0], clamp));
      uint8_t& t = tab[x >> 1];
      t = (x & 1) ? static_cast<uint8_t>((t & 0x0F) | (c << 4)) : static_cast<uint8_t>((t & 0xF0) | c);
    }
  }
}

// Position x's clamped run length from the table (1 outside the reference,
// as bio_str_context gives N).
inline uint8_t bio_ctx_at(const uint8_t* tab, uint64_t ref_len, int64_t x) {
  if (x < 0 || static_cast<uint64_t>(x) >= ref_len) return 1;
  return static_cast<uint8_t>((tab[x >> 1] >> ((x & 1) * 4)) & 0x0F);
}

// A lower bound on the penalty of any alignment of strand st that touches
// every part in `req` (bit j = part j of the read's P parts). A part is
// touched when it does not occur exactly: a base mismatched (or on an N), a
// base inserted, or a deletion between two of its bases. A mismatch or a
// deletion touches one part; a deletion between parts touches none. An
// insertion run touching parts x < y inserts the last base of x, every base
// between and the first base of y. A DP over the parts in read order picks
// the cheapest such events covering `req` (THEORY.md §12, Lemma B).
inline int bio_cover(const BioModel& m, const BioRead& br, int st, int P, uint32_t req) {
  if (!req) return 0;
  const int L = br.L;
  int off[PMAX], len[PMAX], c1[PMAX];
  int pre[LMAX + 1];  // prefix sums of emit
  pre[0] = 0;
  for (int i = 0; i < L; ++i) pre[i + 1] = pre[i] + br.emit[st][i];
  for (int j = 0; j < P; ++j) {
    part_geometry(L, P, j, &off[j], &len[j]);
    int mn = kBioInf, em = kBioInf;
    for (int i = off[j]; i < off[j] + len[j]; ++i) {
      mn = std::min<int>(mn, br.mmin[st][i]);
      em = std::min<int>(em, br.emit[st][i]);
    }
    c1[j] = std::min(mn, m.go_min() + em);
    if (len[j] >= 2) c1[j] = std::min(c1[j], m.go_min());
  }
  int f[PMAX + 1];
  f[0] = 0;
  for (int x = 1; x <= P; ++x) f[x] = kBioInf;
  for (int x = 0; x < P; ++x) {
    if (f[x] >= kBioInf) continue;
    f[x + 1] = std::min(f[x + 1], f[x] + (((req >> x) & 1) ? c1[x] : 0));
    const int last = off[x] + len[x] - 1;
    for (int y = x + 1; y < P; ++y) {
      const int k = off[y] - last + 1;
      const int c = m.go_min() + (k - 1) * m.ge() + pre[off[y] + 1] - pre[last];
      f[y + 1] = std::min(f[y + 1], f[x] + c);
    }
  }
  return f[P];
}

// A candidate band: diagonals [lo, lo + width) of one strand.
struct BioCluster {
  int64_t lo;
  int width;
  uint8_t strand;
  uint16_t pmask;  // parts with a candidate in the band
  int bound;       // every alignment in the band costs at least this
};

// An end position below the horizon: the best alignment ending there.
struct BioEnd {
  int pen;
  uint8_t strand;
  int64_t end;  // reference position of its last aligned base
  int cl, b;    // its band (index into BioScratch::cl) and end diagonal lo + b
};

// One locus a read can be placed at: an end position whose best alignment
// is cheapest among those within L of it on its strand.
struct BioLocus {
  int pen;
  uint8_t strand;
  int64_t end;
  int64_t lo;  // its band, for the traceback
  int width, b;
  int amb = 0;  // end positions within L costing the same span this many bases (tandem shifts)
};

// What --bio found for one read.
struct BioOut {
  bool ok = false;      // p1 < yh: the reported locus is the certified ML placement
  bool tried = false;   // the parts were enumerated on both strands
  bool empty = false;   // tried, and no alignment costs less than yh
  bool rescued = false; // placed by its mate's loci (bio_rescue), not certified alone
  bool wide = false;    // a band too wide to align set the horizon
  bool moved = false;   // unique, and the edit certificate's unique best locus is elsewhere
  int p1 = -1, p2 = -1; // best penalty; second locus's (p1 when tied), -1 if none below yh
  int yh = -1;          // every end position below yh is accounted for
  int n_best = 0;       // loci tied at p1
  int mapq = 0;         // posterior, Phred
  std::vector<BioLocus> loci;  // the cheapest loci (at most max_loci, the chosen one included)
  int list_bound = -1;  // every locus not in `loci` costs at least this
  int chosen = -1;      // index of the reported locus in `loci`
  // Pairs: placed jointly when the pair's optimum is certified.
  bool pair_ok = false;
  int pair_pen = -1, ins_pen = -1, pair_mapq = 0;
  // The work behind it (tools/bio_probe): |S| per strand, the horizon H
  // before any band lowered it, the band reach G, candidates, clusters,
  // bands aligned and their DP cells (L x width).
  int n_s[2] = {0, 0};
  int horizon = -1, reach = -1;
  int n_cand = 0, n_clusters = 0, n_bands = 0;
  int64_t cells = 0;
};

// bio_cover(need) of one read, strand by strand, by required part set (a
// mask below 2^P, P <= PMAX): an array, stamped so that a new read clears
// nothing.
class BioCoverCache {
 public:
  void start(int P) {
    const size_t n = size_t(1) << P;
    for (int st = 0; st < 2; ++st)
      if (val_[st].size() < n) {
        val_[st].resize(n);
        stamp_[st].resize(n, 0);
      }
    if (++gen_ == 0) {  // wrapped: every stamp is stale again
      for (auto& s : stamp_) std::fill(s.begin(), s.end(), 0);
      gen_ = 1;
    }
  }
  template <class F>
  int get(int st, uint32_t need, F compute) {
    if (stamp_[st][need] != gen_) {
      stamp_[st][need] = gen_;
      val_[st][need] = compute();
    }
    return val_[st][need];
  }

 private:
  std::vector<int> val_[2];
  std::vector<uint32_t> stamp_[2];
  uint32_t gen_ = 0;
};

// Pair rescue (bio_rescue): a pair found, and an unseen pair's bound.
struct BioRescueCand {
  int cost, a, v, ins;
  BioLocus f;
};

struct BioScratch {
  BioRead br;
  std::vector<int> rows;               // 6 DP rows
  // The band's reference (from kBioCtxReach before its first base), tandem
  // runs, traceback bits; gap opens: a deletion of the base, an insertion
  // after it.
  std::vector<uint8_t> win, ctx, tb;
  std::vector<int> go_d, go_i;
  std::vector<uint8_t> runs, od, oi;  // bio_opens8: the runs unpacked; bio_window: the opens as bytes
  std::vector<int> endv;               // the last band's end row
  std::vector<BioCluster> cl;
  std::vector<int> order;
  std::vector<BioEnd> ends;
  std::vector<BioLocus> loci;
  std::vector<int64_t> taken[2];  // the loci's end positions so far, per strand, sorted
  BioCoverCache cover;
  // bio_rescue: pairs found, f's loci, the bounds (k's locus or -1), the
  // loci's (k's locus, f's end) so far, sorted.
  std::vector<BioRescueCand> rcand, rloci;
  std::vector<std::pair<int, int>> rbounds;
  std::vector<std::pair<int, int64_t>> rtaken;
  // bio_bands: per window position, 16 bytes (one a lane) of reference codes
  // and gap opens; per diagonal, the lanes outside their band; 6 DP rows; the
  // bands' ends (lane l's from lend[l * lstride]). ln, lq: the bands aligned
  // ahead (bio_place: their positions in `order`; bio_rescue: k's loci).
  // bio_lane_setup writes each band's window as a row of lstage first.
  std::vector<uint8_t> lcode, lgo_d, lgo_i, lmask, lrows, lstage;
  std::vector<int> lend;
  int lstride = 0, ln = 0;
  int lq[16];
};

// Whether a sorted list holds a position within L of `end` (|t - end| < L).
inline bool bio_near(const std::vector<int64_t>& v, int64_t end, int L) {
  const auto it = std::lower_bound(v.begin(), v.end(), end - L + 1);
  return it != v.end() && *it < end + L;
}

namespace bio_detail {

typedef uint8_t u8x16 __attribute__((vector_size(16)));
typedef uint8_t u8x16u __attribute__((vector_size(16), aligned(1), may_alias));
typedef uint64_t u64x2 __attribute__((vector_size(16)));

inline u8x16 splat8(int x) {
  const uint8_t s = static_cast<uint8_t>(x);
  return u8x16{s, s, s, s, s, s, s, s, s, s, s, s, s, s, s, s};
}
inline u8x16 load8(const uint8_t* p) { return *reinterpret_cast<const u8x16u*>(p); }
inline void store8(uint8_t* p, u8x16 v) { *reinterpret_cast<u8x16u*>(p) = v; }
inline u8x16 vmin8(u8x16 a, u8x16 b) { return a < b ? a : b; }
inline u8x16 vmax8(u8x16 a, u8x16 b) { return a > b ? a : b; }
inline bool any_below(u8x16 v, u8x16 cap) {
  const u64x2 w = (u64x2)(v < cap);
  return (w[0] | w[1]) != 0;
}
// t[idx[l]] in every lane (idx < 16).
inline u8x16 lookup8(u8x16 t, u8x16 idx) {
#if defined(__aarch64__)
  return (u8x16)vqtbl1q_u8((uint8x16_t)t, (uint8x16_t)idx);
#elif defined(__SSSE3__)
  return (u8x16)_mm_shuffle_epi8((__m128i)t, (__m128i)idx);
#else
  u8x16 r;
  for (int l = 0; l < 16; ++l) r[l] = t[idx[l] & 15];
  return r;
#endif
}
// The bytes of a and b interleaved (a0 b0 a1 b1 ...), from their low halves
// (zip_lo) or their high halves (zip_hi).
inline u8x16 zip_lo(u8x16 a, u8x16 b) {
  return __builtin_shufflevector(a, b, 0, 16, 1, 17, 2, 18, 3, 19, 4, 20, 5, 21, 6, 22, 7, 23);
}
inline u8x16 zip_hi(u8x16 a, u8x16 b) {
  return __builtin_shufflevector(a, b, 8, 24, 9, 25, 10, 26, 11, 27, 12, 28, 13, 29, 14, 30, 15, 31);
}
// A 16 x 16 block of bytes transposed: v[j][l] becomes the old v[l][j]. A
// round interleaves rows r and r + 8 into rows 2r and 2r + 1, which rotates
// the 8 bits of an element's (row, column) left by one; four rounds swap
// row and column.
inline void transpose16(u8x16* v) {
  for (int round = 0; round < 4; ++round) {
    u8x16 t[16];
    for (int r = 0; r < 8; ++r) {
      t[2 * r] = zip_lo(v[r], v[r + 8]);
      t[2 * r + 1] = zip_hi(v[r], v[r + 8]);
    }
    for (int r = 0; r < 16; ++r) v[r] = t[r];
  }
}

}  // namespace bio_detail

// Whether bio_opens8 can give the gap opens of positions [lo, lo + n): the
// run table is there, every go(c) it can index fits a byte, and the vectors
// it loads lie inside the reference and the table.
inline bool bio_opens8_fit(const IndexView& ix, const BioModel& m, int64_t lo, int n) {
  return m.ctx_table() && m.go_lut() && lo >= 0 && static_cast<uint64_t>(lo) + static_cast<uint64_t>(n) + 64 <= ix.ref_len;
}

// The gap opens of positions lo + x, x < n, from the run table, 16 at a time:
// d[x] = go(run(lo + x)) and in[x] = go(max(run(lo + x), run(lo + x + 1))),
// the costs bio_window's loop computes. The table's nibbles are unpacked 32
// at a time into `runs`, and go is one 16-entry table lookup (a run read from
// the table is at most kBioCtxPacked). Needs bio_opens8_fit(ix, m, lo, n); d
// and in hold n rounded up to 16 (the bytes past n are left unspecified).
inline void bio_opens8(const BioModel& m, int64_t lo, int n, std::vector<uint8_t>& runs, uint8_t* d, uint8_t* in) {
  using namespace bio_detail;
  // runs[k] = the run at 2 (lo / 2) + k, for k < odd + n + 16: the loads of
  // run(lo + x + 1) below reach that far.
  const int odd = static_cast<int>(lo & 1), nr = (odd + n + 16 + 31) / 32 * 32;
  runs.resize(static_cast<size_t>(nr));
  const uint8_t* tab = m.ctx_table() + (lo >> 1);
  const u8x16 low = splat8(0x0F);
  for (int k = 0; k < nr; k += 32) {
    const u8x16 v = load8(tab + k / 2), even = v & low, next = v >> 4;  // x even in the low half
    store8(runs.data() + k, zip_lo(even, next));
    store8(runs.data() + k + 16, zip_hi(even, next));
  }
  const u8x16 lut = load8(m.go_lut());
  const uint8_t* r = runs.data() + odd;
  for (int x = 0; x < n; x += 16) {
    const u8x16 c = load8(r + x), c1 = load8(r + x + 1);
    store8(d + x, lookup8(lut, c));
    store8(in + x, lookup8(lut, vmax8(c, c1)));
  }
}

// The window of a band from diagonal lo, W = L + B + 1 bases: sc.win = the
// reference from lo - kBioCtxReach (W + 2 kBioCtxReach bases); sc.go_d[x],
// sc.go_i[x] = the cost of opening a deletion of reference base lo + x, an
// insertion after it.
inline void bio_window(const IndexView& ix, const BioModel& m, int64_t lo, int W, BioScratch& sc) {
  sc.win.resize(W + 2 * kBioCtxReach);
  sc.go_d.resize(W);
  sc.go_i.resize(W);
  bio_fetch(ix, lo - kBioCtxReach, W + 2 * kBioCtxReach, sc.win.data());
  if (bio_opens8_fit(ix, m, lo, W)) {
    const size_t Wp = static_cast<size_t>(W + 15) / 16 * 16;
    sc.od.resize(Wp);
    sc.oi.resize(Wp);
    bio_opens8(m, lo, W, sc.runs, sc.od.data(), sc.oi.data());
    for (int x = 0; x < W; ++x) {
      sc.go_d[x] = sc.od[x];
      sc.go_i[x] = sc.oi[x];
    }
    sc.go_i[W - 1] = sc.go_d[W - 1];  // the loop below: no base after the window
    return;
  }
  // Gap opens per window base, from its tandem run (and the next base's).
  auto opens = [&](auto run_at) {
    uint8_t c = run_at(0);
    for (int x = 0; x < W; ++x) {
      const uint8_t c1 = x + 1 < W ? run_at(x + 1) : uint8_t(0);
      sc.go_d[x] = m.go(c);
      sc.go_i[x] = m.go(std::max(c, c1));
      c = c1;
    }
  };
  if (const uint8_t* tab = m.ctx_table()) {  // clamped where go no longer falls: the same costs
    opens([&](int x) { return bio_ctx_at(tab, ix.ref_len, lo + x); });
  } else {
    sc.ctx.resize(W);
    bio_str_context_buf(sc.win.data(), W, sc.ctx.data());
    opens([&](int x) { return sc.ctx[x]; });
  }
}

// The read (strand st) on diagonals [lo, lo + B): sc.endv[b] = the least
// penalty of an end-to-end alignment within the band whose last base lies
// on diagonal lo + b (kBioInf if none). Cell (i, b) aligns read base i - 1 to
// reference i - 1 + lo + b; an alignment starts and ends with an aligned
// base. With `trace`, sc.tb[i * B + b] holds: bits 0-1 the M cell's source
// (0 M, 1 I, 2 D; ties to the first, which places gaps leftmost), bit 2 the
// I cell extends, bit 3 the D cell extends (ties extend).
// A cell whose penalty reaches `thr` is dropped (kept at kBioInf). Costs are
// non-negative, so every cell on the best path to a cell below thr is below
// thr too: each row is computed only where the row before it, or a deletion
// run in it, is still below thr, and stops once no cell is. Values below
// thr are exact; an end at or above it reads kBioInf.
inline void bio_band(const IndexView& ix, const BioModel& m, const BioRead& br, int st, int64_t lo, int B,
                     BioScratch& sc, bool trace, int thr = kBioInf) {
  const int L = br.L;
  bio_window(ix, m, lo, L + B + 1, sc);
  // Every row starts at kBioInf: a row's cells are reset where they were set
  // two rows before, not across the whole band.
  sc.rows.assign(6 * (B + 1), kBioInf);
  int* Mp = sc.rows.data();
  int* Ip = Mp + (B + 1);
  int* Dp = Ip + (B + 1);
  int* Mc = Dp + (B + 1);
  int* Ic = Mc + (B + 1);
  int* Dc = Ic + (B + 1);
  for (int b = 0; b < B; ++b) Mp[b] = 0;  // free start
  if (trace) sc.tb.assign(static_cast<size_t>(L + 1) * B, 0);
  const uint8_t* win = sc.win.data() + kBioCtxReach;
  const int* go_d = sc.go_d.data();
  const int* go_i = sc.go_i.data();
  const int ge = m.ge();
  int p0 = 0, p1 = B;  // the live cells of the row before: [p0, p1)
  int set0 = 0, set1 = B;  // the cells the row before set (the rest are kBioInf)
  for (int i = 1; i <= L; ++i) {
    const int16_t* cr = br.cost[st][i - 1];
    const int em = br.emit[st][i - 1];
    uint8_t* tb = trace ? sc.tb.data() + static_cast<size_t>(i) * B : nullptr;
    if (p0 >= p1) {  // no cell below thr: nor is any end
      std::swap(Mp, Mc);
      break;
    }
    // M reads cell b of the row before, I cell b + 1.
    const int b0 = std::max(0, p0 - 1), b1 = p1;
    for (int b = b0; b < b1; ++b) {
      const int x = i - 1 + b;
      int v = Mp[b], src = 0;
      if (Ip[b] < v) { v = Ip[b]; src = 1; }
      if (Dp[b] < v) { v = Dp[b]; src = 2; }
      Mc[b] = v >= kBioInf ? kBioInf : v + cr[win[x]];
      // Read base i - 1 inserted after reference lo + x (cell (i - 1, b + 1)).
      int iv = kBioInf, iext = 0;
      if (b + 1 < B) {
        const int open = i >= 2 && Mp[b + 1] < kBioInf ? Mp[b + 1] + go_i[x] : kBioInf;
        const int ext = Ip[b + 1] < kBioInf ? Ip[b + 1] + ge : kBioInf;
        if (ext <= open && ext < kBioInf) { iv = ext + em; iext = 1; }
        else if (open < kBioInf) iv = open + em;
      }
      Ic[b] = iv >= thr ? kBioInf : iv;
      if (Mc[b] >= thr) Mc[b] = kBioInf;
      if (tb) tb[b] = static_cast<uint8_t>(src | (iext << 2));
    }
    // Reference lo + i - 1 + b deleted after read base i - 1 (cell (i, b - 1)),
    // past the M cells only while a deletion run is still below thr.
    int e = b0 + 1;  // cells from e on stay kBioInf
    for (; e < B && (e <= b1 || Dc[e - 1] < kBioInf); ++e) {
      const int b = e;
      const int x = i - 1 + b;
      const int open = Mc[b - 1] < kBioInf ? Mc[b - 1] + go_d[x] : kBioInf;
      const int ext = Dc[b - 1] < kBioInf ? Dc[b - 1] + ge : kBioInf;
      if (ext <= open && ext < kBioInf) {
        Dc[b] = ext;
        if (tb) tb[b] |= 8;
      } else {
        Dc[b] = open;
      }
      if (Dc[b] >= thr) Dc[b] = kBioInf;
    }
    p0 = B;
    p1 = 0;
    const int end = std::max(e, b1);
    for (int b = b0; b < end; ++b)
      if (Mc[b] < kBioInf || Ic[b] < kBioInf || Dc[b] < kBioInf) {
        p0 = std::min(p0, b);
        p1 = b + 1;
      }
    std::swap(Mp, Mc);
    std::swap(Ip, Ic);
    std::swap(Dp, Dc);
    // The next row is written over the row before this one: reset what it set.
    std::fill(Mc + set0, Mc + set1, kBioInf);
    std::fill(Ic + set0, Ic + set1, kBioInf);
    std::fill(Dc + set0, Dc + set1, kBioInf);
    set0 = b0;
    set1 = end;
  }
  sc.endv.assign(Mp, Mp + B);
}

// Whether bio_bands can align bands pruned at thr: a cell it keeps is at
// most thr, and that plus any one cost must fit a byte.
inline bool bio_lanes_fit(const BioModel& m, int thr) { return thr >= 0 && thr + m.max_add() <= 255; }

// bio_bands' inputs for bands l < n (the other lanes lie past every band).
// Per window position x < L - 1 + Bm (Bm the widest band), 16 bytes, lane l
// in byte l: sc.lcode, the reference code there (+ 5 on strand 1, for the
// read's lane table; 4 past band l's window), and sc.lgo_d, sc.lgo_i, its gap
// opens (bio_window's; 0 past the window). Per diagonal b < Bm, sc.lmask: 0
// inside band l, thr past it. Each band's window is written as a row of
// sc.lstage, 16 positions at a time where bio_opens8 applies (else from
// bio_window), and the rows are transposed into the lanes 16 x 16 bytes at a
// time.
inline void bio_lane_setup(const IndexView& ix, const BioModel& m, int L, int n, const uint8_t* st,
                           const int64_t* lo, const int* B, int thr, BioScratch& sc) {
  using namespace bio_detail;
  constexpr int K = 16;
  int Bm = 0;
  for (int l = 0; l < n; ++l) Bm = std::max(Bm, B[l]);
  const int X = L - 1 + Bm, Xp = (X + 15) / 16 * 16;
  const size_t plane = static_cast<size_t>(K) * Xp;
  sc.lstage.resize(3 * plane);
  for (int l = 0; l < K; ++l) {
    uint8_t* rc = sc.lstage.data() + static_cast<size_t>(l) * Xp;
    uint8_t* rd = rc + plane;
    uint8_t* ri = rd + plane;
    const int e = l < n ? L - 1 + B[l] : 0;  // the window positions band l reads
    if (e > 0 && bio_opens8_fit(ix, m, lo[l], e)) {
      bio_opens8(m, lo[l], e, sc.runs, rd, ri);
      const uint8_t* ref = ix.ref + lo[l];
      const u8x16 s5 = splat8(5 * st[l]);
      for (int x = 0; x < e; x += 16) store8(rc + x, load8(ref + x) + s5);
    } else if (e > 0) {
      bio_window(ix, m, lo[l], L + B[l] + 1, sc);
      const uint8_t* w = sc.win.data() + kBioCtxReach;
      const uint8_t s5 = static_cast<uint8_t>(5 * st[l]);
      for (int x = 0; x < e; ++x) {
        rc[x] = static_cast<uint8_t>(w[x] + s5);
        rd[x] = static_cast<uint8_t>(sc.go_d[x]);
        ri[x] = static_cast<uint8_t>(sc.go_i[x]);
      }
    }
    std::memset(rc + e, 4, static_cast<size_t>(Xp - e));
    std::memset(rd + e, 0, static_cast<size_t>(Xp - e));
    std::memset(ri + e, 0, static_cast<size_t>(Xp - e));
  }
  sc.lcode.resize(plane);
  sc.lgo_d.resize(plane);
  sc.lgo_i.resize(plane);
  uint8_t* const out[3] = {sc.lcode.data(), sc.lgo_d.data(), sc.lgo_i.data()};
  for (int p = 0; p < 3; ++p) {
    const uint8_t* rows = sc.lstage.data() + static_cast<size_t>(p) * plane;
    for (int x = 0; x < Xp; x += 16) {
      u8x16 v[16];
      for (int l = 0; l < K; ++l) v[l] = load8(rows + static_cast<size_t>(l) * Xp + x);
      transpose16(v);
      for (int j = 0; j < 16; ++j) store8(out[p] + static_cast<size_t>(x + j) * K, v[j]);
    }
  }
  sc.lmask.assign(static_cast<size_t>(Bm) * K, static_cast<uint8_t>(thr));
  for (int l = 0; l < n; ++l)
    for (int b = 0; b < B[l]; ++b) sc.lmask[static_cast<size_t>(b) * K + l] = 0;
}

// bio_band without traceback on n <= 16 bands of the read at once, one 8-bit
// lane each: band l on strand st[l], diagonals [lo[l], lo[l] + B[l]). Its
// ends are sc.lend[l * sc.lstride + b], the values bio_band leaves in
// sc.endv. Needs bio_lanes_fit(m, thr).
// - A cell holds min(penalty, thr): every term of the recurrences is
//   non-decreasing in the cells it reads and adds a cost >= 0, so the capped
//   values satisfy the recurrences capped, and a value below thr is exact (a
//   cell at thr is bio_band's kBioInf).
// - A lane's cells past its band (b >= B[l]) stay at thr, as bio_band has
//   none there: row 0 starts them at thr and its D cells are held there, so
//   M (from the same cells a row up) and I (from the next ones) stay there
//   too, and nothing enters the band from them.
// - Each row is computed over the union of the lanes' live cells, as
//   bio_band computes its own; the bands stop together once no lane has one.
inline void bio_bands(const IndexView& ix, const BioModel& m, const BioRead& br, int n, const uint8_t* st,
                      const int64_t* lo, const int* B, int thr, BioScratch& sc) {
  using namespace bio_detail;
  constexpr int K = 16;
  const int L = br.L;
  int Bm = 0;
  for (int l = 0; l < n; ++l) Bm = std::max(Bm, B[l]);
  const uint8_t cap8 = static_cast<uint8_t>(thr);
  bio_lane_setup(ix, m, L, n, st, lo, B, thr, sc);
  uint8_t emi[K];
  for (int l = 0; l < K; ++l) emi[l] = static_cast<uint8_t>(l < n ? 10 + st[l] : 10);
  const size_t R = static_cast<size_t>(Bm + 1) * K;  // a row: Bm + 1 cells (the last stays thr)
  sc.lrows.assign(6 * R, cap8);
  uint8_t* Mp = sc.lrows.data();
  uint8_t* Ip = Mp + R;
  uint8_t* Dp = Ip + R;
  uint8_t* Mc = Dp + R;
  uint8_t* Ic = Mc + R;
  uint8_t* Dc = Ic + R;
  const uint8_t* mask = sc.lmask.data();
  std::memcpy(Mp, mask, static_cast<size_t>(Bm) * K);  // free start: 0 in every band
  const u8x16 cap = splat8(thr), ge = splat8(m.ge()), emix = load8(emi);
  int p0 = 0, p1 = Bm, set0 = 0, set1 = Bm;
  bool live = true;
  for (int i = 1; i <= L; ++i) {
    if (p0 >= p1) {  // no lane has a cell below thr: nor any end
      live = false;
      break;
    }
    const u8x16 T = load8(br.lane[i - 1]), em = lookup8(T, emix);
    const size_t r0 = static_cast<size_t>(i - 1) * K;
    const uint8_t* code = sc.lcode.data() + r0;
    const uint8_t* go_d = sc.lgo_d.data() + r0;
    const uint8_t* go_i = sc.lgo_i.data() + r0;
    const int b0 = std::max(0, p0 - 1), b1 = p1;
    u8x16 mprev = cap, dprev = cap;  // this row's M and D at b - 1
    u8x16 mp = load8(Mp + static_cast<size_t>(b0) * K), ip = load8(Ip + static_cast<size_t>(b0) * K);
    int b = b0;
    for (; b < b1; ++b) {
      const size_t o = static_cast<size_t>(b) * K;
      const u8x16 mpn = load8(Mp + o + K), ipn = load8(Ip + o + K), msk = load8(mask + o);
      const u8x16 v = vmin8(vmin8(mp, ip), load8(Dp + o));
      const u8x16 mc = vmin8(v + lookup8(T, load8(code + o)), cap);
      // Insertions open after an aligned base: none in row 1.
      const u8x16 ic = i >= 2 ? vmin8(vmin8(vmin8(mpn + load8(go_i + o), ipn + ge), cap) + em, cap) : cap;
      const u8x16 dc = vmax8(vmin8(vmin8(mprev + load8(go_d + o), dprev + ge), cap), msk);
      store8(Mc + o, mc);
      store8(Ic + o, ic);
      store8(Dc + o, dc);
      mprev = mc;
      dprev = dc;
      mp = mpn;
      ip = ipn;
    }
    // Past the M cells, while a deletion run is below thr in some lane.
    for (; b < Bm && (b == b1 || any_below(dprev, cap)); ++b) {
      const size_t o = static_cast<size_t>(b) * K;
      const u8x16 dc = vmax8(vmin8(vmin8(mprev + load8(go_d + o), dprev + ge), cap), load8(mask + o));
      store8(Dc + o, dc);
      mprev = cap;
      dprev = dc;
    }
    const int end = b;
    auto alive = [&](int c) {
      const size_t o = static_cast<size_t>(c) * K;
      return any_below(vmin8(vmin8(load8(Mc + o), load8(Ic + o)), load8(Dc + o)), cap);
    };
    p0 = Bm;
    p1 = 0;
    for (int c = b0; c < end; ++c)
      if (alive(c)) {
        p0 = c;
        break;
      }
    for (int c = end; c-- > p0;)
      if (alive(c)) {
        p1 = c + 1;
        break;
      }
    std::swap(Mp, Mc);
    std::swap(Ip, Ic);
    std::swap(Dp, Dc);
    // The next row is written over the row before this one: reset what it set.
    const size_t s0 = static_cast<size_t>(set0) * K, s1 = static_cast<size_t>(set1) * K;
    std::fill(Mc + s0, Mc + s1, cap8);
    std::fill(Ic + s0, Ic + s1, cap8);
    std::fill(Dc + s0, Dc + s1, cap8);
    set0 = b0;
    set1 = end;
  }
  sc.lstride = Bm;
  sc.lend.resize(static_cast<size_t>(n) * Bm);
  for (int l = 0; l < n; ++l) {
    int* e = sc.lend.data() + static_cast<size_t>(l) * Bm;
    for (int b = 0; b < B[l]; ++b) {
      const int v = live ? Mp[static_cast<size_t>(b) * K + l] : thr;
      e[b] = v < thr ? v : kBioInf;
    }
  }
}

// Traceback of the band bio_band(.., trace) last filled, from the alignment
// ending on diagonal lo + b_end: CIGAR, ref_pos, strand, NM and BWA-MEM's
// score (AS) into r. False if the CIGAR would not fit.
inline bool bio_trace(const BioRead& br, int st, int64_t lo, int B, int b_end, const BioScratch& sc, Result& r) {
  const int L = br.L;
  uint32_t rev[MAX_CIGAR];
  int nrev = 0;
  uint32_t run_op = 99, run_len = 0;
  auto emit = [&](uint32_t op) -> bool {
    if (op == run_op) { ++run_len; return true; }
    if (run_len) {
      if (nrev == MAX_CIGAR) return false;
      rev[nrev++] = (run_len << 4) | run_op;
    }
    run_op = op;
    run_len = 1;
    return true;
  };
  int i = L, b = b_end, state = 0, first_b = b_end, nm = 0, score = 0;
  for (;;) {
    const uint8_t bits = sc.tb[static_cast<size_t>(i) * B + b];
    if (state == 0) {
      if (!emit(kOpM)) return false;
      if (sub_cost(br.seq[st][i - 1], sc.win[kBioCtxReach + i - 1 + b])) { ++nm; score -= kMismatch; }
      else score += kMatch;
      first_b = b;
      state = bits & 3;
      if (--i == 0) break;
    } else if (state == 1) {
      if (!emit(kOpI)) return false;
      ++nm;
      state = (bits & 4) ? 1 : 0;
      --i;
      ++b;
    } else {
      if (!emit(kOpD)) return false;
      ++nm;
      state = (bits & 8) ? 2 : 0;
      --b;
    }
  }
  if (run_len) {
    if (nrev == MAX_CIGAR) return false;
    rev[nrev++] = (run_len << 4) | run_op;
  }
  r.n_cigar = 0;
  for (int x = nrev - 1; x >= 0; --x) {
    r.cigar[r.n_cigar++] = rev[x];
    const uint32_t op = rev[x] & 0xF, n = rev[x] >> 4;
    if (op == kOpI || op == kOpD) score -= kGapOpen + static_cast<int>(n) * kGapExt;
  }
  r.ref_pos = lo + first_b;
  r.strand = static_cast<uint8_t>(st);
  r.nm = static_cast<uint8_t>(std::min(nm, 255));
  r.score = static_cast<int16_t>(score);
  return true;
}

// The model's penalty of a reported end-to-end alignment (CIGAR without
// clips, at r.ref_pos on r.strand), -1 if the CIGAR does not span the read.
// The same costs bio_band charges (tests, and certa's own checks).
inline int bio_rescore(const IndexView& ix, const BioModel& m, const BioRead& br, const Result& r) {
  const int L = br.L, st = r.strand;
  int64_t j = r.ref_pos;
  int i = 0, pen = 0;
  for (int c = 0; c < r.n_cigar; ++c) {
    const uint32_t len = r.cigar[c] >> 4, op = r.cigar[c] & 0xF;
    if (op == kOpM) {
      for (uint32_t x = 0; x < len; ++x, ++i, ++j) pen += br.cost[st][i][ref_at(ix, j)];
    } else if (op == kOpI) {
      if (i == 0 || c + 1 == r.n_cigar) return -1;
      uint8_t cx[2];
      bio_str_context(ix, j - 1, j + 1, cx);
      pen += m.go(std::max(cx[0], cx[1])) + static_cast<int>(len - 1) * m.ge();
      for (uint32_t x = 0; x < len; ++x, ++i) pen += br.emit[st][i];
    } else if (op == kOpD) {
      if (i == 0 || c + 1 == r.n_cigar) return -1;
      uint8_t cx[1];
      bio_str_context(ix, j, j + 1, cx);
      pen += m.go(cx[0]) + static_cast<int>(len - 1) * m.ge();
      j += len;
    } else {
      return -1;
    }
  }
  return i == L ? pen : -1;
}

namespace bio_detail {

// 10^(-d / 10).
inline double weight(int d) { return std::pow(10.0, -d / 10.0); }

// Phred of 1 - P, where P = 1 / (1 + others): the posterior of a locus whose
// competitors weigh `others` relative to it. Capped at 60.
inline int posterior_mapq(double others) {
  if (others <= 0) return 60;
  const double q = 10 * std::log10((1 + others) / others);
  return std::min(60, static_cast<int>(q + 0.5));
}

}  // namespace bio_detail

// Re-aligns strand `l.strand` of the read in sc.br at locus l (its band and
// end diagonal) with traceback into r.
inline bool bio_align_locus(const IndexView& ix, const BioModel& m, const BioLocus& l, BioScratch& sc, Result& r) {
  bio_band(ix, m, sc.br, l.strand, l.lo, l.width, sc, true);
  return bio_trace(sc.br, l.strand, l.lo, l.width, l.b, sc, r);
}

// Maps one read with process_read(p) (its edit certificate, or a proven
// repeat, SR), then certifies its maximum-likelihood placement (Theorem B)
// with every part that fits p.budget, on each strand. On success r holds the
// ML alignment (tier BM; d1, d2 and radius stay the edit certificate's when
// it certified the read, else -1) and o the claims; returns true. Otherwise
// r stays SR if it was, or becomes uncertified (kBioHorizon / kBioWide).
inline bool bio_place(const IndexView& ix, const BioModel& m, const Params& p, const uint8_t* codes,
                      const char* qual, int L, uint64_t name_hash, Workspace& ws, BioScratch& sc, Result& r,
                      BioOut& o) {
  o = BioOut();
  process_read(ix, p, codes, L, name_hash, ws, r);
  if (r.reason == kBadLength) return false;
  const Result edit = r;
  const bool edit_cert = edit.certified == 1 && (edit.tier == kTierS0 || edit.tier == kTierS1);
  const int P = r.parts;
  // count_parts ran on both strands; take every part that fits the budget.
  int mst[2];
  uint32_t S[2] = {0, 0};
  for (int st = 0; st < 2; ++st) {
    mst[st] = parts_that_fit(ws, st, P, PMAX, p.budget);
    for (int x = 0; x < mst[st]; ++x) S[st] |= 1u << ws.order[st][x];
  }
  auto fail = [&](uint8_t reason) {
    if (r.certified == 1 && r.tier == kTierSR) return false;  // a proven repeat stays one
    r.certified = 0;
    r.reason = reason;
    return false;
  };
  if (mst[0] < 1 || mst[1] < 1) return fail(r.reason == kOk ? uint8_t(kRadiusNegative) : r.reason);
  o.tried = true;
  for (int st = 0; st < 2; ++st) {
    enumerate_parts(ix, L, P, mst[st], ws, st);
    sort_i64(ws.cand[st], ws.ncand[st]);
  }
  BioRead& br = sc.br;
  bio_read(m, codes, qual, L, br);
  const BioParams& bp = m.params();
  const int Hunw = std::min(bio_cover(m, br, 0, P, S[0]), bio_cover(m, br, 1, P, S[1]));
  int G = 0;
  if (Hunw > m.go_min()) G = std::min(bp.gap_cap, (Hunw - m.go_min() + m.ge() - 1) / m.ge());
  const int H = std::min(Hunw, m.go_min() + G * m.ge());
  o.n_s[0] = mst[0];
  o.n_s[1] = mst[1];
  o.horizon = H;
  o.reach = G;
  o.n_cand = ws.ncand[0] + ws.ncand[1];

  // Clusters: candidates joined while within 2G, so bands (+- G) are disjoint.
  sc.cl.clear();
  sc.cover.start(P);
  for (int st = 0; st < 2; ++st) {
    const int64_t* c = ws.cand[st];
    const int n = ws.ncand[st];
    for (int i = 0; i < n;) {
      int64_t dmin = diag_of(c[i]), dmax = dmin;
      uint16_t pm = static_cast<uint16_t>(1u << part_of(c[i]));
      int j = i + 1;
      for (; j < n && diag_of(c[j]) - dmax <= 2 * G; ++j) {
        dmax = diag_of(c[j]);
        pm |= static_cast<uint16_t>(1u << part_of(c[j]));
      }
      const uint32_t need = S[st] & ~static_cast<uint32_t>(pm);
      const int bound = sc.cover.get(st, need, [&] { return bio_cover(m, br, st, P, need); });
      const int64_t width = dmax - dmin + 2 * G + 1;
      sc.cl.push_back({dmin - G, static_cast<int>(std::min<int64_t>(width, INT32_MAX)), static_cast<uint8_t>(st),
                       pm, bound});
      i = j;
    }
  }
  // Cheapest bound first; the order is fixed, so the output does not depend
  // on how candidates arrived.
  sc.order.resize(sc.cl.size());
  for (size_t x = 0; x < sc.cl.size(); ++x) sc.order[x] = static_cast<int>(x);
  std::sort(sc.order.begin(), sc.order.end(), [&](int a, int b) {
    const BioCluster &x = sc.cl[a], &y = sc.cl[b];
    if (x.bound != y.bound) return x.bound < y.bound;
    if (x.strand != y.strand) return x.strand < y.strand;
    return x.lo < y.lo;
  });
  int p1 = kBioInf, yh = H, wide_bound = kBioInf;
  sc.ends.clear();
  // The band the loop reaches is aligned with the next ones it may reach
  // (bound below the threshold so far, not too wide), up to bp.lanes at once
  // (bio_bands). The loop takes their ends in order, as if each band were
  // aligned when reached; the ends of bands it does not reach are unused.
  const int nl = std::min(16, bp.lanes);
  const bool lanes = nl > 1 && bio_lanes_fit(m, H);
  sc.ln = 0;
  auto ends_of = [&](size_t q, int thr) -> const int* {
    for (int l = 0; l < sc.ln; ++l)
      if (sc.lq[l] == static_cast<int>(q)) return sc.lend.data() + static_cast<size_t>(l) * sc.lstride;
    sc.ln = 0;
    for (size_t z = q; lanes && z < sc.order.size() && sc.ln < nl; ++z) {
      const BioCluster& D = sc.cl[sc.order[z]];
      if (D.bound >= thr) break;
      if (D.width <= bp.max_band) sc.lq[sc.ln++] = static_cast<int>(z);
    }
    if (sc.ln < 2) {
      sc.ln = 0;
      const BioCluster& C = sc.cl[sc.order[q]];
      bio_band(ix, m, br, C.strand, C.lo, C.width, sc, false, H);
      return sc.endv.data();
    }
    uint8_t st[16];
    int64_t lo[16];
    int B[16];
    for (int l = 0; l < sc.ln; ++l) {
      const BioCluster& D = sc.cl[sc.order[sc.lq[l]]];
      st[l] = D.strand;
      lo[l] = D.lo;
      B[l] = D.width;
    }
    bio_bands(ix, m, br, sc.ln, st, lo, B, H, sc);
    return sc.lend.data();
  };
  for (size_t q = 0; q < sc.order.size(); ++q) {
    const int k = sc.order[q];
    const BioCluster& C = sc.cl[k];
    const int thr = std::min(H, p1 >= kBioInf ? kBioInf : p1 + bp.delta);
    if (C.bound >= thr) {  // so is every later one
      yh = std::min(yh, C.bound);
      break;
    }
    if (C.width > bp.max_band) {  // not aligned: bounded only
      wide_bound = std::min(wide_bound, C.bound);
      continue;
    }
    const int* endv = ends_of(q, thr);
    ++o.n_bands;
    o.cells += static_cast<int64_t>(L) * C.width;
    for (int b = 0; b < C.width; ++b) {
      const int v = endv[b];
      if (v >= H) continue;
      sc.ends.push_back({v, C.strand, C.lo + b + L - 1, k, b});
      p1 = std::min(p1, v);
    }
  }
  yh = std::min(yh, wide_bound);
  o.yh = yh;
  o.n_clusters = static_cast<int>(sc.cl.size());
  o.wide = wide_bound <= yh && wide_bound < H;

  // Loci: end positions below yh, cheapest first, each at least L from every
  // cheaper one on its strand.
  std::sort(sc.ends.begin(), sc.ends.end(), [](const BioEnd& a, const BioEnd& b) {
    if (a.pen != b.pen) return a.pen < b.pen;
    if (a.strand != b.strand) return a.strand < b.strand;
    return a.end < b.end;
  });
  sc.loci.clear();
  sc.taken[0].clear();
  sc.taken[1].clear();
  for (const BioEnd& e : sc.ends) {
    if (e.pen >= yh) break;
    std::vector<int64_t>& t = sc.taken[e.strand];
    if (bio_near(t, e.end, L)) continue;
    t.insert(std::lower_bound(t.begin(), t.end(), e.end), e.end);
    const BioCluster& C = sc.cl[e.cl];
    sc.loci.push_back({e.pen, e.strand, e.end, C.lo, C.width, e.b});
  }
  if (sc.loci.empty()) {
    o.empty = true;
    return fail(o.wide ? kBioWide : kBioHorizon);
  }

  const int best = sc.loci[0].pen;
  int n_best = 0;
  while (n_best < static_cast<int>(sc.loci.size()) && sc.loci[n_best].pen == best) ++n_best;
  const int pick = static_cast<int>(name_hash % static_cast<uint64_t>(n_best));
  double others = bio_detail::weight(yh - best);  // one competitor at the horizon
  for (size_t x = 0; x < sc.loci.size(); ++x)
    if (static_cast<int>(x) != pick) others += bio_detail::weight(sc.loci[x].pen - best);
  Result a = r;
  if (!bio_align_locus(ix, m, sc.loci[pick], sc, a)) return fail(kBioHorizon);

  o.ok = true;
  o.p1 = best;
  o.n_best = n_best;
  o.p2 = n_best > 1 ? best : (sc.loci.size() > 1 ? sc.loci[1].pen : -1);
  o.mapq = bio_detail::posterior_mapq(others);
  const size_t keep = std::max<size_t>(static_cast<size_t>(bp.max_loci), static_cast<size_t>(pick) + 1);
  o.loci.assign(sc.loci.begin(), sc.loci.begin() + std::min(keep, sc.loci.size()));
  o.list_bound = o.loci.size() < sc.loci.size() ? sc.loci[o.loci.size()].pen : yh;
  o.chosen = pick;
  o.moved = edit_cert && edit.n_best == 1 && n_best == 1 &&
            (edit.strand != a.strand || std::abs(edit.ref_pos - a.ref_pos) >= L);
  // Positional ambiguity: the end positions as cheap as a locus, within L of
  // it on its strand (a tandem repeat shifted by its period), all exact.
  // sc.ends is sorted by (penalty, strand, end): they are one run of it.
  for (BioLocus& l : o.loci) {
    auto it = std::lower_bound(sc.ends.begin(), sc.ends.end(), l, [&](const BioEnd& e, const BioLocus& k) {
      if (e.pen != k.pen) return e.pen < k.pen;
      if (e.strand != k.strand) return e.strand < k.strand;
      return e.end < k.end - L + 1;
    });
    int64_t a0 = l.end, a1 = l.end;
    for (; it != sc.ends.end() && it->pen == l.pen && it->strand == l.strand && it->end < l.end + L; ++it) {
      a0 = std::min(a0, it->end);
      a1 = std::max(a1, it->end);
    }
    l.amb = static_cast<int>(a1 - a0);
  }

  r = a;
  r.certified = 1;
  r.tier = kTierBM;
  r.reason = kOk;
  r.d1 = edit_cert ? edit.d1 : -1;
  r.d2 = edit_cert ? edit.d2 : -1;
  r.n_best = static_cast<uint16_t>(n_best);
  return true;
}

// The pair objective's insert term for loci a (mate 1, length L1) and b
// (mate 2): a proper pair (opposite strands, one contig, the reverse mate
// starting 0..max_dist after the forward one; starts taken as end - L + 1)
// costs the Gaussian fragment-length penalty, capped at `improper`; any
// other pair costs `improper`.
template <class ContigOf>
inline int bio_insert(const BioParams& bp, const BioLocus& a, int L1, const BioLocus& b, int L2, int max_dist,
                      ContigOf contig_of) {
  if (a.strand == b.strand) return bp.improper;
  const BioLocus& f = a.strand == 0 ? a : b;
  const BioLocus& rv = a.strand == 0 ? b : a;
  const int64_t fs = f.end - (a.strand == 0 ? L1 : L2) + 1, rs = rv.end - (a.strand == 0 ? L2 : L1) + 1;
  if (rs < fs || rs - fs > max_dist) return bp.improper;
  const int cf = contig_of(fs);
  if (cf < 0 || cf != contig_of(rv.end)) return bp.improper;
  const double z = (static_cast<double>(rv.end - fs + 1) - bp.insert_mean) / bp.insert_sd;
  return std::min(bp.improper, static_cast<int>(std::lround(10 / (2 * std::log(10.0)) * z * z)));
}

// The pair's maximum-likelihood placement from both mates' certified loci:
// the least V_a + V_b + insert over listed loci. A pair using an unlisted
// locus of mate 1 costs at least list_bound1 + p1 of mate 2 (and the other
// way round), so the best listed pair is the pair optimum when it is below
// both. x, y: the chosen loci (indices into o1.loci, o2.loci); mapq: each
// mate's posterior over pairs, with one unseen locus at each list bound.
struct BioPair {
  bool ok = false;
  int x = -1, y = -1, pen = -1, ins = -1, n_best = 0;
  int mapq[2] = {0, 0};
};

template <class ContigOf>
inline BioPair bio_pair(const BioParams& bp, const BioOut& o1, int L1, const BioOut& o2, int L2, int max_dist,
                        ContigOf contig_of, uint64_t name_hash) {
  BioPair out;
  const int n1 = static_cast<int>(o1.loci.size()), n2 = static_cast<int>(o2.loci.size());
  if (!o1.ok || !o2.ok || n1 == 0 || n2 == 0) return out;
  std::vector<int> ins(static_cast<size_t>(n1) * n2);
  int best = kBioInf, n_best = 0;
  for (int a = 0; a < n1; ++a)
    for (int b = 0; b < n2; ++b) {
      const int t = bio_insert(bp, o1.loci[a], L1, o2.loci[b], L2, max_dist, contig_of);
      ins[static_cast<size_t>(a) * n2 + b] = t;
      const int v = o1.loci[a].pen + o2.loci[b].pen + t;
      if (v < best) { best = v; n_best = 1; }
      else if (v == best) ++n_best;
    }
  const int horizon = std::min(o1.list_bound + o2.p1, o2.list_bound + o1.p1);
  if (best >= horizon) return out;
  int pick = static_cast<int>(name_hash % static_cast<uint64_t>(n_best));
  for (int a = 0; a < n1 && out.x < 0; ++a)
    for (int b = 0; b < n2; ++b)
      if (o1.loci[a].pen + o2.loci[b].pen + ins[static_cast<size_t>(a) * n2 + b] == best && pick-- == 0) {
        out.x = a;
        out.y = b;
        break;
      }
  // Posteriors, weights relative to the best pair.
  using bio_detail::weight;
  double total = weight(o1.list_bound + o2.list_bound - best), num1 = 0, num2 = 0;
  for (int a = 0; a < n1; ++a) total += weight(o1.loci[a].pen + o2.list_bound - best);
  for (int b = 0; b < n2; ++b) total += weight(o1.list_bound + o2.loci[b].pen - best);
  for (int a = 0; a < n1; ++a)
    for (int b = 0; b < n2; ++b) {
      const double w = weight(o1.loci[a].pen + o2.loci[b].pen + ins[static_cast<size_t>(a) * n2 + b] - best);
      total += w;
      if (a == out.x) num1 += w;
      if (b == out.y) num2 += w;
    }
  out.ok = true;
  out.pen = best;
  out.ins = ins[static_cast<size_t>(out.x) * n2 + out.y];
  out.n_best = n_best;
  out.mapq[0] = bio_detail::posterior_mapq((total - num1) / num1);
  out.mapq[1] = bio_detail::posterior_mapq((total - num2) / num2);
  return out;
}

// Pair rescue (THEORY.md §12): mate k is certified (ok), mate f is not, and
// every alignment of f costs at least lb_f (its horizon when no end lay below
// it; 0 when nothing is known, as for a proven repeat). For each listed locus
// a of k, cheapest first, the end diagonals of f that make a proper pair
// with a, widened by the gap reach G = gap_cap on each side, are aligned
// exactly; f ending at e pairs with a at pen_a + V_f(e) + insert. Every
// other pair costs at least one of
//   lb_k + lb_f                         k at an unlisted locus (lb_k: its list bound);
//   pen_a + lb_f + improper             f outside a's proper window;
//   pen_a + max(lb_f, go_min + G ge)    f inside it, its gaps totalling more than G;
//   pen_a + lb_f                        a locus whose window was not aligned,
// so the best pair found is the optimum of the pair objective (k over its
// loci, f over all its alignments) when it is below all of them. f's end
// positions are grouped into loci per locus of k (cheapest first, L apart);
// MAPQ: each mate's posterior over those pairs, with one unseen pair at each
// bound above. At most rescue_loci windows are aligned.
struct BioRescue {
  bool ok = false;
  int x = -1;      // k's locus (index into its loci)
  BioLocus f{};    // f's locus, with its band for the traceback
  int pen = -1, ins = -1, n_best = 0;
  int mapq[2] = {0, 0};  // k, f
};

// sc.br holds f's read.
template <class ContigOf>
inline BioRescue bio_rescue(const IndexView& ix, const BioModel& m, const BioOut& ok, int Lk, int lb_f, int Lf,
                            int max_dist, ContigOf contig_of, uint64_t name_hash, BioScratch& sc) {
  BioRescue out;
  const BioParams& bp = m.params();
  const int n = static_cast<int>(ok.loci.size());
  if (!ok.ok || n == 0) return out;
  const int G = bp.gap_cap;
  const int thr = lb_f + bp.improper;  // f costs less than this in any pair below the horizon
  using Cand = BioRescueCand;
  std::vector<Cand>& cand = sc.rcand;
  std::vector<std::pair<int, int>>& bounds = sc.rbounds;  // (an unseen pair's bound, k's locus or -1)
  cand.clear();
  bounds.clear();
  bounds.push_back({ok.list_bound + lb_f, -1});
  // f's end diagonals (end - Lf + 1, its start as bio_insert takes it) in a
  // proper pair with k's locus z: after its start when it is forward, before
  // it when not; widened by G.
  const int width = max_dist + 1 + 2 * G;
  auto window = [&](int z, int64_t& lo, uint8_t& sf) {
    const BioLocus& A = ok.loci[z];
    const int64_t sa = A.end - Lk + 1;
    const int64_t d0 = A.strand == 0 ? sa : sa - max_dist;
    lo = d0 - G;
    sf = static_cast<uint8_t>(1 - A.strand);
  };
  // Once the loop reaches a locus it reaches the ones tied with it: the best
  // pair so far costs at least its pen + lb_f, and so does any pair through a
  // tied locus. Their windows are aligned together, up to bp.lanes at once
  // (bio_bands).
  const int nl = std::min(16, bp.lanes);
  const bool lanes = nl > 1 && bio_lanes_fit(m, thr);
  sc.ln = 0;
  auto ends_of = [&](int a) -> const int* {
    for (int l = 0; l < sc.ln; ++l)
      if (sc.lq[l] == a) return sc.lend.data() + static_cast<size_t>(l) * sc.lstride;
    sc.ln = 0;
    for (int z = a; lanes && z < std::min(n, bp.rescue_loci) && sc.ln < nl && ok.loci[z].pen == ok.loci[a].pen; ++z)
      sc.lq[sc.ln++] = z;
    uint8_t st[16];
    int64_t lo[16];
    int B[16];
    if (sc.ln < 2) {
      sc.ln = 0;
      window(a, lo[0], st[0]);
      bio_band(ix, m, sc.br, st[0], lo[0], width, sc, false, thr);
      return sc.endv.data();
    }
    for (int l = 0; l < sc.ln; ++l) {
      window(sc.lq[l], lo[l], st[l]);
      B[l] = width;
    }
    bio_bands(ix, m, sc.br, sc.ln, st, lo, B, thr, sc);
    return sc.lend.data();
  };
  int best = kBioInf;
  for (int a = 0; a < n; ++a) {
    const BioLocus& A = ok.loci[a];
    if (A.pen + lb_f > best || a == bp.rescue_loci) {  // so are the later ones
      for (int z = a; z < n; ++z) bounds.push_back({ok.loci[z].pen + lb_f, z});
      break;
    }
    bounds.push_back({A.pen + lb_f + bp.improper, a});
    bounds.push_back({A.pen + std::max(lb_f, m.go_min() + G * m.ge()), a});
    int64_t lo;
    uint8_t sf;
    window(a, lo, sf);
    const int* endv = ends_of(a);
    for (int b = G; b <= G + max_dist; ++b) {
      const int v = endv[b];
      if (v >= thr) continue;
      const BioLocus F{v, sf, lo + b + Lf - 1, lo, width, b};
      const int ins = bio_insert(bp, A, Lk, F, Lf, max_dist, contig_of);
      cand.push_back({A.pen + v + ins, a, v, ins, F});
      best = std::min(best, A.pen + v + ins);
    }
  }
  int horizon = kBioInf;
  for (const auto& h : bounds) horizon = std::min(horizon, h.first);
  if (cand.empty() || best >= horizon) return out;
  // f's loci per locus of k: cheapest pair first, each L from every cheaper one.
  std::sort(cand.begin(), cand.end(), [](const Cand& x, const Cand& y) {
    if (x.cost != y.cost) return x.cost < y.cost;
    if (x.a != y.a) return x.a < y.a;
    return x.f.end < y.f.end;
  });
  std::vector<Cand>& loci = sc.rloci;
  std::vector<std::pair<int, int64_t>>& taken = sc.rtaken;  // (k's locus, f's end), sorted
  loci.clear();
  taken.clear();
  for (const Cand& c : cand) {
    const std::pair<int, int64_t> key{c.a, c.f.end};
    const auto it = std::lower_bound(taken.begin(), taken.end(), std::make_pair(c.a, c.f.end - Lf + 1));
    if (it != taken.end() && it->first == c.a && it->second < c.f.end + Lf) continue;
    taken.insert(std::lower_bound(taken.begin(), taken.end(), key), key);
    loci.push_back(c);
  }
  int n_best = 0;
  while (n_best < static_cast<int>(loci.size()) && loci[n_best].cost == best) ++n_best;
  const Cand& pick = loci[static_cast<size_t>(name_hash % static_cast<uint64_t>(n_best))];
  using bio_detail::weight;
  double total = 0, num_k = 0, num_f = 0;
  for (const Cand& l : loci) {
    const double w = weight(l.cost - best);
    total += w;
    if (l.a == pick.a) num_k += w;
    if (l.f.strand == pick.f.strand && std::abs(l.f.end - pick.f.end) < Lf) num_f += w;
  }
  for (const auto& h : bounds) {
    const double w = weight(h.first - best);
    total += w;
    if (h.second == pick.a) num_k += w;
  }
  // Positional ambiguity of f's locus: its window's ends as cheap, within Lf.
  out.f = pick.f;
  int64_t a0 = pick.f.end, a1 = pick.f.end;
  for (const Cand& c : cand)
    if (c.a == pick.a && c.v == pick.v && std::abs(c.f.end - pick.f.end) < Lf) {
      a0 = std::min(a0, c.f.end);
      a1 = std::max(a1, c.f.end);
    }
  out.f.amb = static_cast<int>(a1 - a0);
  out.ok = true;
  out.x = pick.a;
  out.pen = best;
  out.ins = pick.ins;
  out.n_best = n_best;
  out.mapq[0] = bio_detail::posterior_mapq((total - num_k) / num_k);
  out.mapq[1] = bio_detail::posterior_mapq((total - num_f) / num_f);
  return out;
}

}  // namespace certa
