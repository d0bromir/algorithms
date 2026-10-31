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
#include <unordered_map>
#include <vector>

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
  }
  const BioParams& params() const { return p_; }
  // Read base b (0-3) at quality q on reference base a (0-3), a != b.
  int mis(int q, int a, int b) const { return mis_[q][a][b]; }
  int mmin(int q) const { return mmin_[q]; }  // the cheapest mismatch at quality q
  int emit(int q) const { return emit_[q]; }  // an inserted base at quality q
  int go(int ctx) const { return go_[ctx]; }  // opening a gap that touches a tandem run of ctx bases
  int ge() const { return ge_; }
  int go_min() const { return go_min_; }

 private:
  BioParams p_;
  int16_t mis_[64][4][4];
  int16_t mmin_[64], emit_[64];
  int go_[kBioCtxCap + 1];
  int ge_, go_min_;
};

// The read's costs on both strands (strand 1: the reverse complement, its
// qualities reversed). cost[st][i][a]: read base i on reference code a
// (4 = N); emit: as an inserted base; mmin: its cheapest mismatch (0 for N).
struct BioRead {
  int L = 0;
  uint8_t seq[2][LMAX];
  int16_t cost[2][LMAX][5];
  int16_t emit[2][LMAX];
  int16_t mmin[2][LMAX];
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
}

// ctx[x - a] for reference positions x in [a, b): the length (capped at
// kBioCtxCap) of the longest tandem run of period 1-4 containing x, at least
// 1. A run of period p is a stretch where every base equals the one p before
// it, plus those p; N ends it. Depends only on the reference within
// kBioCtxCap + 4 of x.
inline void bio_str_context(const IndexView& ix, int64_t a, int64_t b, uint8_t* ctx) {
  for (int64_t x = a; x < b; ++x) ctx[x - a] = 1;
  const int64_t lo = a - kBioCtxCap - 4, hi = b + kBioCtxCap + 4;
  for (int p = 1; p <= 4; ++p) {
    int64_t s = -1;  // start of the current stretch of matches
    for (int64_t x = lo + p; x <= hi; ++x) {
      const uint8_t c = ref_at(ix, x);
      const bool match = x < hi && c < 4 && c == ref_at(ix, x - p);
      if (match) {
        if (s < 0) s = x;
        continue;
      }
      if (s >= 0) {  // stretch [s, x) of matches: run [s - p, x)
        const int64_t r0 = s - p, len = std::min<int64_t>(kBioCtxCap, x - r0);
        for (int64_t y = std::max(r0, a); y < std::min(x, b); ++y)
          if (len > ctx[y - a]) ctx[y - a] = static_cast<uint8_t>(len);
        s = -1;
      }
    }
  }
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
};

struct BioScratch {
  BioRead br;
  std::vector<int> rows;               // 6 DP rows
  std::vector<uint8_t> win, ctx, tb;   // the band's reference, tandem runs, traceback bits
  std::vector<int> endv;               // the last band's end row
  std::vector<BioCluster> cl;
  std::vector<int> order;
  std::vector<BioEnd> ends;
  std::vector<BioLocus> loci;
  std::unordered_map<uint32_t, int> cover[2];
};

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
  const int W = L + B + 1;
  sc.win.resize(W);
  sc.ctx.resize(W);
  for (int x = 0; x < W; ++x) sc.win[x] = ref_at(ix, lo + x);
  bio_str_context(ix, lo, lo + W, sc.ctx.data());
  sc.rows.assign(6 * (B + 1), kBioInf);
  int* Mp = sc.rows.data();
  int* Ip = Mp + (B + 1);
  int* Dp = Ip + (B + 1);
  int* Mc = Dp + (B + 1);
  int* Ic = Mc + (B + 1);
  int* Dc = Ic + (B + 1);
  for (int b = 0; b < B; ++b) Mp[b] = 0;  // free start
  if (trace) sc.tb.assign(static_cast<size_t>(L + 1) * B, 0);
  const uint8_t* win = sc.win.data();
  const uint8_t* ctx = sc.ctx.data();
  const int ge = m.ge();
  int p0 = 0, p1 = B;  // the live cells of the row before: [p0, p1)
  for (int i = 1; i <= L; ++i) {
    const int16_t* cr = br.cost[st][i - 1];
    const int em = br.emit[st][i - 1];
    uint8_t* tb = trace ? sc.tb.data() + static_cast<size_t>(i) * B : nullptr;
    std::fill(Mc, Mc + B + 1, kBioInf);
    std::fill(Ic, Ic + B + 1, kBioInf);
    std::fill(Dc, Dc + B + 1, kBioInf);
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
        const int open = i >= 2 && Mp[b + 1] < kBioInf ? Mp[b + 1] + m.go(std::max(ctx[x], ctx[x + 1])) : kBioInf;
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
    for (int b = b0 + 1; b < B && (b <= b1 || Dc[b - 1] < kBioInf); ++b) {
      const int x = i - 1 + b;
      const int open = Mc[b - 1] < kBioInf ? Mc[b - 1] + m.go(ctx[x]) : kBioInf;
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
    for (int b = b0; b < B; ++b)
      if (Mc[b] < kBioInf || Ic[b] < kBioInf || Dc[b] < kBioInf) {
        p0 = std::min(p0, b);
        p1 = b + 1;
      }
    std::swap(Mp, Mc);
    std::swap(Ip, Ic);
    std::swap(Dp, Dc);
  }
  sc.endv.assign(Mp, Mp + B);
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
      if (sub_cost(br.seq[st][i - 1], sc.win[i - 1 + b])) { ++nm; score -= kMismatch; }
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

  // Clusters: candidates joined while within 2G, so bands (+- G) are disjoint.
  sc.cl.clear();
  sc.cover[0].clear();
  sc.cover[1].clear();
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
      auto it = sc.cover[st].find(need);
      const int bound = it != sc.cover[st].end() ? it->second
                                                 : (sc.cover[st][need] = bio_cover(m, br, st, P, need));
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
  for (int k : sc.order) {
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
    bio_band(ix, m, br, C.strand, C.lo, C.width, sc, false, H);
    for (int b = 0; b < C.width; ++b) {
      const int v = sc.endv[b];
      if (v >= H) continue;
      sc.ends.push_back({v, C.strand, C.lo + b + L - 1, k, b});
      p1 = std::min(p1, v);
    }
  }
  yh = std::min(yh, wide_bound);
  o.yh = yh;
  o.wide = wide_bound <= yh && wide_bound < H;

  // Loci: end positions below yh, cheapest first, each at least L from every
  // cheaper one on its strand.
  std::sort(sc.ends.begin(), sc.ends.end(), [](const BioEnd& a, const BioEnd& b) {
    if (a.pen != b.pen) return a.pen < b.pen;
    if (a.strand != b.strand) return a.strand < b.strand;
    return a.end < b.end;
  });
  sc.loci.clear();
  for (const BioEnd& e : sc.ends) {
    if (e.pen >= yh) break;
    bool near = false;
    for (const BioLocus& l : sc.loci) near |= l.strand == e.strand && std::abs(l.end - e.end) < L;
    if (near) continue;
    const BioCluster& C = sc.cl[e.cl];
    sc.loci.push_back({e.pen, e.strand, e.end, C.lo, C.width, e.b});
  }
  if (sc.loci.empty()) return fail(o.wide ? kBioWide : kBioHorizon);

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
  for (BioLocus& l : o.loci) {
    int64_t a0 = l.end, a1 = l.end;
    for (const BioEnd& e : sc.ends)
      if (e.pen == l.pen && e.strand == l.strand && std::abs(e.end - l.end) < L) {
        a0 = std::min(a0, e.end);
        a1 = std::max(a1, e.end);
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

}  // namespace certa
