// CERTA paired-end certificate (host only).
//
// Objective (tier PR): the pair score score(A1) + score(A2), each the local
// score affine_align maximizes (BWA-MEM defaults, -5 per clipped end), over
// *proper* pairs: the mates on opposite strands, and the reverse mate's first
// aligned base 0..max_dist bases after the forward mate's.
//
// Theorem P. Enumerate the parts S1, S2 of both mates as in the single-end
// certificate (every hit of every chosen part), and let U_x be the Lemma L
// bound of mate x when no part of S_x occurs exactly. Then:
//  - a pair in which no part of S1 or S2 occurs exactly scores <= U1 + U2;
//  - every other pair has an exact part in some mate x, so mate x lies in a
//    chain C of x (an "anchor") and the other mate y lies in the insert
//    window of C. A pair scoring above floor >= U1 + U2 has gap total
//    <= pad = L1 + L2 - floor - 7 in each mate, so the anchor's bands cover
//    mate x;
//  - in the window, mate y is found exhaustively by short parts (Lemma L
//    with parts of ~12 bases, scanned directly in the reference window, no
//    index): an alignment of y in the window scores at most the best found or
//    the window bound.
// Anchors are evaluated best bound first: bound(C) = ub_x(C) + N_y, where N_y
// bounds mate y near C (U_y, or the bound of a chain of y within reach). For
// each evaluated anchor, `up` bounds every pair through it and `r` is a pair
// actually realized (mate x at the best place in C, mate y the best in its
// exact window). The best realized pair R1 is the certified maximum when
// R1 > floor, every anchor of that pair has up == R1, and every other
// anchor's up, every skipped anchor's bound and the floor are <= R1.
// S2 bounds every other proper pair: one with a mate that does not overlap
// R1's (MAPQ). The bounds of the other anchors cover the pairs through
// them; a pair through an anchor of R1's pair either keeps mate x where R1
// has it and moves mate y elsewhere in the window, or moves x elsewhere in
// its chain, and both are bounded by searches that exclude R1's mates.
#pragma once
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <vector>

#include "certa/core.h"

namespace certa {

struct PairParams {
  Params se;           // enumeration (k, budget); min_support must be 1
  int max_dist = 850;  // proper pair: 0 <= start(reverse) - start(forward) <= max_dist
  int wq = 12;         // window part length (exact, scanned)
  int max_anchors = 16;
  int mapq_margin = 10;  // a pair score margin of 10 gives MAPQ 60
  int max_unrealized = 8;  // anchors tried before giving up when none realized a pair
  int max_loss = 1000;     // certify only pairs losing < max_loss points in total
};

struct PairOut {
  bool certified = false;
  int score = 0;  // certified maximum pair score
  int floor = 0;  // threshold it exceeds
  int sub = 0;    // upper bound on every other pair
};

constexpr int kPadMax = (BMAX - 1) / 2;
constexpr int kWinParts = 24;

// Lemma L bound for parts [off[k], end[k]) of a read of length L: the best
// score of an alignment in which no part in `counted` occurs exactly.
// When the counted parts are disjoint, in order, in [0, L) and each >= 5
// long, the full span attains the maximum (Lemma S; lean/Certa/LemmaS.lean,
// lemmaS), so the bound is L - 5 per counted part. Otherwise every span
// candidate is tried.
inline int span_bound(int L, int n, const int* off, const int* end, uint32_t counted) {
  static_assert(kMatch == 1 && kMatch + kMismatch == 5, "Lemma S is stated for +1 / -4");
  {
    int c = 0, prev = 0;
    bool ok = true;
    for (int k = 0; k < n && ok; ++k)
      if ((counted >> k) & 1) {
        ok = off[k] >= prev && end[k] - off[k] >= 5 && end[k] <= L;
        prev = end[k];
        ++c;
      }
    if (ok) return L * kMatch - (kMatch + kMismatch) * c;
  }
  int xs[kWinParts + 1], ys[kWinParts + 1], nx = 0, ny = 0;
  xs[nx++] = 0;
  ys[ny++] = L;
  for (int k = 0; k < n; ++k)
    if ((counted >> k) & 1) {
      xs[nx++] = off[k] + 1;
      ys[ny++] = end[k] - 1;
    }
  int best = -BIG;
  for (int a = 0; a < nx; ++a)
    for (int b = 0; b < ny; ++b) {
      const int x0 = xs[a], y0 = ys[b];
      if (y0 <= x0) continue;
      int c = 0;
      for (int k = 0; k < n; ++k) c += ((counted >> k) & 1) && off[k] >= x0 && end[k] <= y0;
      const int v = (y0 - x0) * kMatch - (kMatch + kMismatch) * c - (x0 > 0 ? kClip : 0) -
                    (y0 < L ? kClip : 0);
      if (v > best) best = v;
    }
  return best;
}

struct WinChain { int ub, beg, end; };

// Scratch for the window search (one per thread).
struct PairScratch {
  std::vector<WinChain> chains;
  std::unique_ptr<Workspace> ws{new Workspace};  // DP rows and traceback
  std::vector<std::pair<uint32_t, int>> keys;    // part key -> part
  std::vector<std::pair<int64_t, int>> hits;     // (diagonal, part)
  // window_local: the window's reference (padded with N), the DP rows
  // (H and E of the previous and current row; 32-bit for the scalar fill,
  // 16-bit for the vector one) and the traceback, a byte per cell.
  std::vector<uint8_t> lref, ltb;
  std::vector<int> lh[2], le[2];
  std::vector<int16_t> lh16[2], le16[2];
};

// Best alignment of `rd` whose first aligned base lies in [smin, smax], found
// exhaustively: every exact occurrence of its ~wq-base parts in the reference
// window is scanned, and chains of them are aligned in bands. Returns the
// score (> *thr) and the band of one best alignment, or kNoSub;
// *thr bounds every alignment not returned. Alignments scoring <= need are
// not looked for (*thr >= need). Returns kNoSub with *thr = L when the window
// is too repetitive to search.
struct Band {
  int64_t lo = 0;
  int B = 0;
};
inline int window_search(const IndexView& ix, const uint8_t* rd, int L, int64_t smin, int64_t smax,
                         int wq, int need, PairScratch& sc, Band* out, int* thr) {
  *thr = L * kMatch;
  int n = L / wq;
  if (n > kWinParts) n = kWinParts;
  if (n < 1) return kNoSub;
  int off[kWinParts], end[kWinParts];
  for (int k = 0; k < n; ++k) {
    int len;
    part_geometry(L, n, k, &off[k], &len);
    end[k] = off[k] + len;
  }
  const uint32_t all = n == 32 ? ~0u : ((1u << n) - 1);
  int t = span_bound(L, n, off, end, all);
  if (t < L * kMatch - 7 - kPadMax) t = L * kMatch - 7 - kPadMax;  // band must fit
  if (t < need) t = need;
  *thr = t;
  if (t >= L * kMatch) return kNoSub;
  // Gap total of an alignment scoring > t; 0 (no gap at all) once t >= L - 7.
  const int pad = std::max(0, L * kMatch - t - 7);
  // Part keys (first wq bases; parts with N can never occur exactly).
  sc.keys.clear();
  for (int k = 0; k < n; ++k) {
    uint32_t key = 0;
    bool ok = true;
    for (int i = 0; i < wq; ++i) {
      if (rd[off[k] + i] > 3) { ok = false; break; }
      key = (key << 2) | rd[off[k] + i];
    }
    if (ok) sc.keys.push_back({key, k});
  }
  std::sort(sc.keys.begin(), sc.keys.end());
  // The first index of each key in sc.keys, by open addressing: most window
  // positions match no part and cost one probe.
  constexpr int kH = 64;
  int16_t hidx[kH];
  uint32_t hkey[kH];
  for (int h = 0; h < kH; ++h) hidx[h] = -1;
  auto slot = [](uint32_t key) { return (int)((key * 2654435761u) >> 26); };
  for (size_t x = 0; x < sc.keys.size(); ++x) {
    if (x > 0 && sc.keys[x].first == sc.keys[x - 1].first) continue;
    int h = slot(sc.keys[x].first);
    while (hidx[h] >= 0) h = (h + 1) & (kH - 1);
    hidx[h] = (int16_t)x;
    hkey[h] = sc.keys[x].first;
  }
  // Scan [smin, smax + L + pad): an alignment starting in [smin, smax] lies there.
  const int64_t lo = smin < 0 ? 0 : smin;
  const int64_t hi = std::min<int64_t>((int64_t)ix.ref_len, smax + L + pad);
  const uint32_t mask = wq >= 16 ? ~0u : ((1u << (2 * wq)) - 1);
  sc.hits.clear();
  uint32_t key = 0;
  int valid = 0;
  for (int64_t j = lo; j < hi; ++j) {
    const uint8_t c = ix.ref[j];
    if (c > 3) { valid = 0; key = 0; continue; }
    key = ((key << 2) | c) & mask;
    if (++valid < wq) continue;
    const int64_t s = j - wq + 1;  // occurrence start
    int h = slot(key);
    while (hidx[h] >= 0 && hkey[h] != key) h = (h + 1) & (kH - 1);
    if (hidx[h] < 0) continue;
    for (auto it = sc.keys.begin() + hidx[h]; it != sc.keys.end() && it->first == key; ++it) {
      const int k = it->second;
      bool ok = s + (end[k] - off[k]) <= hi;
      for (int i = wq; ok && i < end[k] - off[k]; ++i) ok = ix.ref[s + i] == rd[off[k] + i];
      if (ok) sc.hits.push_back({s - off[k], k});
      if (sc.hits.size() > 20000) { *thr = L * kMatch; return kNoSub; }
    }
  }
  std::sort(sc.hits.begin(), sc.hits.end());
  // Chains (consecutive diagonals <= pad apart) with their Lemma L bound;
  // every chain that could beat the threshold is aligned in bands.
  sc.chains.clear();
  size_t i = 0;
  while (i < sc.hits.size()) {
    size_t j = i + 1;
    uint32_t pm = 1u << sc.hits[i].second;
    while (j < sc.hits.size() && sc.hits[j].first - sc.hits[j - 1].first <= pad) {
      pm |= 1u << sc.hits[j].second;
      ++j;
    }
    const int ub = span_bound(L, n, off, end, all & ~pm);
    if (ub > t) sc.chains.push_back({ub, (int)i, (int)j});
    i = j;
  }
  // Best bound first; a chain whose bound cannot beat the best found is
  // skipped (the maximum does not depend on the order).
  std::sort(sc.chains.begin(), sc.chains.end(), [](const WinChain& a, const WinChain& b) {
    return a.ub != b.ub ? a.ub > b.ub : a.beg < b.beg;
  });
  int best = kNoSub;
  for (const WinChain& ch : sc.chains) {
    if (ch.ub <= best) break;
    size_t a = ch.beg;
    while (a < (size_t)ch.end) {
      const int64_t d0 = sc.hits[a].first;
      size_t e = a + 1;
      while (e < (size_t)ch.end && sc.hits[e].first - d0 <= BMAX - 2 * pad - 1) ++e;
      const int64_t blo = d0 - pad;
      const int B = (int)(sc.hits[e - 1].first - d0) + 2 * pad + 1;
      // Only a band beating both the best so far and t can change the result.
      const int v = band_score(ix, rd, L, blo, B, smin, smax, std::max(best, t));
      if (v > best) {
        best = v;
        out->lo = blo;
        out->B = B;
      }
      a = e;
    }
  }
  return best > t ? best : kNoSub;
}

// The DP fill of window_local (below), row by row (read base i = 1..L),
// column by column (window base c = 0..W-1, reference smin + c, in
// sc.lref): cell (i, c) aligns read base i - 1 to window base c. The
// recurrences are affine_fill's, in reference instead of diagonal
// coordinates (a match from (i - 1, c - 1), an insertion from (i - 1, c), a
// deletion from (i, c - 1)); an alignment starts with a match in a column
// <= cmax. Writes affine_fill's traceback bits to sc.ltb[i * Wp + c] and
// returns the best end score (below -2^19: none), with its cell, the first
// in row-major order.
inline int window_fill_scalar(const uint8_t* rd, int L, int W, int Wp, int cmax, PairScratch& sc, int* pbi,
                              int* pbc) {
  const int NEG = -(1 << 20);
  for (int x = 0; x < 2; ++x) {
    sc.lh[x].assign(W, NEG);
    sc.le[x].assign(W, NEG);
  }
  int* hp = sc.lh[0].data();
  int* hc = sc.lh[1].data();
  int* ep = sc.le[0].data();
  int* ec = sc.le[1].data();
  const uint8_t* ref = sc.lref.data();
  int best = NEG, best_i = 0, best_c = 0;
  for (int i = 1; i <= L; ++i) {
    const uint8_t r = rd[i - 1];
    const int start = i == 1 ? 0 : -kClip;
    const int clip_end = i < L ? kClip : 0;
    uint8_t* tb = sc.ltb.data() + (size_t)i * Wp;
    int f = NEG, hl = NEG;  // F and H at (i, c - 1)
    for (int c = 0; c < W; ++c) {
      uint8_t bits = 0;
      int diag = c > 0 ? hp[c - 1] : NEG;
      if (c <= cmax && start >= diag) { diag = start; bits |= 4; }
      const int m = diag + (sub_cost(r, ref[c]) ? -kMismatch : kMatch);
      int e = NEG;
      if (i > 1) {  // no alignment starts with a gap
        const int open = hp[c] - kGapOpen - kGapExt, ext = ep[c] - kGapExt;
        e = open >= ext ? open : ext;
        if (ext > open) bits |= 8;
      }
      int fo = NEG;
      if (c > 0) {
        const int open = hl - kGapOpen - kGapExt, ext = f - kGapExt;
        fo = open >= ext ? open : ext;
        if (ext > open) bits |= 16;
      }
      f = fo;
      int h = m;
      if (e > h) { h = e; bits = (uint8_t)((bits & ~3) | 1); }
      if (f > h) { h = f; bits = (uint8_t)((bits & ~3) | 2); }
      hc[c] = h;
      ec[c] = e;
      hl = h;
      tb[c] = bits;
      // An alignment may end after a match/mismatch.
      const int end = (bits & 3) == 0 ? h - clip_end : NEG;
      if (end > best) { best = end; best_i = i; best_c = c; }
    }
    std::swap(hp, hc);
    std::swap(ep, ec);
  }
  *pbi = best_i;
  *pbc = best_c;
  return best;
}

#ifdef CERTA_AFFINE_VEC
// The same fill, 8 cells at a time with 16-bit lanes, and the same values
// relative to "minus infinity" (NEG -16000 here), as affine_fill_vec16: the
// deletion state from an exclusive prefix maximum of h0[c] + c along the row,
// the best end cell lane by lane, in order. Every compared value (h0[c] + c
// included) is either real, in [-9 LMAX - W - 16, LMAX + W], or NEG plus an
// offset in [-9 LMAX - W - 16, LMAX + W]: with W <= kWinVecMax the two
// ranges stay apart and no lane wraps, so every comparison, traceback bit
// and best cell is the scalar fill's (tested). Bits past column W - 1 are
// not the scalar fill's (it writes none there); nothing reads them.
constexpr int kWinVecMax = 4096;
inline int window_fill_vec16(const uint8_t* rd, int L, int W, int Wp, int cmax, PairScratch& sc, int* pbi,
                             int* pbc) {
  static_assert(kMatch == 1 && kGapExt == 1, "window_fill_vec16 assumes +1 matches and gap extension 1");
  static_assert(__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__, "window_fill_vec16 packs the low byte of each lane");
  static_assert(-16000 - 9 * LMAX - kWinVecMax - 16 > INT16_MIN &&
                    -16000 + LMAX + kWinVecMax < -9 * LMAX - kWinVecMax - 16,
                "16-bit window fill");
  typedef int16_t v8 __attribute__((vector_size(16)));
  typedef uint8_t u8x16 __attribute__((vector_size(16)));
  typedef uint64_t u64x2 __attribute__((vector_size(16)));
  typedef int16_t v8u __attribute__((vector_size(16), aligned(2), may_alias));
  typedef uint8_t u8x8u __attribute__((vector_size(8), aligned(1), may_alias));
  typedef uint64_t u64u __attribute__((aligned(1), may_alias));
  constexpr int V = 8, P = 8;  // lanes; each row's leading pad (column -1 holds NEG)
  const int NEG = -16000;
  auto vmax = [](v8 a, v8 b) { return a > b ? a : b; };
  auto splat = [](int x) {
    const int16_t s = (int16_t)x;
    return v8{s, s, s, s, s, s, s, s};
  };
  const v8 vneg = splat(NEG), lane = {0, 1, 2, 3, 4, 5, 6, 7};
  const v8 ninf = splat(INT16_MIN), vcmax = splat(cmax), vW = splat(W);
  for (int x = 0; x < 2; ++x) {
    sc.lh16[x].assign((size_t)P + Wp + V, (int16_t)NEG);
    sc.le16[x].assign((size_t)P + Wp + V, (int16_t)NEG);
  }
  int16_t* hp = sc.lh16[0].data() + P;
  int16_t* hc = sc.lh16[1].data() + P;
  int16_t* ep = sc.le16[0].data() + P;
  int16_t* ec = sc.le16[1].data() + P;
  const uint8_t* ref = sc.lref.data();
  int best = NEG, best_i = 0, best_c = 0;
  for (int i = 1; i <= L; ++i) {
    const v8 r = splat(rd[i - 1]);
    const v8 start = splat(i == 1 ? 0 : -kClip);
    const int clip_end = i < L ? kClip : 0;
    uint8_t* tb = sc.ltb.data() + (size_t)i * Wp;
    v8 carry = splat(NEG + 6);      // running max of h0[c'] + c'
    v8 fprev = vneg, hprev = vneg;  // the previous vector's F and H (lane 7 used)
    for (int c = 0; c < Wp; c += V) {
      const v8 cv = lane + splat(c);
      const v8 hd = *(const v8u*)(hp + c - 1);  // H(i - 1, c - 1)
      const v8 hu = *(const v8u*)(hp + c), eu = *(const v8u*)(ep + c);
      const v8 st = (start >= hd) & (cv <= vcmax);
      const v8 diag = st ? start : hd;
      const v8 t = __builtin_convertvector(*(const u8x8u*)(ref + c), v8);
      const v8 mis = (t != r) | (r >= splat(4));
      const v8 m = diag + (mis ? splat(-kMismatch) : splat(kMatch));
      v8 e = vneg, ebit = splat(0);
      if (i > 1) {
        const v8 open = hu - splat(kGapOpen + kGapExt), ext = eu - splat(kGapExt);
        e = open >= ext ? open : ext;
        ebit = (ext > open) & splat(8);
      }
      const v8 esel = e > m;
      const v8 h0 = esel ? e : m;
      v8 src = esel & splat(1);
      // F from the exclusive prefix maximum of h0[c] + c.
      const v8 g = h0 + cv;
      v8 inc = vmax(g, __builtin_shufflevector(g, ninf, 8, 0, 1, 2, 3, 4, 5, 6));
      inc = vmax(inc, __builtin_shufflevector(inc, ninf, 8, 9, 0, 1, 2, 3, 4, 5));
      inc = vmax(inc, __builtin_shufflevector(inc, ninf, 8, 9, 10, 11, 0, 1, 2, 3));
      inc = vmax(inc, carry);
      const v8 excl = __builtin_shufflevector(carry, inc, 7, 8, 9, 10, 11, 12, 13, 14);
      carry = __builtin_shufflevector(inc, inc, 7, 7, 7, 7, 7, 7, 7, 7);
      const v8 f = cv > splat(0) ? excl - cv - splat(kGapOpen + kGapExt - 1) : vneg;
      const v8 fsel = f > h0;
      const v8 h = fsel ? f : h0;
      src = fsel ? splat(2) : src;
      const v8 fm1 = __builtin_shufflevector(fprev, f, 7, 8, 9, 10, 11, 12, 13, 14);
      const v8 hm1 = __builtin_shufflevector(hprev, h, 7, 8, 9, 10, 11, 12, 13, 14);
      const v8 fbit = (fm1 + splat(kGapOpen + kGapExt - 1) > hm1) & (cv > splat(0)) & splat(16);
      fprev = f;
      hprev = h;
      *(v8u*)(hc + c) = h;
      *(v8u*)(ec + c) = e;
      const v8 bits = (st & splat(4)) | ebit | fbit | src;
      const u8x16 bytes = (u8x16)bits;
      const u8x16 b8 = __builtin_shufflevector(bytes, bytes, 0, 2, 4, 6, 8, 10, 12, 14, 0, 0, 0, 0, 0, 0, 0, 0);
      *(u64u*)(tb + c) = ((u64x2)b8)[0];
      // An alignment may end after a match/mismatch, in the window.
      const v8 end = (src == splat(0)) & (cv < vW) ? h - splat(clip_end) : vneg;
      const u64x2 over = (u64x2)(end > splat(best));
      if (over[0] | over[1])
        for (int l = 0; l < V; ++l)
          if (end[l] > best) { best = end[l]; best_i = i; best_c = c + l; }
    }
    std::swap(hp, hc);
    std::swap(ep, ec);
  }
  *pbi = best_i;
  *pbc = best_c;
  return best < NEG / 2 ? best - NEG - (1 << 20) : best;
}
#endif

// The best alignment of `rd` (length L) in ref[smin, hi) whose first aligned
// base lies in [smin, smax]: the objective of affine_align (Theorem L's local
// score: +1 / -4, a gap of k 6 + k, 5 per clipped read end) with no band, so
// clips at either end and gaps of any length. Exhaustive: the full DP of the
// window (window_fill_*), ties to the first best end by read base, then
// reference position. This is BWA-MEM's mate rescue (mem_matesw, a
// Smith-Waterman alignment in the mate's window) under CERTA's scoring.
// Returns true and the alignment (ref_pos, CIGAR, nm, score) if it scores
// above `need`; false otherwise, or if the CIGAR does not fit.
inline bool window_local(const IndexView& ix, const uint8_t* rd, int L, int64_t smin, int64_t smax, int64_t hi,
                         int need, PairScratch& sc, Result& out, bool scalar = false) {
  if (smin < 0) smin = 0;
  if (hi > (int64_t)ix.ref_len) hi = (int64_t)ix.ref_len;
  if (smax > hi - 1) smax = hi - 1;
  if (L < 1 || L > LMAX || smin > smax) return false;
  const int W = (int)(hi - smin);  // no alignment starts before smin: no column there
  const int Wp = (W + 7) / 8 * 8;
  sc.lref.assign((size_t)Wp + 8, 4);
  std::copy(ix.ref + smin, ix.ref + hi, sc.lref.begin());
  sc.ltb.resize((size_t)(L + 1) * Wp);
  int best_i = 0, best_c = 0, best;
#ifdef CERTA_AFFINE_VEC
  if (!scalar && W <= kWinVecMax)
    best = window_fill_vec16(rd, L, W, Wp, (int)(smax - smin), sc, &best_i, &best_c);
  else
#endif
    best = window_fill_scalar(rd, L, W, Wp, (int)(smax - smin), sc, &best_i, &best_c);
  (void)scalar;
  if (best < -(1 << 19) || best <= need) return false;
  uint32_t rev[MAX_CIGAR];
  int nrev = 0, nm = 0, i = best_i, state = 0;  // 0 H, 1 insertion, 2 deletion
  int c = best_c;
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
  for (;;) {
    if (i < 1 || c < 0) return false;  // unreachable: every path begins at a start
    const uint8_t bits = sc.ltb[(size_t)i * Wp + c];
    if (state == 0) state = bits & 3;
    if (state == 0) {
      if (!emit(kOpM)) return false;
      nm += sub_cost(rd[i - 1], sc.lref[c]);
      if (bits & 4) break;  // (i, c) is the first aligned pair
      --i;
      --c;
    } else if (state == 1) {
      if (!emit(kOpI)) return false;
      ++nm;
      state = (bits & 8) ? 1 : 0;
      --i;
    } else {
      if (!emit(kOpD)) return false;
      ++nm;
      state = (bits & 16) ? 2 : 0;
      --c;
    }
  }
  if (nrev == MAX_CIGAR) return false;
  rev[nrev++] = (run_len << 4) | run_op;
  const int lead = i - 1, trail = L - best_i;
  if (nrev + (lead > 0) + (trail > 0) > MAX_CIGAR) return false;
  out.n_cigar = 0;
  if (lead > 0) out.cigar[out.n_cigar++] = ((uint32_t)lead << 4) | kOpS;
  for (int x = nrev - 1; x >= 0; --x) out.cigar[out.n_cigar++] = rev[x];
  if (trail > 0) out.cigar[out.n_cigar++] = ((uint32_t)trail << 4) | kOpS;
  out.ref_pos = smin + c;
  out.nm = (uint8_t)std::min(nm, 255);
  out.score = (int16_t)best;
  return true;
}

// One mate: its enumerated parts, chains and bounds.
struct Mate {
  const uint8_t* codes = nullptr;
  int L = 0, P = 0, m = 0, U = 0, B = 0, nch = 0;
  int n0 = 0;  // chains [0, n0) are on the forward strand, [n0, nch) on the reverse
  Workspace* ws = nullptr;
  std::vector<int> ub;
  std::vector<int64_t> slo, shi;  // start range of each chain
};

inline bool prepare_mate(const IndexView& ix, const Params& p, Mate& M) {
  M.P = M.L <= LMAX ? part_count(M.L, ix.q, ix.s) : 0;
  if (M.P < 1) return false;
  Workspace& ws = *M.ws;
  for (int i = 0; i < M.L; ++i) {
    ws.seq[0][i] = M.codes[i];
    const uint8_t c = M.codes[M.L - 1 - i];
    ws.seq[1][i] = c < 4 ? (uint8_t)(3 - c) : (uint8_t)4;
  }
  count_parts(ix, ws.seq[0], M.L, M.P, ws, 0);
  count_parts(ix, ws.seq[1], M.L, M.P, ws, 1);
  const int m0 = parts_that_fit(ws, 0, M.P, p.k + 1, p.budget);
  const int m1 = parts_that_fit(ws, 1, M.P, p.k + 1, p.budget);
  M.m = m0 < m1 ? m0 : m1;
  enumerate_parts(ix, M.L, M.P, M.m, ws, 0);
  enumerate_parts(ix, M.L, M.P, M.m, ws, 1);
  sort_i64(ws.cand[0], ws.ncand[0]);
  sort_i64(ws.cand[1], ws.ncand[1]);
  M.U = -BIG;
  for (int st = 0; st < 2; ++st) M.U = std::max(M.U, local_bound(ws, st, M.L, M.P, M.m, 0, 0));
  return true;
}

inline void chain_mate(Mate& M, int pad) {
  Workspace& ws = *M.ws;
  M.nch = build_chains(ws, pad);
  M.ub.assign(M.nch, 0);
  M.slo.assign(M.nch, 0);
  M.shi.assign(M.nch, 0);
  ubmemo_begin(ws);
  M.B = M.U;
  M.n0 = 0;
  for (int x = 0; x < M.nch; ++x) {
    const Cluster& cl = ws.clusters[x];
    M.n0 += cl.strand == 0;
    const int u = ub_of_mask(ws, cl.strand, M.L, M.P, M.m, cl.pmask);
    M.ub[x] = u;
    M.B = std::max(M.B, (int)u);
    // A member diagonal d puts the first aligned base in [d - pad, d + pad + L - 1].
    M.slo[x] = diag_of(ws.cand[cl.strand][cl.beg]) - pad;
    M.shi[x] = diag_of(ws.cand[cl.strand][cl.end - 1]) + pad + M.L - 1;
  }
}

// Certify the best proper pair. r1/r2 receive the two mates' alignments
// (tier PR) when certified.
inline PairOut certify_pair(const IndexView& ix, const PairParams& pp, const uint8_t* c1, int L1,
                            const uint8_t* c2, int L2, Workspace& w1, Workspace& w2,
                            PairScratch& sc, Result& r1, Result& r2) {
  PairOut po;
  Mate M[2];
  M[0].codes = c1; M[0].L = L1; M[0].ws = &w1;
  M[1].codes = c2; M[1].L = L2; M[1].ws = &w2;
  if (!prepare_mate(ix, pp.se, M[0]) || !prepare_mate(ix, pp.se, M[1])) return po;
  const int Ltot = (L1 + L2) * kMatch;
  int floor = M[0].U + M[1].U;
  if (floor < Ltot - 7 - kPadMax) floor = Ltot - 7 - kPadMax;  // mate bands must fit
  if (floor < Ltot - pp.max_loss) floor = Ltot - pp.max_loss;   // narrower bands, less work
  po.floor = floor;
  if (floor >= Ltot) return po;
  const int pad = std::max(0, Ltot - floor - 7);
  chain_mate(M[0], pad);
  chain_mate(M[1], pad);
  const int64_t D = pp.max_dist;
  // Partner start range of an anchor: the forward mate starts first.
  auto partner_range = [&](int strand, int64_t a_lo, int64_t a_hi, int64_t* plo, int64_t* phi) {
    if (strand == 0) { *plo = a_lo; *phi = a_hi + D; }
    else { *plo = a_lo - D; *phi = a_hi; }
  };
  struct Anchor { int bound, mate, chain, nbeg, nend; };
  std::vector<Anchor> anchors;
  std::vector<int> near;  // per anchor: chains of the other mate within reach
  // build_chains emits each strand's chains in diagonal order, each starting
  // more than pad diagonals after the previous one ends, so within a strand
  // both slo and shi increase: the chains of mate y within reach of an
  // anchor (shi >= plo and slo <= phi) are one index range, found by binary
  // search instead of a scan of every chain (quadratic with --budget3 4096).
  for (int x = 0; x < 2; ++x) {
    const Mate& X = M[x];
    const Mate& Y = M[1 - x];
    for (int c = 0; c < X.nch; ++c) {
      const int st = X.ws->clusters[c].strand;
      int64_t plo, phi;
      partner_range(st, X.slo[c], X.shi[c], &plo, &phi);
      int ny = Y.U;  // mate y near this chain: no exact part, or one of its chains
      const size_t nb = near.size();
      const int yb = st == 1 ? 0 : Y.n0, ye = st == 1 ? Y.n0 : Y.nch;  // chains on strand 1 - st
      const int c2b = (int)(std::lower_bound(Y.shi.begin() + yb, Y.shi.begin() + ye, plo) - Y.shi.begin());
      const int c2e = (int)(std::upper_bound(Y.slo.begin() + yb, Y.slo.begin() + ye, phi) - Y.slo.begin());
      for (int c2 = c2b; c2 < c2e; ++c2) {
        ny = std::max(ny, Y.ub[c2]);
        near.push_back(c2);
      }
      const int b = X.ub[c] + ny;
      if (b > floor) anchors.push_back({b, x, c, (int)nb, (int)near.size()});
      else near.resize(nb);
    }
  }
  // A chain is "done" once every pair whose mate lies in it is accounted for
  // (in upb or S2). Pairs through an anchor whose partner lies in a done
  // chain are then covered, so the anchor's bound drops to ub + (the best
  // partner bound among chains not done, or U).
  std::vector<char> done[2] = {std::vector<char>(M[0].nch, 0), std::vector<char>(M[1].nch, 0)};
  auto dyn_bound = [&](const Anchor& an) {
    const Mate& Y = M[1 - an.mate];
    int ny = Y.U;
    for (int k = an.nbeg; k < an.nend; ++k)
      if (!done[1 - an.mate][near[k]]) ny = std::max(ny, Y.ub[near[k]]);
    return M[an.mate].ub[an.chain] + ny;
  };
  std::sort(anchors.begin(), anchors.end(), [](const Anchor& a, const Anchor& b) {
    if (a.bound != b.bound) return a.bound > b.bound;
    if (a.mate != b.mate) return a.mate < b.mate;
    return a.chain < b.chain;
  });
  int R1 = floor, upb = -BIG, S2 = floor, evaluated = 0;
  // The anchors that realized R1's pair, with what bounding their other
  // pairs needs.
  struct Win {
    int mate, chain, apad, sx, fu, thr, need, relevant;
    int64_t ulo, uhi;
  };
  std::vector<Win> wins;
  bool confirmed = false;  // R1's pair re-found by a full anchor evaluation (its up is known)
  int id_strand = -1;      // mate 1's strand in the best pair
  Result best[2];          // the best pair's two alignments
  Workspace& tw = *sc.ws;  // traceback scratch
  auto same_pair = [&](int s1, const Result& a1, const Result& a2) {
    return id_strand == s1 && std::llabs(best[0].ref_pos - a1.ref_pos) <= 2 * pad + 1 &&
           std::llabs(best[1].ref_pos - a2.ref_pos) <= 2 * pad + 1;
  };
  // One best alignment of the band (traceback); false if its CIGAR does not fit.
  auto locate = [&](const uint8_t* rd, int L, const Band& bd, int score, int64_t smin, int64_t smax,
                    Result* out) {
    return affine_align(ix, rd, L, bd.lo, bd.B, tw, *out, false, smin, smax) && out->score == score;
  };
  // Seed: a pair realized cheaply at the best anchor (narrow band, exact
  // window), so the full evaluations below start with a high R1 and use
  // narrow bands. It counts only once a full evaluation confirms it.
  if (!anchors.empty()) {
    const Anchor& an = anchors[0];
    const Mate& X = M[an.mate];
    const Mate& Y = M[1 - an.mate];
    const Cluster& cl = X.ws->clusters[an.chain];
    const int st = cl.strand;
    Band bx, by;
    Result ax, ay;
    const int sx = chain_best(ix, X.L, cl, std::min(pad, 8), *X.ws, &bx.lo, &bx.B);
    if (locate(X.ws->seq[st], X.L, bx, sx, INT64_MIN, INT64_MAX, &ax)) {
      int64_t elo, ehi;
      partner_range(st, ax.ref_pos, ax.ref_pos, &elo, &ehi);
      int thr;
      const uint8_t* ry = Y.ws->seq[1 - st];
      const int fe = window_search(ix, ry, Y.L, elo, ehi, pp.wq, floor - sx, sc, &by, &thr);
      if (fe != kNoSub && sx + fe > floor && locate(ry, Y.L, by, fe, elo, ehi, &ay)) {
        R1 = sx + fe;
        id_strand = an.mate == 0 ? st : 1 - st;
        best[0] = an.mate == 0 ? ax : ay;
        best[1] = an.mate == 0 ? ay : ax;
      }
    }
  }
  for (const Anchor& an : anchors) {
    if (an.bound <= S2) break;  // sorted: no remaining pair can matter
    // Remaining anchors cannot beat R1; their bounds only refine S2, which
    // already gives the maximum pair MAPQ when R1 - S2 >= mapq_margin.
    // Also stop when S2 is already within the margin (MAPQ is low either way)
    // and no remaining anchor can beat R1.
    if (id_strand >= 0 && an.bound <= R1 &&
        (an.bound <= R1 - pp.mapq_margin || S2 >= R1 - pp.mapq_margin)) {
      S2 = std::max(S2, an.bound);
      break;
    }
    {
      const int db = dyn_bound(an);
      if (db <= S2 || (id_strand >= 0 && db <= R1 &&
                       (db <= R1 - pp.mapq_margin || S2 >= R1 - pp.mapq_margin))) {
        S2 = std::max(S2, db);
        done[an.mate][an.chain] = 1;
        continue;
      }
    }
    if (evaluated == pp.max_anchors) { S2 = std::max(S2, an.bound); break; }
    if (id_strand < 0 && evaluated == pp.max_unrealized) return po;  // give up: no claim
    ++evaluated;
    done[an.mate][an.chain] = 1;  // every path below accounts for its pairs (or gives up)
    const Mate& X = M[an.mate];
    const Mate& Y = M[1 - an.mate];
    const Cluster& cl = X.ws->clusters[an.chain];
    const int st = cl.strand;
    // Only pairs above `relevant` can change R1, or S2 within the margin.
    // Such a pair loses < Ltot - relevant in total, so each mate's gaps fit
    // in a band of +-apad (narrower than the chains' pad once R1 is known).
    const int relevant = id_strand >= 0 ? std::max(S2, R1 - pp.mapq_margin) : S2;
    const int apad = std::max(0, std::min(pad, Ltot - relevant - 7));
    Band bx, bu, be;
    const int sx = chain_best(ix, X.L, cl, apad, *X.ws, &bx.lo, &bx.B);
    const int need = relevant - sx;
    if (need >= Y.L * kMatch) { S2 = std::max(S2, sx + Y.L * kMatch); continue; }  // none relevant here
    // Union window (any placement of mate x in this chain) -> up.
    int64_t ulo, uhi, elo, ehi;
    partner_range(st, X.slo[an.chain], X.shi[an.chain], &ulo, &uhi);
    const uint8_t* ry = Y.ws->seq[1 - st];
    int thr;
    const int fu = window_search(ix, ry, Y.L, ulo, uhi, pp.wq, need, sc, &bu, &thr);
    const int up = sx + std::max(fu, thr);
    if (fu == kNoSub) { S2 = std::max(S2, up); continue; }  // nothing in the union: nothing exact
    // A pair actually realized: mate x at one best place in the chain, mate y
    // the best in that place's exact window.
    Result ax, ay;
    if (!locate(X.ws->seq[st], X.L, bx, sx, INT64_MIN, INT64_MAX, &ax)) return po;
    partner_range(st, ax.ref_pos, ax.ref_pos, &elo, &ehi);
    const int fe = window_search(ix, ry, Y.L, elo, ehi, pp.wq, need, sc, &be, &thr);
    if (fe == kNoSub) { S2 = std::max(S2, up); continue; }
    if (!locate(ry, Y.L, be, fe, elo, ehi, &ay)) return po;
    const int r = sx + fe;
    const Result& a1 = an.mate == 0 ? ax : ay;
    const Result& a2 = an.mate == 0 ? ay : ax;
    const int s1 = an.mate == 0 ? st : 1 - st;  // mate 1's strand
    const Win win{an.mate, an.chain, apad, sx, fu, thr, need, relevant, ulo, uhi};
    if (id_strand >= 0 && same_pair(s1, a1, a2)) {
      upb = confirmed ? std::max(upb, up) : up;
      confirmed = true;
      wins.push_back(win);
      if (r > R1) {  // the full band found a better alignment of the same pair
        R1 = r;
        best[0] = a1;
        best[1] = a2;
      }
    } else if (r > R1) {
      if (id_strand >= 0) S2 = std::max(S2, confirmed ? upb : R1);
      R1 = r;
      upb = up;
      confirmed = true;
      id_strand = s1;
      best[0] = a1;
      best[1] = a2;
      wins.assign(1, win);
    } else {
      S2 = std::max(S2, up);
    }
  }
  if (id_strand < 0 || !confirmed || R1 <= floor || upb > R1 || S2 > R1) return po;
  // The other pairs through R1's anchors. Pairs at or below `relevant` (the
  // anchor's bands were narrowed to it) are bounded by it.
  for (const Win& w : wins) {
    const Mate& X = M[w.mate];
    const Mate& Y = M[1 - w.mate];
    const Cluster& cl = X.ws->clusters[w.chain];
    const Result& ax = best[w.mate];
    const Result& ay = best[1 - w.mate];
    auto span = [](const Result& a) {
      int64_t n = 0;
      for (int c = 0; c < a.n_cigar; ++c)
        if ((a.cigar[c] & 0xF) == kOpM || (a.cigar[c] & 0xF) == kOpD) n += a.cigar[c] >> 4;
      return n;
    };
    S2 = std::max(S2, w.relevant);
    // x where R1 has it (<= sx), y starting before R1's y or after its end.
    const uint8_t* ry = Y.ws->seq[1 - cl.strand];
    auto ywin = [&](int64_t a, int64_t b) {
      if (b < a) return;
      Band bd;
      int t;
      const int v = window_search(ix, ry, Y.L, a, b, pp.wq, w.need, sc, &bd, &t);
      S2 = std::max(S2, w.sx + std::max(v, t));
    };
    ywin(w.ulo, std::min(w.uhi, ay.ref_pos - 1));
    ywin(std::max(w.ulo, ay.ref_pos + span(ay)), w.uhi);
    // x elsewhere in its chain (starting before R1's x or after its end), y
    // anywhere in the window (<= max(fu, thr)).
    int64_t blo;
    int bB;
    const int xo = std::max(chain_best(ix, X.L, cl, w.apad, *X.ws, &blo, &bB, INT64_MIN, ax.ref_pos - 1),
                            chain_best(ix, X.L, cl, w.apad, *X.ws, &blo, &bB, ax.ref_pos + span(ax), INT64_MAX));
    if (xo != kNoSub) S2 = std::max(S2, xo + std::max(w.fu, w.thr));
  }
  S2 = std::min(S2, R1);  // every pair scores <= R1 = OPT
  po.certified = true;
  po.score = R1;
  po.sub = S2;
  Result* rr[2] = {&r1, &r2};
  for (int x = 0; x < 2; ++x) {
    Result& o = *rr[x];
    const Result& a = best[x];
    o.ref_pos = a.ref_pos;
    o.n_cigar = a.n_cigar;
    for (int c = 0; c < a.n_cigar; ++c) o.cigar[c] = a.cigar[c];
    o.nm = a.nm;
    o.score = a.score;
    o.strand = (uint8_t)(x == 0 ? id_strand : 1 - id_strand);
    o.certified = 1;
    o.tier = kTierPR;
    o.reason = kOk;
    o.radius = -1;
    o.d1 = -1;
    o.d2 = -1;
    o.d2x = -1;
    o.n_best = (uint16_t)(S2 == R1 ? 2 : 1);
    o.floor = (int16_t)floor;
    o.sub_score = (int16_t)S2;
    o.parts = (uint8_t)M[x].P;
    o.used = (uint8_t)M[x].m;
    // Witness (Theorem V): the parts this mate enumerated per strand, with
    // support 1. Lemma P0's floor and S2 rest on them.
    for (int st = 0; st < 2; ++st) {
      o.probes[st] = 0;
      for (int k = 0; k < M[x].m; ++k) o.probes[st] |= (uint16_t)(1u << M[x].ws->order[st][k]);
    }
    o.support = 1;
  }
  return po;
}

}  // namespace certa
