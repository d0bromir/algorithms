// CERTA tier R: repeats (host only). EXACT_ALIGNMENT_PROPOSAL.md Aim 2.
//
// Objective: the local score of Theorem L (affine_align, BWA-MEM defaults,
// -5 per clipped end). For one read, tier R returns an interval [lo, hi] that
// provably contains OPT, the best score over the whole reference, both
// strands. lo is the score of an alignment actually computed; the read is
// *exact* when lo == hi.
//
// Parts (Theorem G with a hash index, s = 1). From position s, a part
// [s, b) is grown until it occurs at most tau times, or not at all:
//  - absent: [s, b) contains a q-mer with count 0 (or an N); then the part
//    cannot occur exactly anywhere and needs no enumeration;
//  - rare: the occurrences of [s, b) are exactly the hits of one anchor
//    q-mer inside it whose exact match extends left to s and right to b.
//    Extending every hit of the anchor once gives occ(b) for all b, so the
//    shortest b with occ(b) <= tau is found exactly.
// A part that cannot be made rare within the anchor limits ends the
// partition (that stretch of the read occurs > tau times, or its anchors
// are too frequent to extend). Parts are disjoint and >= q long, so
// Lemma L applies with t = 1:
//  - an alignment in which no part occurs exactly scores <= U (span bound);
//  - an alignment in which some rare part occurs exactly lies, if it scores
//    above U_eff = max(U, L - 7 - kPadMax), within +-pad of that
//    occurrence's diagonal (pad = L - U_eff - 7), in a chain C whose parts
//    in pmask(C) are the only ones it can match exactly, so it scores
//    <= ub(C) = span bound with pmask(C) free.
// Chains are evaluated best bound first and skipped once ub <= max(lo,
// U_eff). Hence OPT <= hi = max(lo, U_eff), and OPT = lo when lo >= U_eff.
//
// Lower bound. Besides the chains, the read is aligned in a band around the
// copies whose exact match with the read is longest (from the anchor
// extensions). In a near-identical repeat these copies are the best ones,
// and every read error or SNP makes the parts that cover it absent, so lo
// reaches U and the read is exact although it occurs hundreds of times.
//
// Ties (Theorem R). If the reference string covered by the reported
// alignment occurs at another position (non-overlapping, or on the other
// strand), the same alignment exists there: two loci score lo. With lo ==
// OPT that proves a tie; otherwise it proves SUB >= lo, so OPT - SUB <= hi - lo.
//
// SUB interval (Aim 3). SUB is the best score at a locus not overlapping
// the reference interval [p, p + span) of A* on its strand. sub_lo is the
// best score of an evaluated band lying entirely at another locus (or lo,
// with a copy). sub_hi bounds SUB: SUB <= OPT <= hi always; for an exact
// read, every other alignment scores <= U_eff or lies in a chain, and the
// alignments of a chain at other loci are bounded by its best score (when
// the whole chain lies away from A*), else by its best alignments starting
// before p or at/after p + span (a superset); unevaluated chains by ub.
#pragma once
#include <algorithm>
#include <cstdint>
#include <vector>

#include "certa/core.h"
#include "certa/longread.h"
#include "certa/pair.h"

namespace certa {

struct RepeatParams {
  int tau = 16;     // max exact occurrences of a rare part
  int cmax = 4096;  // max hits of an anchor q-mer that are extended
  int window = 40;  // anchors are searched in [s, s + window]
  int samples = 8;  // copies per strand aligned for the lower bound
  int spad = 8;     // band half-width around a sampled copy
  // SUB bound of exact reads: chains are evaluated down to lo - sub_margin
  // (0: sub_hi = lo, no extra work).
  int sub_margin = 0;
  // When no q-mer of the read has 1..cmax hits (satellites), this many hits
  // per strand of its rarest q-mer are aligned anyway, for lo only (0 = off).
  int fallback_samples = 0;
  // Lower bound only: this many more chains, best bound first, are aligned
  // although their bound cannot beat U_eff (diverged repeats, where every
  // part but a few is absent and the sampled copies align poorly).
  int lb_chains = 0;
};

struct RepeatOut {
  bool claim = false;  // false: read not handled (bad length, index with s > 1)
  bool exact = false;
  int lo = kNoSub;     // score of A* (kNoSub if no alignment was computed)
  int hi = 0;          // OPT <= hi
  int copies = 0;      // exact copies of A*'s reference string found (A* included, max 2)
  int64_t copy_pos[2] = {-1, -1};
  int copy_strand[2] = {0, 0};
  int nparts = 0, nabsent = 0, nocc = 0;
  int np[2] = {0, 0};    // the parts per strand [poff, pend) (strand 1: in the
  int16_t poff[2][PMAX], pend[2][PMAX];  // reverse complement), the witness
  bool has_aln = false;  // aln holds A* (its CIGAR fits)
  Result aln;            // A*
  int sub_lo = kNoSub;   // SUB >= sub_lo (an alignment at another locus)
  int sub_hi = 0;        // SUB <= sub_hi
};

struct RepeatScratch {
  std::vector<uint64_t> clo, chi;            // q-mer index ranges per read position
  std::vector<uint64_t> keys;                // q-mer keys per read position
  std::vector<uint8_t> kok;                  // ... and whether each is N-free
  std::vector<uint32_t> kcnt;                // ... and its occurrence count
  std::vector<std::pair<int, int64_t>> elig;  // (right extension, hit position)
  struct Samp { int len; int64_t diag; };
  std::vector<Samp> samp[2];
  std::vector<uint8_t> win;                  // A*'s reference string
  struct Eval { int score, strand; int64_t lo; int B; };
  std::vector<Eval> evals;                   // bands evaluated (lower bounds)
  std::vector<int> cval;                     // per chain: best score, or kUnset
  std::vector<uint8_t> tmp;                  // find_copies: the string in one orientation
};

constexpr int kUnset = -BIG;

namespace repeat_detail {

inline bool same_base(uint8_t a, uint8_t b) { return a < 4 && a == b; }

// How many leading positions i < n have same_base(rd[i], ref[j + i]), as the
// scalar loop counts them; 8 positions per step where the window lies inside
// the reference (the first mismatching byte of a word is its lowest set byte
// on little-endian hosts).
inline int match_run(const IndexView& ix, const uint8_t* rd, int64_t j, int n) {
  int i = 0;
#if defined(__GNUC__) && !defined(__CUDACC__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
  if (j >= 0 && (uint64_t)(j + n) <= ix.ref_len) {
    const uint8_t* rf = ix.ref + j;
    const uint64_t high = 0x8080808080808080ull;
    for (; i + 8 <= n; i += 8) {
      uint64_t a, b;
      std::memcpy(&a, rd + i, 8);
      std::memcpy(&b, rf + i, 8);
      // Per byte (bases are 0..4, so no sum carries into the next byte):
      // (a ^ b) + 0x7f has bit 7 set when the bases differ, a + 0x7c when
      // the read base is N.
      const uint64_t bad = (((a ^ b) + 0x7f7f7f7f7f7f7f7full) | (a + 0x7c7c7c7c7c7c7c7cull)) & high;
      if (bad) return i + (__builtin_ctzll(bad) >> 3);
    }
  }
#endif
  for (; i < n && same_base(rd[i], ref_at(ix, j + i)); ++i) {}
  return i;
}

// Prefetches the reference bytes [a, a + n), at most 256 of them: a copy whose
// comparison with the read is a few hits away (hits are scattered over the
// reference, so each comparison would otherwise wait on memory).
inline void prefetch_ref(const IndexView& ix, int64_t a, int n) {
  if (a < 0) return;
  const int64_t e = a + std::min(n, 256);
  for (int64_t x = a & ~int64_t(63); x < e && (uint64_t)x < ix.ref_len; x += 64) __builtin_prefetch(ix.ref + x);
}
constexpr uint64_t kPrefetchAhead = 16;  // hits

// Exact non-overlapping copies of w (length n) in the reference, beside
// [p, p + n) on its own strand; stops once `need` are found. Copies on the
// other strand are searched with the reverse complement. Returns the number
// found and their positions/strands. Uses sc.tmp, keys, kok, clo and kcnt.
inline int find_copies(const IndexView& ix, const std::vector<uint8_t>& w, int64_t p, int strand_of_w,
                       int need, int cmax, int64_t* pos, int* strand, RepeatScratch& sc) {
  const int n = (int)w.size(), q = ix.q;
  if (n < q) return 0;
  for (uint8_t c : w) if (c > 3) return 0;  // N never matches: no exact copy
  const int nq = n - q + 1;
  std::vector<uint8_t>& tmp = sc.tmp;
  int found = 0;
  for (int orient = 0; orient < 2 && found < need; ++orient) {
    tmp.resize(n);
    for (int i = 0; i < n; ++i) tmp[i] = orient == 0 ? w[i] : (uint8_t)(3 - w[n - 1 - i]);
    // Rarest q-mer of the string (the first of the rarest), the lookups
    // pipelined (lookup_batch gives lookup's ranges).
    sc.keys.resize(nq);
    sc.kok.resize(nq);
    sc.clo.resize(nq);
    sc.kcnt.resize(nq);
    qmer_keys(tmp.data(), n, q, sc.keys.data(), sc.kok.data());
    lookup_batch(ix, sc.keys.data(), sc.kok.data(), nq, sc.clo.data(), sc.kcnt.data());
    int best = 0;
    for (int i = 1; i < nq; ++i)
      if (sc.kcnt[i] < sc.kcnt[best]) best = i;
    const uint64_t blo = sc.clo[best], bhi = blo + sc.kcnt[best];
    if (bhi - blo > (uint64_t)cmax) continue;
    for (uint64_t e = blo; e < std::min(bhi, blo + kPrefetchAhead); ++e) prefetch_ref(ix, (int64_t)ix.pos[e] - best, n);
    for (uint64_t e = blo; e < bhi && found < need; ++e) {
      if (e + kPrefetchAhead < bhi) prefetch_ref(ix, (int64_t)ix.pos[e + kPrefetchAhead] - best, n);
      const int64_t start = (int64_t)ix.pos[e] - best;
      if (orient == 0 && start + n > p && p + n > start) continue;  // overlaps A*'s own interval
      if (start < 0 || match_run(ix, tmp.data(), start, n) < n) continue;
      bool dup = false;
      for (int k = 0; k < found; ++k)
        dup |= strand[k] == (orient == 0 ? strand_of_w : 1 - strand_of_w) && start + n > pos[k] &&
               pos[k] + n > start;
      if (dup) continue;
      pos[found] = start;
      strand[found] = orient == 0 ? strand_of_w : 1 - strand_of_w;
      ++found;
    }
  }
  return found;
}

}  // namespace repeat_detail

inline RepeatOut repeat_tier(const IndexView& ix, const RepeatParams& rp, const uint8_t* codes, int L,
                             Workspace& ws, RepeatScratch& sc) {
  using namespace repeat_detail;
  RepeatOut out;
  const int q = ix.q;
  if (ix.s != 1 || L > LMAX || L < q) return out;
  out.claim = true;
  for (int i = 0; i < L; ++i) {
    ws.seq[0][i] = codes[i];
    const uint8_t c = codes[L - 1 - i];
    ws.seq[1][i] = c < 4 ? (uint8_t)(3 - c) : (uint8_t)4;
  }
  int off[2][PMAX], end[2][PMAX], np[2] = {0, 0};
  int fb_i[2] = {-1, -1};  // rarest q-mer that occurs, per strand (fallback samples)
  uint64_t fb_n[2] = {0, 0}, fb_lo[2] = {0, 0};
  const int nq = L - q + 1;
  sc.clo.resize(nq);
  sc.chi.resize(nq);
  for (int st = 0; st < 2; ++st) {
    const uint8_t* rd = ws.seq[st];
    sc.samp[st].clear();
    ws.ncand[st] = 0;
    // Every q-mer's index range, the lookups pipelined (lookup_batch gives
    // lookup's ranges; a q-mer with N: count 0, it never occurs exactly).
    sc.keys.resize(nq);
    sc.kok.resize(nq);
    sc.kcnt.resize(nq);
    qmer_keys(rd, L, q, sc.keys.data(), sc.kok.data());
    lookup_batch(ix, sc.keys.data(), sc.kok.data(), nq, sc.clo.data(), sc.kcnt.data());
    for (int i = 0; i < nq; ++i) sc.chi[i] = sc.clo[i] + sc.kcnt[i];
    auto cnt = [&](int i) { return sc.chi[i] - sc.clo[i]; };
    for (int i = 0; i < nq; ++i)
      if (cnt(i) > 0 && (fb_i[st] < 0 || cnt(i) < fb_n[st])) { fb_i[st] = i; fb_n[st] = cnt(i); fb_lo[st] = sc.clo[i]; }
    int s = 0;
    while (s < nq && np[st] < PMAX) {
      int babs = BIG;
      for (int i = s; i < nq; ++i)
        if (cnt(i) == 0) { babs = i + q; break; }
      int c = -1;
      for (int i = s; i < nq && i <= s + rp.window; ++i)
        if (cnt(i) > 0 && cnt(i) <= (uint64_t)rp.cmax && (c < 0 || cnt(i) < cnt(c))) c = i;
      int banc = BIG;
      sc.elig.clear();
      if (c >= 0) {
        // Each hit's copy is compared with rd from s on: its reference bytes
        // are prefetched kPrefetchAhead hits ahead.
        const uint64_t e0 = sc.clo[c], e1 = sc.chi[c];
        for (uint64_t e = e0; e < std::min(e1, e0 + kPrefetchAhead); ++e)
          prefetch_ref(ix, (int64_t)ix.pos[e] - (c - s), L - s);
        for (uint64_t e = e0; e < e1; ++e) {
          if (e + kPrefetchAhead < e1) prefetch_ref(ix, (int64_t)ix.pos[e + kPrefetchAhead] - (c - s), L - s);
          const int64_t h = ix.pos[e];
          // The copy must match rd[s, c) exactly to its left; then extend right.
          if (c > s && match_run(ix, rd + s, h - (c - s), c - s) < c - s) continue;
          const int R = c + q < L ? match_run(ix, rd + c + q, h + q, L - c - q) : 0;
          sc.elig.push_back({R, h});
        }
        std::sort(sc.elig.begin(), sc.elig.end(),
                  [](const std::pair<int, int64_t>& a, const std::pair<int, int64_t>& b) {
                    return a.first != b.first ? a.first > b.first : a.second < b.second;
                  });
        if ((int)sc.elig.size() <= rp.tau) banc = c + q;
        else if (c + q + sc.elig[rp.tau].first < L) banc = c + q + sc.elig[rp.tau].first + 1;
        // The copies matching the read longest here are the lower-bound samples.
        for (int k = 0; k < (int)sc.elig.size() && k < rp.samples; ++k)
          sc.samp[st].push_back({(c - s) + q + sc.elig[k].first, sc.elig[k].second - c});
      }
      const int b = std::min(babs, banc);
      if (b > L) break;  // the rest occurs > tau times (or has no affordable anchor)
      const bool rare = b == banc && b < babs;  // otherwise absent: contains a count-0 q-mer
      int occ = 0;
      if (rare)
        for (const auto& x : sc.elig) occ += c + q + x.first >= b;
      // Every occurrence must be enumerated for Lemma L: if they do not fit,
      // end the partition before this part (fewer parts, still sound).
      if (ws.ncand[st] + occ > MAX_CAND) break;
      const int j = np[st]++;
      off[st][j] = s;
      end[st][j] = b;
      if (rare)
        for (const auto& x : sc.elig)
          if (c + q + x.first >= b) ws.cand[st][ws.ncand[st]++] = (x.second - c) * 16 + j;
      out.nocc += occ;
      out.nabsent += occ == 0;
      s = b;
    }
    out.nparts += np[st];
    out.np[st] = np[st];
    for (int j = 0; j < np[st]; ++j) {
      out.poff[st][j] = (int16_t)off[st][j];
      out.pend[st][j] = (int16_t)end[st][j];
    }
  }

  // Upper bound for alignments with no exact part, and the band they force.
  int U = -BIG;
  for (int st = 0; st < 2; ++st) {
    const uint32_t all = np[st] ? ((1u << np[st]) - 1) : 0;
    U = std::max(U, span_bound(L, np[st], off[st], end[st], all));
  }
  const int Ueff = std::max(U, L * kMatch - 7 - kPadMax);
  const int pad = std::max(0, L * kMatch - Ueff - 7);

  // Satellites: no anchor was affordable on either strand.
  if (rp.fallback_samples > 0 && sc.samp[0].empty() && sc.samp[1].empty())
    for (int st = 0; st < 2; ++st)
      for (uint64_t e = fb_lo[st]; fb_i[st] >= 0 && e < fb_lo[st] + fb_n[st] && e < fb_lo[st] + rp.fallback_samples; ++e)
        sc.samp[st].push_back({q, (int64_t)ix.pos[e] - fb_i[st]});

  // Lower bound from the sampled copies (longest exact matches first).
  int lo = kNoSub, lo_st = 0, lo_B = 0;
  int64_t lo_lo = 0;
  sc.evals.clear();
  const int nsamp = std::max(rp.samples, rp.fallback_samples);
  for (int st = 0; st < 2; ++st) {
    auto& v = sc.samp[st];
    std::sort(v.begin(), v.end(), [](const RepeatScratch::Samp& a, const RepeatScratch::Samp& b) {
      return a.len != b.len ? a.len > b.len : a.diag < b.diag;
    });
    for (int k = 0; k < (int)v.size() && k < nsamp; ++k) {
      if (k > 0 && v[k].diag == v[k - 1].diag) continue;
      const int64_t blo = v[k].diag - rp.spad;
      const int B = 2 * rp.spad + 1;
      const int sc_ = band_score(ix, ws.seq[st], L, blo, B);
      if (sc_ != kNoSub) sc.evals.push_back({sc_, st, blo, B});
      if (sc_ > lo) { lo = sc_; lo_st = st; lo_lo = blo; lo_B = B; }
    }
  }

  // Chains of rare-part occurrences, best bound first (branch and bound).
  for (int st = 0; st < 2; ++st) sort_i64(ws.cand[st], ws.ncand[st]);
  const int nch = build_chains(ws, pad);
  std::vector<std::pair<int, int>> order;
  order.reserve(nch);
  for (int x = 0; x < nch; ++x) {
    const Cluster& cl = ws.clusters[x];
    const int st = cl.strand;
    const uint32_t all = np[st] ? ((1u << np[st]) - 1) : 0;
    order.push_back({span_bound(L, np[st], off[st], end[st], all & ~(uint32_t)cl.pmask), x});
  }
  std::sort(order.begin(), order.end(), [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
    return a.first != b.first ? a.first > b.first : a.second < b.second;
  });
  sc.cval.assign(nch, kUnset);
  for (const auto& o : order) {
    if (o.first <= std::max(lo, Ueff)) break;
    int64_t blo;
    int B;
    const Cluster& cl = ws.clusters[o.second];
    const int v = chain_best(ix, L, cl, pad, ws, &blo, &B);
    sc.cval[o.second] = v;
    if (v != kNoSub) sc.evals.push_back({v, cl.strand, blo, B});
    if (v > lo) { lo = v; lo_st = cl.strand; lo_lo = blo; lo_B = B; }
  }
  int extra = rp.lb_chains;
  for (const auto& o : order) {
    if (extra == 0 || o.first <= lo) break;
    if (sc.cval[o.second] != kUnset) continue;
    --extra;
    int64_t blo;
    int B;
    const Cluster& cl = ws.clusters[o.second];
    // A wider band than pad only enlarges the set searched, so cval stays
    // an upper bound for the chain's alignments within pad.
    const int v = chain_best(ix, L, cl, std::max(pad, rp.spad), ws, &blo, &B);
    sc.cval[o.second] = v;
    if (v != kNoSub) sc.evals.push_back({v, cl.strand, blo, B});
    if (v > lo) { lo = v; lo_st = cl.strand; lo_lo = blo; lo_B = B; }
  }
  out.hi = std::max(lo, Ueff);
  out.lo = lo;
  out.exact = lo != kNoSub && lo >= Ueff;
  out.sub_hi = out.hi;  // SUB <= OPT <= hi
  if (lo == kNoSub) return out;

  // A* and its copies (Theorem R).
  Result a;
  if (!affine_align(ix, ws.seq[lo_st], L, lo_lo, lo_B, ws, a) || a.score != lo) return out;
  a.strand = (uint8_t)lo_st;
  out.aln = a;
  out.has_aln = true;
  int64_t span = 0;
  for (int i = 0; i < a.n_cigar; ++i) {
    const uint32_t op = a.cigar[i] & 0xF;
    if (op == kOpM || op == kOpD) span += a.cigar[i] >> 4;
  }
  sc.win.resize((size_t)span);
  for (int64_t i = 0; i < span; ++i) sc.win[(size_t)i] = ref_at(ix, a.ref_pos + i);
  out.copy_pos[0] = a.ref_pos;
  out.copy_strand[0] = lo_st;
  out.copies = 1 + find_copies(ix, sc.win, a.ref_pos, lo_st, 1, rp.cmax, &out.copy_pos[1], &out.copy_strand[1], sc);

  // SUB interval. A band [blo, blo + B) holds alignments starting at or
  // after blo and ending at or before blo + B + L - 2.
  const int64_t p0 = a.ref_pos;
  auto away = [&](int st, int64_t blo, int B) {
    return st != lo_st || blo >= p0 + span || blo + B + L - 2 < p0;
  };
  for (const auto& e : sc.evals)
    if (away(e.strand, e.lo, e.B) && e.score > out.sub_lo) out.sub_lo = e.score;
  if (out.copies >= 2) out.sub_lo = lo;
  if (!out.exact) return out;  // sub_hi = hi
  if (out.copies >= 2 || rp.sub_margin <= 0) {
    out.sub_hi = lo;
    return out;
  }
  int sh = Ueff, budget = 32;
  for (const auto& o : order) {
    if (o.first <= std::max(sh, lo - rp.sub_margin)) {  // the rest: ub <= o.first
      sh = std::max(sh, o.first);
      break;
    }
    if (sh >= lo || budget-- == 0) {  // SUB <= OPT = lo
      sh = lo;
      break;
    }
    const Cluster& cl = ws.clusters[o.second];
    const int64_t* c = ws.cand[cl.strand];
    const int64_t dmin = diag_of(c[cl.beg]), dmax = diag_of(c[cl.end - 1]);
    int64_t blo;
    int B;
    int b;
    if (away(cl.strand, dmin - pad, (int)(dmax - dmin) + 2 * pad + 1)) {
      b = sc.cval[o.second] != kUnset ? sc.cval[o.second] : chain_best(ix, L, cl, pad, ws, &blo, &B);
      if (b > out.sub_lo) out.sub_lo = b;
    } else {
      const int left = chain_best(ix, L, cl, pad, ws, &blo, &B, INT64_MIN, p0 - 1);
      const int right = chain_best(ix, L, cl, pad, ws, &blo, &B, p0 + span, INT64_MAX);
      if (right > out.sub_lo) out.sub_lo = right;  // starts after the end of A*: another locus
      b = std::max(left, right);
    }
    sh = std::max(sh, b);
  }
  out.sub_hi = std::min(sh, lo);
  return out;
}

}  // namespace certa
