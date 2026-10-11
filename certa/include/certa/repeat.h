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
#pragma once
#include <algorithm>
#include <cstdint>
#include <vector>

#include "certa/core.h"
#include "certa/pair.h"

namespace certa {

struct RepeatParams {
  int tau = 16;     // max exact occurrences of a rare part
  int cmax = 4096;  // max hits of an anchor q-mer that are extended
  int window = 40;  // anchors are searched in [s, s + window]
  int samples = 8;  // copies per strand aligned for the lower bound
  int spad = 8;     // band half-width around a sampled copy
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
  Result aln;          // A*, when lo != kNoSub
};

struct RepeatScratch {
  std::vector<uint64_t> clo, chi;            // q-mer index ranges per read position
  std::vector<std::pair<int, int64_t>> elig;  // (right extension, hit position)
  struct Samp { int len; int64_t diag; };
  std::vector<Samp> samp[2];
  std::vector<uint8_t> win;                  // A*'s reference string
};

namespace repeat_detail {

inline bool same_base(uint8_t a, uint8_t b) { return a < 4 && a == b; }

// Exact non-overlapping copies of w (length n) in the reference, beside
// [p, p + n) on its own strand; stops once `need` are found. Copies on the
// other strand are searched with the reverse complement. Returns the number
// found and their positions/strands.
inline int find_copies(const IndexView& ix, const std::vector<uint8_t>& w, int64_t p, int strand_of_w,
                       int need, int cmax, int64_t* pos, int* strand, std::vector<uint8_t>& tmp) {
  const int n = (int)w.size(), q = ix.q;
  if (n < q) return 0;
  for (uint8_t c : w) if (c > 3) return 0;  // N never matches: no exact copy
  int found = 0;
  for (int orient = 0; orient < 2 && found < need; ++orient) {
    tmp.resize(n);
    for (int i = 0; i < n; ++i) tmp[i] = orient == 0 ? w[i] : (uint8_t)(3 - w[n - 1 - i]);
    // Rarest q-mer of the string.
    int best = -1;
    uint64_t blo = 0, bhi = 0;
    for (int i = 0; i + q <= n; ++i) {
      uint64_t key, lo, hi;
      kmer_key(tmp.data() + i, q, &key);
      lookup(ix, key, &lo, &hi);
      if (best < 0 || hi - lo < bhi - blo) { best = i; blo = lo; bhi = hi; }
    }
    if (bhi - blo > (uint64_t)cmax) continue;
    for (uint64_t e = blo; e < bhi && found < need; ++e) {
      const int64_t start = (int64_t)ix.pos[e] - best;
      if (orient == 0 && start + n > p && p + n > start) continue;  // overlaps A*'s own interval
      bool ok = start >= 0;
      for (int i = 0; ok && i < n; ++i) ok = same_base(tmp[i], ref_at(ix, start + i));
      if (!ok) continue;
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
  const int nq = L - q + 1;
  sc.clo.resize(nq);
  sc.chi.resize(nq);
  for (int st = 0; st < 2; ++st) {
    const uint8_t* rd = ws.seq[st];
    sc.samp[st].clear();
    ws.ncand[st] = 0;
    for (int i = 0; i < nq; ++i) {
      uint64_t key;
      if (kmer_key(rd + i, q, &key)) lookup(ix, key, &sc.clo[i], &sc.chi[i]);
      else sc.clo[i] = sc.chi[i] = 0;  // contains N: count 0, never occurs exactly
    }
    auto cnt = [&](int i) { return sc.chi[i] - sc.clo[i]; };
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
        for (uint64_t e = sc.clo[c]; e < sc.chi[c]; ++e) {
          const int64_t h = ix.pos[e];
          bool ok = true;
          for (int k = c - 1; k >= s && ok; --k) ok = same_base(rd[k], ref_at(ix, h - (c - k)));
          if (!ok) continue;
          int R = 0;
          while (c + q + R < L && same_base(rd[c + q + R], ref_at(ix, h + q + R))) ++R;
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
  }

  // Upper bound for alignments with no exact part, and the band they force.
  int U = -BIG;
  for (int st = 0; st < 2; ++st) {
    const uint32_t all = np[st] ? ((1u << np[st]) - 1) : 0;
    U = std::max(U, span_bound(L, np[st], off[st], end[st], all));
  }
  const int Ueff = std::max(U, L * kMatch - 7 - kPadMax);
  const int pad = std::max(0, L * kMatch - Ueff - 7);

  // Lower bound from the sampled copies (longest exact matches first).
  int lo = kNoSub, lo_st = 0, lo_B = 0;
  int64_t lo_lo = 0;
  for (int st = 0; st < 2; ++st) {
    auto& v = sc.samp[st];
    std::sort(v.begin(), v.end(), [](const RepeatScratch::Samp& a, const RepeatScratch::Samp& b) {
      return a.len != b.len ? a.len > b.len : a.diag < b.diag;
    });
    for (int k = 0; k < (int)v.size() && k < rp.samples; ++k) {
      if (k > 0 && v[k].diag == v[k - 1].diag) continue;
      const int64_t blo = v[k].diag - rp.spad;
      const int B = 2 * rp.spad + 1;
      const int sc_ = band_score(ix, ws.seq[st], L, blo, B);
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
  for (const auto& o : order) {
    if (o.first <= std::max(lo, Ueff)) break;
    int64_t blo;
    int B;
    const Cluster& cl = ws.clusters[o.second];
    const int v = chain_best(ix, L, cl, pad, ws, &blo, &B);
    if (v > lo) { lo = v; lo_st = cl.strand; lo_lo = blo; lo_B = B; }
  }
  out.hi = std::max(lo, Ueff);
  out.lo = lo;
  out.exact = lo != kNoSub && lo >= Ueff;
  if (lo == kNoSub) return out;

  // A* and its copies (Theorem R).
  Result a;
  if (!affine_align(ix, ws.seq[lo_st], L, lo_lo, lo_B, ws, a) || a.score != lo) return out;
  a.strand = (uint8_t)lo_st;
  out.aln = a;
  int64_t span = 0;
  for (int i = 0; i < a.n_cigar; ++i) {
    const uint32_t op = a.cigar[i] & 0xF;
    if (op == kOpM || op == kOpD) span += a.cigar[i] >> 4;
  }
  sc.win.resize((size_t)span);
  for (int64_t i = 0; i < span; ++i) sc.win[(size_t)i] = ref_at(ix, a.ref_pos + i);
  out.copy_pos[0] = a.ref_pos;
  out.copy_strand[0] = lo_st;
  std::vector<uint8_t> tmp;
  out.copies = 1 + find_copies(ix, sc.win, a.ref_pos, lo_st, 1, rp.cmax, &out.copy_pos[1], &out.copy_strand[1], tmp);
  return out;
}

}  // namespace certa
