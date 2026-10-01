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
// anchor's up, every skipped anchor's bound and the floor are <= R1. The
// largest of these other bounds, S2, bounds every other pair (MAPQ).
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
  int max_anchors = 64;
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
inline int span_bound(int L, int n, const int* off, const int* end, uint32_t counted) {
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

// Scratch for the window search (one per thread).
struct PairScratch {
  std::unique_ptr<Workspace> ws{new Workspace};  // DP rows and traceback
  std::vector<std::pair<uint32_t, int>> keys;    // part key -> part
  std::vector<std::pair<int64_t, int>> hits;     // (diagonal, part)
};

// Best alignment of `rd` whose first aligned base lies in [smin, smax], found
// exhaustively: every exact occurrence of its ~wq-base parts in the reference
// window is scanned, and chains of them are aligned in bands. Returns the
// score (> *thr) and the alignment, or kNoSub; *thr bounds every alignment
// not returned. Returns kNoSub with *thr = L when the window is too
// repetitive to search.
inline int window_search(const IndexView& ix, const uint8_t* rd, int L, int64_t smin, int64_t smax,
                         int wq, PairScratch& sc, Result* out, int* thr) {
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
  *thr = t;
  if (t >= L * kMatch) return kNoSub;
  const int pad = L * kMatch - t - 7;
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
    auto it = std::lower_bound(sc.keys.begin(), sc.keys.end(), std::make_pair(key, -1));
    for (; it != sc.keys.end() && it->first == key; ++it) {
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
  int best = kNoSub;
  int64_t best_lo = 0;
  int best_B = 0;
  size_t i = 0;
  while (i < sc.hits.size()) {
    size_t j = i + 1;
    uint32_t pm = 1u << sc.hits[i].second;
    while (j < sc.hits.size() && sc.hits[j].first - sc.hits[j - 1].first <= pad) {
      pm |= 1u << sc.hits[j].second;
      ++j;
    }
    if (span_bound(L, n, off, end, all & ~pm) > t) {
      size_t a = i;
      while (a < j) {
        const int64_t d0 = sc.hits[a].first;
        size_t e = a + 1;
        while (e < j && sc.hits[e].first - d0 <= BMAX - 2 * pad - 1) ++e;
        const int64_t blo = d0 - pad;
        const int B = (int)(sc.hits[e - 1].first - d0) + 2 * pad + 1;
        Result tmp;
        if (affine_align(ix, rd, L, blo, B, *sc.ws, tmp, true, smin, smax) && tmp.score > best) {
          best = tmp.score;
          best_lo = blo;
          best_B = B;
        }
        a = e;
      }
    }
    i = j;
  }
  if (best <= t) return kNoSub;
  if (out) {
    Result tmp;
    if (!affine_align(ix, rd, L, best_lo, best_B, *sc.ws, tmp, false, smin, smax) || tmp.score != best) {
      *thr = L * kMatch;  // CIGAR does not fit: no claim
      return kNoSub;
    }
    *out = tmp;
  }
  return best;
}

// One mate: its enumerated parts, chains and bounds.
struct Mate {
  const uint8_t* codes = nullptr;
  int L = 0, P = 0, m = 0, U = 0, B = 0, nch = 0;
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
  for (int st = 0; st < 2; ++st)
    for (int x = 0; x < 256; ++x) ws.ubmemo[st][x] = kNoSub;
  M.B = M.U;
  for (int x = 0; x < M.nch; ++x) {
    const Cluster& cl = ws.clusters[x];
    int16_t& u = ws.ubmemo[cl.strand][cl.pmask];
    if (u == kNoSub) u = (int16_t)local_bound(ws, cl.strand, M.L, M.P, M.m, cl.pmask, 0);
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
  po.floor = floor;
  if (floor >= Ltot) return po;
  const int pad = Ltot - floor - 7;
  chain_mate(M[0], pad);
  chain_mate(M[1], pad);
  const int64_t D = pp.max_dist;
  // Partner start range of an anchor: the forward mate starts first.
  auto partner_range = [&](int strand, int64_t a_lo, int64_t a_hi, int64_t* plo, int64_t* phi) {
    if (strand == 0) { *plo = a_lo; *phi = a_hi + D; }
    else { *plo = a_lo - D; *phi = a_hi; }
  };
  struct Anchor { int bound, mate, chain; };
  std::vector<Anchor> anchors;
  for (int x = 0; x < 2; ++x) {
    const Mate& X = M[x];
    const Mate& Y = M[1 - x];
    for (int c = 0; c < X.nch; ++c) {
      const int st = X.ws->clusters[c].strand;
      int64_t plo, phi;
      partner_range(st, X.slo[c], X.shi[c], &plo, &phi);
      int ny = Y.U;  // mate y near this chain: no exact part, or one of its chains
      for (int c2 = 0; c2 < Y.nch; ++c2)
        if (Y.ws->clusters[c2].strand == 1 - st && Y.shi[c2] >= plo && Y.slo[c2] <= phi)
          ny = std::max(ny, Y.ub[c2]);
      const int b = X.ub[c] + ny;
      if (b > floor) anchors.push_back({b, x, c});
    }
  }
  std::sort(anchors.begin(), anchors.end(), [](const Anchor& a, const Anchor& b) {
    if (a.bound != b.bound) return a.bound > b.bound;
    if (a.mate != b.mate) return a.mate < b.mate;
    return a.chain < b.chain;
  });
  int R1 = floor, upb = -BIG, S2 = floor, evaluated = 0;
  int64_t id_pos[2] = {0, 0};
  int id_strand = -1;
  Result best[2];
  auto same_pair = [&](int s1, const Result& a1, const Result& a2) {
    return id_strand == s1 && std::llabs(id_pos[0] - a1.ref_pos) <= 2 * pad + 1 &&
           std::llabs(id_pos[1] - a2.ref_pos) <= 2 * pad + 1;
  };
  for (const Anchor& an : anchors) {
    if (an.bound <= S2) break;  // sorted: no remaining pair can matter
    if (evaluated == pp.max_anchors) { S2 = std::max(S2, an.bound); break; }
    ++evaluated;
    const Mate& X = M[an.mate];
    const Mate& Y = M[1 - an.mate];
    const Cluster& cl = X.ws->clusters[an.chain];
    const int st = cl.strand;
    int64_t blo = 0;
    int bB = 0;
    chain_best(ix, X.L, cl, pad, *X.ws, &blo, &bB);
    Result ax;
    if (!affine_align(ix, X.ws->seq[st], X.L, blo, bB, *X.ws, ax)) return po;  // CIGAR overflow
    // Union window (any placement of mate x in this chain) -> up.
    int64_t ulo, uhi, elo, ehi;
    partner_range(st, X.slo[an.chain], X.shi[an.chain], &ulo, &uhi);
    partner_range(st, ax.ref_pos, ax.ref_pos, &elo, &ehi);
    const uint8_t* ry = Y.ws->seq[1 - st];
    Result ay;
    int thr_u;
    const int fu = window_search(ix, ry, Y.L, ulo, uhi, pp.wq, sc, &ay, &thr_u);
    const int up = ax.score + std::max(fu, thr_u);
    int fe = kNoSub;
    if (fu != kNoSub && ay.ref_pos >= elo && ay.ref_pos <= ehi) {
      fe = fu;  // the union's best lies in the exact window: it is the exact best
    } else {
      int thr_e;
      fe = window_search(ix, ry, Y.L, elo, ehi, pp.wq, sc, &ay, &thr_e);
    }
    if (fe == kNoSub) { S2 = std::max(S2, up); continue; }
    const int r = ax.score + fe;
    const Result& a1 = an.mate == 0 ? ax : ay;
    const Result& a2 = an.mate == 0 ? ay : ax;
    const int s1 = an.mate == 0 ? st : 1 - st;  // mate 1's strand
    if (id_strand >= 0 && same_pair(s1, a1, a2)) {
      upb = std::max(upb, up);
    } else if (r > R1) {
      if (id_strand >= 0) S2 = std::max(S2, upb);
      R1 = r;
      upb = up;
      id_strand = s1;
      id_pos[0] = a1.ref_pos;
      id_pos[1] = a2.ref_pos;
      best[0] = a1;
      best[1] = a2;
    } else {
      S2 = std::max(S2, up);
    }
  }
  if (id_strand < 0 || R1 <= floor || upb > R1 || S2 > R1) return po;
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
  }
  return po;
}

}  // namespace certa
