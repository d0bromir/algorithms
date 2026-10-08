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

struct WinChain { int ub, beg, end; };

// Scratch for the window search (one per thread).
struct PairScratch {
  std::vector<WinChain> chains;
  std::unique_ptr<Workspace> ws{new Workspace};  // DP rows and traceback
  std::vector<std::pair<uint32_t, int>> keys;    // part key -> part
  std::vector<std::pair<int64_t, int>> hits;     // (diagonal, part)
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
      const int v = band_score(ix, rd, L, blo, B, smin, smax);
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
  ubmemo_begin(ws);
  M.B = M.U;
  for (int x = 0; x < M.nch; ++x) {
    const Cluster& cl = ws.clusters[x];
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
  const int pad = Ltot - floor - 7;
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
  for (int x = 0; x < 2; ++x) {
    const Mate& X = M[x];
    const Mate& Y = M[1 - x];
    for (int c = 0; c < X.nch; ++c) {
      const int st = X.ws->clusters[c].strand;
      int64_t plo, phi;
      partner_range(st, X.slo[c], X.shi[c], &plo, &phi);
      int ny = Y.U;  // mate y near this chain: no exact part, or one of its chains
      const size_t nb = near.size();
      for (int c2 = 0; c2 < Y.nch; ++c2)
        if (Y.ws->clusters[c2].strand == 1 - st && Y.shi[c2] >= plo && Y.slo[c2] <= phi) {
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
    if (id_strand >= 0 && same_pair(s1, a1, a2)) {
      upb = confirmed ? std::max(upb, up) : up;
      confirmed = true;
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
    } else {
      S2 = std::max(S2, up);
    }
  }
  if (id_strand < 0 || !confirmed || R1 <= floor || upb > R1 || S2 > R1) return po;
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
