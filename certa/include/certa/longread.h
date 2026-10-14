// Long reads (EXACT_ALIGNMENT_PROPOSAL.md Aim 6; tools/lr_ceiling,
// tests/test_longread.cpp): exact parts, Lemma L bounds, and the certificates
// of the local optimum and of the split (jump) optimum. BWA-MEM's scores:
// +1, -4, gaps -(6 + k), -5 per clipped read end; an inexact part inside an
// alignment costs >= 5 (Lemma L). The index must have s = 1.
#pragma once
#include <algorithm>
#include <chrono>
#include <climits>
#include <cstdint>
#include <vector>

#include "certa/core.h"

namespace certa {

struct ReadPart {
  int off, end;
  uint32_t occ;  // exact occurrences (0: absent)
};

inline bool lr_same(uint8_t a, uint8_t b) { return a < 4 && a == b; }
inline uint8_t lr_comp(uint8_t c) { return c < 4 ? (uint8_t)(3 - c) : (uint8_t)4; }

// ---- Lemma L bounds -------------------------------------------------------

// Lemma L, t = 1, with the parts flagged in `free` uncharged (null: none):
// every local alignment with no exact occurrence of a charged part scores
// <= max over spans [x, y) of (y - x) - 5 (charged parts inside) - 5
// (clipped ends). Span candidates x in {0, off + 1}, y in {L, end - 1}
// (Lemma S). O(P).
inline long lemma_l_span(int L, const std::vector<ReadPart>& p, const std::vector<char>* free = nullptr) {
  const int P = (int)p.size();
  const long c5 = kMismatch + kMatch;
  long G = 0, best = LONG_MIN, C = 0;  // C: charged parts among 0..j-1
  for (int j = 0; j <= P; ++j) {
    const long y = j == P ? L : p[j].end - 1;
    best = std::max(best, y - c5 * C - (y < L ? kClip : 0) + G);
    if (j < P) {
      C += !(free && (*free)[j]);
      G = std::max(G, -(long)(p[j].off + 1) + c5 * C - kClip);
    }
  }
  return best;
}

// Jump alignments (split alignments of a read across structural variants):
// pieces A_1 .. A_k, each a gapped local alignment of a read interval
// [a_j, b_j) on its own strand and locus, in read order (b_j <= a_{j+1});
// score sum sigma(A_j) - P (k - 1) - 5 [a_1 > 0] - 5 [b_k < L], read bases
// outside the pieces scoring 0. k = 1 is the local objective.
//
// A piece type is a strand and a set of free parts (null: none); a piece of
// that type is bounded by Lemma L with those parts uncharged: (b - a) - 5
// (charged parts of its strand inside [a, b)). jump_bound returns the best
// decomposition into pieces of the given types (with req >= 0: using type
// req at least once), exactly, in one pass over the read.
//
// part[s]: the parts of strand s in forward read coordinates (a strand-1
// part [o, e) of the reverse complement is [L - e, L - o) here), sorted;
// `free` is indexed in that order. States: inside a piece of type t, with
// c = 0 when the part holding the current base has been covered from its
// first base (charged at its last base unless free), and u = whether type
// req has been used; or between pieces.
struct PieceType {
  int strand;
  const std::vector<char>* free;
};
inline long jump_bound(int L, const std::vector<ReadPart>* part, const std::vector<PieceType>& types, long P,
                       int req = -1) {
  const long NEG = LONG_MIN / 4, c5 = kMismatch + kMatch;
  std::vector<int> first[2], last[2], owner[2], owner_off[2];
  for (int s = 0; s < 2; ++s) {
    first[s].assign(L, 0);
    last[s].assign(L, 0);
    owner[s].assign(L, -1);
    owner_off[s].assign(L, -1);
    for (size_t k = 0; k < part[s].size(); ++k) {
      const ReadPart& q = part[s][k];
      first[s][q.off] = 1;
      last[s][q.end - 1] = 1;
      for (int i = q.off; i < q.end; ++i) {
        owner[s][i] = (int)k;
        owner_off[s][i] = q.off;
      }
    }
  }
  const int T = (int)types.size();
  auto at = [](int t, int c, int u) { return (t * 2 + c) * 2 + u; };
  std::vector<long> in(T * 4, NEG), nx(T * 4);
  long mid[2] = {NEG, NEG};
  for (int i = 0; i < L; ++i) {
    std::fill(nx.begin(), nx.end(), NEG);
    const long pre = i > 0 ? -(long)kClip : 0L;  // the first piece
    for (int t = 0; t < T; ++t) {
      const int s = types[t].strand, ut = t == req;
      for (int u = 0; u < 2; ++u) {
        // Continue a piece of type t through base i.
        const long cont = std::max(in[at(t, 0, u)], in[at(t, 1, u)]);
        if (first[s][i]) {
          nx[at(t, 0, u)] = std::max(nx[at(t, 0, u)], cont);
        } else if (owner[s][i] >= 0) {
          nx[at(t, 0, u)] = std::max(nx[at(t, 0, u)], in[at(t, 0, u)]);
          nx[at(t, 1, u)] = std::max(nx[at(t, 1, u)], in[at(t, 1, u)]);
        } else {
          nx[at(t, 1, u)] = std::max(nx[at(t, 1, u)], cont);
        }
      }
      // Start a piece of type t at base i (inside a part past its first base:
      // that part is not charged).
      const int c = owner_off[s][i] == i ? 0 : 1;
      nx[at(t, c, ut)] = std::max(nx[at(t, c, ut)], pre);
      nx[at(t, c, ut)] = std::max(nx[at(t, c, ut)], mid[0] - P);
      nx[at(t, c, 1)] = std::max(nx[at(t, c, 1)], mid[1] - P);
      const bool charge = last[s][i] && !(types[t].free && (*types[t].free)[owner[s][i]]);
      for (int cc = 0; cc < 2; ++cc)
        for (int u = 0; u < 2; ++u) {
          long& v = nx[at(t, cc, u)];
          if (v <= NEG) continue;
          v += kMatch;
          if (cc == 0 && charge) v -= c5;
        }
    }
    in.swap(nx);
    for (int t = 0; t < T; ++t)
      for (int cc = 0; cc < 2; ++cc)
        for (int u = 0; u < 2; ++u) mid[u] = std::max(mid[u], in[at(t, cc, u)]);  // a piece may end here
  }
  long best = NEG;
  for (int u = 0; u < 2; ++u) {
    if (req >= 0 && u == 0) continue;
    best = std::max(best, mid[u] - kClip);  // the last piece ends before L
    for (int t = 0; t < T; ++t)
      for (int cc = 0; cc < 2; ++cc) best = std::max(best, in[at(t, cc, u)]);  // or at L
  }
  return best;
}

// Lemma L-J: a jump alignment in which no piece holds an exact occurrence of
// a part of its strand scores <= U_J, the best decomposition into pieces each
// bounded by Lemma L (each piece is a local alignment of its read interval).
inline long lemma_l_jump(int L, const std::vector<ReadPart>* part, long P) {
  return jump_bound(L, part, {{0, nullptr}, {1, nullptr}}, P);
}

// ---- Parts of a long read (Theorem G at tau, Lemma R) ----------------------

struct PartitionParams {
  int tau = 16, cmax = 4096, window = 40, step = 8;
  bool keep_occ = false;  // keep every occurrence (for the certificates)
};

struct LrScratch {
  std::vector<uint64_t> lo;
  std::vector<uint32_t> cnt;
  std::vector<int> nz;
  std::vector<int64_t> elig;
  std::vector<int> ext;
};

struct StrandParts {
  std::vector<ReadPart> parts;
  long skipped = 0;      // read bases where no part could start
  long occ_sum = 0;      // enumerated occurrences
  int absent = 0;
  int longest_free = 0;  // longest stretch covered by no part
  // keep_occ: every exact occurrence of a rare part, as (diagonal = the
  // reference position of read base 0, part index), on this strand.
  std::vector<std::pair<int64_t, int>> occ;
};

// The greedy partition of one strand of a read (rd: that strand's sequence):
// at each position s the shortest prefix occurring <= tau times, found by
// anchor extension; parts containing a count-0 q-mer are absent. A stretch
// where no such part starts is skipped (Lemma L needs disjoint parts, not a
// tiling).
inline void partition_strand(const IndexView& ix, const uint8_t* rd, int L, const PartitionParams& o, LrScratch& sc,
                             StrandParts& out) {
  const int q = ix.q, kBig = INT_MAX / 2;
  out = StrandParts();
  if (L < q) {
    out.skipped = L;
    out.longest_free = L;
    return;
  }
  const int nq = L - q + 1;
  sc.lo.resize(nq);
  sc.cnt.resize(nq);
  sc.nz.resize(nq + 1);
  for (int i = 0; i < nq; ++i) {
    uint64_t key, lo = 0, hi = 0;
    if (kmer_key(rd + i, q, &key)) lookup(ix, key, &lo, &hi);
    sc.lo[i] = lo;
    sc.cnt[i] = (uint32_t)std::min<uint64_t>(hi - lo, UINT32_MAX);
  }
  sc.nz[nq] = nq;
  for (int i = nq - 1; i >= 0; --i) sc.nz[i] = sc.cnt[i] == 0 ? i : sc.nz[i + 1];
  int s = 0;
  while (s < nq) {
    const int babs = sc.nz[s] < nq ? sc.nz[s] + q : kBig;
    int c = -1;
    for (int i = s; i < nq && i <= s + o.window; ++i)
      if (sc.cnt[i] > 0 && (int)sc.cnt[i] <= o.cmax && (c < 0 || sc.cnt[i] < sc.cnt[c])) c = i;
    int banc = kBig;
    sc.elig.clear();
    sc.ext.clear();
    if (c >= 0) {
      for (uint64_t e = sc.lo[c]; e < sc.lo[c] + sc.cnt[c]; ++e) {
        const int64_t h = ix.pos[e];
        bool ok = true;
        for (int k = c - 1; k >= s && ok; --k) ok = lr_same(rd[k], ref_at(ix, h - (c - k)));
        if (ok) sc.elig.push_back(h);
      }
      if ((int)sc.elig.size() <= o.tau) {
        banc = c + q;
      } else {
        for (int64_t h : sc.elig) {
          int R = 0;
          while (c + q + R < L && lr_same(rd[c + q + R], ref_at(ix, h + q + R))) ++R;
          sc.ext.push_back(R);
        }
        std::vector<int> e2(sc.ext);
        std::nth_element(e2.begin(), e2.begin() + o.tau, e2.end(), std::greater<int>());
        const int Rt = e2[o.tau];  // the (tau+1)-th longest extension
        if (c + q + Rt < L) banc = c + q + Rt + 1;
      }
    }
    const int b = std::min(babs, banc);
    if (b > L) {  // no part of <= tau occurrences starts at s
      const int nx = std::min(s + o.step, L);
      out.skipped += nx - s;
      s = nx;
      continue;
    }
    const bool rare = b == banc && b < babs;
    uint32_t occ = 0;
    if (rare) {
      const int part = (int)out.parts.size();
      for (size_t k = 0; k < sc.elig.size(); ++k)
        if (sc.ext.empty() || c + q + sc.ext[k] >= b) {
          ++occ;
          if (o.keep_occ) out.occ.push_back({sc.elig[k] - c, part});
        }
    }
    out.parts.push_back({s, b, occ});
    out.occ_sum += occ;
    out.absent += !rare;
    s = b;
  }
  if (s < L) out.skipped += L - s;  // a tail shorter than q
  int prev = 0;
  for (const ReadPart& p : out.parts) {
    out.longest_free = std::max(out.longest_free, p.off - prev);
    prev = p.end;
  }
  out.longest_free = std::max(out.longest_free, L - prev);
}

// Strand 1's parts (reverse-complement coordinates) in forward coordinates.
inline void forward_parts(int L, const StrandParts* so, std::vector<ReadPart>* out) {
  out[0] = so[0].parts;
  out[1].clear();
  for (auto q = so[1].parts.rbegin(); q != so[1].parts.rend(); ++q) out[1].push_back({L - q->end, L - q->off, q->occ});
}

// ---- Banded local DP ---------------------------------------------------------

struct BandDp {
  std::vector<int> H, I, Hn, In;
};

// The best local score (-5 per clipped end) of rd over alignments whose cells
// all lie on diagonals [dlo, dhi] (reference position minus read position).
// Gotoh, one row of the band at a time.
inline long band_best(const IndexView& ix, const uint8_t* rd, int L, int64_t dlo, int64_t dhi, BandDp& dp,
                      double* cells) {
  const int B = (int)(dhi - dlo + 1);
  const int NEG = INT_MIN / 4;
  dp.H.assign(B + 1, NEG);
  dp.I.assign(B + 1, NEG);
  dp.Hn.assign(B + 1, NEG);
  dp.In.assign(B + 1, NEG);
  int best = NEG;
  for (int i = 1; i <= L; ++i) {
    const uint8_t c = rd[i - 1];
    const int start = i == 1 ? 0 : -kClip;
    const int end = i < L ? kClip : 0;
    int d = NEG, left = NEG;  // deletion state, and H at (i, b - 1)
    const int64_t j0 = i - 1 + dlo;
    for (int b = 0; b < B; ++b) {
      const int m = std::max(dp.H[b], start) + (lr_same(c, ref_at(ix, j0 + b)) ? kMatch : -kMismatch);
      const int ins = i > 1 ? std::max(dp.H[b + 1] - kGapOpen - kGapExt, dp.I[b + 1] - kGapExt) : NEG;
      d = std::max(left - kGapOpen - kGapExt, d - kGapExt);
      const int h = std::max(m, std::max(ins, d));
      dp.Hn[b] = h;
      dp.In[b] = ins;
      left = h;
      best = std::max(best, m - end);
    }
    dp.Hn[B] = dp.In[B] = NEG;
    std::swap(dp.H, dp.Hn);
    std::swap(dp.I, dp.In);
  }
  *cells += (double)L * B;
  return best;
}

// ---- The certificate of the local optimum (Corollary L-long) ---------------

// Given a fitting read (lb = a known alignment's score > U1): every alignment
// above lb has gap total <= pad = L - lb - 7 (Lemma W) and, having an exact
// part (it beats U1), lies in a chain of part occurrences on its strand
// (split where diagonals are more than pad apart), within its diagonals +-
// pad (Lemma C), scoring at most the chain's Lemma L bound with its parts
// free. So OPT = max(lb, the band maxima of the chains whose bound exceeds
// the best so far). rd[s]: strand s's sequence; so[s] with keep_occ.
struct LrCert {
  int status = 0;  // 0: not attempted (no fit); 1: OPT proven; 2: over the cell cap; 3: an interval [lo, hi]
  long opt = 0;    // status 1: OPT; status 3: the lower end
  long hi = 0;     // status 3: the upper end
  int chains = 0, searched = 0;
  double cells = 0, us = 0;
};
inline LrCert certify_single(const IndexView& ix, const uint8_t* const* rd, int L, StrandParts* so, long lb, long U,
                             double max_cells, BandDp& dp) {
  LrCert c;
  if (lb <= U) return c;
  const auto t0 = std::chrono::steady_clock::now();
  const long pad = std::max(0L, (long)L - lb - 7);
  long best = lb;
  std::vector<char> free;
  for (int st = 0; st < 2; ++st) {
    auto& occ = so[st].occ;
    std::sort(occ.begin(), occ.end());
    const std::vector<ReadPart>& P = so[st].parts;
    free.assign(P.size(), 0);
    for (size_t i = 0, j; i < occ.size(); i = j) {
      j = i + 1;
      while (j < occ.size() && occ[j].first - occ[j - 1].first <= pad) ++j;
      ++c.chains;
      for (size_t k = i; k < j; ++k) free[occ[k].second] = 1;
      const long ub = lemma_l_span(L, P, &free);
      for (size_t k = i; k < j; ++k) free[occ[k].second] = 0;
      if (ub <= best) continue;
      const int64_t dlo = occ[i].first - pad, dhi = occ[j - 1].first + pad;
      if (c.cells + (double)L * (double)(dhi - dlo + 1) > max_cells) {
        c.status = 2;
        return c;
      }
      ++c.searched;
      best = std::max(best, band_best(ix, rd[st], L, dlo, dhi, dp, &c.cells));
    }
  }
  c.status = 1;
  c.opt = best;
  c.us = std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - t0).count();
  return c;
}

// ---- The certificate of the split (jump) optimum ----------------------------

// Given lb, the score of a known jump alignment, and lb > U_J:
//  - every jump alignment above lb has gap total <= pad = L - lb - 7 in each
//    piece (the piece's loss is part of the whole's, Lemma W);
//  - a piece holding an exact part of its strand lies in a chain of part
//    occurrences (pad), on the chain's diagonals +- pad (Lemma C);
//  - pot(C), the best decomposition using a piece of chain C (bounded by
//    Lemma L with C's parts free, every other piece with every occurring
//    part free), bounds all alignments using C: chains with pot <= lb are
//    dropped;
//  - one pass over the read then combines, for the kept chains, the exact
//    best piece ending at each read position (banded DP, a piece starting
//    from the best prefix value minus P), and, in an upper layer only, the
//    pieces holding no exact part (Lemma L).
// OPT_J lies in [max(lb, lower), max(lb, upper)]; equal ends prove it.
// Strand-1 chains run on the forward read against the reverse-complemented
// reference, so every chain advances with the read.
inline LrCert certify_jump(const IndexView& ix, const uint8_t* fw, int L, StrandParts* so, long lb, long P,
                           double max_cells) {
  LrCert c;
  std::vector<ReadPart> fp[2];
  forward_parts(L, so, fp);
  const long UJ = lemma_l_jump(L, fp, P);
  if (lb <= UJ) return c;
  const auto t0 = std::chrono::steady_clock::now();
  const long NEGL = LONG_MIN / 4;
  const long pad = std::max(0L, (long)L - lb - 7);
  const int64_t N = (int64_t)ix.ref_len;
  const int n1 = (int)so[1].parts.size();
  // Chains, with their free parts in forward order.
  struct Chain {
    int strand;
    int64_t elo, ehi;  // forward-coordinate diagonals (target position minus read position)
    std::vector<char> free;
  };
  std::vector<Chain> chains;
  std::vector<char> any[2];
  for (int s = 0; s < 2; ++s) any[s].assign(fp[s].size(), 0);
  for (int s = 0; s < 2; ++s) {
    auto& occ = so[s].occ;
    std::sort(occ.begin(), occ.end());
    for (size_t i = 0, j; i < occ.size(); i = j) {
      j = i + 1;
      while (j < occ.size() && occ[j].first - occ[j - 1].first <= pad) ++j;
      Chain ch;
      ch.strand = s;
      const int64_t dlo = occ[i].first - pad, dhi = occ[j - 1].first + pad;
      if (s == 0) {
        ch.elo = dlo;
        ch.ehi = dhi;
      } else {  // strand-1 diagonal d (ref - rc read position) is N - L - d here
        ch.elo = N - L - dhi;
        ch.ehi = N - L - dlo;
      }
      ch.free.assign(fp[s].size(), 0);
      for (size_t k = i; k < j; ++k) {
        const int idx = s == 0 ? occ[k].second : n1 - 1 - occ[k].second;
        ch.free[idx] = 1;
        any[s][idx] = 1;
      }
      chains.push_back(std::move(ch));
    }
  }
  c.chains = (int)chains.size();
  // pot(C): drop the chains no alignment above lb can use.
  std::vector<int> keep;
  for (size_t k = 0; k < chains.size(); ++k) {
    const std::vector<PieceType> types = {{chains[k].strand, &chains[k].free}, {0, &any[0]}, {1, &any[1]}};
    if (jump_bound(L, fp, types, P, 0) > lb) keep.push_back((int)k);
  }
  c.searched = (int)keep.size();
  double need = 0;
  for (int k : keep) need += 2.0 * L * (double)(chains[k].ehi - chains[k].elo + 1);
  if (need > max_cells) {
    c.status = 2;
    return c;
  }
  // The pass over the read: layer 0 upper (real pieces of kept chains and
  // unwitnessed pieces), layer 1 lower (real pieces only).
  const int NEG = INT_MIN / 4;
  struct Band {
    std::vector<int> H, I, Hn, In;
  };
  const int K = (int)keep.size();
  std::vector<Band> band(2 * K);
  for (int k = 0; k < K; ++k)
    for (int l = 0; l < 2; ++l) {
      const int B = (int)(chains[keep[k]].ehi - chains[keep[k]].elo + 1);
      Band& bd = band[2 * k + l];
      bd.H.assign(B + 1, NEG);
      bd.I.assign(B + 1, NEG);
      bd.Hn.assign(B + 1, NEG);
      bd.In.assign(B + 1, NEG);
    }
  // Unwitnessed pieces of the upper layer: Lemma L states per strand.
  std::vector<int> first[2], last[2], owner_off[2];
  for (int s = 0; s < 2; ++s) {
    first[s].assign(L, 0);
    last[s].assign(L, 0);
    owner_off[s].assign(L, -1);
    for (const ReadPart& q : fp[s]) {
      first[s][q.off] = 1;
      last[s][q.end - 1] = 1;
      for (int i = q.off; i < q.end; ++i) owner_off[s][i] = q.off;
    }
  }
  long vin[2][2] = {{NEGL, NEGL}, {NEGL, NEGL}};
  long M[2] = {NEGL, NEGL};  // best decomposition ending at or before the current base, per layer
  long Fend[2] = {NEGL, NEGL};  // best decomposition ending exactly at L
  for (int f = 0; f < L; ++f) {
    long G[2];
    for (int l = 0; l < 2; ++l) G[l] = std::max(f > 0 ? -(long)kClip : 0L, M[l] - P);
    long F[2] = {NEGL, NEGL};  // pieces ending after base f
    const uint8_t rb = fw[f];
    for (int k = 0; k < K; ++k) {
      const Chain& ch = chains[keep[k]];
      const int B = (int)(ch.ehi - ch.elo + 1);
      const int64_t u0 = f + ch.elo;
      for (int l = 0; l < 2; ++l) {
        Band& bd = band[2 * k + l];
        const int start = G[l] > NEG ? (int)std::max<long>(G[l], NEG) : NEG;
        int d = NEG, left = NEG, rowbest = NEG;
        for (int b = 0; b < B; ++b) {
          const int64_t u = u0 + b;
          const uint8_t t = ch.strand == 0 ? ref_at(ix, u) : lr_comp(ref_at(ix, N - 1 - u));
          const int m = std::max(bd.H[b], start) + (lr_same(rb, t) ? kMatch : -kMismatch);
          const int ins = f > 0 ? std::max(bd.H[b + 1] - kGapOpen - kGapExt, bd.I[b + 1] - kGapExt) : NEG;
          d = std::max(left - kGapOpen - kGapExt, d - kGapExt);
          const int h = std::max(m, std::max(ins, d));
          bd.Hn[b] = h;
          bd.In[b] = ins;
          left = h;
          rowbest = std::max(rowbest, m);
        }
        bd.Hn[B] = bd.In[B] = NEG;
        std::swap(bd.H, bd.Hn);
        std::swap(bd.I, bd.In);
        F[l] = std::max(F[l], (long)rowbest);
      }
      c.cells += 2.0 * B;
    }
    // Unwitnessed pieces (upper layer), as in jump_bound.
    long nv[2][2] = {{NEGL, NEGL}, {NEGL, NEGL}};
    for (int s = 0; s < 2; ++s) {
      const long cont = std::max(vin[s][0], vin[s][1]);
      if (first[s][f]) nv[s][0] = std::max(nv[s][0], cont);
      else if (owner_off[s][f] >= 0) {
        nv[s][0] = std::max(nv[s][0], vin[s][0]);
        nv[s][1] = std::max(nv[s][1], vin[s][1]);
      } else {
        nv[s][1] = std::max(nv[s][1], cont);
      }
      const int cc = owner_off[s][f] == f ? 0 : 1;
      nv[s][cc] = std::max(nv[s][cc], G[0]);
      for (int k = 0; k < 2; ++k) {
        if (nv[s][k] <= NEGL) continue;
        nv[s][k] += kMatch;
        if (k == 0 && last[s][f]) nv[s][k] -= kMismatch + kMatch;
      }
    }
    for (int s = 0; s < 2; ++s)
      for (int k = 0; k < 2; ++k) {
        vin[s][k] = nv[s][k];
        F[0] = std::max(F[0], vin[s][k]);
      }
    for (int l = 0; l < 2; ++l) {
      M[l] = std::max(M[l], F[l]);
      if (f == L - 1) Fend[l] = F[l];
    }
  }
  long ends[2];
  for (int l = 0; l < 2; ++l) ends[l] = std::max(M[l] - kClip, Fend[l]);
  const long lo = std::max(lb, ends[1]), hi = std::max(lo, ends[0]);
  c.status = hi == lo ? 1 : 3;
  c.opt = lo;
  c.hi = hi;
  c.us = std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - t0).count();
  return c;
}

}  // namespace certa
