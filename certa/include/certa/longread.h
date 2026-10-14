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
#include <functional>
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

// The targets of a pass over the read: 0 the reference, 1 its reverse
// complement (strand-1 pieces of the forward read); for the pass over the
// reversed read, 2 the reference reversed, 3 the reference complemented.
inline uint8_t jump_target(const IndexView& ix, int t, int64_t u) {
  const int64_t N = (int64_t)ix.ref_len;
  switch (t) {
    case 0: return ref_at(ix, u);
    case 1: return lr_comp(ref_at(ix, N - 1 - u));
    case 2: return ref_at(ix, N - 1 - u);
    default: return lr_comp(ref_at(ix, u));
  }
}

// A chain as a band of a pass: target t, diagonals [elo, ehi] (target
// position minus read position), read rows [r0, r1) of the pass.
struct JumpBand {
  int t;
  int64_t elo, ehi;
  int r0, r1;
};

// Per-strand tables of the parts (forward order of the pass).
struct PartTables {
  std::vector<int> first, last, owner, owner_off;
  void build(int L, const std::vector<ReadPart>& p) {
    first.assign(L, 0);
    last.assign(L, 0);
    owner.assign(L, -1);
    owner_off.assign(L, -1);
    for (size_t k = 0; k < p.size(); ++k) {
      first[p[k].off] = 1;
      last[p[k].end - 1] = 1;
      for (int i = p[k].off; i < p[k].end; ++i) {
        owner[i] = (int)k;
        owner_off[i] = p[k].off;
      }
    }
  }
};

// One step of a Lemma L-bounded piece (as in jump_bound) at base i: in[c]
// continued, a start with value `start`, then +1 and the charge of a part
// ending here that was covered from its first base and is not free.
inline void lemma_step(const PartTables& tb, const std::vector<char>* free, int i, long start, long* in) {
  const long NEGL = LONG_MIN / 4;
  long nx[2] = {NEGL, NEGL};
  const long cont = std::max(in[0], in[1]);
  if (tb.first[i]) nx[0] = cont;
  else if (tb.owner[i] >= 0) {
    nx[0] = in[0];
    nx[1] = in[1];
  } else {
    nx[1] = cont;
  }
  const int c = tb.owner_off[i] == i ? 0 : 1;
  nx[c] = std::max(nx[c], start);
  const bool charge = tb.last[i] && !(free && (*free)[tb.owner[i]]);
  for (int k = 0; k < 2; ++k) {
    if (nx[k] <= NEGL) continue;
    nx[k] += kMatch;
    if (k == 0 && charge) nx[k] -= kMismatch + kMatch;
  }
  in[0] = nx[0];
  in[1] = nx[1];
}

// One pass over a read (rd, in the pass's direction): the best decomposition
// into pieces of the given bands (exact: banded Gotoh on the band's rows, a
// piece starting from the best prefix value minus P, or 0 / -5 for the first
// piece) and of the given virtual types (Lemma L bounds). M (if not null) receives, for f =
// 0..L, the best decomposition whose pieces all end at or before f. Returns
// the best complete one (-5 if the last piece ends before L).
inline long jump_pass(const IndexView& ix, const uint8_t* rd, int L, const std::vector<JumpBand>& bands,
                      const PartTables* tb, const std::vector<PieceType>& virt, long P, std::vector<long>* M,
                      double* cells) {
  const long NEGL = LONG_MIN / 4;
  const int NEG = INT_MIN / 4;
  const int K = (int)bands.size();
  std::vector<std::vector<int>> H(K), I(K), Hn(K), In(K);
  for (int k = 0; k < K; ++k) {
    const int B = (int)(bands[k].ehi - bands[k].elo + 1);
    H[k].assign(B + 1, NEG);
    I[k].assign(B + 1, NEG);
    Hn[k].assign(B + 1, NEG);
    In[k].assign(B + 1, NEG);
  }
  // The targets of each band, read row f and band cell b at index f - r0 + b.
  std::vector<std::vector<uint8_t>> tg(K);
  for (int k = 0; k < K; ++k) {
    const JumpBand& bd = bands[k];
    const int64_t n = std::max<int64_t>(0, bd.r1 - bd.r0 + bd.ehi - bd.elo);
    tg[k].resize(n);
    for (int64_t x = 0; x < n; ++x) tg[k][x] = jump_target(ix, bd.t, bd.elo + bd.r0 + x);
  }
  const int T = (int)virt.size();
  std::vector<long> vin(2 * T, NEGL);
  long prefix = NEGL, end_at_L = NEGL;
  if (M) M->assign(L + 1, NEGL);
  for (int f = 0; f < L; ++f) {
    const long G = std::max(f > 0 ? -(long)kClip : 0L, prefix - P);
    const int g = (int)G;  // >= -5
    long F = NEGL;
    const uint8_t rb = rd[f];
    for (int k = 0; k < K; ++k) {
      const JumpBand& bd = bands[k];
      if (f < bd.r0 || f >= bd.r1) continue;
      const int B = (int)(bd.ehi - bd.elo + 1);
      const uint8_t* tr = tg[k].data() + (f - bd.r0);
      int d = NEG, left = NEG, rowbest = NEG;
      for (int b = 0; b < B; ++b) {
        const int m = std::max(H[k][b], g) + (lr_same(rb, tr[b]) ? kMatch : -kMismatch);
        const int ins = f > bd.r0 ? std::max(H[k][b + 1] - kGapOpen - kGapExt, I[k][b + 1] - kGapExt) : NEG;
        d = std::max(left - kGapOpen - kGapExt, d - kGapExt);
        const int h = std::max(m, std::max(ins, d));
        Hn[k][b] = h;
        In[k][b] = ins;
        left = h;
        rowbest = std::max(rowbest, m);
      }
      Hn[k][B] = In[k][B] = NEG;
      std::swap(H[k], Hn[k]);
      std::swap(I[k], In[k]);
      F = std::max(F, (long)rowbest);
      *cells += B;
    }
    for (int t = 0; t < T; ++t) {
      lemma_step(tb[virt[t].strand], virt[t].free, f, G, &vin[2 * t]);
      F = std::max(F, std::max(vin[2 * t], vin[2 * t + 1]));
    }
    prefix = std::max(prefix, F);
    if (M) (*M)[f + 1] = prefix;
    if (f == L - 1) end_at_L = F;
  }
  return std::max(prefix - kClip, end_at_L);
}

// The parts of the reversed read: part k of p, [o, e), is part n - 1 - k
// here, [L - e, L - o).
inline void mirror_parts(int L, const std::vector<ReadPart>& p, std::vector<ReadPart>& out) {
  out.clear();
  for (auto q = p.rbegin(); q != p.rend(); ++q) out.push_back({L - q->end, L - q->off, q->occ});
}

// pot: the best decomposition holding a piece of a chain (tables tb, the
// chain's parts free) over some [a, b) containing one of the chain's parts
// (a piece of the chain holds one exactly; a piece of its band holding none
// is unwitnessed), every other piece bounded by pre[a] (pieces ending at or
// before a) or suf[b] (pieces starting at or after b). The piece: Lemma L
// states as in lemma_step, u = a free part covered. With last_end: the
// largest b of such a decomposition above lb (-1: none).
inline long jump_pot(int L, const PartTables& tb, const std::vector<char>& free, const std::vector<long>& pre,
                     const std::vector<long>& suf, long P, long lb = LONG_MAX, int* last_end = nullptr) {
  const long NEGL = LONG_MIN / 4;
  long in[4] = {NEGL, NEGL, NEGL, NEGL}, best = NEGL;  // in[c * 2 + u]
  if (last_end) *last_end = -1;
  for (int i = 0; i < L; ++i) {
    long nx[4] = {NEGL, NEGL, NEGL, NEGL};
    for (int u = 0; u < 2; ++u) {
      const long cont = std::max(in[u], in[2 + u]);
      if (tb.first[i]) nx[u] = cont;
      else if (tb.owner[i] >= 0) {
        nx[u] = in[u];
        nx[2 + u] = in[2 + u];
      } else {
        nx[2 + u] = cont;
      }
    }
    const int c = tb.owner_off[i] == i ? 0 : 1;
    nx[c * 2] = std::max(nx[c * 2], std::max(i > 0 ? -(long)kClip : 0L, pre[i] - P));
    const bool is_free = tb.last[i] && free[tb.owner[i]];
    for (int k = 0; k < 4; ++k) {
      if (nx[k] <= NEGL) continue;
      nx[k] += kMatch;
      if (k < 2 && tb.last[i] && !is_free) nx[k] -= kMismatch + kMatch;
    }
    if (is_free) nx[1] = std::max(nx[1], nx[0]);  // the part, covered from its first base, is free: u = 1
    std::copy(nx, nx + 4, in);
    const long piece = std::max(in[1], in[3]);
    if (piece <= NEGL) continue;
    const long v = piece + std::max(i + 1 < L ? -(long)kClip : 0L, suf[i + 1] - P);
    best = std::max(best, v);
    if (last_end && v > lb) *last_end = i + 1;
  }
  return best;
}

// The read window of a chain (tables tb / tbr of the read and of the
// reversed read, free parts free / free_rev): the smallest [w0, w1) holding
// every [a, b) that contains one of the chain's parts and whose pot value
// (pre and suf as in jump_pot) exceeds lb. The first piece start comes from
// the same scan over the reversed read.
inline void jump_window(int L, const PartTables& tb, const PartTables& tbr, const std::vector<char>& free,
                        const std::vector<char>& free_rev, const std::vector<long>& pre, const std::vector<long>& suf,
                        long P, long lb, int* w0, int* w1) {
  std::vector<long> pre_r(L + 1), suf_r(L + 1);
  for (int x = 0; x <= L; ++x) {
    pre_r[x] = suf[L - x];
    suf_r[x] = pre[L - x];
  }
  int b = -1, y = -1;
  jump_pot(L, tb, free, pre, suf, P, lb, &b);
  jump_pot(L, tbr, free_rev, pre_r, suf_r, P, lb, &y);
  if (b < 0 || y < 0) {  // none above lb: the whole read (cannot happen for a chain with pot > lb)
    *w0 = 0;
    *w1 = L;
  } else {
    *w0 = L - y;
    *w1 = b;
  }
}

// Given lb <= OPT_J and lb > U_J:
//  - lbt: lb, raised each round to the best alignment of the evaluated
//    chains (a real one, so lbt <= OPT_J; every statement below holds with
//    the current lbt in place of lb);
//  - every jump alignment above lbt has gap total <= pad = L - lbt - 7 in
//    each piece (the piece's loss is part of the whole's, Lemma W);
//  - a piece holding an exact part of its strand lies in a chain of part
//    occurrences (grouped with the first pad, a coarsening of every later
//    one), on the chain's diagonals +- pad (Lemma C);
//  - chains are evaluated (banded DP on a read window) iteratively, best
//    single-piece bound first. pre: a pass over the read, the evaluated
//    chains exact on their windows, unevaluated chains' pieces by Lemma L
//    with their parts free, unwitnessed pieces by Lemma L; suf2: the same
//    pass over the reversed read; suf1: that pass without the unevaluated
//    chains;
//  - an alignment using unevaluated chains is bounded at its last piece of
//    one, C (the pieces after it are of evaluated chains or unwitnessed):
//    by pot(C) = the best decomposition holding a piece of C (Lemma L, C's
//    parts free), pre before it and suf1 after it. A chain with pot > lbt
//    is evaluated (the `batch` best per round; 0: at least 4, or a quarter
//    of them; the first round: the best single-piece bound), until none is
//    left: every alignment using an unevaluated chain then scores <= lbt;
//  - the window of a chain evaluated in a round holds every [a, b) whose
//    value with pre before and suf2 after exceeds lbt. By induction over the
//    rounds, every piece of an evaluated chain in an alignment above lbt
//    lies in its window (the pieces around it are bounded by pre and suf2,
//    those of chains evaluated before on their windows);
//  - the static path: the first round's pre and suf2 bound every
//    alignment (no chain is evaluated yet), so pot from them prunes at once,
//    and their windows hold: the chains with that pot above lbt are
//    evaluated together. Taken when the rounds would cost more than half of
//    it (batch 0), or at once (batch < 0);
//  - a final pass combines the evaluated chains' exact pieces, with
//    (upper) or without (lower) the unwitnessed pieces' Lemma L bounds.
// OPT_J lies in [max(lbt, lower), max(lbt, upper)]; equal ends prove it.
// Over the cell cap (status 2), OPT_J lies in [lbt, max(lbt, the last pre
// pass's value)]: that pass bounds every alignment above lbt.
// (A piece of a chain holds one of its parts exactly, so its read interval
// contains that part; a piece of the chain's band holding none is
// unwitnessed.)
// Strand-1 chains run on the forward read against the reverse-complemented
// reference (diagonal N - L - d), so every chain advances with the read.
inline LrCert certify_jump(const IndexView& ix, const uint8_t* fw, int L, StrandParts* so, long lb, long P,
                           double max_cells, int batch = 0) {
  LrCert c;
  std::vector<ReadPart> fp[2];
  forward_parts(L, so, fp);
  const long UJ = lemma_l_jump(L, fp, P);
  if (lb <= UJ) return c;
  const auto t0 = std::chrono::steady_clock::now();
  long lbt = lb;
  auto pad_of = [&]() { return std::max(0L, (long)L - lbt - 7); };
  const long pad0 = pad_of();
  const int64_t N = (int64_t)ix.ref_len;
  const int n1 = (int)so[1].parts.size();
  // The reversed read and its parts (a pass over it gives the suffix values).
  std::vector<uint8_t> rv(fw, fw + L);
  std::reverse(rv.begin(), rv.end());
  std::vector<ReadPart> rp[2];
  for (int s = 0; s < 2; ++s) mirror_parts(L, fp[s], rp[s]);
  PartTables tb[2], tbr[2];
  for (int s = 0; s < 2; ++s) {
    tb[s].build(L, fp[s]);
    tbr[s].build(L, rp[s]);
  }
  struct Chain {
    int strand;
    int64_t clo, chi;  // the occurrences' diagonals in the forward pass
    std::vector<char> free, free_rev;
    long ub;
    int w0, w1;  // the read window, once evaluated
  };
  // The bands of a chain in the pass over the read and over the reversed read.
  auto fwd_band = [&](const Chain& ch, long pad) {
    return JumpBand{ch.strand, ch.clo - pad, ch.chi + pad, ch.w0, ch.w1};
  };
  auto rev_band = [&](const Chain& ch, long pad) {
    return JumpBand{ch.strand + 2, N - L - ch.chi - pad, N - L - ch.clo + pad, L - ch.w1, L - ch.w0};
  };
  std::vector<Chain> chains;
  for (int s = 0; s < 2; ++s) {
    auto& occ = so[s].occ;
    std::sort(occ.begin(), occ.end());
    for (size_t i = 0, j; i < occ.size(); i = j) {
      j = i + 1;
      while (j < occ.size() && occ[j].first - occ[j - 1].first <= pad0) ++j;
      Chain ch;
      ch.strand = s;
      ch.clo = s == 0 ? occ[i].first : N - L - occ[j - 1].first;  // strand-1 diagonal d is N - L - d here
      ch.chi = s == 0 ? occ[j - 1].first : N - L - occ[i].first;
      ch.free.assign(fp[s].size(), 0);
      for (size_t k = i; k < j; ++k) ch.free[s == 0 ? occ[k].second : n1 - 1 - occ[k].second] = 1;
      ch.free_rev.assign(ch.free.rbegin(), ch.free.rend());
      ch.ub = lemma_l_span(L, fp[s], &ch.free);
      ch.w0 = ch.w1 = 0;
      chains.push_back(std::move(ch));
    }
  }
  c.chains = (int)chains.size();
  const size_t K = chains.size();
  std::vector<char> evaluated(K, 0);
  std::vector<JumpBand> bf, br;
  std::vector<long> pre, mrev, suf1(L + 1), suf2(L + 1), pre0, suf20, pot0(K, LONG_MIN / 4);
  std::vector<char> any[2], any_rev[2];
  auto band_cost = [&](long pad) {  // the cells of one pass over the evaluated chains
    double x = 0;
    for (size_t k = 0; k < K; ++k)
      if (evaluated[k])
        x += (double)(chains[k].w1 - chains[k].w0) * (double)(chains[k].chi - chains[k].clo + 1 + 2 * pad);
    return x;
  };
  auto static_cost = [&](long pad) {  // the final passes over the static set, whole-read windows
    double x = 0;
    for (size_t k = 0; k < K; ++k)
      if (pot0[k] > lbt) x += 2.0 * L * (double)(chains[k].chi - chains[k].clo + 1 + 2 * pad);
    return x;
  };
  long ub_all = L;  // every jump alignment scores <= L
  auto capped = [&]() {
    c.status = 2;
    c.opt = lbt;
    c.hi = std::max(lbt, ub_all);
    c.searched = (int)std::count(evaluated.begin(), evaluated.end(), 1);
    c.us = std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - t0).count();
    return c;
  };
  for (int iter = 0;; ++iter) {
    const long pad = pad_of();
    if (c.cells + 6.0 * band_cost(pad) > max_cells || iter > 64) return capped();
    if (iter > 0 && (batch < 0 || (batch == 0 && c.cells + 6.0 * band_cost(pad) > 0.5 * static_cost(pad)))) {
      // The static path.
      for (size_t k = 0; k < K; ++k) {
        Chain& ch = chains[k];
        evaluated[k] = pot0[k] > lbt;
        if (evaluated[k])
          jump_window(L, tb[ch.strand], tbr[ch.strand], ch.free, ch.free_rev, pre0, suf20, P, lbt, &ch.w0, &ch.w1);
      }
      if (c.cells + 2.0 * band_cost(pad) > max_cells) return capped();
      break;
    }
    bf.clear();
    br.clear();
    for (int s = 0; s < 2; ++s) {
      any[s].assign(fp[s].size(), 0);
      any_rev[s].assign(fp[s].size(), 0);
    }
    for (size_t k = 0; k < K; ++k) {
      const Chain& ch = chains[k];
      if (evaluated[k]) {
        bf.push_back(fwd_band(ch, pad));
        br.push_back(rev_band(ch, pad));
      } else {
        for (size_t x = 0; x < ch.free.size(); ++x) {
          any[ch.strand][x] |= ch.free[x];
          any_rev[ch.strand][x] |= ch.free_rev[x];
        }
      }
    }
    if (!bf.empty()) lbt = std::max(lbt, jump_pass(ix, fw, L, bf, tb, {}, P, nullptr, &c.cells));
    const std::vector<PieceType> vf = {{0, nullptr}, {1, nullptr}, {0, &any[0]}, {1, &any[1]}};
    ub_all = jump_pass(ix, fw, L, bf, tb, vf, P, &pre, &c.cells);
    jump_pass(ix, rv.data(), L, br, tbr, {{0, nullptr}, {1, nullptr}}, P, &mrev, &c.cells);
    for (int b = 0; b <= L; ++b) suf1[b] = mrev[L - b];
    std::vector<std::pair<long, int>> add;
    for (size_t k = 0; k < K; ++k) {
      if (evaluated[k]) continue;
      const long pot = jump_pot(L, tb[chains[k].strand], chains[k].free, pre, suf1, P);
      if (pot > lbt) add.push_back({iter == 0 ? chains[k].ub : pot, (int)k});
    }
    if (add.empty()) break;
    std::sort(add.begin(), add.end(), std::greater<std::pair<long, int>>());
    const size_t n_add = iter == 0 ? 1 : batch > 0 ? (size_t)batch : std::max<size_t>(4, add.size() / 4);
    // The windows of the chains evaluated now: suf2 bounds unevaluated
    // chains' pieces too.
    const std::vector<PieceType> vr = {{0, nullptr}, {1, nullptr}, {0, &any_rev[0]}, {1, &any_rev[1]}};
    jump_pass(ix, rv.data(), L, br, tbr, vr, P, &mrev, &c.cells);
    for (int b = 0; b <= L; ++b) suf2[b] = mrev[L - b];
    if (iter == 0) {  // the first round's values, for the static path
      pre0 = pre;
      suf20 = suf2;
      for (size_t k = 0; k < K; ++k) pot0[k] = jump_pot(L, tb[chains[k].strand], chains[k].free, pre0, suf20, P);
    }
    for (size_t x = 0; x < add.size() && x < n_add; ++x) {
      Chain& ch = chains[add[x].second];
      jump_window(L, tb[ch.strand], tbr[ch.strand], ch.free, ch.free_rev, pre, suf2, P, lbt, &ch.w0, &ch.w1);
      evaluated[add[x].second] = 1;
    }
  }
  const long pad = pad_of();
  bf.clear();
  for (size_t k = 0; k < K; ++k)
    if (evaluated[k]) bf.push_back(fwd_band(chains[k], pad));
  c.searched = (int)bf.size();
  const long upper = jump_pass(ix, fw, L, bf, tb, {{0, nullptr}, {1, nullptr}}, P, nullptr, &c.cells);
  const long lower = jump_pass(ix, fw, L, bf, tb, {}, P, nullptr, &c.cells);
  const long lo = std::max(lbt, lower), hi = std::max(lo, upper);
  c.status = hi == lo ? 1 : 3;
  c.opt = lo;
  c.hi = hi;
  c.us = std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - t0).count();
  return c;
}

}  // namespace certa
