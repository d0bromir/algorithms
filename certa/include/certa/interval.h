// CERTA interval output (host only). EXACT_ALIGNMENT_PROPOSAL.md Aim 3.
//
// With `certa map --intervals` no read goes to a heuristic mapper. Every
// read leaves with:
//  - an alignment, or none (unmapped);
//  - [XL, XU] containing OPT, the best local score (the objective of
//    Theorem L) over the whole reference, both strands;
//  - XQ, a certified lower bound of BWA-MEM's MAPQ formula. The formula
//    rises with the score and falls with SUB, the best score at a locus that
//    does not overlap the reported one; evaluated at a proven SUB_hi >= SUB,
//    it bounds the formula at the true SUB from below.
//
// Bounds by tier:
//  - S0/S1 certify end-to-end edit distance, not the local score. Every
//    other locus has distance >= d* (the exact second best, or R + 1), so by
//    Lemma E OPT <= max(score, edit_bound(L, d*)); SUB_hi from process_read;
//  - SR: two exact full-length copies, OPT = SUB = L;
//  - SL: OPT = score (Theorem L); SUB_hi from certify_local;
//  - PR mates: the pair certificate; a mate's single-end bounds come from its
//    single-end result if that was certified, else XU = L;
//  - every other read: tier R, [lo, hi] and SUB in [sub_lo, sub_hi].
//
// Pairs with an uncertified mate are placed the way BWA-MEM places pairs,
// from what the tiers found: each mate's single-end alignment (and its proven
// copy, Theorem R), and the best alignment of one mate in the proper-pair
// window of the other's (the window search of Theorem P, exhaustive in the
// window). The best pair is taken if it beats the sum of the mates' best
// single-end scores minus BWA-MEM's unpaired penalty; MAPQ as in mem_sam_pe,
// except that the pair is trusted no more than its more confident mate (on
// simulated pairs this cut the misplaced records with MAPQ >= 20 from 6,203
// to 539 of 1.91 M). This placement is a choice among proven alignments, not
// a certificate: the record says so (XT:Z:MR, XQ 0, and [XL, XU] stays the
// single-end interval).
//
// XQ bounds the formula, not the error rate: it makes the formula's inputs
// exact. In duplications a one-mismatch margin gives MAPQ 30 by the formula,
// while one read error or SNV flips the optimum far more often than 1/1000.
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>

#include "certa/core.h"
#include "certa/pair.h"
#include "certa/repeat.h"

namespace certa {

constexpr int kPenUnpaired = 17;  // BWA-MEM -U
constexpr int kMinScore = 30;     // BWA-MEM -T: alignments below are not reported

// mem_approx_mapq_se (mapQ_coef_len = 50) on a read of length L, before
// clamping and the sub_n correction.
inline int mapq_se_raw(int score, int sub, int L) {
  const double l = L, identity = 1.0 - (l * kMatch - score) / (kMatch + kMismatch) / l;
  double t = l < 50 ? 1.0 : std::log(50.0) / std::log(l);
  t *= identity * identity;
  return static_cast<int>(6.02 * (score - sub) / kMatch * t * t + 0.499);
}

// The same with BWA-MEM's floor of 19 for sub (min_seed_len x match), 0 when
// sub >= score, clamped to [0, 60]. Non-increasing in sub.
inline int mapq_se(int score, int sub, int L) {
  if (sub < 19) sub = 19;
  if (score <= 0 || sub >= score) return 0;
  const int q = mapq_se_raw(score, sub, L);
  return q < 0 ? 0 : (q > 60 ? 60 : q);
}

// BWA-MEM's raw_mapq of a pair score margin, clamped to [0, 60].
inline int mapq_margin(int d) {
  const int q = static_cast<int>(6.02 * d / kMatch + 0.499);
  return q < 0 ? 0 : (q > 60 ? 60 : q);
}

inline int64_t ref_span_of(const Result& r) {
  int64_t span = 0;
  for (int i = 0; i < r.n_cigar; ++i)
    if ((r.cigar[i] & 0xF) == kOpM || (r.cigar[i] & 0xF) == kOpD) span += r.cigar[i] >> 4;
  return span;
}

// Same locus: same strand and overlapping reference intervals.
inline bool same_locus(const Result& a, const Result& b) {
  return a.strand == b.strand && a.ref_pos < b.ref_pos + ref_span_of(b) && b.ref_pos < a.ref_pos + ref_span_of(a);
}

// A result holding only an alignment (affine_align fills ref_pos, CIGAR, nm
// and score; everything else is reset).
inline Result alignment_only(const Result& a, int strand, uint8_t tier) {
  Result o{};
  o.ref_pos = a.ref_pos;
  o.n_cigar = a.n_cigar;
  for (int i = 0; i < a.n_cigar; ++i) o.cigar[i] = a.cigar[i];
  o.nm = a.nm;
  o.score = a.score;
  o.strand = (uint8_t)strand;
  o.tier = tier;
  o.radius = -1;
  o.d1 = -1;
  o.d2 = -1;
  o.d2x = -1;
  o.sub_score = kNoSub;
  o.sub_hi = kNoSub;
  return o;
}

// The witness of a record's bounds (Theorem V, checked by certa check): the
// disjoint read intervals ("probes") whose exact occurrences the bound was
// proven from, per strand (strand 1 in reverse-complement coordinates), and
// the support t: alignments with fewer than t exact probes are bounded by
// Lemma L alone. t = 0: no witness (the record claims nothing that needs one).
struct Witness {
  int t = 0;
  int n[2] = {0, 0};
  int16_t off[2][PMAX], end[2][PMAX];
};

// The parts a certified result enumerated (tiers S0, S1, SL; for PR, the
// parts of the pair certificate, see pair.h).
inline Witness witness_of(const Result& r, int L) {
  Witness w;
  if (r.tier != kTierS0 && r.tier != kTierS1 && r.tier != kTierSL && r.tier != kTierPR) return w;
  w.t = r.support;
  for (int st = 0; st < 2; ++st)
    for (int j = 0; j < r.parts; ++j)
      if ((r.probes[st] >> j) & 1) {
        int off, len;
        part_geometry(L, r.parts, j, &off, &len);
        w.off[st][w.n[st]] = (int16_t)off;
        w.end[st][w.n[st]] = (int16_t)(off + len);
        ++w.n[st];
      }
  return w;
}

// Tier R's parts (t = 1).
inline Witness witness_of(const RepeatOut& iv) {
  Witness w;
  w.t = 1;
  for (int st = 0; st < 2; ++st) {
    w.n[st] = iv.np[st];
    for (int j = 0; j < iv.np[st]; ++j) {
      w.off[st][j] = iv.poff[st][j];
      w.end[st][j] = iv.pend[st][j];
    }
  }
  return w;
}

// One SAM record of interval output.
struct Placement {
  Result r{};                // the alignment, when mapped
  bool mapped = false;
  int mapq = 0, xq = 0;      // MAPQ and its certified lower bound
  int xl = kNoSub, xu = 0;   // OPT in [xl, xu]; xl = kNoSub: no alignment found
  int sub_lo = kNoSub, sub_hi = 0;
  Witness w;                 // proves xu and sub_hi
  bool xl_sep = false;       // xl is attained by xl_r, not by r (PR: xl_r is the single-end alignment)
  Result xl_r{};
  // PR mates: the pair term of XQ (q_pe, from the certified margin R1 - S2;
  // -1 for other records), and whether sub_hi is the single-end SUB bound
  // that the single-end term of XQ rests on (relative to xl_r if xl_sep,
  // else to r).
  int q_pe = -1;
  bool pr_xh = false;
  Witness pw;  // PR mates: the pair certificate's parts, which prove S2
};

// Single-end bounds of a certified result (tiers S0, S1, SR, SL).
inline void certified_bounds(const Result& r, int L, int* xl, int* xu, int* sub_hi) {
  *xl = r.score;
  if (r.tier == kTierSL) {
    *xu = r.score;
    *sub_hi = r.sub_hi;
  } else if (r.tier == kTierSR) {
    *xu = *sub_hi = L * kMatch;
  } else {
    const int dstar = r.n_best > 1 ? r.d1 : (r.d2 >= 0 ? r.d2 : r.radius + 1);
    *xu = std::max<int>(r.score, edit_bound(L, dstar));
    *sub_hi = r.sub_hi;
  }
}

// Single-end placement: a certified result (cert_mapq is its MAPQ), else
// tier R's alignment (unmapped below kMinScore, as in BWA-MEM).
inline Placement place_single(const Result& r, int cert_mapq, const RepeatOut* iv, int L) {
  Placement o;
  if (r.certified == 1 && r.tier != kTierPR) {
    o.r = r;
    o.mapped = true;
    o.mapq = cert_mapq;
    certified_bounds(r, L, &o.xl, &o.xu, &o.sub_hi);
    o.sub_lo = r.sub_score;
    o.xq = std::min(o.mapq, mapq_se(r.score, o.sub_hi, L));
    o.w = witness_of(r, L);
    return o;
  }
  o.xu = L * kMatch;
  if (!iv || !iv->claim) return o;
  o.xl = iv->lo;
  o.xu = iv->hi;
  o.sub_lo = iv->sub_lo;
  o.sub_hi = iv->sub_hi;
  o.w = witness_of(*iv);
  if (iv->has_aln && iv->lo < kMinScore) {  // too low to report: XL's alignment goes in the witness
    o.xl_sep = true;
    o.xl_r = alignment_only(iv->aln, iv->aln.strand, kTierIV);
  }
  if (!iv->has_aln || iv->lo < kMinScore) {
    if (!iv->has_aln) o.xl = kNoSub;  // no alignment to show for it
    return o;
  }
  o.r = alignment_only(iv->aln, iv->aln.strand, iv->exact ? kTierR : kTierIV);
  o.r.certified = iv->exact ? 1 : 0;
  o.r.n_best = (uint16_t)(iv->copies >= 2 && iv->exact ? 2 : 1);
  o.r.sub_score = (int16_t)std::max<int>(iv->sub_lo, kNoSub);
  o.mapped = true;
  o.mapq = iv->copies >= 2 ? 0 : mapq_se(iv->lo, iv->sub_lo, L);
  o.xq = std::min(o.mapq, mapq_se(iv->lo, iv->sub_hi, L));
  return o;
}

inline bool same_alignment(const Result& a, const Result& b) {
  if (a.ref_pos != b.ref_pos || a.strand != b.strand || a.n_cigar != b.n_cigar) return false;
  for (int i = 0; i < a.n_cigar; ++i)
    if (a.cigar[i] != b.cigar[i]) return false;
  return true;
}

// A PR mate r (mate: the other mate's alignment; se: its single-end result,
// whose MAPQ is se_mapq). The single-end bounds come from se when it is
// certified, else XU = L. XQ follows BWA-MEM's pairing rule,
// max(q_se, min(q_pe, q_se + 40)):
//  - q_pe from the certified pair margin R1 - S2 (0 for a tie);
//  - q_se from se and its proven SUB bound, which is about loci not
//    overlapping se. So se is shown (xl_r) whenever it is not r itself.
inline Placement pr_placement(const Result& r, const Result& mate, const Result& se, int L, int mapq,
                              int se_mapq) {
  Placement o;
  o.r = r;
  o.mapped = true;
  o.mapq = mapq;
  int q_se = 0;
  o.xl = r.score;
  o.xu = L * kMatch;
  if (se.certified == 1) {
    int sh;
    certified_bounds(se, L, &o.xl, &o.xu, &sh);
    if (!same_alignment(se, r)) {
      o.xl_sep = true;
      o.xl_r = se;
    }
    o.xl = std::max<int>(o.xl, r.score);
    o.sub_hi = sh;
    o.pr_xh = true;
    q_se = std::min(se_mapq, mapq_se(se.score, sh, L));
    o.w = witness_of(se, L);
  }
  o.q_pe = r.n_best == 1 ? mapq_margin(r.score + mate.score - r.sub_score) : 0;
  o.xq = std::min(mapq, std::max(q_se, std::min(o.q_pe, q_se + 40)));
  o.pw = witness_of(r, L);
  return o;
}

// A* moved to its proven copy (Theorem R): the same CIGAR there, or, for a
// copy of the reverse complement, the other strand with the CIGAR reversed.
// Same score.
inline Result copy_placement(const RepeatOut& iv) {
  Result c = alignment_only(iv.aln, iv.aln.strand, kTierIV);
  c.ref_pos = iv.copy_pos[1];
  if (iv.copy_strand[1] != iv.copy_strand[0]) {
    c.strand = (uint8_t)(1 - c.strand);
    std::reverse(c.cigar, c.cigar + c.n_cigar);
  }
  return c;
}

// The best alignment of a mate (codes forward, length L) that makes a proper
// pair with alignment x of the other mate: opposite strand, the reverse
// one starting 0..max_dist after the forward one. Exhaustive in that window
// above max(need, the window threshold); false if nothing scores above it.
inline bool rescue_mate(const IndexView& ix, int max_dist, int wq, const Result& x, const uint8_t* codes, int L,
                        int need, PairScratch& sc, Result& out) {
  if (L > LMAX) return false;
  const int st = 1 - x.strand;
  uint8_t rd[LMAX];
  for (int i = 0; i < L; ++i) {
    const uint8_t c = codes[st ? L - 1 - i : i];
    rd[i] = st ? (c < 4 ? (uint8_t)(3 - c) : (uint8_t)4) : c;
  }
  int64_t smin, smax;
  if (x.strand == 0) {
    smin = x.ref_pos;
    smax = x.ref_pos + max_dist;
  } else {
    smin = x.ref_pos - max_dist;
    smax = x.ref_pos;
  }
  Band b;
  int thr;
  const int v = window_search(ix, rd, L, smin, smax, wq, need, sc, &b, &thr);
  if (v == kNoSub) return false;
  Result a;
  if (!affine_align(ix, rd, L, b.lo, b.B, *sc.ws, a, false, smin, smax) || a.score != v) return false;
  out = alignment_only(a, st, kTierMR);
  return true;
}

// Places a pair in which at least one mate is not certified. se[x]: the
// single-end placements; iv[x]: tier R's output for mate x, or null.
// contig_of(pos) gives a reference position's contig. Returns true for a
// proper pair (out = the pair), false otherwise (out = se).
template <class ContigOf>
inline bool place_pair(const IndexView& ix, int max_dist, int wq, const Placement* se, const RepeatOut* const* iv,
                       const uint8_t* const* codes, const int* L, PairScratch& sc, ContigOf contig_of,
                       Placement* out) {
  out[0] = se[0];
  out[1] = se[1];
  struct Cand { Result a[2]; int score; };
  Cand cand[8];
  int nc = 0;
  auto add = [&](const Result& a0, const Result& a1) {
    if (nc == 8 || a0.strand == a1.strand) return;
    const Result& f = a0.strand ? a1 : a0;
    const Result& r = a0.strand ? a0 : a1;
    if (r.ref_pos < f.ref_pos || r.ref_pos - f.ref_pos > max_dist) return;
    const int c = contig_of(a0.ref_pos);
    if (c < 0 || contig_of(a0.ref_pos + ref_span_of(a0) - 1) != c || contig_of(a1.ref_pos) != c ||
        contig_of(a1.ref_pos + ref_span_of(a1) - 1) != c)
      return;
    cand[nc].a[0] = a0;
    cand[nc].a[1] = a1;
    cand[nc].score = a0.score + a1.score;
    ++nc;
  };
  if (se[0].mapped && se[1].mapped) add(se[0].r, se[1].r);
  for (int x = 0; x < 2; ++x) {
    const int y = 1 - x;
    Result anc[2];
    int na = 0;
    if (se[x].mapped) anc[na++] = se[x].r;
    if (se[x].mapped && iv[x] && iv[x]->has_aln && iv[x]->copies >= 2 && se[x].r.tier != kTierMR &&
        (se[x].r.tier == kTierR || se[x].r.tier == kTierIV))
      anc[na++] = copy_placement(*iv[x]);
    for (int k = 0; k < na; ++k) {
      const int need = se[y].mapped ? se[y].r.score - kPenUnpaired : kMinScore;
      Result w;
      if (rescue_mate(ix, max_dist, wq, anc[k], codes[y], L[y], need, sc, w)) {
        if (x == 0) add(anc[k], w);
        else add(w, anc[k]);
      }
      if (k > 0 && se[y].mapped) {
        if (x == 0) add(anc[k], se[y].r);
        else add(se[y].r, anc[k]);
      }
    }
  }
  if (nc == 0) return false;
  int b = 0;
  for (int c = 1; c < nc; ++c)
    if (cand[c].score > cand[b].score) b = c;
  int best_se[2];
  for (int x = 0; x < 2; ++x) {
    best_se[x] = se[x].mapped ? se[x].r.score : -BIG;
    for (int c = 0; c < nc; ++c) best_se[x] = std::max(best_se[x], (int)cand[c].a[x].score);
  }
  const int un = best_se[0] + best_se[1] - kPenUnpaired;
  if (cand[b].score <= un) return false;
  int subo = un;
  for (int c = 0; c < nc; ++c)
    if (c != b && (!same_locus(cand[c].a[0], cand[b].a[0]) || !same_locus(cand[c].a[1], cand[b].a[1])))
      subo = std::max(subo, cand[c].score);
  // BWA-MEM pairs every alignment it found of one mate with every one of the
  // other; here a repeat mate contributes one alignment (and one copy), so
  // the other pairs of a duplicated region are mostly unseen. The pair is
  // therefore trusted no more than its more confident mate: a unique mate
  // vouches for its partner's window, two repeat mates vouch for nothing.
  bool at_se[2];
  int q_se[2];
  for (int x = 0; x < 2; ++x) {
    const Result& p = cand[b].a[x];
    at_se[x] = se[x].mapped && same_locus(p, se[x].r) && se[x].r.score >= p.score;
    q_se[x] = at_se[x] ? se[x].mapq : 0;
  }
  const int q_pe = std::min(mapq_margin(cand[b].score - subo), std::max(q_se[0], q_se[1]));
  for (int x = 0; x < 2; ++x) {
    const Result& p = cand[b].a[x];
    Placement& o = out[x];
    if (!at_se[x]) {
      // XL stays the single-end best when that beats the pair placement.
      const Result best_se = se[x].mapped ? se[x].r : se[x].xl_r;
      const bool se_best = (se[x].mapped || se[x].xl_sep) && best_se.score > p.score;
      o.r = p;
      o.r.tier = kTierMR;
      o.r.certified = 0;
      o.mapped = true;
      o.xl = se_best ? (int)best_se.score : (int)p.score;
      o.xl_sep = se_best;
      if (se_best) o.xl_r = best_se;
      o.xq = 0;
    }
    o.mapq = std::max(q_se[x], std::min(q_pe, q_se[x] + 40));
  }
  return true;
}

}  // namespace certa
