// --mapq pair (paired interval mode): MAPQ from boosted trees over both
// mates' records, fitted at the locus on HG002 (bench/truth_mapq_pair.py;
// bench/VALIDATION_RESULTS.md, "A MAPQ from both mates"). With --t2t the
// trees also read each mate's placement on T2T-CHM13, lifted to the
// reference (certa/mapq_model_t2t.h).
//
// --mapq cal ranks a record by its own tier, edits, gap to the second-best
// locus and clusters. At the operating points where most mates count as
// confident, the mates that decide the ranking are in repeats: whether the
// other mate is unique, and how far the pair score stands from the next
// pair, then say more than the record itself. minibwa, as BWA-MEM, raises a
// repeat mate's MAPQ by its unique mate's; the trees learn how far that
// holds on real reads. On T2T-CHM13, a mate from a copy GRCh38 lacks tends
// to align better elsewhere or not lift.
//
// The placement, XQ and every certificate are unchanged: only MAPQ is
// replaced. The model is a fixed set of trees, so MAPQ is deterministic.
#pragma once
#include <cstdint>

#include "certa/interval.h"
#include "certa/mapq_model.h"
#include "certa/mapq_model_t2t.h"

namespace certa {

// A mate's placement on T2T-CHM13 (--t2t) as the model reads it.
struct T2tView {
  // 0: unmapped there; 1: lifted within 150 bp of the record's unclipped
  // leftmost base (same contig and strand); 2: lifted elsewhere; 3: mapped
  // but not lifted. A record that is unmapped on the reference keeps 0.
  int cls = 0;
  int tier = -1, mapq = 0, xe = -1, xd = -1, xo = -1;  // its tier, MAPQ, XE, XD and XO there
};

constexpr int kPairOwn = 12, kPairPart = 8;          // features of the mate and of the other mate
constexpr int kPairFeatures = kPairOwn + kPairPart;  // == kModelFeatures
constexpr int kT2tFeatures = kPairFeatures + 7;      // == kT2tModelFeatures
static_assert(kModelFeatures == kPairFeatures, "mapq_model.h was fitted on other features");
static_assert(kT2tModelFeatures == kT2tFeatures, "mapq_model_t2t.h was fitted on other features");

inline int clampi(int v, int lo, int hi) { return v < lo ? lo : (v > hi ? hi : v); }

inline int floor_log2(uint32_t v) {  // v >= 1
  int k = -1;
  while (v) {
    v >>= 1;
    ++k;
  }
  return k;
}

inline int soft_clips(const Result& r) {
  int n = 0;
  for (int i = 0; i < r.n_cigar; ++i)
    if ((r.cigar[i] & 0xF) == kOpS) n += static_cast<int>(r.cigar[i] >> 4);
  return n;
}

// One record's own inputs (truth_mapq_pair.py features: tier, mq, gap, xe,
// xn, xb, xq, nm, clip; margin, pp and tlen need the other mate).
inline void record_inputs(const Placement& p, int* f) {
  if (!p.mapped) {
    for (int k = 0; k < 9; ++k) f[k] = 0;
    f[2] = 11;
    return;
  }
  const Result& r = p.r;
  const int second = r.d2 >= 0 ? r.d2 : r.d2x;
  f[0] = r.tier + 1;
  f[1] = clampi(p.mapq, 0, 60);
  f[2] = second >= 0 ? clampi(second - r.d1, 0, 10) : 11;
  f[3] = clampi(r.d1, -1, 10) + 1;
  f[4] = std::min(floor_log2(static_cast<uint32_t>(r.n_clusters) + 1), 15);
  f[5] = clampi(r.n_best, 0, 3);
  f[6] = clampi(p.xq, 0, 60);
  f[7] = clampi(r.nm, 0, 15);
  f[8] = std::min(soft_clips(r) / 5, 15);
}

// A PR mate's pair margin: its AS plus the other mate's minus XS, the bound
// on every other proper pair (0: not a PR mate with a mapped partner).
inline int pair_margin(const Placement& a, const Placement& b) {
  if (!a.mapped || !b.mapped || a.r.tier != kTierPR || a.r.sub_score <= -1000) return 0;
  return clampi(a.r.score + b.r.score - a.r.sub_score, 0, 40) + 1;
}

// The inputs of mate a's MAPQ (its partner b; proper: the pair's proper flag;
// tlen: |TLEN| of the records, 0 if not on one contig).
inline void pair_model_features(const Placement& a, const Placement& b, bool proper, int64_t tlen, int* f) {
  int own[9], part[9];
  record_inputs(a, own);
  record_inputs(b, part);
  for (int k = 0; k < 9; ++k) f[k] = own[k];
  f[9] = pair_margin(a, b);
  f[10] = (proper ? 2 : 0) + (b.mapped ? 0 : 1);
  static const int64_t kTlen[] = {1, 200, 300, 400, 500, 700, 1000};
  int t = 0;
  if (a.mapped)
    for (int64_t v : kTlen) t += tlen >= v;
  f[11] = t;
  // The other mate: tier, mq, gap, xe, xq, nm, clip, margin.
  const int pk[] = {0, 1, 2, 3, 6, 7, 8};
  for (int k = 0; k < 7; ++k) f[12 + k] = part[pk[k]];
  f[19] = pair_margin(b, a);
}

// The T2T inputs of mate a (ta, tb: both mates' T2T views), after
// pair_model_features: t_cls, t_gap, t_dx, t_mq, t_tier, and the other
// mate's t_cls and t_gap.
inline void t2t_model_features(const Placement& a, const T2tView& ta, const T2tView& tb, int* f) {
  auto gap = [](const T2tView& t) {
    if (t.cls == 0) return 12;
    const int second = t.xd >= 0 ? t.xd : t.xo;
    return second >= 0 ? clampi(second - t.xe, 0, 10) : 11;
  };
  f[0] = ta.cls;
  f[1] = gap(ta);
  f[2] = ta.cls > 0 && ta.xe >= 0 && a.mapped && a.r.d1 >= 0 ? clampi(a.r.d1 - ta.xe, -3, 5) + 3 : 9;
  f[3] = ta.cls > 0 ? std::min(ta.mapq, 60) / 5 : 0;
  f[4] = ta.cls > 0 ? ta.tier + 1 : 0;
  f[5] = tb.cls;
  f[6] = gap(tb);
}

// Sum of the trees' leaves over x, in tree order (the score is the log odds
// of a misplacement), and the MAPQ it gives: -10 log10 of the probability,
// rounded, at most 60, by the fitted thresholds.
template <int NT, int NN, int NL>
double tree_score(const int* x, double base, const int8_t (&feat)[NT][NN], const int8_t (&thr)[NT][NN],
                  const double (&leaf)[NT][NL]) {
  double f = base;
  for (int t = 0; t < NT; ++t) {
    int node = 0;
    for (int d = 0; (1 << d) < NL; ++d) {
      const int k = (1 << d) - 1 + node;
      const int j = feat[t][k];
      node = node * 2 + (j >= 0 && x[j] > thr[t][k] ? 1 : 0);
    }
    f += leaf[t][node];
  }
  return f;
}

inline int score_mapq(double s, const double (&q)[60]) {
  int m = 0;
  for (int k = 0; k < 60; ++k) m += s <= q[k];
  return m;
}

// --mapq pair: both mates' MAPQ (pl: their placements after --iv-mapq; t2t:
// their T2T views, or null without --t2t). Unmapped mates keep theirs.
inline void pair_model_mapq(Placement* pl, bool proper, int64_t tlen, const T2tView* t2t) {
  int q[2] = {pl[0].mapq, pl[1].mapq};
  for (int x = 0; x < 2; ++x) {
    if (!pl[x].mapped) continue;
    int f[kT2tFeatures];
    pair_model_features(pl[x], pl[1 - x], proper, tlen, f);
    if (t2t) {
      t2t_model_features(pl[x], t2t[x], t2t[1 - x], f + kPairFeatures);
      q[x] = score_mapq(tree_score(f, kT2tModelBase, kT2tModelFeat, kT2tModelThr, kT2tModelLeaf), kT2tModelQ);
    } else {
      q[x] = score_mapq(tree_score(f, kModelBase, kModelFeat, kModelThr, kModelLeaf), kModelQ);
    }
  }
  pl[0].mapq = q[0];
  pl[1].mapq = q[1];
}

}  // namespace certa
