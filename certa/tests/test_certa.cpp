// Brute-force validation of the certificate.
//
// The oracle is Sellers' semi-global edit-distance DP over the *entire*
// reference on both strands, with Ukkonen's cutoff: it lists every
// reference end position where the read aligns with <= KMAX edits.
// For every read and every (k, cap) setting it checks:
//   1. soundness      - the reported CIGAR re-scores to exactly d1 edits;
//   2. optimality     - a certified d1 equals the oracle's global minimum;
//   3. losslessness   - if the oracle has any locus within the radius R, the
//                       read is not rejected as "not found";
//   4. uniqueness     - if CERTA reports a unique best with no second locus
//                       within R, the oracle has no other locus within R;
//   5. determinism    - results are identical with 1 and 4 threads.
// Tier SL (certified local alignment) is checked against a second oracle:
// the best affine local score (bwa-mem defaults, clip penalty per clipped
// end) over the entire reference, both strands, by full Gotoh DP.
// Interval output (certa/interval.h): every [XL, XU] must contain that
// optimum, and every SUB bound must hold against the best score at a locus
// not overlapping the reported one (the same DP on the other strand and on
// the reference left and right of the reported interval).
// The independent checker (check.h, Theorem V) must accept every honest
// claim and reject every claim made false against the oracles, and index
// validation must reject corrupted indexes and FASTAs.
#include <algorithm>
#include <climits>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <map>
#include <memory>
#include <random>
#include <string>
#include <functional>
#include <tuple>
#include <vector>

#include <unistd.h>

#include "index.h"
#include "certa/interval.h"
#include "check.h"
#include "certa/pair.h"
#include "certa/repeat.h"
#include "mapper.h"

using namespace certa;

namespace {

int failures = 0;
#define CHECK(cond, ...)                                  \
  do {                                                    \
    if (!(cond)) {                                        \
      if (++failures <= 20) {                             \
        std::fprintf(stderr, "FAIL %s:%d: ", __FILE__, __LINE__); \
        std::fprintf(stderr, __VA_ARGS__);                \
        std::fprintf(stderr, "\n");                       \
      }                                                   \
    }                                                     \
  } while (0)

struct Hit { int64_t end; int dist; int strand; };

// All end positions (inclusive) with distance <= K, both strands.
std::vector<Hit> oracle(const std::vector<uint8_t>& ref, const uint8_t* fwd, int L, int K) {
  std::vector<uint8_t> rc(L);
  for (int i = 0; i < L; ++i) rc[i] = fwd[L - 1 - i] < 4 ? 3 - fwd[L - 1 - i] : 4;
  std::vector<Hit> hits;
  std::vector<int> C(L + 1);
  for (int st = 0; st < 2; ++st) {
    const uint8_t* rd = st ? rc.data() : fwd;
    for (int i = 0; i <= L; ++i) C[i] = std::min(i, K + 1);
    int lact = std::min(K, L);
    for (size_t j = 0; j < ref.size(); ++j) {
      int diag = 0;  // D[0][j-1] = 0 (free start)
      int top = std::min(lact + 1, L);
      for (int i = 1; i <= top; ++i) {
        int old = C[i];
        int v = diag + sub_cost(rd[i - 1], ref[j]);
        if (C[i - 1] + 1 < v) v = C[i - 1] + 1;
        if (old + 1 < v) v = old + 1;
        C[i] = std::min(v, K + 1);
        diag = old;
      }
      lact = top;
      while (lact > 0 && C[lact] > K) --lact;
      if (lact == L) hits.push_back({static_cast<int64_t>(j), C[L], st});
    }
  }
  return hits;
}

// Best affine local score over the whole reference (the objective of
// affine_align with an unbounded band): alignments end with a match/mismatch
// and pay kClip per clipped read end.
int oracle_local(const std::vector<uint8_t>& ref, const uint8_t* fwd, int L) {
  const int NEG = -(1 << 20);
  std::vector<uint8_t> rc(L);
  for (int i = 0; i < L; ++i) rc[i] = fwd[L - 1 - i] < 4 ? 3 - fwd[L - 1 - i] : 4;
  std::vector<int> hp(L + 1), hc(L + 1), yp(L + 1), yc(L + 1);
  int best = NEG;
  for (int st = 0; st < 2; ++st) {
    const uint8_t* rd = st ? rc.data() : fwd;
    hp.assign(L + 1, NEG);
    yp.assign(L + 1, NEG);
    hp[0] = 0;
    for (size_t j = 0; j < ref.size(); ++j) {
      hc[0] = 0;
      yc[0] = NEG;
      int x = NEG;  // insertion state at (i - 1, j)
      for (int i = 1; i <= L; ++i) {
        const int start = i == 1 ? 0 : -kClip;
        const int m = std::max(hp[i - 1], start) +
                      (sub_cost(rd[i - 1], ref[j]) ? -kMismatch : kMatch);
        x = i > 1 ? std::max(hc[i - 1] - kGapOpen - kGapExt, x - kGapExt) : NEG;
        yc[i] = std::max(hp[i] - kGapOpen - kGapExt, yp[i] - kGapExt);
        hc[i] = std::max(m, std::max(x, yc[i]));
        best = std::max(best, m - (i < L ? kClip : 0));
      }
      std::swap(hp, hc);
      std::swap(yp, yc);
    }
  }
  return best;
}

// Best affine local score of `rd` (one orientation) against ref[j0, j1),
// the same recurrences as oracle_local; kNoSub-like -2^20 if the range is empty.
int local_one(const std::vector<uint8_t>& ref, const uint8_t* rd, int L, int64_t j0, int64_t j1) {
  const int NEG = -(1 << 20);
  std::vector<int> hp(L + 1, NEG), hc(L + 1), yp(L + 1, NEG), yc(L + 1);
  hp[0] = 0;
  int best = NEG;
  for (int64_t j = std::max<int64_t>(0, j0); j < j1 && j < static_cast<int64_t>(ref.size()); ++j) {
    hc[0] = 0;
    yc[0] = NEG;
    int x = NEG;
    for (int i = 1; i <= L; ++i) {
      const int start = i == 1 ? 0 : -kClip;
      const int m = std::max(hp[i - 1], start) + (sub_cost(rd[i - 1], ref[j]) ? -kMismatch : kMatch);
      x = i > 1 ? std::max(hc[i - 1] - kGapOpen - kGapExt, x - kGapExt) : NEG;
      yc[i] = std::max(hp[i] - kGapOpen - kGapExt, yp[i] - kGapExt);
      hc[i] = std::max(m, std::max(x, yc[i]));
      best = std::max(best, m - (i < L ? kClip : 0));
    }
    std::swap(hp, hc);
    std::swap(yp, yc);
  }
  return best;
}

// SUB: the best local score at a locus not overlapping [p, p + span) on
// `strand` (anything on the other strand; on the same strand, alignments
// lying wholly left or right of the interval).
int oracle_sub(const std::vector<uint8_t>& ref, const std::vector<uint8_t>& fwd, int strand, int64_t p,
               int64_t span) {
  const int L = static_cast<int>(fwd.size());
  std::vector<uint8_t> rc(L);
  for (int i = 0; i < L; ++i) rc[i] = fwd[L - 1 - i] < 4 ? 3 - fwd[L - 1 - i] : 4;
  const uint8_t* same = strand ? rc.data() : fwd.data();
  const uint8_t* other = strand ? fwd.data() : rc.data();
  const int64_t n = static_cast<int64_t>(ref.size());
  return std::max(local_one(ref, other, L, 0, n),
                  std::max(local_one(ref, same, L, 0, p), local_one(ref, same, L, p + span, n)));
}

// For every reference position j: the best affine local score of `rd`
// (aligned forward) among alignments whose first aligned base is j. Local DP
// on the reversed read against the reference read backwards (its end column
// is the start in forward coordinates; clip penalties are symmetric).
std::vector<int> oracle_starts(const std::vector<uint8_t>& ref, const uint8_t* rd, int L) {
  const int NEG = -(1 << 20);
  std::vector<int> best(ref.size(), NEG), hp(L + 1), hc(L + 1), yp(L + 1), yc(L + 1);
  hp.assign(L + 1, NEG);
  yp.assign(L + 1, NEG);
  hp[0] = 0;
  for (size_t jj = ref.size(); jj-- > 0;) {
    hc[0] = 0;
    yc[0] = NEG;
    int x = NEG;
    for (int i = 1; i <= L; ++i) {
      const uint8_t c = rd[L - i];
      const int start = i == 1 ? 0 : -kClip;
      const int m = std::max(hp[i - 1], start) + (sub_cost(c, ref[jj]) ? -kMismatch : kMatch);
      x = i > 1 ? std::max(hc[i - 1] - kGapOpen - kGapExt, x - kGapExt) : NEG;
      yc[i] = std::max(hp[i] - kGapOpen - kGapExt, yp[i] - kGapExt);
      hc[i] = std::max(m, std::max(x, yc[i]));
      best[jj] = std::max(best[jj], m - (i < L ? kClip : 0));
    }
    std::swap(hp, hc);
    std::swap(yp, yc);
  }
  return best;
}

// Best proper pair: opposite strands, reverse mate starting 0..D after the
// forward mate (sliding-window maximum).
int oracle_pair(const std::vector<uint8_t>& ref, const std::vector<uint8_t>& m1,
                const std::vector<uint8_t>& m2, int D) {
  auto rc = [](const std::vector<uint8_t>& v) {
    std::vector<uint8_t> o(v.size());
    for (size_t i = 0; i < v.size(); ++i) o[i] = v[v.size() - 1 - i] < 4 ? 3 - v[v.size() - 1 - i] : 4;
    return o;
  };
  const std::vector<uint8_t> r1 = rc(m1), r2 = rc(m2);
  const auto F1 = oracle_starts(ref, m1.data(), static_cast<int>(m1.size()));
  const auto R1 = oracle_starts(ref, r1.data(), static_cast<int>(r1.size()));
  const auto F2 = oracle_starts(ref, m2.data(), static_cast<int>(m2.size()));
  const auto R2 = oracle_starts(ref, r2.data(), static_cast<int>(r2.size()));
  const int n = static_cast<int>(ref.size());
  int best = INT32_MIN;
  for (int pass = 0; pass < 2; ++pass) {
    const std::vector<int>& F = pass ? F2 : F1;
    const std::vector<int>& R = pass ? R1 : R2;
    std::vector<int> dq;  // indices into R, decreasing values, window [f, f + D]
    size_t head = 0;
    int next = 0;
    for (int f = 0; f < n; ++f) {
      while (next < n && next <= f + D) {
        while (dq.size() > head && R[dq.back()] <= R[next]) dq.pop_back();
        dq.push_back(next++);
      }
      while (head < dq.size() && dq[head] < f) ++head;
      if (head < dq.size()) best = std::max(best, F[f] + R[dq[head]]);
    }
  }
  return best;
}

// The claims of an interval-output record, as certa check reads them.
check::Claim claim_of(const Placement& p, const std::vector<uint8_t>& read, int pair_dist = -1) {
  check::Claim c;
  c.read = read;
  c.mapped = p.mapped;
  auto aln = [](const Result& r) {
    check::Aln a;
    a.pos = r.ref_pos;
    a.strand = r.strand;
    a.cigar.assign(r.cigar, r.cigar + r.n_cigar);
    return a;
  };
  if (p.mapped) {
    c.aln = aln(p.r);
    c.as = p.r.score;
  }
  c.mapq = p.mapq;
  c.xl = p.xl == kNoSub ? check::kNone : p.xl;
  c.xu = p.xu;
  c.xq = p.xq;
  c.pair_xq = p.mapped && p.r.tier == kTierPR;
  if (p.xq > 0 && (!c.pair_xq || p.pr_xh)) c.xh = p.sub_hi;
  if (c.pair_xq) c.xs = p.r.sub_score;
  if (p.xl_sep) {
    c.has_alt = true;
    c.alt = aln(p.xl_r);
  }
  c.t = p.w.t;
  for (int st = 0; st < 2; ++st)
    for (int j = 0; j < p.w.n[st]; ++j) c.probes[st].push_back({p.w.off[st][j], p.w.end[st][j]});
  if (p.mapped && p.pw.t > 0) {
    c.pd = pair_dist;
    c.pt = p.pw.t;
    for (int st = 0; st < 2; ++st)
      for (int j = 0; j < p.pw.n[st]; ++j) c.pprobes[st].push_back({p.pw.off[st][j], p.pw.end[st][j]});
  }
  return c;
}

// A certified pair (r: both mates, tier PR; se: their single-end results)
// as certa map --intervals --witness writes it, through the checker's pair
// check. Honest claims must be accepted. Claims made false must not be:
// XS one below sub (the best proper pair differing in some mate, from the
// oracle) with XQ raised to what that margin gives, the pair witness
// dropped, D too small for the pair to be proper, XQ one above the pairing
// rule, a mate's alignment moved.
struct PrCheck {
  int pairs = 0, needed = 0, proven = 0, tampered = 0, caught = 0;
};
void check_pr(const check::Checker& ck, const Result* r, const Result* se, const std::vector<uint8_t>* m, int D,
              int sub, const char* what, size_t id, PrCheck& st) {
  check::Claim c[2];
  for (int x = 0; x < 2; ++x) {
    const Placement pl = pr_placement(r[x], r[1 - x], se[x], static_cast<int>(m[x].size()), 60, 60);
    c[x] = claim_of(pl, m[x], D);
  }
  check::Work w[2];
  check::Verdict v[2];
  ck.check_pair(c[0], c[1], &w[0], &w[1], &v[0], &v[1]);
  ++st.pairs;
  for (int x = 0; x < 2; ++x) {
    st.needed += w[x].needs_pair;
    st.proven += w[x].xs_checked;
    CHECK(v[x] == check::kAccepted, "%s %zu mate %d: checker says %s (XQ %d, q_se %d, XS %d, AS %d + %d)", what, id,
          x + 1, check::kVerdictNames[v[x]], c[x].xq, w[x].q_se, c[x].xs, c[0].as, c[1].as);
  }
  auto expect_reject = [&](check::Claim a, check::Claim b, const char* how) {
    check::Work wa, wb;
    check::Verdict va, vb;
    ck.check_pair(a, b, &wa, &wb, &va, &vb);
    ++st.tampered;
    const bool caught = va != check::kAccepted || vb != check::kAccepted;
    st.caught += caught;
    CHECK(caught, "%s %zu: false pair claim accepted (%s)", what, id, how);
  };
  const int R1 = c[0].as + c[1].as;
  if (sub > INT32_MIN / 2 && sub - 1 < c[0].xs) {
    check::Claim t[2] = {c[0], c[1]};
    int raised = 0;
    for (int x = 0; x < 2; ++x) {
      t[x].xs = sub - 1;
      const int q = std::max(w[x].q_se, std::min(mapq_margin(R1 - (sub - 1)), w[x].q_se + 40));
      if (q > w[x].q_se) {
        t[x].xq = q;
        t[x].mapq = 60;
        ++raised;
      }
    }
    if (raised) expect_reject(t[0], t[1], "XS = the best other pair - 1");  // also ties (XQ 0 honestly)
  }
  if (!w[0].needs_pair && !w[1].needs_pair) return;  // XS carries no other claim here
  {
    check::Claim t[2] = {c[0], c[1]};
    t[0].pt = t[1].pt = 0;
    expect_reject(t[0], t[1], "no pair witness");
  }
  {
    check::Claim t[2] = {c[0], c[1]};
    const check::Aln& f = t[0].aln.strand ? t[1].aln : t[0].aln;
    const check::Aln& rv = t[0].aln.strand ? t[0].aln : t[1].aln;
    t[0].pd = t[1].pd = static_cast<int>(rv.pos - f.pos) - 1;
    expect_reject(t[0], t[1], "D one below the pair's distance");
  }
  for (int x = 0; x < 2; ++x) {
    if (!w[x].needs_pair) continue;
    const int q = std::max(w[x].q_se, std::min(mapq_margin(R1 - c[x].xs), w[x].q_se + 40));
    if (q >= 60) continue;
    check::Claim t[2] = {c[0], c[1]};
    t[x].xq = q + 1;
    t[x].mapq = 60;
    expect_reject(t[0], t[1], "XQ one above the pairing rule");
  }
  {
    check::Claim t[2] = {c[0], c[1]};
    t[0].aln.pos += 1;
    expect_reject(t[0], t[1], "mate 1 moved by one base");
  }
}

// The best pair from per-start bests: forward mate starting at f, reverse
// mate starting at r in [f, f + D] (sliding-window maximum).
int pair_from_starts(const std::vector<int>& F, const std::vector<int>& R, int D) {
  const int n = static_cast<int>(F.size());
  int best = INT32_MIN;
  std::vector<int> dq;
  size_t head = 0;
  int next = 0;
  for (int f = 0; f < n; ++f) {
    while (next < n && next <= f + D) {
      while (dq.size() > head && R[dq.back()] <= R[next]) dq.pop_back();
      dq.push_back(next++);
    }
    while (head < dq.size() && dq[head] < f) ++head;
    if (head < dq.size()) best = std::max(best, F[f] + R[dq[head]]);
  }
  return best;
}

// Per-start bests of `rd` over alignments that do not overlap [p, p + span):
// the oracle on the reference left of it and right of it.
std::vector<int> starts_outside(const std::vector<uint8_t>& ref, const uint8_t* rd, int L, int64_t p, int64_t span) {
  std::vector<int> out(ref.size(), -(1 << 20));
  const std::vector<uint8_t> left(ref.begin(), ref.begin() + std::max<int64_t>(0, p));
  const std::vector<uint8_t> right(ref.begin() + std::min<int64_t>(static_cast<int64_t>(ref.size()), p + span), ref.end());
  const std::vector<int> a = oracle_starts(left, rd, L), b = oracle_starts(right, rd, L);
  for (size_t j = 0; j < a.size(); ++j) out[j] = a[j];
  for (size_t j = 0; j < b.size(); ++j) out[ref.size() - right.size() + j] = b[j];
  return out;
}

// SUB of a pair: the best proper pair in which some mate does not overlap
// the reported one (a1 for m1, a2 for m2; strand 1 = reverse complement).
int oracle_pair_sub(const std::vector<uint8_t>& ref, const std::vector<uint8_t>& m1, const std::vector<uint8_t>& m2,
                    int D, const Result& a1, const Result& a2) {
  auto rc = [](const std::vector<uint8_t>& v) {
    std::vector<uint8_t> o(v.size());
    for (size_t i = 0; i < v.size(); ++i) o[i] = v[v.size() - 1 - i] < 4 ? 3 - v[v.size() - 1 - i] : 4;
    return o;
  };
  auto span = [](const Result& a) {
    int64_t n = 0;
    for (int c = 0; c < a.n_cigar; ++c)
      if ((a.cigar[c] & 0xF) == kOpM || (a.cigar[c] & 0xF) == kOpD) n += a.cigar[c] >> 4;
    return n;
  };
  const std::vector<uint8_t> seq[2][2] = {{m1, rc(m1)}, {m2, rc(m2)}};  // [mate][strand]
  const Result* rep[2] = {&a1, &a2};
  std::vector<int> full[2][2], other[2][2];
  for (int m = 0; m < 2; ++m)
    for (int st = 0; st < 2; ++st) {
      const int L = static_cast<int>(seq[m][st].size());
      full[m][st] = oracle_starts(ref, seq[m][st].data(), L);
      other[m][st] = st == rep[m]->strand
                         ? starts_outside(ref, seq[m][st].data(), L, rep[m]->ref_pos, span(*rep[m]))
                         : full[m][st];
    }
  int best = INT32_MIN;
  for (int f = 0; f < 2; ++f) {  // f: the forward mate; 1 - f the reverse one
    const int r = 1 - f;
    best = std::max(best, pair_from_starts(other[f][0], full[r][1], D));
    best = std::max(best, pair_from_starts(full[f][0], other[r][1], D));
  }
  return best;
}

// Affine score of a reported alignment, clip penalties included.
int rescore_affine(const IndexView& ix, const Result& r, const std::vector<uint8_t>& fwd) {
  const int L = static_cast<int>(fwd.size());
  std::vector<uint8_t> rd(fwd);
  if (r.strand)
    for (int i = 0; i < L; ++i) rd[i] = fwd[L - 1 - i] < 4 ? 3 - fwd[L - 1 - i] : 4;
  int64_t j = r.ref_pos;
  int i = 0, sc = 0;
  for (int c = 0; c < r.n_cigar; ++c) {
    const int len = static_cast<int>(r.cigar[c] >> 4);
    switch (r.cigar[c] & 0xF) {
      case kOpS: sc -= kClip; i += len; break;
      case kOpI: sc -= kGapOpen + kGapExt * len; i += len; break;
      case kOpD: sc -= kGapOpen + kGapExt * len; j += len; break;
      default:
        for (int x = 0; x < len; ++x) sc += sub_cost(rd[i++], ref_at(ix, j++)) ? -kMismatch : kMatch;
    }
  }
  return i == L ? sc : -100000;
}

uint8_t rand_base(std::mt19937_64& g) { return static_cast<uint8_t>(g() % 4); }

std::vector<uint8_t> random_seq(std::mt19937_64& g, size_t n) {
  std::vector<uint8_t> s(n);
  for (auto& c : s) c = rand_base(g);
  return s;
}

void mutate(std::mt19937_64& g, std::vector<uint8_t>& s, int edits) {
  for (int e = 0; e < edits && !s.empty(); ++e) {
    size_t p = g() % s.size();
    switch (g() % 3) {
      case 0: s[p] = static_cast<uint8_t>((s[p] + 1 + g() % 3) % 4); break;
      case 1: s.insert(s.begin() + p, rand_base(g)); break;
      default: s.erase(s.begin() + p); break;
    }
  }
}

// Reference with the structures that stress a certificate: diverged copies
// of a segment, an exact duplicate across contigs, a tandem repeat, N runs.
Reference make_reference(std::mt19937_64& g) {
  std::vector<uint8_t> c1 = random_seq(g, 90000), c2 = random_seq(g, 40000);
  std::vector<uint8_t> seg = random_seq(g, 2000);
  for (int copy = 0; copy < 3; ++copy) {
    std::vector<uint8_t> m = seg;
    mutate(g, m, 3 + copy * 2);  // near-identical copies
    std::copy(m.begin(), m.end(), c1.begin() + 10000 + copy * 20000);
  }
  std::copy(seg.begin(), seg.end(), c2.begin() + 5000);  // exact duplicate
  std::vector<uint8_t> unit = random_seq(g, 37);
  for (int r = 0; r < 30; ++r)
    std::copy(unit.begin(), unit.end(), c1.begin() + 75000 + r * 37);
  std::fill(c1.begin() + 60000, c1.begin() + 60100, 4);  // N run
  // High-copy exact repeat (Alu-like): 40 copies in chr1, 10 in chr2, with a
  // few diverged copies, so budgets overflow and tier SR is exercised.
  std::vector<uint8_t> elem = random_seq(g, 300);
  for (int r = 0; r < 40; ++r) {
    std::vector<uint8_t> e = elem;
    if (r % 7 == 3) mutate(g, e, 2);
    e.resize(300, 0);  // a deletion may have shortened it
    std::copy(e.begin(), e.end(), c1.begin() + 60200 + r * 350);
  }
  for (int r = 0; r < 10; ++r)
    std::copy(elem.begin(), elem.end(), c2.begin() + 20000 + r * 350);
  Reference ref;
  ref.seq.owned.assign(kContigPad, 4);
  const std::vector<uint8_t>* cs[2] = {&c1, &c2};
  for (int c = 0; c < 2; ++c) {
    ref.names.push_back("chr" + std::to_string(c + 1));
    ref.offsets.push_back(ref.seq.owned.size());
    ref.lengths.push_back(cs[c]->size());
    ref.seq.owned.insert(ref.seq.owned.end(), cs[c]->begin(), cs[c]->end());
    ref.seq.owned.insert(ref.seq.owned.end(), kContigPad, 4);
  }
  return ref;
}

int rescore(const IndexView& ix, const Result& r, const std::vector<uint8_t>& fwd) {
  const int L = static_cast<int>(fwd.size());
  std::vector<uint8_t> rd(fwd);
  if (r.strand)
    for (int i = 0; i < L; ++i) rd[i] = fwd[L - 1 - i] < 4 ? 3 - fwd[L - 1 - i] : 4;
  int64_t j = r.ref_pos;
  int i = 0, d = 0;
  for (int c = 0; c < r.n_cigar; ++c) {
    uint32_t len = r.cigar[c] >> 4, op = r.cigar[c] & 0xF;
    for (uint32_t x = 0; x < len; ++x) {
      if (op == kOpM) d += sub_cost(rd[i++], ref_at(ix, j++));
      else if (op == kOpI) { ++i; ++d; }
      else if (op == kOpS) { ++i; }  // soft clip: read base, no cost
      else { ++j; ++d; }
    }
  }
  return i == L ? d : -1000;
}

// MappedFastq's parallel boundary scan must give the serial scan's batches.
// Random FASTQ files with LF or CRLF line ends, rare blank lines and stray
// '\r' between records, rare empty reads, varying lengths and sometimes no
// final newline are read in batches of several sizes on 1 thread (serial
// scan) and on several; batches and records must match the written records.
void test_fastq_scan(std::mt19937_64& g) {
  const std::string path = "test_certa." + std::to_string(getpid()) + ".tmp.fq";
  int files = 0, batches = 0;
  for (int shape = 0; shape < 6; ++shape) {
    std::vector<FastqRecord> truth;
    std::string text;
    const bool crlf = shape == 1 || shape == 4;
    const char* eol = crlf ? "\r\n" : "\n";
    const int nrec = 20000 + static_cast<int>(g() % 5000);
    for (int i = 0; i < nrec; ++i) {
      if (shape >= 2 && g() % 700 == 0) text += g() % 2 ? "\n" : "\r\n";  // blank line
      if (shape >= 3 && g() % 900 == 0) text += "\r";                    // stray '\r' before '@'
      FastqRecord r;
      r.name = "r" + std::to_string(i) + (g() % 3 == 0 ? "/1" : "");
      int len = shape == 0 ? 150 : 50 + static_cast<int>(g() % 260);
      if (shape >= 3 && g() % 1000 == 0) len = 0;  // an empty read: empty lines inside a record
      for (int j = 0; j < len; ++j) r.seq += "ACGTN"[g() % 5];
      for (int j = 0; j < len; ++j) r.qual += static_cast<char>('!' + g() % 40);
      text += "@" + r.name + (g() % 2 ? " comment x" : "") + eol + r.seq + eol + "+" + eol + r.qual;
      if (i + 1 < nrec || shape != 5) text += eol;  // shape 5: no final newline
      truth.push_back(r);
    }
    std::FILE* f = std::fopen(path.c_str(), "wb");
    std::fwrite(text.data(), 1, text.size(), f);
    std::fclose(f);
    ++files;
    for (size_t max_records : {size_t(1), size_t(7), size_t(1000), size_t(3333), size_t(9000)}) {
      if (max_records == 1 && shape > 1) continue;  // slow and covered by the others
      for (int threads : {1, 3, 8}) {
        MappedFastq mf(path);
        std::vector<FastqRecord> out;
        size_t got = 0;
        bool ok = true;
        std::string err;
        try {
          for (size_t n; ok && (n = mf.next_batch(out, max_records, threads)) > 0;) {
            ++batches;
            ok = n == out.size() && (n == max_records || got + n == truth.size());
            for (size_t i = 0; ok && i < n; ++i) {
              const FastqRecord& t = truth[got + i];
              ok = out[i].name == t.name && out[i].seq == t.seq && out[i].qual == t.qual;
            }
            got += n;
          }
        } catch (const std::exception& e) {
          ok = false;
          err = e.what();
        }
        CHECK(ok && got == truth.size(), "FASTQ scan shape %d, batch %zu, %d threads: records differ (%zu of %zu read) %s",
              shape, max_records, threads, got, truth.size(), err.c_str());
      }
    }
  }
  std::remove(path.c_str());
  std::fprintf(stderr, "FASTQ batches: %d files, %d batches identical on 1, 3 and 8 threads\n", files, batches);
}

}  // namespace

int main(int argc, char** argv) {
  const uint64_t seed = argc > 1 ? std::strtoull(argv[1], nullptr, 10) : 20260927;
  std::fprintf(stderr, "seed %llu\n", static_cast<unsigned long long>(seed));
  std::mt19937_64 g(seed);
  {
    std::mt19937_64 gq(seed ^ 0x5eedf00dULL);  // its own stream: the other tests' data stay as they were
    test_fastq_scan(gq);
  }
  Reference ref = make_reference(g);
  // Adversarial reads for tier SL, planted before the index is built:
  //  (a) the read has a mismatch in most parts at its true locus (so few or
  //      no enumerated parts occur exactly there), and an exact copy of a
  //      piece of it is planted elsewhere as a slightly worse decoy, so the
  //      certificate threshold decides between the two;
  //  (b) a long deletion or insertion near a read end, so the optimum needs
  //      a wide band and competes with clipping.
  struct Planted { std::vector<uint8_t> seq; };
  std::vector<Planted> planted;
  {
    const int L = 150, P = part_count(L, 14, 4);
    auto& R = ref.seq.owned;
    for (int n = 0; n < 100; ++n) {
      const int64_t a = static_cast<int64_t>(ref.offsets[0]) + 1000 + static_cast<int64_t>(g() % 8000);
      std::vector<uint8_t> s(R.begin() + a, R.begin() + a + L);
      // One mismatch in each of `mut` random parts: the true alignment then
      // has exactly P - mut exact parts and scores L - 5 mut.
      const int mut = 2 + static_cast<int>(g() % (P - 1));
      int order[PMAX];
      for (int j = 0; j < P; ++j) order[j] = j;
      for (int j = P - 1; j > 0; --j) std::swap(order[j], order[g() % (j + 1)]);
      for (int y = 0; y < mut; ++y) {
        int off, len;
        part_geometry(L, P, order[y], &off, &len);
        const int x = off + 2 + static_cast<int>(g() % (len - 4));
        s[x] = static_cast<uint8_t>((s[x] + 1 + g() % 3) % 4);
      }
      // Decoy: an exact copy of a piece, scoring (y0 - x0) - 5 per clip.
      const int x0 = g() % 2 ? 0 : static_cast<int>(g() % 30), y0 = L - 5 - static_cast<int>(g() % 40);
      // Free stretches: chr1 32000-50000 and chr2 24000-40000.
      const int64_t d = n < 60 ? static_cast<int64_t>(ref.offsets[0]) + 32000 + n * 240
                               : static_cast<int64_t>(ref.offsets[1]) + 24000 + (n - 60) * 240;
      std::copy(s.begin() + x0, s.begin() + y0, R.begin() + d);
      planted.push_back({s});
    }
    for (int n = 0; n < 120; ++n) {
      const int c = static_cast<int>(g() % 2);
      const int64_t a = static_cast<int64_t>(ref.offsets[c]) + 1000 + static_cast<int64_t>(g() % 3000);
      // Gap lengths around the band limit, near a read end, where clipping
      // the short side competes with the gapped alignment.
      const int G = 8 + static_cast<int>(g() % 20);
      const int at = n % 2 ? 12 + static_cast<int>(g() % 16) : L - 12 - static_cast<int>(g() % 16);
      std::vector<uint8_t> s(R.begin() + a, R.begin() + a + at);
      if (n % 4 < 2) {  // deletion of G reference bases
        s.insert(s.end(), R.begin() + a + at + G, R.begin() + a + L + G);
      } else {          // insertion of G random bases
        std::vector<uint8_t> ins = random_seq(g, G);
        s.insert(s.end(), ins.begin(), ins.end());
        s.insert(s.end(), R.begin() + a + at, R.begin() + a + at + std::max(0, L - at - G));
      }
      s.resize(L);
      planted.push_back({s});
    }
  }
  // Adversarial pairs for tier PR: at the true locus mate 1 has k1
  // mismatches (in distinct parts) and mate 2 is exact; at a planted decoy
  // mate 1 is exact and mate 2 has kd mismatches, at a proper distance. The
  // two pair scores are close, either way round, so a loose bound or an early
  // stop picks the wrong pair.
  std::vector<std::pair<std::vector<uint8_t>, std::vector<uint8_t>>> decoy_pairs;
  {
    const int L = 150, P = part_count(L, 14, 4);
    auto& R = ref.seq.owned;
    auto rcv = [](std::vector<uint8_t> v) {
      std::reverse(v.begin(), v.end());
      for (auto& c : v) c = c < 4 ? static_cast<uint8_t>(3 - c) : c;
      return v;
    };
    for (int n = 0; n < 40; ++n) {
      const int F = 300 + static_cast<int>(g() % 150), Fd = 300 + static_cast<int>(g() % 100);
      const int64_t T = static_cast<int64_t>(ref.offsets[1]) + 1000 + static_cast<int64_t>(g() % 3000);
      std::vector<uint8_t> m1(R.begin() + T, R.begin() + T + L);
      const std::vector<uint8_t> fwd2(R.begin() + T + F - L, R.begin() + T + F);
      const int k1 = 1 + static_cast<int>(g() % 4), kd = 1 + static_cast<int>(g() % 6);
      int order[PMAX];
      for (int j = 0; j < P; ++j) order[j] = j;
      for (int j = P - 1; j > 0; --j) std::swap(order[j], order[g() % (j + 1)]);
      for (int y = 0; y < k1; ++y) {
        int off, len;
        part_geometry(L, P, order[y], &off, &len);
        const int x = off + 2 + static_cast<int>(g() % (len - 4));
        m1[x] = static_cast<uint8_t>((m1[x] + 1 + g() % 3) % 4);
      }
      // Decoy in chr1 12500-29500 (free): mate 1 exact, mate 2 with kd mismatches.
      const int64_t D = static_cast<int64_t>(ref.offsets[0]) + 12500 + 420 * n;
      std::copy(m1.begin(), m1.end(), R.begin() + D);
      std::vector<uint8_t> d2 = fwd2;
      for (int y = 0; y < kd; ++y) {
        const int x = static_cast<int>((y * L) / kd + g() % (L / kd));
        d2[x] = static_cast<uint8_t>((d2[x] + 1 + g() % 3) % 4);
      }
      std::copy(d2.begin(), d2.end(), R.begin() + D + Fd - L);
      decoy_pairs.push_back({m1, rcv(fwd2)});
    }
  }
  // Adversarial reads for the interval bounds (Aim 3), planted before the
  // index is built:
  //  (c) Lemma E tight: the read has 2 mismatches at its true locus (chr2
  //      1000-4000), and a decoy (chr2 34000-39850) equals the read except
  //      for its last 3 bases. The decoy's clipped alignment scores L - 8,
  //      above the true locus's L - 10, yet lies 3 edits away: outside
  //      radius 2, so only edit_bound(L, 3) = L - 8 covers it;
  //  (d) the own chain: a tandem array (16 copies of a 20-base unit, a
  //      period below the chain pad) in unique sequence (chr1 77000-83300).
  //      The read covers one end of the array (145 bases) and 5 flank
  //      bases; the copy of its tandem part 160 bases inside the array lies
  //      in the same chain, does not overlap the read's alignment, and
  //      scores L - 10 (above lo - 20, so tier R's SUB search must find
  //      it). Even n: the array's right end (copy to the left); odd n: its
  //      left end (copy to the right);
  //  (e) the SL floor (Aim 5): the read occurs exactly at its true locus
  //      (chr1 1000-9000), and a decoy (chr1 84000-89500) has one exact part
  //      and a mismatch in each of the other 7. With support 2 the decoy
  //      scores exactly the floor L - 35 and lies in no chain above it, so
  //      only the floor bounds SUB there;
  //  (f-i) tier SL at its other thresholds (local-only, 8 parts; support 2:
  //      floor L - 35, pad 28; support 3: floor L - 30, pad 23), on chr2
  //      7000-19000:
  //    f  the band pad: the read has a 27-bp deletion at a part boundary
  //       (117, gap 27 <= 28); a decoy clips 29 bases (116);
  //    g  the chain bound: the read is exact (150); one locus has a 1-bp
  //       deletion between parts (143, bound 150), another one mismatch (145,
  //       bound 145): SUB = 145 rests on evaluating a chain whose bound is
  //       2 above s2;
  //    h  chain splitting (support 3): a 12-bp deletion at a part boundary
  //       (132) against two decoys clipping 15 bases (130 each);
  //    i  the skip test: the read has one mismatch at its locus (145, bound
  //       145) against two loci with a 1-bp deletion (143 each, bound 150);
  //  (j) tier R, Theorem G at the read end (Aim 5): a 108-base element in 20
  //      copies (chr1 52000-58000), each after its own flank. The read is the
  //      last copy's 42-base flank, with mismatches at offsets 0, 20 and 34,
  //      then the element: three absent parts, then a stretch occurring 20
  //      times (> tau) up to the read end, which ends the partition. OPT is
  //      135, exactly U; the sampled copies score ~103. Charging the element
  //      as a fourth part would give U = 130.
  std::vector<std::vector<uint8_t>> bound_reads;
  {
    auto& R = ref.seq.owned;
    const int L = 150;
    for (int n = 0; n < 20; ++n) {
      const int64_t a = static_cast<int64_t>(ref.offsets[1]) + 1000 + static_cast<int64_t>(g() % 3000);
      std::vector<uint8_t> s(R.begin() + a, R.begin() + a + L);
      for (int y = 0; y < 2; ++y) {
        const int x = 10 + 60 * y + static_cast<int>(g() % 40);
        s[x] = static_cast<uint8_t>((s[x] + 1 + g() % 3) % 4);
      }
      std::vector<uint8_t> d = s;
      for (int x = L - 3; x < L; ++x) d[x] = static_cast<uint8_t>((s[x] + 1 + g() % 3) % 4);
      std::copy(d.begin(), d.end(), R.begin() + static_cast<int64_t>(ref.offsets[1]) + 34000 + 300 * n);
      bound_reads.push_back(s);
    }
    for (int n = 0; n < 10; ++n) {
      const int64_t at = static_cast<int64_t>(ref.offsets[0]) + 77000 + 700 * n;
      const std::vector<uint8_t> unit = random_seq(g, 20);
      for (int c = 0; c < 16; ++c) std::copy(unit.begin(), unit.end(), R.begin() + at + 20 * c);
      const int64_t from = n % 2 ? at - 5 : at + 175;
      bound_reads.push_back(std::vector<uint8_t>(R.begin() + from, R.begin() + from + L));
    }
    for (int n = 0; n < 10; ++n) {
      const int64_t a = static_cast<int64_t>(ref.offsets[0]) + 1000 + 800 * n;
      std::vector<uint8_t> s(R.begin() + a, R.begin() + a + L);
      std::vector<uint8_t> d = s;
      const int P = part_count(L, 14, 4), keep = static_cast<int>(g() % P);
      for (int j = 0; j < P; ++j) {
        if (j == keep) continue;
        int off, len;
        part_geometry(L, P, j, &off, &len);
        const int x = off + len / 2;
        d[x] = static_cast<uint8_t>((d[x] + 1 + g() % 3) % 4);
      }
      std::copy(d.begin(), d.end(), R.begin() + static_cast<int64_t>(ref.offsets[0]) + 84000 + 550 * n);
      bound_reads.push_back(s);
    }
    const int64_t c2 = static_cast<int64_t>(ref.offsets[1]);
    auto from = [&](int64_t a, int64_t b) { return std::vector<uint8_t>(R.begin() + a, R.begin() + b); };
    auto put = [&](int64_t at, const std::vector<uint8_t>& v) { std::copy(v.begin(), v.end(), R.begin() + at); };
    auto with_extra = [&](const std::vector<uint8_t>& v, int at) {  // one more reference base at read offset `at`
      std::vector<uint8_t> o(v.begin(), v.begin() + at);
      o.push_back(rand_base(g));
      o.insert(o.end(), v.begin() + at, v.end());
      return o;
    };
    auto changed = [&](std::vector<uint8_t> v, int at) {
      v[at] = static_cast<uint8_t>((v[at] + 1 + g() % 3) % 4);
      return v;
    };
    for (int n = 0; n < 6; ++n) {  // f
      const int64_t a = c2 + 7000 + 300 * n;
      std::vector<uint8_t> s = from(a, a + 76), tail = from(a + 103, a + 177);
      s.insert(s.end(), tail.begin(), tail.end());
      std::vector<uint8_t> d(s.begin(), s.begin() + 121);
      const std::vector<uint8_t> rnd = random_seq(g, 29);
      d.insert(d.end(), rnd.begin(), rnd.end());
      put(c2 + 8800 + 160 * n, d);
      bound_reads.push_back(s);
    }
    for (int n = 0; n < 6; ++n) {  // g
      const std::vector<uint8_t> s = from(c2 + 9800 + 160 * n, c2 + 9950 + 160 * n);
      put(c2 + 10800 + 320 * n, with_extra(s, 76));
      put(c2 + 10960 + 320 * n, changed(s, 47));
      bound_reads.push_back(s);
    }
    for (int n = 0; n < 6; ++n) {  // h
      const int64_t a = c2 + 12800 + 180 * n;
      std::vector<uint8_t> s = from(a, a + 76), tail = from(a + 88, a + 162);
      s.insert(s.end(), tail.begin(), tail.end());
      for (int k = 0; k < 2; ++k) {
        std::vector<uint8_t> d(s.begin(), s.begin() + 135);
        const std::vector<uint8_t> rnd = random_seq(g, 15);
        d.insert(d.end(), rnd.begin(), rnd.end());
        put(c2 + 13900 + 320 * n + 160 * k, d);
      }
      bound_reads.push_back(s);
    }
    for (int n = 0; n < 6; ++n) {  // i
      const std::vector<uint8_t> loc = from(c2 + 15900 + 160 * n, c2 + 16050 + 160 * n);
      const std::vector<uint8_t> s = changed(loc, 47);
      put(c2 + 16900 + 320 * n, with_extra(s, 76));
      put(c2 + 17060 + 320 * n, with_extra(s, 76));
      bound_reads.push_back(s);
    }
    for (int n = 0; n < 2; ++n) {  // j
      const std::vector<uint8_t> elem = random_seq(g, 108);
      int64_t last = 0;
      for (int c = 0; c < 20; ++c) {
        last = static_cast<int64_t>(ref.offsets[0]) + 52000 + 300 * c + 150 * n;
        put(last, elem);
      }
      std::vector<uint8_t> s = from(last - 42, last);
      for (int x : {0, 20, 34}) s = changed(s, x);
      s.insert(s.end(), elem.begin(), elem.end());
      bound_reads.push_back(s);
    }
  }
  Index ix = Index::build(ref, 14, 4, 4);
  const IndexView view = ix.view(ref);
  // Round trip through the on-disk format: with q = 14 the file stores only
  // 16-bit key suffixes and is memory-mapped; results must not change.
  // Per process: parallel runs in one directory must not share the file.
  const std::string tmp = "test_certa." + std::to_string(getpid()) + ".tmp.cidx";
  save_index(tmp, ref, ix);
  Reference ref_m;
  Index ix_m;
  load_index(tmp, ref_m, ix_m);
  std::remove(tmp.c_str());  // the mapping stays valid until ix_m is gone
  const IndexView view_m = ix_m.view(ref_m);
  CHECK(view_m.keys16 != nullptr, "saved index should use 16-bit keys (q=14)");
  // The file stores the directory as 64 bits, memory holds it as 32.
  CHECK(ix_m.dir.size() == ix.dir.size() &&
            std::equal(ix.dir.data(), ix.dir.data() + ix.dir.size(), ix_m.dir.data()),
        "directory changed by the round trip");
  {
    // A directory offset that 32 bits cannot hold is rejected on load.
    save_index(tmp, ref, ix);
    {
      std::fstream f(tmp, std::ios::in | std::ios::out | std::ios::binary);
      f.seekp(-8, std::ios::end);  // the directory's last entry ends the file
      const uint64_t big = (uint64_t(1) << 32) + ix.pos.size();
      f.write(reinterpret_cast<const char*>(&big), sizeof big);
    }
    Reference r;
    Index x;
    bool thrown = false;
    try {
      load_index(tmp, r, x);
    } catch (const std::runtime_error&) {
      thrown = true;
    }
    std::remove(tmp.c_str());
    CHECK(thrown, "a directory offset beyond 2^32 - 1 accepted");
  }

  // Reads: from unique sequence, repeat copies, the exact duplicate, the
  // tandem repeat and near the N run; with 0..7 edits; plus random reads.
  struct Sim { std::vector<uint8_t> seq; int edits; };
  std::vector<Sim> reads;
  const int64_t anchors[] = {-1, -1, -1, 10000, 30000, 50000, 75000, 59900, 60250, 66000};
  for (int n = 0; n < 480; ++n) {
    int L = (n % 4 == 0) ? 110 : 150;
    std::vector<uint8_t> s;
    int edits = static_cast<int>(g() % 8);
    if (n % 25 == 0) {
      s = random_seq(g, L);
      edits = -1;
    } else {
      int c = static_cast<int>(g() % 2);
      int64_t a = anchors[g() % 10];
      int64_t start = a >= 0 && c == 0 ? a + static_cast<int64_t>(g() % 1500)
                                       : static_cast<int64_t>(g() % (ref.lengths[c] - L - 20));
      start += static_cast<int64_t>(ref.offsets[c]);
      s.assign(ref.seq.owned.begin() + start, ref.seq.owned.begin() + start + L + 10);
      mutate(g, s, edits);
      s.resize(L);
      if (g() % 2) {  // reverse strand
        std::vector<uint8_t> rc(L);
        for (int i = 0; i < L; ++i) rc[i] = s[L - 1 - i] < 4 ? 3 - s[L - 1 - i] : 4;
        s.swap(rc);
      }
    }
    if (n % 40 == 7) s[g() % L] = 4;  // an N in the read
    reads.push_back({s, edits});
  }
  for (auto& pl : planted) {
    std::vector<uint8_t> s = pl.seq;
    if (g() % 2) {
      const int L = static_cast<int>(s.size());
      std::vector<uint8_t> rc(L);
      for (int i = 0; i < L; ++i) rc[i] = s[L - 1 - i] < 4 ? 3 - s[L - 1 - i] : 4;
      s.swap(rc);
    }
    reads.push_back({s, -1});
  }
  // Reads for tier SL: an adapter-like random tail or head, or a chimera of
  // two loci, with 0..3 further edits.
  for (int n = 0; n < 120; ++n) {
    const int L = 150, cut = 2 + static_cast<int>(g() % 40);
    auto take = [&](int len) {
      int c = static_cast<int>(g() % 2);
      int64_t start = static_cast<int64_t>(g() % (ref.lengths[c] - len - 20)) + static_cast<int64_t>(ref.offsets[c]);
      if (n % 6 == 5) start = static_cast<int64_t>(ref.offsets[0]) + 10000 + static_cast<int64_t>(g() % 1500);
      return std::vector<uint8_t>(ref.seq.owned.begin() + start, ref.seq.owned.begin() + start + len);
    };
    std::vector<uint8_t> s = take(L + 10);
    mutate(g, s, static_cast<int>(g() % 4));
    s.resize(L);
    std::vector<uint8_t> other = n % 3 == 2 ? take(cut) : random_seq(g, cut);
    if (n % 2) std::copy(other.begin(), other.end(), s.begin());   // head
    else std::copy(other.begin(), other.end(), s.end() - cut);     // tail
    if (g() % 2) {
      std::vector<uint8_t> rc(L);
      for (int i = 0; i < L; ++i) rc[i] = s[L - 1 - i] < 4 ? 3 - s[L - 1 - i] : 4;
      s.swap(rc);
    }
    reads.push_back({s, -1});
  }
  // Reads for tier R: from inside copies of the 50-copy element (some copies
  // diverged), with 0-3 edits; the repeat tier must bound them exactly.
  const size_t repeat_first = reads.size();
  for (int n = 0; n < 60; ++n) {
    const int L = 150;
    const int64_t start = static_cast<int64_t>(ref.offsets[0]) + 60200 + 350 * static_cast<int64_t>(g() % 40) +
                          static_cast<int64_t>(g() % 150);
    std::vector<uint8_t> s(ref.seq.owned.begin() + start, ref.seq.owned.begin() + start + L + 10);
    mutate(g, s, static_cast<int>(g() % 4));
    s.resize(L);
    if (g() % 2) {
      std::vector<uint8_t> rc(L);
      for (int i = 0; i < L; ++i) rc[i] = s[L - 1 - i] < 4 ? 3 - s[L - 1 - i] : 4;
      s.swap(rc);
    }
    reads.push_back({s, -1});
  }

  // The planted interval-bound reads, on both strands; always checked.
  const size_t bounds_first = reads.size();
  for (const auto& b : bound_reads) {
    std::vector<uint8_t> s = b;
    if (g() % 2) {
      const int L = static_cast<int>(s.size());
      std::vector<uint8_t> rc(L);
      for (int i = 0; i < L; ++i) rc[i] = s[L - 1 - i] < 4 ? 3 - s[L - 1 - i] : 4;
      s.swap(rc);
    }
    reads.push_back({s, -1});
  }
  auto planted_bound = [&](size_t i) { return i >= bounds_first; };

  std::fprintf(stderr, "computing oracle for %zu reads over %zu bases...\n", reads.size(), ref.seq.size());
  std::vector<std::vector<Hit>> truth;
  for (auto& r : reads) truth.push_back(oracle(ref.seq.owned, r.seq.data(), static_cast<int>(r.seq.size()), KMAX));

  ReadBatch batch;
  for (size_t i = 0; i < reads.size(); ++i) {
    batch.offs.push_back(batch.codes.size());
    batch.codes.insert(batch.codes.end(), reads[i].seq.begin(), reads[i].seq.end());
    batch.lens.push_back(static_cast<uint16_t>(reads[i].seq.size()));
    batch.hashes.push_back(i * 0x9E3779B97F4A7C15ULL);
  }

  // {k, budget, reverse_order, s2_limit}
  const Params configs[] = {{2, 256}, {3, 256}, {5, 256}, {4, 64}, {4, 16}, {3, 4}, {2, 1},
                            {3, 16, 0, 8}, {5, 256, 0, 10}, {2, 4, 0, 6},
                            {5, 2048}, {4, 4096},   // host-only budgets (CPU pass 3)
                            // {k, budget, reverse, s2, mapq_limit, min_support}: q-gram filter
                            {3, 256, 0, 0, 0, 2}, {4, 4096, 0, 0, 0, 2}, {2, 16, 0, 0, 0, 2}, {2, 64, 0, 0, 0, 3},
                            // {.., min_support, local, local_only}: tier SL
                            {5, 256, 0, 0, 8, 1, 2}, {5, 4096, 0, 0, 8, 1, 1}, {3, 16, 0, 0, 0, 1, 2},
                            {5, 256, 0, 0, 8, 1, 2, 1}, {5, 4096, 0, 0, 0, 1, 3, 1}};
  std::vector<int> local_truth(reads.size(), INT32_MIN);
  int local_checked = 0;
  // Interval output of certified reads: [XL, XU] and SUB_hi against the
  // oracles, for every third read (SUB oracles cached by reported locus).
  std::map<std::tuple<size_t, int, int64_t, int64_t>, int> sub_cache;
  int bounds_checked = 0, sub_checked = 0, sub_nontrivial = 0;
  auto oracle_opt = [&](size_t i) {
    if (local_truth[i] == INT32_MIN) {
      local_truth[i] = oracle_local(ref.seq.owned, reads[i].seq.data(), static_cast<int>(reads[i].seq.size()));
      ++local_checked;
    }
    return local_truth[i];
  };
  auto oracle_sub_at = [&](size_t i, const Result& r) {
    const int64_t span = ref_span_of(r);
    const auto key = std::make_tuple(i, static_cast<int>(r.strand), r.ref_pos, span);
    auto it = sub_cache.find(key);
    if (it != sub_cache.end()) return it->second;
    const int v = oracle_sub(ref.seq.owned, reads[i].seq, r.strand, r.ref_pos, span);
    sub_cache[key] = v;
    return v;
  };
  // Checker on the certified results (index q = 14, s = 4).
  const check::Checker checker4(ref, ix);
  int claims_accepted = 0, claims_tampered = 0, tamper_caught = 0;
  int caught_by[check::kNumVerdicts] = {};  // rejections of falsified claims, by verdict
  auto check_claims = [&](size_t i, const Result& r, const Params& p) {
    if (r.tier != kTierS0 && r.tier != kTierS1 && r.tier != kTierSL) return;
    const int L = static_cast<int>(reads[i].seq.size());
    const Placement pl = place_single(r, 60, nullptr, L);
    check::Claim c = claim_of(pl, reads[i].seq);
    check::Work w;
    const check::Verdict v = checker4.check(c, &w);
    CHECK(v == check::kAccepted, "k=%d budget=%d local=%d read %zu tier %d: checker says %s (XL %d XU %d XH %d XQ %d)",
          p.k, p.budget, p.local, i, r.tier, check::kVerdictNames[v], c.xl, c.xu, c.xh, c.xq);
    claims_accepted += v == check::kAccepted;
    // Falsified claims must be rejected: XU below OPT, XH below SUB, AS off
    // by one, the alignment moved, the witness dropped.
    const int opt = oracle_opt(i);
    auto expect_reject = [&](check::Claim t, const char* what) {
      check::Work w2;
      const check::Verdict v2 = checker4.check(t, &w2);
      ++claims_tampered;
      tamper_caught += v2 != check::kAccepted;
      ++caught_by[v2];
      CHECK(v2 != check::kAccepted, "read %zu tier %d: false claim accepted (%s)", i, r.tier, what);
    };
    { check::Claim t = c; t.xu = opt - 1; expect_reject(t, "XU = OPT - 1"); }
    {  // the same with no alignment shown: only the bound search can refute it
      check::Claim t = c;
      t.mapped = false;
      t.xl = t.xh = check::kNone;
      t.xq = 0;
      t.has_alt = false;
      t.xu = opt - 1;
      expect_reject(t, "XU = OPT - 1, unmapped");
    }
    if (c.xh != check::kNone) { check::Claim t = c; t.xq = c.xq + 5; t.mapq = 99; expect_reject(t, "XQ + 5"); }
    { check::Claim t = c; t.as = c.as + 1; expect_reject(t, "AS + 1"); }
    { check::Claim t = c; t.aln.pos += 1; expect_reject(t, "alignment moved by one base"); }
    if (c.xu < L * kMatch) { check::Claim t = c; t.t = 0; expect_reject(t, "no witness"); }
    if (i % 6 == 0 || planted_bound(i)) {
      const int sub = oracle_sub_at(i, r);
      if (sub >= 19 && sub < c.as) {  // a claim XH = SUB - 1 with a positive XQ is false
        check::Claim t = c;
        t.xh = sub - 1;
        t.xq = 1;
        t.mapq = 60;
        expect_reject(t, "XH = SUB - 1");
      }
    }
  };
  auto check_certified_bounds = [&](size_t i, const Result& r, const Params& p) {
    const int L = static_cast<int>(reads[i].seq.size());
    int xl, xu, sh;
    certified_bounds(r, L, &xl, &xu, &sh);
    const int opt = oracle_opt(i);
    ++bounds_checked;
    CHECK(xl <= opt && opt <= xu, "k=%d budget=%d local=%d read %zu tier %d: OPT %d outside [%d, %d]", p.k,
          p.budget, p.local, i, r.tier, opt, xl, xu);
    CHECK(sh != kNoSub, "read %zu tier %d: certified without a SUB bound", i, r.tier);
    if (sh >= opt) return;  // SUB <= OPT: nothing to check
    ++sub_nontrivial;
    if (i % 6 != 0 && !planted_bound(i)) return;  // oracle cost: half of them
    const int sub = oracle_sub_at(i, r);
    ++sub_checked;
    CHECK(sub <= sh, "k=%d budget=%d local=%d read %zu tier %d: SUB %d above SUB_hi %d (score %d, d1 %d d2 %d R %d)",
          p.k, p.budget, p.local, i, r.tier, sub, sh, r.score, r.d1, r.d2, r.radius);
  };
  for (const Params& p : configs) {
    std::vector<Result> res, res1, resrev, resm;
    map_cpu(view, p, batch, res, 4);
    map_cpu(view, p, batch, res1, 1);
    map_cpu(view_m, p, batch, resm, 4);
    for (size_t i = 0; i < reads.size(); ++i)
      CHECK(std::memcmp(&res[i], &resm[i], sizeof(Result)) == 0,
            "read %zu differs between in-memory and memory-mapped (16-bit key) index", i);
    Params prev = p;
    prev.reverse_order = 1;
    map_cpu(view, prev, batch, resrev, 4);
    int certified = 0, easy = 0, easy_cert = 0, sr = 0, s2 = 0, sl = 0;
    int reasons[kNumReasons] = {};
    for (size_t i = 0; i < reads.size(); ++i) {
      const Result& r = res[i];
      CHECK(std::memcmp(&r, &res1[i], sizeof r) == 0, "read %zu differs between 1 and 4 threads", i);
      // Metamorphic: verifying clusters in reverse order must give the same
      // certificate (d1; d2 whenever the best locus is unique).
      const Result& rv = resrev[i];
      CHECK(r.tier != kTierSL || (rv.tier == kTierSL && r.score == rv.score),
            "local=%d read %zu: order-dependent SL score (%d vs %d)", p.local, i, r.score, rv.score);
      CHECK(r.certified == rv.certified && r.d1 == rv.d1 && r.radius == rv.radius,
            "k=%d budget=%d read %zu: order-dependent d1 (%d vs %d)", p.k, p.budget, i, r.d1, rv.d1);
      if (r.certified && r.n_best == 1 && rv.n_best == 1)
        CHECK(r.d2 == rv.d2, "k=%d budget=%d read %zu: order-dependent d2 (%d vs %d)",
              p.k, p.budget, i, r.d2, rv.d2);
      CHECK((r.n_best > 1) == (rv.n_best > 1), "k=%d budget=%d read %zu: order-dependent tie",
            p.k, p.budget, i);
      ++reasons[r.reason];
      int omin = 99;
      for (auto& h : truth[i]) omin = std::min(omin, h.dist);
      if (reads[i].edits >= 0 && reads[i].edits <= p.k) {
        ++easy;
        easy_cert += r.certified == 1;
      }
      if (r.reason == kNotFound && !p.local_only)
        CHECK(omin > r.radius, "k=%d budget=%d read %zu: lossless violated, oracle %d <= R %d",
              p.k, p.budget, i, omin, r.radius);
      if (!r.certified) continue;
      if (r.certified == 2) {  // S2: sound, not below the optimum, nothing within R
        ++s2;
        CHECK(r.tier == kTierS2 && r.d1 <= p.s2_limit, "read %zu: bad S2 result", i);
        CHECK(rescore(view, r, reads[i].seq) == r.nm, "S2 read %zu: CIGAR rescores to %d, NM %d",
              i, rescore(view, r, reads[i].seq), r.nm);
        CHECK(omin > r.radius, "S2 read %zu: a locus within R=%d exists (oracle %d)", i, r.radius, omin);
        CHECK(r.d1 >= (omin <= KMAX ? omin : KMAX + 1), "S2 read %zu: d1 %d below oracle %d", i, r.d1, omin);
        continue;
      }
      ++certified;
      if (i % 3 == 0 || planted_bound(i)) {
        check_certified_bounds(i, r, p);
        if (i % 6 == 0 || planted_bound(i)) check_claims(i, r, p);
      }
      if (r.tier == kTierSL) {  // the certified maximum local score
        ++sl;
        CHECK(rescore(view, r, reads[i].seq) == r.nm, "SL read %zu: CIGAR rescores to %d edits, NM %d",
              i, rescore(view, r, reads[i].seq), r.nm);
        CHECK(rescore_affine(view, r, reads[i].seq) == r.score, "SL read %zu: CIGAR scores %d, AS %d",
              i, rescore_affine(view, r, reads[i].seq), r.score);
        CHECK(r.score > r.floor, "SL read %zu: score %d not above floor %d", i, r.score, r.floor);
        if (local_truth[i] == INT32_MIN) {
          local_truth[i] = oracle_local(ref.seq.owned, reads[i].seq.data(), static_cast<int>(reads[i].seq.size()));
          ++local_checked;
        }
        CHECK(r.score == local_truth[i], "local=%d budget=%d read %zu: SL score %d != oracle %d",
              p.local, p.budget, i, r.score, local_truth[i]);
        continue;
      }
      // The reported alignment (affine, possibly clipped) must be what NM says;
      // unclipped, it cannot have fewer edits than the certified minimum d1.
      CHECK(rescore(view, r, reads[i].seq) == r.nm,
            "k=%d read %zu: CIGAR rescores to %d, NM %d", p.k, i,
            rescore(view, r, reads[i].seq), r.nm);
      bool clipped = false;
      for (int c = 0; c < r.n_cigar; ++c) clipped |= (r.cigar[c] & 0xF) == kOpS;
      CHECK(clipped || r.nm >= r.d1, "read %zu: unclipped NM %d below certified d1 %d", i, r.nm, r.d1);
      CHECK(r.d1 == omin, "k=%d budget=%d read %zu: d1 %d != oracle %d", p.k, p.budget, i, r.d1, omin);
      if (r.tier == kTierSR) {  // must be backed by >= 2 distinct exact loci
        ++sr;
        int exact = 0;
        for (auto& h : truth[i]) exact += h.dist == 0;
        CHECK(exact >= 2, "k=%d budget=%d read %zu: SR claims a repeat, oracle has %d exact loci",
              p.k, p.budget, i, exact);
      }
      if (r.n_best == 1 && r.d2 < 0) {
        int64_t span = 0;
        for (int c = 0; c < r.n_cigar; ++c)
          if ((r.cigar[c] & 0xF) == kOpM || (r.cigar[c] & 0xF) == kOpD) span += r.cigar[c] >> 4;
        const int64_t end = r.ref_pos + span - 1;
        for (auto& h : truth[i]) {
          if (h.dist > r.radius) continue;
          bool near = h.strand == r.strand && (h.end - end <= 2 * BMAX && end - h.end <= 2 * BMAX);
          CHECK(near, "k=%d read %zu: oracle locus (strand %d end %lld d %d) missed; chosen end %lld",
                p.k, i, h.strand, static_cast<long long>(h.end), h.dist, static_cast<long long>(end));
        }
      }
    }
    std::fprintf(stderr,
                 "k=%d budget=%-3d s2=%-2d local=%d%s certified %3d/%zu (SR %d, SL %d, S2 placed %d); reads with <=k edits certified %d/%d; "
                 "uncertified: length %d, radius<0 %d, not-found %d, wide %d\n",
                 p.k, p.budget, p.s2_limit, p.local, p.local_only ? "o" : " ", certified, reads.size(), sr, sl, s2,
                 easy_cert, easy, reasons[kBadLength],
                 reasons[kRadiusNegative], reasons[kNotFound], reasons[kClusterTooWide]);
  }
  std::fprintf(stderr, "SL scores checked against the local oracle for %d reads\n", local_checked);
  std::fprintf(stderr, "certified [XL, XU] checked %d times; SUB bounds below OPT %d, checked against the SUB oracle %d\n",
               bounds_checked, sub_nontrivial, sub_checked);
  std::fprintf(stderr, "checker: %d certified records accepted; %d of %d falsified claims rejected (", claims_accepted,
               tamper_caught, claims_tampered);
  for (int v = 1; v < check::kNumVerdicts; ++v)
    std::fprintf(stderr, "%s%s %d", v > 1 ? ", " : "", check::kVerdictNames[v], caught_by[v]);
  std::fprintf(stderr, ")\n");

  // ---- band_del16 (vectorized where the compiler allows) must equal its
  // scalar reference on random rows, sentinels included.
  {
    std::mt19937_64 gd(seed * 3 + 1);
    int rows = 0, bad = 0;
    for (int it = 0; it < 20000; ++it) {
      const int B = 1 + static_cast<int>(gd() % BMAX);
      std::vector<int16_t> hq(B), h1(B), h2(B);
      for (auto& x : hq) x = gd() % 6 == 0 ? (int16_t)-16000 : (int16_t)((int)(gd() % 400) - 200);
      band_del16(B, hq.data(), h1.data(), (int16_t)-16000);
      band_del16_scalar(B, hq.data(), h2.data(), (int16_t)-16000);
      ++rows;
      bad += h1 != h2;
    }
    CHECK(bad == 0, "band_del16 differs from its scalar reference on %d of %d rows", bad, rows);
    std::fprintf(stderr, "band_del16 = band_del16_scalar on %d random rows: %s\n", rows, bad ? "no" : "yes");
  }
#ifdef CERTA_AFFINE_VEC
  // ---- affine_fill_vec must equal the scalar fill cell for cell: every
  // traceback byte (sentinel cells included) and the best end cell, on
  // random bands, start ranges and reads with N.
  {
    std::unique_ptr<Workspace> w1(new Workspace), w2(new Workspace);
    std::mt19937_64 ga(seed * 5 + 3);
    int bands = 0, bad = 0;
    for (int it = 0; it < 4000; ++it) {
      const int B = 1 + static_cast<int>(ga() % BMAX);
      const int64_t lo = static_cast<int64_t>(ga() % (ref.seq.size() - 600));
      std::vector<uint8_t> rd = reads[ga() % reads.size()].seq;
      if (it % 4 >= 2) {  // a read from this very place, so scores are high
        rd.assign(ref.seq.owned.begin() + lo + B / 2, ref.seq.owned.begin() + lo + B / 2 + 160);
        mutate(ga, rd, static_cast<int>(ga() % 8));
        rd.resize(150, 0);
      }
      if (it % 8 == 1) rd[ga() % rd.size()] = 4;  // an N in the read
      const int L = static_cast<int>(rd.size());
      int64_t smin = INT64_MIN, smax = INT64_MAX;
      if (it % 2) {
        smin = lo + static_cast<int64_t>(ga() % 200) - 50;
        smax = smin + static_cast<int64_t>(ga() % 120);
      }
      if (lo + L - 1 + (B + 3) / 4 * 4 > static_cast<int64_t>(view.ref_len)) continue;
      int i1 = -1, b1 = -1, i2 = -2, b2 = -2;
      const int s1 = affine_fill_vec(view, rd.data(), L, lo, B, *w1, smin, smax, &i1, &b1);
      const int s2 = affine_fill_scalar(view, rd.data(), L, lo, B, *w2, smin, smax, &i2, &b2);
      bool same = s1 == s2 && i1 == i2 && b1 == b2 &&
                  std::memcmp(w1->tb + B, w2->tb + B, static_cast<size_t>(L) * B) == 0;
      if (lo + L - 1 + (B + 7) / 8 * 8 <= static_cast<int64_t>(view.ref_len)) {  // the 16-bit fill too
        int i3 = -3, b3 = -3;
        const int s3 = affine_fill_vec16(view, rd.data(), L, lo, B, *w1, smin, smax, &i3, &b3);
        same = same && s3 == s2 && i3 == i2 && b3 == b2 &&
               std::memcmp(w1->tb + B, w2->tb + B, static_cast<size_t>(L) * B) == 0;
      }
      ++bands;
      bad += !same;
    }
    CHECK(bad == 0, "affine_fill_vec differs from affine_fill_scalar on %d of %d bands", bad, bands);
    std::fprintf(stderr, "affine_fill_vec (32- and 16-bit) = affine_fill_scalar (every traceback byte) on %d random bands: %s\n",
                 bands, bad ? "no" : "yes");
  }
#endif
  // ---- repeat_detail::match_run must count the leading same_base positions
  // as the scalar loop does: reads copied from the reference with a few
  // mismatches and N, windows reaching past either end.
  {
    std::mt19937_64 gm(seed * 11 + 7);
    const int64_t n = static_cast<int64_t>(view.ref_len);
    int runs = 0, bad = 0;
    for (int it = 0; it < 20000; ++it) {
      const int len = 1 + static_cast<int>(gm() % 200);
      int64_t j = static_cast<int64_t>(gm() % static_cast<uint64_t>(n));
      if (it % 10 == 0) j = n - static_cast<int64_t>(gm() % 60);
      if (it % 10 == 1) j = -static_cast<int64_t>(gm() % 20);
      std::vector<uint8_t> rd(len);
      for (int i = 0; i < len; ++i) rd[i] = ref_at(view, j + i) < 4 ? ref_at(view, j + i) : static_cast<uint8_t>(gm() % 4);
      if (it % 3) rd[gm() % len] = static_cast<uint8_t>(gm() % 4);  // a mismatch, or not
      if (it % 7 == 2) rd[gm() % len] = 4;                          // an N in the read
      int want = 0;
      while (want < len && repeat_detail::same_base(rd[want], ref_at(view, j + want))) ++want;
      ++runs;
      bad += repeat_detail::match_run(view, rd.data(), j, len) != want;
    }
    CHECK(bad == 0, "match_run differs from the scalar count on %d of %d windows", bad, runs);
    std::fprintf(stderr, "match_run = scalar same_base count on %d windows: %s\n", runs, bad ? "no" : "yes");
  }
#ifdef CERTA_EDIT_VEC
  // ---- band_distance_vec and hamming_vec must equal their scalar versions on
  // any band, limit and start: reads with N, bands reaching past either end of
  // the reference, early stops (limit + 1) included.
  {
    std::mt19937_64 ge(seed * 7 + 5);
    std::vector<int> pa(BMAX + 2), pb(BMAX + 2);
    const int64_t n = static_cast<int64_t>(view.ref_len);
    int bands = 0, bad_b = 0, windows = 0, bad_h = 0;
    for (int it = 0; it < 20000; ++it) {
      const int B = 1 + static_cast<int>(ge() % BAND_E2E);
      std::vector<uint8_t> rd = reads[ge() % reads.size()].seq;
      int64_t lo;
      if (it % 10 == 0) lo = -static_cast<int64_t>(ge() % 60);                            // before the start
      else if (it % 10 == 1) lo = n - static_cast<int64_t>(rd.size()) + static_cast<int64_t>(ge() % 60) - 30;  // past the end
      else lo = static_cast<int64_t>(ge() % static_cast<uint64_t>(n - 400));
      if (it % 4 >= 2 && lo >= 0 && lo + 200 < n) {  // a read from this very place: small distances
        rd.assign(ref.seq.owned.begin() + lo + B / 2, ref.seq.owned.begin() + lo + B / 2 + 160);
        mutate(ge, rd, static_cast<int>(ge() % 8));
        rd.resize(150, 0);
      }
      if (it % 8 == 1) rd[ge() % rd.size()] = 4;  // an N in the read
      if (it % 16 == 3 && ref.offsets.size() > 1) {
        // The read copies a window across a contig's N padding: N facing N
        // must still cost 1.
        lo = static_cast<int64_t>(ref.offsets[1]) - 40 - static_cast<int64_t>(ge() % 40);
        rd.assign(ref.seq.owned.begin() + lo + B / 2, ref.seq.owned.begin() + lo + B / 2 + 150);
      }
      int limit = static_cast<int>(ge() % 12) - 1;
      int Bt = B;
      if (it % 16 == 7) {
        // A long deletion: the optimum needs the deletion chain across many
        // lanes of one vector (every step of the prefix minimum), so the read
        // starts on diagonal 1 or 2 and the deletion keeps it in the first.
        const int dl = 5 + static_cast<int>(ge() % 9), b0 = 1 + static_cast<int>(ge() % 2);
        Bt = std::min(BAND_E2E, b0 + dl + 2 + static_cast<int>(ge() % 8));
        lo = static_cast<int64_t>(ge() % static_cast<uint64_t>(n - 400));
        const int64_t s0 = lo + b0, cut = 40 + static_cast<int64_t>(ge() % 70);
        rd.assign(ref.seq.owned.begin() + s0, ref.seq.owned.begin() + s0 + cut);
        rd.insert(rd.end(), ref.seq.owned.begin() + s0 + cut + dl, ref.seq.owned.begin() + s0 + 150 + dl);
        limit = dl + static_cast<int>(ge() % 3);
      }
      const int L = static_cast<int>(rd.size());
      const int d1 = band_distance_vec(view, rd.data(), L, lo, Bt, limit);      // the 8-bit form where it applies
      const int d3 = band_distance_vec16(view, rd.data(), L, lo, Bt, limit);
      const int d2 = band_distance_scalar(view, rd.data(), L, lo, Bt, limit, pa.data(), pb.data());
      ++bands;
      bad_b += d1 != d2 || d3 != d2;
      const int64_t st = it % 16 == 3 ? lo + B / 2 : lo + static_cast<int64_t>(ge() % static_cast<uint64_t>(B));
      const int h1 = hamming_vec(view, rd.data(), L, st, limit), h2 = hamming_scalar(view, rd.data(), L, st, limit);
      ++windows;
      bad_h += h1 != h2;
    }
    CHECK(bad_b == 0, "band_distance_vec differs from band_distance_scalar on %d of %d bands", bad_b, bands);
    CHECK(bad_h == 0, "hamming_vec differs from hamming_scalar on %d of %d windows", bad_h, windows);
    std::fprintf(stderr, "band_distance_vec (8- and 16-bit) = band_distance_scalar on %d random bands, hamming_vec = hamming_scalar on %d: %s\n",
                 bands, windows, bad_b || bad_h ? "no" : "yes");
  }
#endif
  // ---- band_score must equal affine_align's score on any band and start range.
  {
    std::unique_ptr<Workspace> ws(new Workspace);
    int compared = 0;
    for (int n = 0; n < 4000; ++n) {
      const int B = 1 + static_cast<int>(g() % BMAX);
      int64_t lo = static_cast<int64_t>(g() % (ref.seq.size() - 600));
      std::vector<uint8_t> rd = reads[g() % reads.size()].seq;
      if (n % 4 >= 2) {  // a read from this very place, so scores are high
        rd.assign(ref.seq.owned.begin() + lo + B / 2, ref.seq.owned.begin() + lo + B / 2 + 160);
        mutate(g, rd, static_cast<int>(g() % 8));
        rd.resize(150, 0);
      }
      const int L = static_cast<int>(rd.size());
      int64_t smin = INT64_MIN, smax = INT64_MAX;
      if (n % 2) {
        smin = lo + static_cast<int64_t>(g() % 200) - 50;
        smax = smin + static_cast<int64_t>(g() % 120);
      }
      if (n % 8 == 5) rd[g() % rd.size()] = 4;  // an N in the read
      Result tmp;
      const bool ok = affine_align(view, rd.data(), L, lo, B, *ws, tmp, true, smin, smax);
      const int v = band_score(view, rd.data(), L, lo, B, smin, smax);
      CHECK(ok ? v == tmp.score : v == kNoSub, "band_score %d vs affine_align %d (ok %d) lo %lld B %d", v,
            ok ? tmp.score : 0, ok, static_cast<long long>(lo), B);
      const int vs = band_score_scalar(view, rd.data(), L, lo, B, smin, smax);
      CHECK(v == vs, "band_score %d vs band_score_scalar %d lo %lld B %d", v, vs, static_cast<long long>(lo), B);
      // With a threshold: exact above it, at or below it otherwise.
      const int theta = (vs == kNoSub ? L - 40 : vs) + static_cast<int>(g() % 21) - 10;
      const int vt = band_score(view, rd.data(), L, lo, B, smin, smax, theta);
      CHECK(vs != kNoSub && vs > theta ? vt == vs : (vt == kNoSub || vt <= theta),
            "band_score with theta %d gives %d, exact %d (lo %lld B %d)", theta, vt, vs, static_cast<long long>(lo), B);

      ++compared;
    }
    std::fprintf(stderr, "band_score compared with affine_align on %d random bands\n", compared);
  }

  // ---- Paired-end certificate (tier PR) against the proper-pair oracle.
  {
    const int D = 500, L = 150;
    struct SimPair { std::vector<uint8_t> a, b; };
    std::vector<SimPair> pairs;
    auto& R = ref.seq.owned;
    auto rcv = [](std::vector<uint8_t> v) {
      std::reverse(v.begin(), v.end());
      for (auto& c : v) c = c < 4 ? static_cast<uint8_t>(3 - c) : c;
      return v;
    };
    for (int n = 0; n < 90; ++n) {
      const int F = L + static_cast<int>(g() % (D + 1));  // fragment length
      int64_t f;
      switch (n % 6) {
        case 0: case 1:  // unique sequence
          f = static_cast<int64_t>(ref.offsets[g() % 2]) + 1000 + static_cast<int64_t>(g() % 3000); break;
        case 2: case 3:  // one mate inside a copy of the 40/10-copy element
          f = static_cast<int64_t>(ref.offsets[0]) + 60200 + 350 * static_cast<int64_t>(g() % 40) -
              static_cast<int64_t>(F) + 100 + static_cast<int64_t>(g() % 150); break;
        case 4:          // diverged segment copies
          f = static_cast<int64_t>(ref.offsets[0]) + 10000 + 20000 * static_cast<int64_t>(g() % 3) +
              static_cast<int64_t>(g() % 1500); break;
        default:         // tandem repeat
          f = static_cast<int64_t>(ref.offsets[0]) + 74800 + static_cast<int64_t>(g() % 600); break;
      }
      std::vector<uint8_t> frag(R.begin() + f, R.begin() + f + F + 10);
      std::vector<uint8_t> a(frag.begin(), frag.begin() + L);
      std::vector<uint8_t> b(frag.begin() + F - L, frag.begin() + F);
      mutate(g, a, static_cast<int>(g() % 5));
      mutate(g, b, static_cast<int>(g() % 9));
      a.resize(L, 0);
      b.resize(L, 0);
      if (n % 7 == 3) {  // adapter-like tail on one mate
        const int cut = 5 + static_cast<int>(g() % 30);
        for (int i = L - cut; i < L; ++i) b[i] = rand_base(g);
      }
      b = rcv(b);
      if (g() % 2) std::swap(a, b);
      pairs.push_back({a, b});
    }
    for (auto& dp : decoy_pairs) {
      if (g() % 2) pairs.push_back({dp.first, dp.second});
      else pairs.push_back({dp.second, dp.first});
    }
    PairParams pp;
    pp.max_dist = D;
    const Params cfgs[] = {{5, 256}, {5, 16}};
    std::vector<int> truth_pair(pairs.size(), INT32_MIN);
    int checked = 0, sub_checked = 0;
    // The mates single-end (the single-end term of a PR mate's XQ).
    ReadBatch mb;
    for (const auto& pr : pairs)
      for (const std::vector<uint8_t>* s : {&pr.a, &pr.b}) {
        mb.offs.push_back(mb.codes.size());
        mb.codes.insert(mb.codes.end(), s->begin(), s->end());
        mb.lens.push_back(static_cast<uint16_t>(s->size()));
        mb.hashes.push_back(mb.lens.size() * 0x9E3779B97F4A7C15ULL);
      }
    PrCheck prc;
    for (const Params& pc : cfgs) {
      pp.se = pc;
      std::vector<Result> se_res;
      map_cpu(view, pc, mb, se_res, 4);
      int cert = 0;
      std::unique_ptr<Workspace> w1(new Workspace), w2(new Workspace);
      PairScratch sc;
      for (size_t i = 0; i < pairs.size(); ++i) {
        Result r1{}, r2{};
        const PairOut po = certify_pair(view, pp, pairs[i].a.data(), L, pairs[i].b.data(), L, *w1, *w2, sc, r1, r2);
        if (!po.certified) continue;
        ++cert;
        CHECK(r1.tier == kTierPR && r2.tier == kTierPR && r1.strand != r2.strand, "pair %zu: bad PR result", i);
        CHECK(rescore_affine(view, r1, pairs[i].a) == r1.score && rescore_affine(view, r2, pairs[i].b) == r2.score,
              "pair %zu: CIGARs score %d/%d, reported %d/%d", i, rescore_affine(view, r1, pairs[i].a),
              rescore_affine(view, r2, pairs[i].b), r1.score, r2.score);
        const Result& fw = r1.strand ? r2 : r1;
        const Result& rv = r1.strand ? r1 : r2;
        CHECK(rv.ref_pos >= fw.ref_pos && rv.ref_pos - fw.ref_pos <= D, "pair %zu: improper (%lld, %lld)", i,
              static_cast<long long>(fw.ref_pos), static_cast<long long>(rv.ref_pos));
        CHECK(po.score == r1.score + r2.score && po.score > po.floor && po.sub <= po.score,
              "pair %zu: inconsistent pair score", i);
        if (truth_pair[i] == INT32_MIN) {
          truth_pair[i] = oracle_pair(R, pairs[i].a, pairs[i].b, D);
          ++checked;
        }
        CHECK(po.score == truth_pair[i], "budget=%d pair %zu: PR score %d != oracle %d", pc.budget, i, po.score,
              truth_pair[i]);
        // S2 bounds every other proper pair (a mate not overlapping R1's).
        const int sub = oracle_pair_sub(R, pairs[i].a, pairs[i].b, D, r1, r2);
        ++sub_checked;
        CHECK(po.sub >= sub, "budget=%d pair %zu: S2 %d below the best other pair %d (R1 %d)", pc.budget, i, po.sub,
              sub, po.score);
        const Result rr[2] = {r1, r2}, ss[2] = {se_res[2 * i], se_res[2 * i + 1]};
        const std::vector<uint8_t> mm[2] = {pairs[i].a, pairs[i].b};
        check_pr(checker4, rr, ss, mm, D, sub, "PR pair", i, prc);
      }
      std::fprintf(stderr, "pairs budget=%-3d certified %d/%zu\n", pc.budget, cert, pairs.size());
    }
    std::fprintf(stderr, "PR scores checked against the proper-pair oracle for %d pairs; S2 against the other-pair oracle %d times\n",
                 checked, sub_checked);
    std::fprintf(stderr, "checker on PR pairs: %d pairs, %d mates' XQ resting on XS, %d proven; %d of %d false pair claims rejected\n",
                 prc.pairs, prc.needed, prc.proven, prc.caught, prc.tampered);
  }
  // ---- Tier PR at its thresholds (Aim 5): pairs built so that one bound of
  // Theorem P is tight, on a reference of their own. A mutant that moves the
  // bound by a few points certifies a wrong pair score or a too-low S2.
  //  A (floor, k = 3, so the floor is U1 + U2 = 260): at the true locus each
  //    mate has one mismatch in each of its 4 searched parts (no exact part:
  //    unseen, 260); a decoy pair scores 259.
  //  B (band pad, k = 5, floor 246, pad 47): mate 1 has a 46-bp deletion at a
  //    part boundary (98 + 150 = 248); a decoy pair scores 247.
  //  D (window threshold, 110-bp mates): the best pair scores 180; another
  //    has a partner with a mismatch in every window part, scoring exactly
  //    the window bound (65, never returned), so S2 must reach 175.
  //  E (window chain skip): two placements of each mate in one window, one
  //    with a 1-bp deletion between window parts (bound 150, score 143), one
  //    with a mismatch (bound 145, score 145); the best pair is 290.
  //  F (S2 through R1's own chain): a tandem mate with a non-overlapping
  //    exact copy and a partner without exact parts: a tie at 270.
  // For the checker's pair proof (check.h, P1-P3), an exact pair (300) and:
  //  G (k = 3): another pair with a mismatch in each searched part of both
  //    mates, scoring exactly U1 + U2 = 260: only (P1) refutes XS = 259;
  //  H: another pair, 900 bp on, whose mate 1 has a 3-bp deletion at a part
  //    boundary (141 + 150 = 291). XS = 290 gives pad 3, which joins its
  //    parts into one chain (bound 150). A pad of 1 splits them (bounds 140
  //    and 130): with the partner's 150 neither exceeds 290, from either
  //    mate, and the false XS would pass;
  //  I: F with a mismatch in x's second copy (265 against 270): the other
  //    pair keeps R1's y (no exact part), so only "x elsewhere in R1's
  //    chain" finds it;
  //  J (110-bp mates): R1's x with two partners without exact parts (180
  //    and 175): only "y elsewhere in R1's window" finds the second.
  {
    std::mt19937_64 h(seed ^ 0xA5A5A5A5ULL);
    Reference tr;
    const int64_t len = 50000;
    tr.seq.owned.assign(kContigPad, 4);
    tr.names.push_back("thr");
    tr.offsets.push_back(kContigPad);
    tr.lengths.push_back(len);
    const std::vector<uint8_t> body = random_seq(h, len);
    tr.seq.owned.insert(tr.seq.owned.end(), body.begin(), body.end());
    tr.seq.owned.insert(tr.seq.owned.end(), kContigPad, 4);
    auto& T = tr.seq.owned;
    auto at = [&](int64_t x) { return static_cast<int64_t>(kContigPad) + x; };
    auto plant = [&](int64_t x, const std::vector<uint8_t>& v) { std::copy(v.begin(), v.end(), T.begin() + at(x)); };
    auto change = [&](std::vector<uint8_t> v, std::initializer_list<int> where) {
      for (int i : where) v[i] = static_cast<uint8_t>((v[i] + 1 + h() % 3) % 4);
      return v;
    };
    auto rcv = [](std::vector<uint8_t> v) {
      std::reverse(v.begin(), v.end());
      for (auto& c : v) c = c < 4 ? static_cast<uint8_t>(3 - c) : c;
      return v;
    };
    auto cat = [](std::vector<uint8_t> a, const std::vector<uint8_t>& b) {
      a.insert(a.end(), b.begin(), b.end());
      return a;
    };
    auto sub = [](const std::vector<uint8_t>& v, int a, int b) { return std::vector<uint8_t>(v.begin() + a, v.begin() + b); };
    struct Fam { char kind; std::vector<uint8_t> m1, m2; int k; };
    std::vector<Fam> fam;
    int64_t base = 200;
    for (int n = 0; n < 4; ++n, base += 1300) {  // A: floor
      const std::vector<uint8_t> frag = random_seq(h, 600);
      plant(base, frag);
      const std::vector<uint8_t> a = change(sub(frag, 0, 150), {9, 28, 47, 66});
      const std::vector<uint8_t> f2 = change(sub(frag, 250, 400), {9, 28, 47, 66});
      plant(base + 700, a);
      plant(base + 950, cat(sub(f2, 0, 114), random_seq(h, 36)));
      plant(base + 1150, sub(f2, 114, 150));  // parts 6-7 of mate 2 twice: parts 0-3 are the rarest
      fam.push_back({'A', a, rcv(f2), 3});
    }
    for (int n = 0; n < 4; ++n, base += 1300) {  // B: band pad
      const std::vector<uint8_t> frag = random_seq(h, 700);
      plant(base, frag);
      const std::vector<uint8_t> a = cat(sub(frag, 0, 76), sub(frag, 122, 196));  // 76M46D74M
      const std::vector<uint8_t> f2 = sub(frag, 300, 450);
      // The decoy's tail differs from mate 1 at every base, so it scores
      // exactly 102 - 5 (with a random tail it could extend by chance).
      std::vector<uint8_t> tail = sub(a, 102, 150);
      for (auto& c : tail) c = static_cast<uint8_t>((c + 1 + h() % 3) % 4);
      plant(base + 750, cat(sub(a, 0, 102), tail));
      plant(base + 1050, f2);
      fam.push_back({'B', a, rcv(f2), 5});
    }
    for (int n = 0; n < 4; ++n, base += 1300) {  // D: window threshold, 110-bp mates
      const std::vector<uint8_t> x = random_seq(h, 110), yf = random_seq(h, 110);
      plant(base, x);
      plant(base + 290, change(yf, {6, 19, 31, 43, 55, 67, 79, 91}));  // 8 window parts: 70
      plant(base + 650, x);
      plant(base + 940, change(yf, {6, 19, 31, 43, 55, 67, 79, 91, 103}));  // all 9: 65
      fam.push_back({'D', x, rcv(yf), 5});
    }
    for (int n = 0; n < 4; ++n, base += 1300) {  // F: S2, mate x elsewhere in R1's own chain
      // x lies in a tandem array (period 20, 320 bp): its placements at
      // shifts 0 and 160 are exact and do not overlap. y has a mismatch in 6
      // parts (no exact part, 120), 400 bp on: both pairs are proper and
      // score 270, so S2 must reach R1 through x's own chain.
      const std::vector<uint8_t> unit = random_seq(h, 20);
      std::vector<uint8_t> arr;
      for (int c = 0; c < 16; ++c) arr.insert(arr.end(), unit.begin(), unit.end());
      plant(base, arr);
      const std::vector<uint8_t> yref = random_seq(h, 150);
      plant(base + 400, yref);
      fam.push_back({'F', sub(arr, 0, 150), rcv(change(yref, {9, 28, 47, 66, 85, 104})), 5});
    }
    for (int n = 0; n < 4; ++n, base += 1300) {  // E: window chain skip
      const std::vector<uint8_t> x = random_seq(h, 150), yf = random_seq(h, 150);
      auto extra = [&](const std::vector<uint8_t>& v) {  // one more reference base at window boundary 78
        return cat(cat(sub(v, 0, 78), random_seq(h, 1)), sub(v, 78, 150));
      };
      plant(base, change(x, {20}));
      plant(base + 160, extra(x));
      plant(base + 330, extra(yf));
      plant(base + 490, change(yf, {20}));
      fam.push_back({'E', x, rcv(yf), 5});
    }
    for (int n = 0; n < 3; ++n, base += 1300) {  // G: the checker's (P1), U1 + U2
      const std::vector<uint8_t> x = random_seq(h, 150), yf = random_seq(h, 150);
      plant(base, x);
      plant(base + 400, yf);
      plant(base + 700, change(x, {9, 28, 47, 66}));  // a mismatch in each searched part
      plant(base + 1100, change(yf, {9, 28, 47, 66}));
      fam.push_back({'G', x, rcv(yf), 3});
    }
    for (int n = 0; n < 3; ++n, base += 1600) {  // H: the checker's pad
      const std::vector<uint8_t> x = random_seq(h, 150), yf = random_seq(h, 150);
      plant(base, x);
      plant(base + 400, yf);
      plant(base + 900, cat(cat(sub(x, 0, 76), random_seq(h, 3)), sub(x, 76, 150)));  // 76M3D74M: 141
      plant(base + 1200, yf);
      fam.push_back({'H', x, rcv(yf), 5});
    }
    for (int n = 0; n < 3; ++n, base += 1300) {  // I: the checker's x elsewhere in R1's chain
      const std::vector<uint8_t> unit = random_seq(h, 20);
      std::vector<uint8_t> arr;
      for (int c = 0; c < 16; ++c) arr.insert(arr.end(), unit.begin(), unit.end());
      plant(base, change(arr, {300}));  // x's copy at shift 160 now has a mismatch: 145
      const std::vector<uint8_t> yref = random_seq(h, 150);
      plant(base + 400, yref);
      fam.push_back({'I', sub(arr, 0, 150), rcv(change(yref, {9, 28, 47, 66, 85, 104})), 5});
    }
    for (int n = 0; n < 3; ++n, base += 1300) {  // J: the checker's other partner of R1's x
      const std::vector<uint8_t> x = random_seq(h, 110), yf = random_seq(h, 110);
      plant(base, x);
      plant(base + 290, change(yf, {6, 19, 31, 43, 55, 67, 79, 91}));       // 70
      plant(base + 450, change(yf, {6, 19, 31, 43, 55, 67, 79, 91, 103}));  // 65
      fam.push_back({'J', x, rcv(yf), 5});
    }
    Index tix = Index::build(tr, 14, 4, 4);
    const IndexView tv = tix.view(tr);
    PairParams tp;
    tp.max_dist = 500;
    std::unique_ptr<Workspace> w1(new Workspace), w2(new Workspace);
    PairScratch tsc;
    int certified = 0, as_designed = 0;
    const check::Checker tck(tr, tix);
    PrCheck tprc;
    for (size_t i = 0; i < fam.size(); ++i) {
      const Fam& f = fam[i];
      tp.se = Params{f.k, 256};
      tp.se.min_support = 1;
      Result r1{}, r2{};
      const int L1 = static_cast<int>(f.m1.size()), L2 = static_cast<int>(f.m2.size());
      const PairOut po = certify_pair(tv, tp, f.m1.data(), L1, f.m2.data(), L2, *w1, *w2, tsc, r1, r2);
      const int opt = oracle_pair(T, f.m1, f.m2, tp.max_dist);
      const int expect = f.kind == 'A' ? 260 : f.kind == 'B' ? 248 : f.kind == 'D' ? 180 : f.kind == 'F' ? 270
                         : f.kind == 'G' || f.kind == 'H' ? 300 : f.kind == 'I' ? 270 : f.kind == 'J' ? 180 : 290;
      const bool designed = opt == expect && (f.kind == 'A' ? !po.certified : po.certified);
      as_designed += designed;
      if (!designed)
        std::fprintf(stderr, "threshold pair %c%zu: oracle %d (designed %d), certified %d\n", f.kind, i, opt, expect,
                     static_cast<int>(po.certified));
      if (!po.certified) continue;
      ++certified;
      CHECK(po.score == opt, "threshold pair %c%zu: PR score %d != oracle %d", f.kind, i, po.score, opt);
      const int sub2 = oracle_pair_sub(T, f.m1, f.m2, tp.max_dist, r1, r2);
      CHECK(po.sub >= sub2, "threshold pair %c%zu: S2 %d below the best other pair %d", f.kind, i, po.sub, sub2);
      // The checker, with no single-end certificate: XQ rests on XS alone.
      const Result rr[2] = {r1, r2}, ss[2] = {Result{}, Result{}};
      const std::vector<uint8_t> mm[2] = {f.m1, f.m2};
      check_pr(tck, rr, ss, mm, tp.max_dist, sub2, "threshold pair", i, tprc);
    }
    std::fprintf(stderr, "threshold pairs: %zu built, %d behave as designed, %d certified and checked\n", fam.size(),
                 as_designed, certified);
    std::fprintf(stderr, "checker on threshold pairs: %d pairs, %d mates' XQ resting on XS, %d proven; %d of %d false "
                         "pair claims rejected\n",
                 tprc.pairs, tprc.needed, tprc.proven, tprc.caught, tprc.tampered);
    // Lemma WS itself: for each mate (reference orientation) and a range of
    // thresholds up to need = L - 2, in windows around its best start
    // (centred, ending there, starting there), window_search returns the
    // best alignment starting in the window, or that one scores <= *thr.
    int ws_checked = 0;
    for (size_t i = 0; i < fam.size(); ++i) {
      for (int x = 0; x < 2; ++x) {
        const std::vector<uint8_t> m = x ? rcv(fam[i].m2) : fam[i].m1;
        const int L = static_cast<int>(m.size());
        const std::vector<int> st = oracle_starts(T, m.data(), L);
        const int64_t p = std::max_element(st.begin(), st.end()) - st.begin();
        const int64_t win[3][2] = {{p - 150, p + 150}, {p - 300, p}, {p, p + 300}};
        for (const auto& w : win) {
          int opt = INT32_MIN;
          for (int64_t j = std::max<int64_t>(0, w[0]); j <= w[1] && j < static_cast<int64_t>(st.size()); ++j)
            opt = std::max(opt, st[j]);
          for (int need : {0, L - 30, L - 10, L - 7, L - 6, L - 5, L - 2}) {
            Band bd;
            int thr;
            const int v = window_search(tv, m.data(), L, w[0], w[1], tp.wq, need, tsc, &bd, &thr);
            ++ws_checked;
            CHECK(v == kNoSub ? opt <= thr : v == opt,
                  "Lemma WS %c%zu mate %d, starts [%lld, %lld], need %d: returned %d (threshold %d), best %d",
                  fam[i].kind, i, x + 1, static_cast<long long>(w[0]), static_cast<long long>(w[1]), need, v, thr,
                  opt);
          }
        }
      }
    }
    std::fprintf(stderr, "Lemma WS checked against the start oracle %d times\n", ws_checked);
  }


  // ---- Tier R (repeats): lo <= OPT <= hi against the whole-reference local
  // oracle; exact claims equal OPT; copies are verified string by string.
  {
    Index ix1 = Index::build(ref, 14, 1, 4);  // tier R needs s = 1
    const IndexView view1 = ix1.view(ref);
    const auto& R = ref.seq.owned;
    auto rcs = [](std::vector<uint8_t> v) {
      std::reverse(v.begin(), v.end());
      for (auto& c : v) c = c < 4 ? static_cast<uint8_t>(3 - c) : c;
      return v;
    };
    RepeatParams cfgs[5];
    cfgs[1].tau = 2;
    cfgs[1].cmax = 64;
    cfgs[2].tau = 64;
    cfgs[2].samples = 2;
    for (int c = 3; c < 5; ++c) {  // interval mode: SUB bounds, extra chains, satellite fallback
      cfgs[c].tau = 256;
      cfgs[c].sub_margin = 20;
      cfgs[c].fallback_samples = 8;
      cfgs[c].lb_chains = 8;
    }
    cfgs[4].cmax = 8;  // most anchors unaffordable: the fallback samples
    std::unique_ptr<Workspace> ws(new Workspace);
    RepeatScratch sc;
    int checked = 0;
    int sub_r_checked = 0;
    for (int c = 0; c < 5; ++c) {
      int exact = 0, exact_rep = 0, ties = 0;
      for (size_t i = 0; i < reads.size(); ++i) {
        const bool rep = i >= repeat_first && i < bounds_first;
        if (!rep && !planted_bound(i) && i % 3 != 0) continue;  // oracle cost: all repeat and planted reads, a third of the others
        const auto& rd = reads[i].seq;
        const int L = static_cast<int>(rd.size());
        const RepeatOut o = repeat_tier(view1, cfgs[c], rd.data(), L, *ws, sc);
        if (!o.claim || o.lo == kNoSub) continue;
        if (local_truth[i] == INT32_MIN) local_truth[i] = oracle_local(R, rd.data(), L);
        ++checked;
        const int opt = local_truth[i];
        CHECK(o.lo <= opt && opt <= o.hi, "tier R cfg %d read %zu: OPT %d outside [%d, %d]", c, i, opt, o.lo, o.hi);
        CHECK(!o.exact || o.lo == opt, "tier R cfg %d read %zu: exact %d but OPT %d", c, i, o.lo, opt);
        CHECK(rescore_affine(view1, o.aln, rd) == o.lo, "tier R read %zu: CIGAR scores %d, lo %d", i,
              rescore_affine(view1, o.aln, rd), o.lo);
        CHECK(o.has_aln, "tier R read %zu: lo %d without an alignment", i, o.lo);
        // SUB interval: sub_lo <= SUB <= sub_hi, sub_hi <= hi.
        CHECK(o.sub_hi <= o.hi && (o.sub_lo == kNoSub || o.sub_lo <= o.sub_hi),
              "tier R cfg %d read %zu: SUB interval [%d, %d], hi %d", c, i, o.sub_lo, o.sub_hi, o.hi);
        if (o.has_aln && (o.sub_hi < opt || o.sub_lo != kNoSub) &&
            (planted_bound(i) || ((c == 3 || c == 4) && i % 2 == 0))) {
          const int sub = oracle_sub(R, rd, o.aln.strand, o.aln.ref_pos, ref_span_of(o.aln));
          ++sub_r_checked;
          CHECK(sub <= o.sub_hi, "tier R cfg %d read %zu: SUB %d above sub_hi %d (lo %d, exact %d)", c, i, sub,
                o.sub_hi, o.lo, o.exact);
          CHECK(o.sub_lo == kNoSub || o.sub_lo <= sub, "tier R cfg %d read %zu: sub_lo %d above SUB %d", c, i,
                o.sub_lo, sub);
        }
        exact += o.exact;
        exact_rep += o.exact && rep;
        if (o.copies >= 2) {
          ++ties;
          int64_t span = 0;
          for (int k = 0; k < o.aln.n_cigar; ++k)
            if ((o.aln.cigar[k] & 0xF) == kOpM || (o.aln.cigar[k] & 0xF) == kOpD) span += o.aln.cigar[k] >> 4;
          std::vector<uint8_t> w0(R.begin() + o.copy_pos[0], R.begin() + o.copy_pos[0] + span);
          std::vector<uint8_t> w1(R.begin() + o.copy_pos[1], R.begin() + o.copy_pos[1] + span);
          const bool same = o.copy_strand[1] == o.copy_strand[0];
          CHECK((same ? w1 == w0 : w1 == rcs(w0)), "tier R read %zu: claimed copy differs", i);
          CHECK(!same || o.copy_pos[1] + span <= o.copy_pos[0] || o.copy_pos[0] + span <= o.copy_pos[1],
                "tier R read %zu: claimed copy overlaps", i);
          for (uint8_t b : w0) CHECK(b < 4, "tier R read %zu: copy contains N", i);
        }
      }
      std::fprintf(stderr, "tier R cfg %d (tau %d, cmax %d): exact %d (repeat reads %d of %zu), copies>=2 %d\n", c,
                   cfgs[c].tau, cfgs[c].cmax, exact, exact_rep, bounds_first - repeat_first, ties);
    }
    std::fprintf(stderr, "tier R intervals checked against the local oracle: %d; SUB intervals: %d\n", checked,
                 sub_r_checked);

    // ---- Interval-mode placement of pairs (certa/interval.h): every record
    // rescores to its AS, [XL, XU] contains the mate's optimum, a proper pair
    // is proper, XQ <= MAPQ.
    {
      const int D = 500, Lp = 150;
      RepeatParams ip;
      ip.tau = 256;
      ip.sub_margin = 20;
      ip.fallback_samples = 8;
      ip.lb_chains = 8;
      PairScratch psc;
      int placed = 0, proper_n = 0, rescued = 0, pair_checked = 0;
      auto contig_of = [&](int64_t pos) { return ref.contig_of(pos); };
      const check::Checker checker_p(ref, ix1);
      for (int n = 0; n < 60; ++n) {
        // Fragments from unique sequence, the 50-copy element and diverged
        // copies; mate 2 sometimes heavily mutated (only rescue can place it).
        const int F = Lp + static_cast<int>(g() % (D - Lp));
        int64_t f;
        switch (n % 3) {
          case 0: f = static_cast<int64_t>(ref.offsets[g() % 2]) + 1000 + static_cast<int64_t>(g() % 3000); break;
          case 1: f = static_cast<int64_t>(ref.offsets[0]) + 60200 + 350 * static_cast<int64_t>(g() % 40) -
                      static_cast<int64_t>(F) + 100 + static_cast<int64_t>(g() % 150); break;
          default: f = static_cast<int64_t>(ref.offsets[0]) + 10000 + 20000 * static_cast<int64_t>(g() % 3) +
                       static_cast<int64_t>(g() % 1500); break;
        }
        std::vector<uint8_t> a(R.begin() + f, R.begin() + f + Lp + 10);
        std::vector<uint8_t> b(R.begin() + f + F - Lp, R.begin() + f + F + 10);
        mutate(g, a, static_cast<int>(g() % 4));
        mutate(g, b, n % 4 == 3 ? 12 + static_cast<int>(g() % 8) : static_cast<int>(g() % 4));
        a.resize(Lp, 0);
        b.resize(Lp, 0);
        b = rcs(b);
        if (g() % 2) std::swap(a, b);
        const std::vector<uint8_t>* m[2] = {&a, &b};
        RepeatOut iv[2];
        Placement se[2], out[2];
        const RepeatOut* ivp[2] = {&iv[0], &iv[1]};
        const uint8_t* cp[2] = {a.data(), b.data()};
        const int Ls[2] = {Lp, Lp};
        for (int x = 0; x < 2; ++x) {
          iv[x] = repeat_tier(view1, ip, m[x]->data(), Lp, *ws, sc);
          Result none{};
          se[x] = place_single(none, 0, &iv[x], Lp);
        }
        const bool proper = place_pair(view1, D, 12, se, ivp, cp, Ls, psc, contig_of, out);
        proper_n += proper;
        for (int x = 0; x < 2; ++x) {
          const Placement& o = out[x];
          const int opt = local_one(R, m[x]->data(), Lp, 0, static_cast<int64_t>(R.size()));
          const std::vector<uint8_t> rcm = rcs(*m[x]);
          const int opt2 = std::max(opt, local_one(R, rcm.data(), Lp, 0, static_cast<int64_t>(R.size())));
          CHECK(o.xu >= opt2 && (o.xl == kNoSub || o.xl <= opt2), "pair %d mate %d: OPT %d outside [%d, %d]", n, x,
                opt2, o.xl, o.xu);
          CHECK(o.xq <= o.mapq, "pair %d mate %d: XQ %d above MAPQ %d", n, x, o.xq, o.mapq);
          if (!o.mapped) continue;
          ++placed;
          rescued += o.r.tier == kTierMR;
          CHECK(rescore_affine(view1, o.r, *m[x]) == o.r.score, "pair %d mate %d: CIGAR scores %d, AS %d", n, x,
                rescore_affine(view1, o.r, *m[x]), o.r.score);
          CHECK(o.xl != kNoSub && o.r.score <= o.xl, "pair %d mate %d: AS %d above XL %d", n, x, o.r.score, o.xl);
        }
        for (int x = 0; x < 2; ++x) {
          check::Work w;
          const check::Verdict v = checker_p.check(claim_of(out[x], *m[x]), &w);
          ++pair_checked;
          CHECK(v == check::kAccepted, "pair %d mate %d: checker says %s", n, x, check::kVerdictNames[v]);
        }
        if (proper) {
          const Result& fw = out[0].r.strand ? out[1].r : out[0].r;
          const Result& rv = out[0].r.strand ? out[0].r : out[1].r;
          CHECK(out[0].mapped && out[1].mapped && out[0].r.strand != out[1].r.strand && rv.ref_pos >= fw.ref_pos &&
                    rv.ref_pos - fw.ref_pos <= D,
                "pair %d: claimed proper but not", n);
        }
      }
      std::fprintf(stderr, "interval pairs: %d records placed (%d by the pair window), %d proper pairs of 60; %d checked\n",
                   placed, rescued, proper_n, pair_checked);
    }

    // ---- The checker on tier R and pair records (index q = 14, s = 1), and
    // index validation.
    {
      const check::Checker checker1(ref, ix1);
      RepeatParams ip;
      ip.tau = 256;
      ip.sub_margin = 20;
      ip.fallback_samples = 8;
      ip.lb_chains = 8;
      int acc = 0, total = 0, tampered = 0, caught = 0, by_proof = 0;
      for (size_t i = 0; i < reads.size(); ++i) {
        if (i % 2 && !planted_bound(i)) continue;
        const auto& rd = reads[i].seq;
        const int L = static_cast<int>(rd.size());
        const RepeatOut o = repeat_tier(view1, ip, rd.data(), L, *ws, sc);
        if (!o.claim) continue;
        Result none{};
        const Placement pl = place_single(none, 0, &o, L);
        const check::Claim c = claim_of(pl, rd);
        check::Work w;
        const check::Verdict v = checker1.check(c, &w);
        ++total;
        acc += v == check::kAccepted;
        CHECK(v == check::kAccepted, "tier R read %zu: checker says %s (XL %d XU %d XH %d XQ %d)", i,
              check::kVerdictNames[v], c.xl, c.xu, c.xh, c.xq);
        if (local_truth[i] == INT32_MIN) local_truth[i] = oracle_local(R, rd.data(), L);
        const int opt = local_truth[i];
        if (opt - 1 >= 0) {  // with no alignment shown, only the bound search can refute it
          check::Claim t = c;
          t.mapped = false;
          t.xl = t.xh = check::kNone;
          t.xq = 0;
          t.has_alt = false;
          t.xu = opt - 1;
          check::Work w2;
          ++tampered;
          const check::Verdict v2 = checker1.check(t, &w2);
          const bool rej = v2 != check::kAccepted;
          caught += rej;
          by_proof += v2 == check::kXUUnproven;
          CHECK(rej, "tier R read %zu: false XU = OPT - 1 = %d accepted", i, opt - 1);
        }
        if (pl.mapped && pl.xq > 0 && i % 4 == 0) {
          const int sub = oracle_sub(R, rd, pl.r.strand, pl.r.ref_pos, ref_span_of(pl.r));
          if (sub >= 19 && sub < c.as) {
            check::Claim t = c;
            t.xh = sub - 1;
            check::Work w2;
            ++tampered;
            const check::Verdict v2 = checker1.check(t, &w2);
            const bool rej = v2 != check::kAccepted;
            caught += rej;
            by_proof += v2 == check::kXHUnproven;
            CHECK(rej, "tier R read %zu: false XH = SUB - 1 = %d accepted", i, sub - 1);
          }
        }
      }
      std::fprintf(stderr,
                   "checker: tier R records accepted %d of %d; %d of %d falsified claims rejected (%d by the bound "
                   "search, the rest as XL > XU or wrong XQ)\n",
                   acc, total, caught, tampered, by_proof);

      // Index validation: the in-memory, the memory-mapped (16-bit keys) and
      // the s = 1 index are valid, with the reference written as FASTA;
      // corruptions are not.
      const std::string fa = "test_certa." + std::to_string(getpid()) + ".tmp.fa";
      {
        std::ofstream out(fa);
        for (size_t c = 0; c < ref.names.size(); ++c) {
          out << '>' << ref.names[c] << " test\n";
          for (uint64_t k = 0; k < ref.lengths[c]; ++k) {
            out << "ACGTN"[R[ref.offsets[c] + k]];
            if (k % 60 == 59 || k + 1 == ref.lengths[c]) out << '\n';
          }
        }
      }
      std::string why;
      CHECK(check::validate_index(ref, ix, fa, 4, &why), "index q14 s4 invalid: %s", why.c_str());
      CHECK(check::validate_index(ref, ix1, fa, 4, &why), "index q14 s1 invalid: %s", why.c_str());
      CHECK(check::validate_index(ref_m, ix_m, fa, 4, &why), "memory-mapped index invalid: %s", why.c_str());
      auto corrupt = [&](const char* what, const std::function<void(Index&)>& f) {
        Index bad = ix;
        f(bad);
        std::string w;
        CHECK(!check::validate_index(ref, bad, "", 4, &w), "corrupted index accepted (%s)", what);
      };
      corrupt("position moved", [](Index& b) { b.pos.owned[b.pos.owned.size() / 2] += 4; });
      corrupt("key changed", [](Index& b) { b.keys.owned[b.keys.owned.size() / 3] ^= 1; });
      corrupt("entry dropped", [](Index& b) {
        b.keys.owned.pop_back();
        b.pos.owned.pop_back();
        b.dir.owned.back() -= 1;
      });
      corrupt("entry duplicated", [](Index& b) { b.pos.owned[1] = b.pos.owned[0]; b.keys.owned[1] = b.keys.owned[0]; });
      {
        Reference r2 = ref;
        r2.seq.owned[ref.offsets[1] + 77] = static_cast<uint8_t>((r2.seq.owned[ref.offsets[1] + 77] + 1) % 4);
        std::string w;
        CHECK(!check::validate_index(r2, ix, fa, 4, &w), "a reference differing from the FASTA accepted");
      }
      std::remove(fa.c_str());
      std::fprintf(stderr, "checker: index validation (3 valid, 5 corruptions rejected) done\n");
    }
  }

  // ---- Pass V (adaptive parts, certify_vparts; THEORY.md §4b) against the
  // edit-distance oracle with radii up to VMAX: a certified d1 is the global
  // minimum, nothing within the radius is missed, a unique best with no second
  // locus within R has none, the result does not depend on threads or on the
  // order clusters are verified in, and [XL, XU] holds OPT. Lemma N and the
  // DP over disjoint q-mers are checked directly.
  {
    std::mt19937_64 gv(seed ^ 0xadab7ee5ULL);  // its own stream: the other tests' data stay as they were
    const auto& R = ref.seq.owned;
    {  // Lemma N: the reference q-mer at the aligned start of a stretch within one edit is a neighbour
      uint64_t keys[NBMAX];
      int checked = 0, missing = 0;
      for (int it = 0; it < 20000; ++it) {
        const int qq = 8 + static_cast<int>(gv() % 25);
        const std::vector<uint8_t> s = random_seq(gv, qq);
        std::vector<uint8_t> t = s;
        const int kind = static_cast<int>(gv() % 4), at = static_cast<int>(gv() % qq);
        if (kind == 1) t[at] = static_cast<uint8_t>((t[at] + 1 + gv() % 3) % 4);  // a substitution
        if (kind == 2) t.erase(t.begin() + at);                                    // a read base missing from the reference
        if (kind == 3) t.insert(t.begin() + at, rand_base(gv));                    // a reference base before read base `at`
        while (static_cast<int>(t.size()) < qq) t.push_back(rand_base(gv));       // the reference continues
        t.resize(qq);
        const int n = vparts_neighbours(s.data(), qq, keys);
        uint64_t key;
        kmer_key(t.data(), qq, &key);
        ++checked;
        missing += !std::binary_search(keys, keys + n, key);
      }
      CHECK(missing == 0, "Lemma N: %d of %d one-edit reference q-mers not among the neighbours", missing, checked);
      std::fprintf(stderr, "Lemma N: %d one-edit reference q-mers, all among the neighbours: %s\n", checked,
                   missing ? "no" : "yes");
    }
    {  // vparts_dp / vparts_pick against brute force: the most disjoint q-mers within the budget, least cost
      std::unique_ptr<Workspace> ws(new Workspace);
      int checked = 0;
      for (int it = 0; it < 600; ++it) {
        const int L = 20 + static_cast<int>(gv() % 17), qq = 6 + static_cast<int>(gv() % 5), n = L - qq + 1;
        uint64_t cost[LMAX];
        for (int a = 0; a < n; ++a) cost[a] = gv() % 7 == 0 ? kVInf : gv() % 12;
        const uint64_t budget = it % 10 == 0 ? kVInf - 1 : gv() % 40;
        const int m = vparts_dp(L, qq, cost, budget, *ws, 0);
        int bm = 0;
        uint64_t bc = 0;
        std::function<void(int, int, uint64_t)> rec = [&](int from, int k, uint64_t c) {
          if (k > bm || (k == bm && c < bc)) { bm = k; bc = c; }
          for (int a = from; a < n; ++a)
            if (cost[a] < kVInf && c + cost[a] <= budget) rec(a + qq, k + 1, c + cost[a]);
        };
        rec(0, 0, 0);
        ++checked;
        CHECK(m == bm, "vparts_dp: %d parts, brute force %d (L %d q %d budget %llu)", m, bm, L, qq,
              static_cast<unsigned long long>(budget));
        if (m != bm || m == 0) continue;
        vparts_pick(L, m, *ws, 0);
        uint64_t c = 0;
        bool ok = true;
        for (int x = 0; x < m; ++x) {
          const int a = ws->voff[0][x];
          ok &= a >= 0 && a + qq <= L && (x == 0 || a >= ws->voff[0][x - 1] + qq) && cost[a] < kVInf;
          if (ok) c += cost[a];
        }
        CHECK(ok && c == bc, "vparts_pick: parts %s, cost %llu, brute force %llu", ok ? "disjoint" : "NOT disjoint",
              static_cast<unsigned long long>(c), static_cast<unsigned long long>(bc));
      }
      std::fprintf(stderr, "pass V DP checked against brute force %d times\n", checked);
    }
    const int q = 14;
    Index ixv = Index::build(ref, q, 1, 4);  // pass V needs s = 1
    const IndexView viewv = ixv.view(ref);
    // Reads: every second read above, and reads of 60-150 bases with 0..13
    // edits from unique sequence, the diverged copies, the 50-copy element and
    // the tandem repeat, on both strands (a few with an N).
    std::vector<std::vector<uint8_t>> vreads;
    for (size_t i = 0; i < reads.size(); i += 2) vreads.push_back(reads[i].seq);
    const int64_t vanchor[] = {-1, -1, -1, 10000, 30000, 50000, 60250, 66000, 75000};
    for (int n = 0; n < 400; ++n) {
      const int L = 60 + 20 * (n % 4) + (n % 4 == 3 ? 30 : 0);  // 60, 80, 100, 150
      const int c = static_cast<int>(gv() % 2);
      const int64_t an = vanchor[gv() % 9];
      int64_t start = an >= 0 && c == 0 ? an + static_cast<int64_t>(gv() % 1500)
                                        : static_cast<int64_t>(gv() % (ref.lengths[c] - L - 40));
      start += static_cast<int64_t>(ref.offsets[c]);
      std::vector<uint8_t> s(R.begin() + start, R.begin() + start + L + 20);
      mutate(gv, s, static_cast<int>(gv() % 14));
      s.resize(L);
      if (n % 37 == 11) s[gv() % L] = 4;
      if (gv() % 2) {
        std::vector<uint8_t> rc(L);
        for (int i = 0; i < L; ++i) rc[i] = s[L - 1 - i] < 4 ? 3 - s[L - 1 - i] : 4;
        s.swap(rc);
      }
      vreads.push_back(s);
    }
    ReadBatch vbatch;
    std::vector<int> vmin(vreads.size(), 99);
    std::vector<std::vector<Hit>> vtruth(vreads.size());
    for (size_t i = 0; i < vreads.size(); ++i) {
      vbatch.offs.push_back(vbatch.codes.size());
      vbatch.codes.insert(vbatch.codes.end(), vreads[i].begin(), vreads[i].end());
      vbatch.lens.push_back(static_cast<uint16_t>(vreads[i].size()));
      vbatch.hashes.push_back(i * 0x9E3779B97F4A7C15ULL + 7);
      vtruth[i] = oracle(R, vreads[i].data(), static_cast<int>(vreads[i].size()), VMAX);
      for (auto& h : vtruth[i]) vmin[i] = std::min(vmin[i], h.dist);
    }
    std::vector<int> vopt(vreads.size(), INT32_MIN);
    auto vp = [](int vparts, int budget, int k) {
      Params p{k, budget};
      p.vparts = vparts;
      p.mapq_limit = 8;
      return p;
    };
    const Params vcfgs[] = {vp(1, 4096, VMAX), vp(2, 4096, VMAX), vp(2, 256, VMAX), vp(2, 64, VMAX), vp(2, 16, VMAX),
                            vp(2, 4096, 7)};
    int beyond_kmax = 0, bounds = 0, subs = 0;
    for (const Params& p : vcfgs) {
      std::vector<Result> res, res1, resrev;
      map_cpu(viewv, p, vbatch, res, 4);
      map_cpu(viewv, p, vbatch, res1, 1);
      Params prev = p;
      prev.reverse_order = 1;
      map_cpu(viewv, prev, vbatch, resrev, 4);
      int cert = 0, deep = 0, nf = 0;
      for (size_t i = 0; i < vreads.size(); ++i) {
        const Result& r = res[i];
        const int L = static_cast<int>(vreads[i].size());
        CHECK(std::memcmp(&r, &res1[i], sizeof r) == 0, "pass V read %zu differs between 1 and 4 threads", i);
        const Result& rv = resrev[i];
        CHECK(r.certified == rv.certified && r.d1 == rv.d1 && r.radius == rv.radius,
              "pass V vparts=%d budget=%d read %zu: order-dependent result", p.vparts, p.budget, i);
        if (r.certified == 1 && r.n_best == 1 && rv.n_best == 1)
          CHECK(r.d2 == rv.d2, "pass V read %zu: order-dependent d2 (%d vs %d)", i, r.d2, rv.d2);
        CHECK((r.n_best > 1) == (rv.n_best > 1), "pass V read %zu: order-dependent tie", i);
        CHECK(r.radius <= std::min(p.k, VMAX), "pass V read %zu: radius %d above the cap", i, r.radius);
        // The budget bounds the work (and ws.cand): at most `budget` hits per strand.
        CHECK(r.ncand <= 2u * static_cast<uint32_t>(p.budget), "pass V budget=%d read %zu: %u hits enumerated",
              p.budget, i, r.ncand);
        if (r.certified != 1) {
          // Lossless: an uncertified read has nothing within the radius.
          if (r.reason == kNotFound) {
            ++nf;
            CHECK(vmin[i] > r.radius, "pass V vparts=%d budget=%d read %zu: oracle %d within R %d, not found",
                  p.vparts, p.budget, i, vmin[i], r.radius);
          }
          continue;
        }
        ++cert;
        deep += r.d1 > KMAX;
        CHECK(r.tier == kTierS0 || r.tier == kTierS1, "pass V read %zu: tier %d", i, r.tier);
        CHECK(r.d1 <= r.radius && r.support == 0 && r.probes[0] == 0 && r.probes[1] == 0,
              "pass V read %zu: d1 %d R %d support %d", i, r.d1, r.radius, r.support);
        CHECK(rescore(viewv, r, vreads[i]) == r.nm, "pass V read %zu: CIGAR rescores to %d, NM %d", i,
              rescore(viewv, r, vreads[i]), r.nm);
        bool clipped = false;
        for (int c = 0; c < r.n_cigar; ++c) clipped |= (r.cigar[c] & 0xF) == kOpS;
        CHECK(clipped || r.nm >= r.d1, "pass V read %zu: unclipped NM %d below d1 %d", i, r.nm, r.d1);
        CHECK(r.d1 == vmin[i], "pass V vparts=%d budget=%d read %zu: d1 %d != oracle %d", p.vparts, p.budget, i,
              r.d1, vmin[i]);
        if (r.n_best == 1 && r.d2 < 0) {  // no second locus within R
          int64_t span = 0;
          for (int c = 0; c < r.n_cigar; ++c)
            if ((r.cigar[c] & 0xF) == kOpM || (r.cigar[c] & 0xF) == kOpD) span += r.cigar[c] >> 4;
          const int64_t end = r.ref_pos + span - 1;
          for (auto& h : vtruth[i]) {
            if (h.dist > r.radius) continue;
            const bool near = h.strand == r.strand && h.end - end <= 2 * BMAX && end - h.end <= 2 * BMAX;
            CHECK(near, "pass V read %zu: oracle locus (strand %d end %lld d %d) within R %d missed", i, h.strand,
                  static_cast<long long>(h.end), h.dist, r.radius);
          }
        }
        if (i % 4 == 0) {  // [XL, XU] and SUB_hi against the local oracles
          if (vopt[i] == INT32_MIN) vopt[i] = oracle_local(R, vreads[i].data(), L);
          int xl, xu, sh;
          certified_bounds(r, L, &xl, &xu, &sh);
          ++bounds;
          CHECK(xl <= vopt[i] && vopt[i] <= xu, "pass V read %zu: OPT %d outside [%d, %d]", i, vopt[i], xl, xu);
          if (i % 8 == 0 && sh < vopt[i]) {
            const int sub = oracle_sub(R, vreads[i], r.strand, r.ref_pos, ref_span_of(r));
            ++subs;
            CHECK(sub <= sh, "pass V read %zu: SUB %d above SUB_hi %d (d1 %d d2 %d R %d)", i, sub, sh, r.d1, r.d2,
                  r.radius);
          }
        }
      }
      beyond_kmax += deep;
      std::fprintf(stderr, "pass V vparts=%d budget=%-4d k=%-2d certified %3d/%zu (%d beyond %d edits), not found %d\n",
                   p.vparts, p.budget, p.k, cert, vreads.size(), deep, KMAX, nf);
    }
    CHECK(beyond_kmax > 0, "pass V certified no read beyond %d edits: the test does not exercise it", KMAX);
    std::fprintf(stderr, "pass V: [XL, XU] checked %d times, SUB bounds %d times\n", bounds, subs);
  }
  if (failures) {
    std::fprintf(stderr, "%d check(s) FAILED\n", failures);
    return 1;
  }
  std::fprintf(stderr, "all certificate checks passed\n");
  return 0;
}
