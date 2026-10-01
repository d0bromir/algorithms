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
#include <algorithm>
#include <climits>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <vector>

#include "index.h"
#include "certa/pair.h"
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

}  // namespace

int main(int argc, char** argv) {
  const uint64_t seed = argc > 1 ? std::strtoull(argv[1], nullptr, 10) : 20260927;
  std::fprintf(stderr, "seed %llu\n", static_cast<unsigned long long>(seed));
  std::mt19937_64 g(seed);
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
  Index ix = Index::build(ref, 14, 4, 4);
  const IndexView view = ix.view(ref);
  // Round trip through the on-disk format: with q = 14 the file stores only
  // 16-bit key suffixes and is memory-mapped; results must not change.
  const std::string tmp = "test_certa.tmp.cidx";
  save_index(tmp, ref, ix);
  Reference ref_m;
  Index ix_m;
  load_index(tmp, ref_m, ix_m);
  std::remove(tmp.c_str());  // the mapping stays valid until ix_m is gone
  const IndexView view_m = ix_m.view(ref_m);
  CHECK(view_m.keys16 != nullptr, "saved index should use 16-bit keys (q=14)");

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
    PairParams pp;
    pp.max_dist = D;
    const Params cfgs[] = {{5, 256}, {5, 16}};
    std::vector<int> truth_pair(pairs.size(), INT32_MIN);
    int checked = 0;
    for (const Params& pc : cfgs) {
      pp.se = pc;
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
      }
      std::fprintf(stderr, "pairs budget=%-3d certified %d/%zu\n", pc.budget, cert, pairs.size());
    }
    std::fprintf(stderr, "PR scores checked against the proper-pair oracle for %d pairs\n", checked);
  }
  if (failures) {
    std::fprintf(stderr, "%d check(s) FAILED\n", failures);
    return 1;
  }
  std::fprintf(stderr, "all certificate checks passed\n");
  return 0;
}
