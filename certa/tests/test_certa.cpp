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
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <vector>

#include "index.h"
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
  Index ix = Index::build(ref, 14, 4, 4);
  const IndexView view = ix.view(ref);

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

  const Params configs[] = {{2, 256}, {3, 256}, {5, 256}, {4, 64}, {4, 16}, {3, 4}, {2, 1}};
  for (const Params& p : configs) {
    std::vector<Result> res, res1, resrev;
    map_cpu(view, p, batch, res, 4);
    map_cpu(view, p, batch, res1, 1);
    Params prev = p;
    prev.reverse_order = 1;
    map_cpu(view, prev, batch, resrev, 4);
    int certified = 0, easy = 0, easy_cert = 0, sr = 0;
    int reasons[kNumReasons] = {};
    for (size_t i = 0; i < reads.size(); ++i) {
      const Result& r = res[i];
      CHECK(std::memcmp(&r, &res1[i], sizeof r) == 0, "read %zu differs between 1 and 4 threads", i);
      // Metamorphic: verifying clusters in reverse order must give the same
      // certificate (d1; d2 whenever the best locus is unique).
      const Result& rv = resrev[i];
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
        easy_cert += r.certified;
      }
      if (r.reason == kNotFound)
        CHECK(omin > r.radius, "k=%d budget=%d read %zu: lossless violated, oracle %d <= R %d",
              p.k, p.budget, i, omin, r.radius);
      if (!r.certified) continue;
      ++certified;
      CHECK(rescore(view, r, reads[i].seq) == r.d1,
            "k=%d read %zu: CIGAR rescores to %d, reported %d", p.k, i,
            rescore(view, r, reads[i].seq), r.d1);
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
          if ((r.cigar[c] & 0xF) != kOpI) span += r.cigar[c] >> 4;
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
                 "k=%d budget=%-3d certified %3d/%zu (SR %d); reads with <=k edits certified %d/%d; "
                 "uncertified: length %d, radius<0 %d, not-found %d, wide %d\n",
                 p.k, p.budget, certified, reads.size(), sr, easy_cert, easy, reasons[kBadLength],
                 reasons[kRadiusNegative], reasons[kNotFound], reasons[kClusterTooWide]);
  }
  if (failures) {
    std::fprintf(stderr, "%d check(s) FAILED\n", failures);
    return 1;
  }
  std::fprintf(stderr, "all certificate checks passed\n");
  return 0;
}
