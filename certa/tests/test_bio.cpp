// Brute-force validation of --bio (certa/bio.h, Theorem B).
//
// The oracle is the model's DP over the *entire* reference, both strands
// (Gotoh, end to end, the same costs): for every end position, the least
// penalty of an alignment ending there. For every read and every setting
// (two indexes, default and tight parameters) it checks:
//   1. soundness    - the reported alignment rescores (bio_rescore) to YP;
//   2. optimality   - a certified YP equals the oracle's global minimum;
//   3. completeness - the loci below the horizon (cheapest first, L apart)
//                     are the oracle's, value for value, so the second best
//                     and the ties are exact too;
//   4. canonical    - no gap of the reported alignment moves one base left at
//                     the same penalty (indels left-aligned);
//   5. determinism  - the same results in reverse order, on another workspace;
//   6. pairs        - a jointly placed pair is the oracle's best pair of loci.
// Lemma B (bio_cover) is checked against an exact DP over transcripts on
// random small reads. Placement against the simulated truth is reported, next
// to the fewest-edits placement, and must improve on reads where the base
// qualities decide.
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <random>
#include <string>
#include <thread>
#include <vector>

#include "certa/bio.h"
#include "index.h"

using namespace certa;

namespace {

std::atomic<int> failures{0};
std::mutex print_mu;
#define CHECK(cond, ...)                                          \
  do {                                                            \
    if (!(cond)) {                                                \
      if (++failures <= 25) {                                     \
        std::lock_guard<std::mutex> lk(print_mu);                 \
        std::fprintf(stderr, "FAIL %s:%d: ", __FILE__, __LINE__); \
        std::fprintf(stderr, __VA_ARGS__);                        \
        std::fprintf(stderr, "\n");                               \
      }                                                           \
    }                                                             \
  } while (0)

template <class F>
void parallel(size_t n, F&& f) {
  const unsigned T = std::max(1u, std::min(8u, std::thread::hardware_concurrency()));
  std::atomic<size_t> next{0};
  std::vector<std::thread> pool;
  for (unsigned t = 0; t < T; ++t)
    pool.emplace_back([&] {
      for (size_t i; (i = next.fetch_add(1)) < n;) f(i);
    });
  for (auto& th : pool) th.join();
}

uint8_t rand_base(std::mt19937_64& g) { return static_cast<uint8_t>(g() % 4); }

std::vector<uint8_t> random_seq(std::mt19937_64& g, size_t n) {
  std::vector<uint8_t> s(n);
  for (auto& c : s) c = rand_base(g);
  return s;
}

uint8_t transition(uint8_t b) { return static_cast<uint8_t>(b ^ 2); }

std::vector<uint8_t> revcomp(const std::vector<uint8_t>& s) {
  std::vector<uint8_t> r(s.rbegin(), s.rend());
  for (auto& c : r) c = c < 4 ? static_cast<uint8_t>(3 - c) : c;
  return r;
}

// ---------------------------------------------------------------------------
// The oracle: V[st][e] = the least penalty of an end-to-end alignment of
// strand st of the read whose last base lies on reference position e. Column
// by column over the whole reference, the recurrences of bio_band.
void oracle(const IndexView& ix, const BioModel& m, const BioRead& br, const std::vector<uint8_t>& ctx,
            std::vector<int>* V) {
  const int L = br.L;
  const int64_t N = static_cast<int64_t>(ix.ref_len);
  std::vector<int> Mp(L + 1), Ip(L + 1), Dp(L + 1), Mc(L + 1), Ic(L + 1), Dc(L + 1);
  for (int st = 0; st < 2; ++st) {
    V[st].assign(N, kBioInf);
    std::fill(Mp.begin(), Mp.end(), kBioInf);
    std::fill(Ip.begin(), Ip.end(), kBioInf);
    std::fill(Dp.begin(), Dp.end(), kBioInf);
    for (int64_t j = 0; j < N; ++j) {
      const uint8_t a = ix.ref[j];
      Mc[0] = 0;  // free start: row 0 before every column
      Ic[0] = Dc[0] = kBioInf;
      for (int i = 1; i <= L; ++i) {
        // M: read base i - 1 on reference j, after (i - 1, j - 1).
        const int prev = j == 0 ? (i == 1 ? 0 : kBioInf) : std::min({i == 1 ? 0 : Mp[i - 1], Ip[i - 1], Dp[i - 1]});
        Mc[i] = prev >= kBioInf ? kBioInf : prev + br.cost[st][i - 1][a];
        // I: read base i - 1 inserted after reference j (read base i - 2 on j).
        int iv = kBioInf;
        const int open = i >= 2 && Mc[i - 1] < kBioInf
                             ? Mc[i - 1] + m.go(std::max(ctx[j], j + 1 < N ? ctx[j + 1] : uint8_t(1)))
                             : kBioInf;
        const int ext = Ic[i - 1] < kBioInf ? Ic[i - 1] + m.ge() : kBioInf;
        iv = std::min(open, ext);
        Ic[i] = iv >= kBioInf ? kBioInf : iv + br.emit[st][i - 1];
        // D: reference j deleted after read base i - 1 on j - 1.
        const int dop = j > 0 && Mp[i] < kBioInf ? Mp[i] + m.go(ctx[j]) : kBioInf;
        const int dex = j > 0 && Dp[i] < kBioInf ? Dp[i] + m.ge() : kBioInf;
        Dc[i] = std::min(dop, dex);
      }
      V[st][j] = Mc[L];
      std::swap(Mp, Mc);
      std::swap(Ip, Ic);
      std::swap(Dp, Dc);
    }
  }
}

struct OLocus {
  int pen;
  int strand;
  int64_t end;
};

// Loci below `bound` from V, as bio_place picks them: end positions cheapest
// first (then strand, then position), each at least L from every cheaper one
// on its strand.
std::vector<OLocus> oracle_loci(const std::vector<int>* V, int L, int bound) {
  std::vector<OLocus> ends;
  for (int st = 0; st < 2; ++st)
    for (int64_t e = 0; e < static_cast<int64_t>(V[st].size()); ++e)
      if (V[st][e] < bound) ends.push_back({V[st][e], st, e});
  std::sort(ends.begin(), ends.end(), [](const OLocus& a, const OLocus& b) {
    if (a.pen != b.pen) return a.pen < b.pen;
    if (a.strand != b.strand) return a.strand < b.strand;
    return a.end < b.end;
  });
  std::vector<OLocus> loci;
  for (const OLocus& e : ends) {
    bool near = false;
    for (const OLocus& l : loci) near |= l.strand == e.strand && std::llabs(l.end - e.end) < L;
    if (!near) loci.push_back(e);
  }
  return loci;
}

int global_min(const std::vector<int>* V) {
  int m = kBioInf;
  for (int st = 0; st < 2; ++st)
    for (int v : V[st]) m = std::min(m, v);
  return m;
}

// ---------------------------------------------------------------------------
// The reference: random contigs with diverged and exact copies of a segment,
// homopolymers, short tandem repeats, a tandem array and N runs.
struct Genome {
  Reference ref;
  std::vector<int64_t> hp;   // starts of homopolymer runs (concatenated coordinates)
  std::vector<int64_t> str;  // starts of tandem repeats
  int64_t seg[4];            // the segment's copies: three diverged (chr1), one exact (chr2)
  int64_t free1, free2, free3, free4, free5;  // stretches kept free for planted reads and decoys
  int64_t array;             // the tandem array (25 copies of 37 bases)
};

Genome make_genome(std::mt19937_64& g) {
  std::vector<uint8_t> c1 = random_seq(g, 70000), c2 = random_seq(g, 80000);
  std::vector<int64_t> hp1, str1;
  for (int n = 0; n < 60; ++n) {  // homopolymers of 5-16 bases
    const int64_t at = 2000 + n * 700 + static_cast<int64_t>(g() % 300);
    const int len = 5 + static_cast<int>(g() % 12);
    const uint8_t b = rand_base(g);
    for (int x = 0; x < len; ++x) c1[at + x] = b;
    hp1.push_back(at);
  }
  for (int n = 0; n < 40; ++n) {  // tandem repeats, period 2-4, 4-14 copies
    const int64_t at = 2350 + n * 1000 + static_cast<int64_t>(g() % 200);
    const int p = 2 + static_cast<int>(g() % 3), copies = 4 + static_cast<int>(g() % 11);
    const std::vector<uint8_t> unit = random_seq(g, p);
    for (int x = 0; x < p * copies; ++x) c1[at + x] = unit[x % p];
    str1.push_back(at);
  }
  // A 1,500-base segment: three diverged copies in chr1 (0.4 %, 1 %, 2 %
  // substitutions), one exact copy in chr2.
  const std::vector<uint8_t> seg = random_seq(g, 1500);
  const int64_t sat[3] = {44000, 50000, 56000};
  const double div[3] = {0.004, 0.01, 0.02};
  for (int k = 0; k < 3; ++k) {
    std::vector<uint8_t> s = seg;
    for (auto& c : s)
      if (std::uniform_real_distribution<double>(0, 1)(g) < div[k]) c = static_cast<uint8_t>((c + 1 + g() % 3) % 4);
    std::copy(s.begin(), s.end(), c1.begin() + sat[k]);
  }
  std::copy(seg.begin(), seg.end(), c2.begin() + 3000);
  const std::vector<uint8_t> unit = random_seq(g, 37);  // tandem array
  for (int r = 0; r < 25; ++r) std::copy(unit.begin(), unit.end(), c1.begin() + 62000 + r * 37);
  std::fill(c1.begin() + 64000, c1.begin() + 64100, 4);  // N run
  for (int n = 0; n < 5; ++n) c1[66000 + n * 400] = 4;  // single Ns
  Genome G;
  Reference& ref = G.ref;
  ref.seq.owned.assign(kContigPad, 4);
  const std::vector<uint8_t>* cs[2] = {&c1, &c2};
  for (int c = 0; c < 2; ++c) {
    ref.names.push_back("chr" + std::to_string(c + 1));
    ref.offsets.push_back(ref.seq.owned.size());
    ref.lengths.push_back(cs[c]->size());
    ref.seq.owned.insert(ref.seq.owned.end(), cs[c]->begin(), cs[c]->end());
    ref.seq.owned.insert(ref.seq.owned.end(), kContigPad, 4);
  }
  const int64_t o1 = static_cast<int64_t>(ref.offsets[0]), o2 = static_cast<int64_t>(ref.offsets[1]);
  for (int64_t x : hp1) G.hp.push_back(o1 + x);
  for (int64_t x : str1) G.str.push_back(o1 + x);
  for (int k = 0; k < 3; ++k) G.seg[k] = o1 + sat[k];
  G.seg[3] = o2 + 3000;
  G.free1 = o2 + 8000;   // chr2 8,000-20,000: decoys
  G.free2 = o2 + 20000;  // chr2 20,000-32,000: decoys
  G.free3 = o2 + 34000;  // chr2 34,000-50,000: reads at the band's reach
  G.free4 = o2 + 50000;  // chr2 50,000-62,000: their decoys, and the tandem array's
  G.free5 = o2 + 62000;  // chr2 62,000-79,000: pairs at the list bound
  G.array = o1 + 62000;
  return G;
}

// ---------------------------------------------------------------------------
// Reads.
struct TestRead {
  std::vector<uint8_t> codes;
  std::string qual;
  int64_t truth;    // reference position of the first read base on the forward strand, or -1
  int strand;
  int family;       // 0 sampled, 1 quality decides, 2 unwitnessed, 3 indel in a run, 4 repeat
};

char qchar(int q) { return static_cast<char>(33 + q); }

// NovaSeq-like binned qualities: mostly 37, more 23/12/2 towards the 3' end,
// and sometimes a Q2 tail.
std::string novaseq_quals(std::mt19937_64& g, int L) {
  std::string q(L, qchar(37));
  std::uniform_real_distribution<double> U(0, 1);
  for (int i = 0; i < L; ++i) {
    const double t = static_cast<double>(i) / L, low = 0.03 + 0.15 * t * t * t, u = U(g);
    if (u < low * 0.2) q[i] = qchar(2);
    else if (u < low * 0.5) q[i] = qchar(12);
    else if (u < low) q[i] = qchar(23);
  }
  if (U(g) < 0.15) {
    const int tail = 5 + static_cast<int>(g() % 25);
    for (int i = L - tail; i < L; ++i) q[i] = qchar(2);
  }
  return q;
}

// A read sampled from the sample genome: the reference at `at` with SNVs
// (two thirds transitions) and indels (ten times likelier next to a run),
// then sequencing errors drawn from its qualities; strand 1 is reverse
// complemented (qualities reversed with it).
TestRead sample_read(std::mt19937_64& g, const IndexView& ix, int64_t at, int strand, int L) {
  std::uniform_real_distribution<double> U(0, 1);
  std::vector<uint8_t> s;
  for (int64_t j = at; static_cast<int>(s.size()) < L + 40 && j < static_cast<int64_t>(ix.ref_len); ++j) {
    const uint8_t c = ix.ref[j];
    const bool run = c < 4 && ref_at(ix, j - 1) == c && ref_at(ix, j + 1) == c;
    const double ind = run ? 0.01 : 0.0008;
    const double u = U(g);
    if (u < ind / 2) continue;                                        // deletion
    if (u < ind) s.push_back(run ? c : rand_base(g));                 // insertion
    if (c < 4 && U(g) < 0.004) s.push_back(U(g) < 2.0 / 3 ? transition(c) : static_cast<uint8_t>((c + 1 + 2 * (g() % 2)) % 4));
    else s.push_back(c);
  }
  s.resize(L, 0);
  TestRead r;
  r.qual = novaseq_quals(g, L);
  r.codes = strand ? revcomp(s) : s;
  for (int i = 0; i < L; ++i) {
    const double e = std::pow(10.0, -(r.qual[i] - 33) / 10.0);
    if (r.codes[i] < 4 && U(g) < e) r.codes[i] = static_cast<uint8_t>((r.codes[i] + 1 + g() % 3) % 4);
  }
  r.truth = at;
  r.strand = strand;
  r.family = 0;
  return r;
}

// ---------------------------------------------------------------------------
// Lemma B: bio_cover against the exact least cost of a transcript touching
// every required part. A transcript of the read (any reference) is a
// sequence of aligned bases (match 0, mismatch mmin), insertions (open or
// extend, plus the base's emission) and deletions (open or extend), starting
// and ending with an aligned base, no insertion next to a deletion.
int exact_cover(const BioModel& m, const BioRead& br, int P, uint32_t req) {
  const int L = br.L;
  std::vector<int> part(L);
  for (int j = 0; j < P; ++j) {
    int off, len;
    part_geometry(L, P, j, &off, &len);
    for (int i = off; i < off + len; ++i) part[i] = j;
  }
  const int NM = 1 << P;
  // best[i][mode][mask]: i bases consumed; mode 0 aligned, 1 insertion, 2 deletion.
  std::vector<int> best(static_cast<size_t>(L + 1) * 3 * NM, kBioInf);
  auto at = [&](int i, int mode, int mask) -> int& { return best[(static_cast<size_t>(i) * 3 + mode) * NM + mask]; };
  auto relax = [&](int i, int mode, int mask, int v) { if (v < at(i, mode, mask)) at(i, mode, mask) = v; };
  for (int i = 0; i < L; ++i) {
    // Deletions before base i (between i - 1 and i), from an aligned base.
    if (i >= 1)
      for (int mask = 0; mask < NM; ++mask) {
        const int v = at(i, 0, mask);
        if (v >= kBioInf) continue;
        const int t = part[i - 1] == part[i] ? 1 << part[i] : 0;
        relax(i, 2, mask | t, v + m.go_min());
      }
    for (int mode = 0; mode < 3; ++mode)
      for (int mask = 0; mask < NM; ++mask) {
        const int v = i == 0 ? (mode == 0 && mask == 0 ? 0 : kBioInf) : at(i, mode, mask);
        if (v >= kBioInf) continue;
        const int t = 1 << part[i];
        relax(i + 1, 0, mask, v);                         // match
        relax(i + 1, 0, mask | t, v + br.mmin[0][i]);     // mismatch
        if (i >= 1 && i + 1 < L && mode != 2)             // insertion of base i
          relax(i + 1, 1, mask | t, v + (mode == 1 ? m.ge() : m.go_min()) + br.emit[0][i]);
      }
  }
  int r = kBioInf;
  for (int mask = 0; mask < NM; ++mask)
    if ((mask & req) == req) r = std::min(r, at(L, 0, mask));
  return r;
}

void test_lemma_b(std::mt19937_64& g) {
  int cases = 0, tight = 0;
  for (int n = 0; n < 4000; ++n) {
    BioParams bp;
    bp.indel = std::pow(10.0, -1 - static_cast<double>(g() % 40) / 10);
    bp.gap_ext = 0.1 + 0.8 * (g() % 100) / 100.0;
    bp.str_cap = static_cast<int>(g() % 40);
    const BioModel m(bp);
    const int P = 2 + static_cast<int>(g() % 4);
    const int L = P * 3 + static_cast<int>(g() % (P * 3));
    std::vector<uint8_t> codes = random_seq(g, L);
    std::string qual(L, ' ');
    for (auto& c : qual) c = qchar(static_cast<int>(g() % 2 ? 2 + g() % 10 : 2 + g() % 40));
    BioRead br;
    bio_read(m, codes.data(), qual.data(), L, br);
    for (uint32_t req = 1; req < (1u << P); ++req) {
      const int c = bio_cover(m, br, 0, P, req), e = exact_cover(m, br, P, req);
      ++cases;
      tight += c == e;
      CHECK(c <= e, "Lemma B: cover %d above the exact %d (L %d, P %d, req %u)", c, e, L, P, req);
    }
  }
  std::fprintf(stderr, "Lemma B: %d cases, cover <= exact on all, equal on %d\n", cases, tight);
  CHECK(tight > cases / 2, "Lemma B: the bound is rarely tight (%d of %d)", tight, cases);
}

// ---------------------------------------------------------------------------
// One placement: what the checks compare across settings and orders.
struct Placed {
  Result r;
  BioOut o;
  bool ok;
};

bool same(const Placed& a, const Placed& b) {
  if (a.ok != b.ok || a.r.certified != b.r.certified || a.r.tier != b.r.tier || a.r.reason != b.r.reason) return false;
  if (!a.ok) return true;
  if (a.r.ref_pos != b.r.ref_pos || a.r.strand != b.r.strand || a.r.n_cigar != b.r.n_cigar) return false;
  for (int x = 0; x < a.r.n_cigar; ++x)
    if (a.r.cigar[x] != b.r.cigar[x]) return false;
  return a.o.p1 == b.o.p1 && a.o.p2 == b.o.p2 && a.o.yh == b.o.yh && a.o.n_best == b.o.n_best &&
         a.o.mapq == b.o.mapq && a.o.loci.size() == b.o.loci.size();
}

// The reported alignment with gap c moved one base left, if that is still an
// alignment: the base before the gap moves after it.
bool shift_left(const Result& r, int c, Result& out) {
  if (c == 0 || (r.cigar[c - 1] & 0xF) != kOpM) return false;
  out = r;
  out.cigar[c - 1] -= 1u << 4;
  if (c + 1 < r.n_cigar && (r.cigar[c + 1] & 0xF) == kOpM) out.cigar[c + 1] += 1u << 4;
  else return false;
  if ((out.cigar[c - 1] >> 4) == 0) return false;  // the gap would start the alignment
  return true;
}

struct Setting {
  const char* name;
  int q, s;
  BioParams bp;
};

}  // namespace

int main(int argc, char** argv) {
  const uint64_t seed = argc > 1 ? std::strtoull(argv[1], nullptr, 10) : 20261029;
  std::fprintf(stderr, "seed %llu\n", static_cast<unsigned long long>(seed));
  std::mt19937_64 g(seed);
  {
    std::mt19937_64 gl(seed ^ 0xb10b10ULL);
    test_lemma_b(gl);
  }
  Genome G = make_genome(g);
  const BioModel m0;  // the default model: reads are planted against its costs
  const int L = 150;
  std::vector<TestRead> reads;

  // Family 1, the qualities decide: at the true locus the read has three
  // mismatches at Q2 bases (in three parts); a planted decoy differs from the
  // read at one Q37 base only. Fewest edits picks the decoy, the likelihood
  // the truth.
  {
    auto& R = G.ref.seq.owned;
    for (int n = 0; n < 30; ++n) {
      const int64_t a = static_cast<int64_t>(G.ref.offsets[0]) + 1000 + static_cast<int64_t>(g() % 40000);
      std::vector<uint8_t> s(R.begin() + a, R.begin() + a + L);
      if (std::count(s.begin(), s.end(), 4)) { --n; continue; }
      std::string q(L, qchar(37));
      for (int y = 0; y < 3; ++y) {
        const int x = 10 + 50 * y + static_cast<int>(g() % 30);
        s[x] = transition(s[x]);
        q[x] = qchar(2);
      }
      std::vector<uint8_t> d = s;
      const int x = 5 + static_cast<int>(g() % 140);
      if (q[x] == qchar(2)) { --n; continue; }
      d[x] = static_cast<uint8_t>((d[x] + 1) % 4);
      std::copy(d.begin(), d.end(), R.begin() + G.free1 + 400 * n);
      TestRead r;
      r.strand = static_cast<int>(g() % 2);
      r.codes = r.strand ? revcomp(s) : s;
      r.qual = r.strand ? std::string(q.rbegin(), q.rend()) : q;
      r.truth = a;
      r.family = 1;
      reads.push_back(r);
    }
  }
  // Family 2, at the horizon: at the true locus every part has one mismatch
  // at its Q2 base, inside every q-mer the index looks up for that part (so
  // no part has a candidate there), and the alignment touches every part and
  // costs exactly cover(S). A planted decoy equals the read except at j
  // extra Q2 bases of the first part, j = P - 1, P or P + 1 (index q = 14,
  // s = 4: P = 8): just below, at and above the horizon.
  {
    auto& R = G.ref.seq.owned;
    const int P = part_count(L, 14, 4);
    for (int n = 0; n < 24; ++n) {
      const int64_t a = static_cast<int64_t>(G.ref.offsets[0]) + 1000 + static_cast<int64_t>(g() % 40000);
      std::vector<uint8_t> s(R.begin() + a, R.begin() + a + L);
      if (std::count(s.begin(), s.end(), 4)) { --n; continue; }
      std::string q(L, qchar(37));
      int off0, len0;
      part_geometry(L, P, 0, &off0, &len0);
      for (int j = 0; j < P; ++j) {
        int off, len;
        part_geometry(L, P, j, &off, &len);
        const int x = off + len / 2;  // in every q-mer of the part (q = 14, s = 4): mismatched at the truth
        s[x] = transition(s[x]);
        q[x] = qchar(2);
      }
      const int extra = P - 1 + n % 3;  // Q2 bases the decoy changes (in part 0 and on)
      std::vector<uint8_t> d = s;
      for (int y = 0; y < extra; ++y) {
        const int x = 2 * y;  // never the middle of part 0 (9)
        q[x] = qchar(2);
        d[x] = transition(d[x]);
      }
      std::copy(d.begin(), d.end(), R.begin() + G.free2 + 400 * n);
      TestRead r;
      r.strand = static_cast<int>(g() % 2);
      r.codes = r.strand ? revcomp(s) : s;
      r.qual = r.strand ? std::string(q.rbegin(), q.rend()) : q;
      r.truth = a;
      r.family = 2;
      reads.push_back(r);
    }
  }
  // Family 6, the band's reach: a deletion of 6 bases (the tight setting's
  // reach, gap_cap 6) opening in a 12-base homopolymer (go_min) near one end
  // of the read, where no part has a candidate at that end's diagonal. The
  // alignment lies 6 diagonals beyond every candidate and costs
  // 22 + 5 * 3 = 37, below the tight horizon 22 + 6 * 3 = 40. A planted decoy
  // differs from the read by one Q60 transversion (38): a band one diagonal
  // short on that side certifies the decoy. Even n: the deletion after the
  // first 12 bases (6 before the run, 6 of it), the head 6 diagonals left of
  // every candidate; odd n: before the last 14 (the run's last 6, 8 after
  // it), the tail 6 diagonals right.
  {
    auto& R = G.ref.seq.owned;
    for (int n = 0; n < 16; ++n) {
      const int64_t r = G.free3 + 1000 * n + 300;
      const uint8_t b = rand_base(g);
      for (int x = 0; x < 12; ++x) R[r + x] = b;
      auto not_b = [&](int64_t at) {
        if (R[at] == b) R[at] = static_cast<uint8_t>((b + 1 + g() % 3) % 4);
      };
      std::vector<uint8_t> s;
      if (n % 2 == 0) {
        for (int64_t x = r - 6; x < r; ++x) not_b(x);
        not_b(r + 12);
        s.assign(R.begin() + r - 6, R.begin() + r + 6);
        s.insert(s.end(), R.begin() + r + 12, R.begin() + r + 12 + (L - 12));
      } else {
        not_b(r - 1);
        not_b(r + 12);
        s.assign(R.begin() + r - 136, R.begin() + r);
        s.insert(s.end(), R.begin() + r + 6, R.begin() + r + 20);
      }
      std::string q(L, qchar(37));
      q[70] = qchar(60);
      std::vector<uint8_t> d = s;
      d[70] = static_cast<uint8_t>((s[70] + 1 + 2 * (g() % 2)) % 4);  // a transversion
      const int dc = m0.mis(60, d[70], s[70]);
      CHECK(dc > 37 && dc < 40, "reach family: the decoy costs %d, not between the truth (37) and the horizon (40)", dc);
      std::copy(d.begin(), d.end(), R.begin() + G.free4 + 300 * n);
      TestRead t;
      t.codes = s;
      t.qual = q;
      t.strand = 0;
      t.truth = n % 2 == 0 ? r - 6 : r - 136;
      t.family = 6;
      reads.push_back(t);
    }
  }
  // Family 7, bands too wide to align: reads from the tandem array, whose
  // candidates span ~1 kb; in the tight setting (max_band 40) that band is
  // bounded, not aligned, and the horizon falls to its bound. A planted decoy
  // differs from the read by one Q37 transition (31): not lowering the
  // horizon certifies the decoy.
  {
    auto& R = G.ref.seq.owned;
    for (int n = 0; n < 8; ++n) {
      const int64_t a = G.array + static_cast<int64_t>(g() % (25 * 37 - L));
      std::vector<uint8_t> s(R.begin() + a, R.begin() + a + L), d = s;
      d[75] = transition(d[75]);
      std::copy(d.begin(), d.end(), R.begin() + G.free4 + 6000 + 300 * n);
      TestRead t;
      t.codes = s;
      t.qual = std::string(L, qchar(37));
      t.strand = 0;
      t.truth = a;
      t.family = 7;
      reads.push_back(t);
    }
  }
  // Pairs at the list bound (the tight setting keeps 3 loci per read): mate 1
  // occurs exactly at four copies c0 < c1 < c2 < c3 of a 150-base sequence,
  // the truth at c3; mate 2 is unique sequence ending 400 bases after c3's
  // start (a proper pair at the mean fragment length: insert term 0), and a
  // copy of it with one base changed sits at the same distance from c2 (mate 2
  // costs 2 there, a Q2 transition). Mate 1's list keeps three of the four
  // tied copies, so the pair optimum (0, at c3) is not provable; bounding the
  // left-out copy too high, or the pair horizon by 3, certifies c2's pair at 2.
  std::vector<TestRead> list_pairs;
  {
    auto& R = G.ref.seq.owned;
    for (int n = 0; n < 8; ++n) {
      const int64_t base = G.free5 + 2000 * n;
      int64_t c[4];
      for (int k = 0; k < 4; ++k) c[k] = base + 450 * k;
      for (int k = 0; k < 3; ++k) std::copy(R.begin() + c[3], R.begin() + c[3] + L, R.begin() + c[k]);
      std::vector<uint8_t> f4(R.begin() + c[3] + 250, R.begin() + c[3] + 400), f3 = f4;
      f3[75] = transition(f3[75]);
      std::copy(f3.begin(), f3.end(), R.begin() + c[2] + 250);
      TestRead m1, m2;
      m1.codes.assign(R.begin() + c[3], R.begin() + c[3] + L);
      m1.qual = std::string(L, qchar(37));
      m1.strand = 0;
      m1.truth = c[3];
      m2.codes = revcomp(f4);
      m2.qual = std::string(L, qchar(37));
      m2.qual[L - 1 - 75] = qchar(2);
      m2.strand = 1;
      m2.truth = c[3] + 250;
      m1.family = m2.family = 8;
      list_pairs.push_back(m1);
      list_pairs.push_back(m2);
    }
  }
  // Family 9, clusters left unaligned: the read occurs exactly at its locus
  // (Q37) except for one Q12 base in the middle of part 3 (inside every q-mer
  // looked up for it); a planted decoy has a transition there (16, a part's
  // whole bound). In the tight setting (delta 12) the decoy's cluster, bound
  // 16 >= 0 + 12, is not aligned and the horizon is its bound, 16; a horizon
  // above it would claim a list without the decoy.
  {
    auto& R = G.ref.seq.owned;
    const int P = part_count(L, 14, 4);
    int off3, len3;
    part_geometry(L, P, 3, &off3, &len3);
    const int x = off3 + len3 / 2;
    for (int n = 0; n < 8; ++n) {
      const int64_t a = static_cast<int64_t>(G.ref.offsets[0]) + 1000 + static_cast<int64_t>(g() % 40000);
      std::vector<uint8_t> s(R.begin() + a, R.begin() + a + L);
      if (std::count(s.begin(), s.end(), 4)) { --n; continue; }
      std::string q(L, qchar(37));
      q[x] = qchar(12);
      std::vector<uint8_t> d = s;
      d[x] = transition(d[x]);
      CHECK(m0.mis(12, d[x], s[x]) >= 12 && m0.mis(12, d[x], s[x]) < m0.go_min(),
            "skipped-cluster family: the decoy costs %d", m0.mis(12, d[x], s[x]));
      std::copy(d.begin(), d.end(), R.begin() + G.free4 + 9000 + 300 * n);
      TestRead t;
      t.codes = s;
      t.qual = q;
      t.strand = 0;
      t.truth = a;
      t.family = 9;
      reads.push_back(t);
    }
  }
  // The index (built after the decoys are planted).
  const std::vector<Setting> settings = [] {
    std::vector<Setting> v;
    v.push_back({"q14 s4", 14, 4, BioParams()});
    BioParams tight;
    tight.gap_cap = 6;     // reach shorter than the horizon allows
    tight.delta = 12;      // most clusters bounded, not aligned
    tight.max_band = 40;   // and the wide ones too
    tight.max_loci = 3;
    v.push_back({"q14 s4 tight", 14, 4, tight});
    v.push_back({"q12 s1", 12, 1, BioParams()});
    return v;
  }();
  std::vector<Index> idx;
  for (const Setting& s : settings) idx.push_back(Index::build(G.ref, s.q, s.s, 4));
  std::vector<IndexView> views;
  for (size_t k = 0; k < settings.size(); ++k) views.push_back(idx[k].view(G.ref));
  const IndexView& ix = views[0];

  // Family 0: sampled reads, anywhere (repeats, runs, tandem array, N).
  for (int n = 0; n < 160; ++n) {
    const int c = n % 5 == 4 ? 1 : 0;
    const int64_t a = static_cast<int64_t>(G.ref.offsets[c]) +
                      static_cast<int64_t>(g() % (G.ref.lengths[c] - L - 60));
    reads.push_back(sample_read(g, ix, a, static_cast<int>(g() % 2), L));
  }
  // Family 3: an indel in a homopolymer or tandem repeat (left-alignment).
  for (int n = 0; n < 50; ++n) {
    const int64_t run = n % 2 ? G.hp[g() % G.hp.size()] : G.str[g() % G.str.size()];
    const int64_t a = run - 40 - static_cast<int64_t>(g() % 60);
    std::vector<uint8_t> s(G.ref.seq.owned.begin() + a, G.ref.seq.owned.begin() + a + L + 10);
    const int x = static_cast<int>(run - a) + 2;
    const int k = 1 + static_cast<int>(g() % 3);
    if (n % 4 < 2) s.erase(s.begin() + x, s.begin() + x + k);
    else s.insert(s.begin() + x, s.begin() + x, s.begin() + x + k);  // a copy: an expansion
    s.resize(L);
    TestRead r;
    r.qual = std::string(L, qchar(37));
    r.strand = static_cast<int>(g() % 2);
    r.codes = r.strand ? revcomp(s) : s;
    r.truth = a;
    r.family = 3;
    reads.push_back(r);
  }
  // Family 4: the segment's copies (ties at the exact copy, near-ties at the
  // diverged ones).
  for (int n = 0; n < 40; ++n) {
    const int64_t a = G.seg[n % 4] + static_cast<int64_t>(g() % (1500 - L));
    TestRead r = sample_read(g, ix, a, static_cast<int>(g() % 2), L);
    r.family = 4;
    reads.push_back(r);
  }
  // Pairs: fragments of 250-550 bases, mate 1 forward from f, mate 2 the
  // reverse complement of the fragment's end; some in the segment's copies,
  // where the unique mate decides.
  const size_t first_pair = reads.size();
  for (int n = 0; n < 60; ++n) {
    const int F = 250 + static_cast<int>(g() % 300);
    int64_t f;
    if (n % 3 == 0) f = G.seg[n % 4] + 1500 - L - static_cast<int64_t>(g() % 100);  // mate 1 in a copy
    else f = static_cast<int64_t>(G.ref.offsets[0]) + 500 + static_cast<int64_t>(g() % 60000);
    TestRead m1 = sample_read(g, ix, f, 0, L), m2 = sample_read(g, ix, f + F - L, 1, L);
    m1.family = m2.family = 5;
    reads.push_back(m1);
    reads.push_back(m2);
  }
  reads.insert(reads.end(), list_pairs.begin(), list_pairs.end());
  const int npairs = static_cast<int>((reads.size() - first_pair) / 2);
  std::fprintf(stderr, "reads: %zu (%d pairs); genome %llu bases\n", reads.size(), npairs,
               static_cast<unsigned long long>(ix.ref_len));

  // Map every read in every setting, forward and in reverse order.
  const size_t nr = reads.size(), ns = settings.size();
  std::vector<std::vector<Placed>> res(ns, std::vector<Placed>(nr));
  std::vector<BioModel> models;
  for (const Setting& s : settings) models.push_back(BioModel(s.bp));
  for (size_t k = 0; k < ns; ++k) {
    Params p;
    p.k = KMAX;
    p.budget = HOST_BUDGET_MAX;
    std::unique_ptr<Workspace> ws(new Workspace), ws2(new Workspace);
    std::unique_ptr<BioScratch> sc(new BioScratch), sc2(new BioScratch);
    for (size_t i = 0; i < nr; ++i) {
      Placed& P = res[k][i];
      const TestRead& t = reads[i];
      P.ok = bio_place(views[k], models[k], p, t.codes.data(), t.qual.data(), L, 0x9e3779b97f4a7c15ULL * (i + 1), *ws,
                       *sc, P.r, P.o);
    }
    int diff = 0;
    for (size_t i = nr; i-- > 0;) {
      Placed Q;
      const TestRead& t = reads[i];
      Q.ok = bio_place(views[k], models[k], p, t.codes.data(), t.qual.data(), L, 0x9e3779b97f4a7c15ULL * (i + 1), *ws2,
                       *sc2, Q.r, Q.o);
      diff += !same(Q, res[k][i]);
    }
    CHECK(diff == 0, "%s: %d reads differ in reverse order on another workspace", settings[k].name, diff);
  }

  // The oracle, read by read (pairs: both mates at once), and every check.
  std::vector<std::vector<uint8_t>> ctx(1);
  ctx[0].resize(ix.ref_len);
  bio_str_context(ix, 0, static_cast<int64_t>(ix.ref_len), ctx[0].data());
  std::atomic<int> certified[3] = {{0}, {0}, {0}}, completeness[3] = {{0}, {0}, {0}};
  std::atomic<int> n_pairs_ok{0}, n_moved{0}, n_ties{0}, n_shift_checked{0};
  // Placement against the truth (default setting): bio, fewest edits.
  std::atomic<int> at_truth_bio[10] = {}, at_truth_edit[10] = {}, n_fam[10] = {};
  auto at_truth = [&](const TestRead& t, const Result& r) {
    if (t.truth < 0) return false;
    return r.strand == t.strand && std::llabs(r.ref_pos - t.truth) <= 20;
  };
  const size_t units = first_pair + npairs;  // single reads, then pairs
  parallel(units, [&](size_t u) {
    const bool pair = u >= first_pair;
    const size_t i0 = pair ? first_pair + 2 * (u - first_pair) : u;
    const int nm = pair ? 2 : 1;
    std::vector<int> V[2][2];
    BioRead br[2];
    for (int x = 0; x < nm; ++x) {
      bio_read(m0, reads[i0 + x].codes.data(), reads[i0 + x].qual.data(), L, br[x]);
      oracle(ix, m0, br[x], ctx[0], V[x]);
    }
    for (int x = 0; x < nm; ++x) {
      const size_t i = i0 + x;
      const TestRead& t = reads[i];
      const int gmin = global_min(V[x]);
      for (size_t k = 0; k < ns; ++k) {
        const Placed& P = res[k][i];
        if (settings[k].bp.snv != m0.params().snv) continue;  // every setting shares the model's costs
        if (!P.ok) continue;
        ++certified[k];
        const BioOut& o = P.o;
        const int pen = bio_rescore(views[k], models[k], br[x], P.r);
        CHECK(pen == o.loci[o.chosen].pen, "%s read %zu: alignment rescores to %d, reported %d", settings[k].name, i,
              pen, o.loci[o.chosen].pen);
        CHECK(o.p1 == gmin, "%s read %zu (family %d): certified %d, oracle minimum %d (YH %d)", settings[k].name, i,
              t.family, o.p1, gmin, o.yh);
        // Completeness: the oracle's loci below YH, value for value.
        const std::vector<OLocus> ol = oracle_loci(V[x], L, o.yh);
        bool eq = ol.size() >= o.loci.size();
        for (size_t y = 0; eq && y < o.loci.size(); ++y)
          eq = ol[y].pen == o.loci[y].pen && ol[y].strand == o.loci[y].strand && ol[y].end == o.loci[y].end;
        // YA: the span of the end positions as cheap as the locus, within L.
        for (size_t y = 0; eq && y < o.loci.size(); ++y) {
          const BioLocus& l = o.loci[y];
          const std::vector<int>& Vs = V[x][l.strand];
          int64_t a0 = l.end, a1 = l.end;
          for (int64_t e = std::max<int64_t>(0, l.end - L + 1); e < std::min<int64_t>(Vs.size(), l.end + L); ++e)
            if (Vs[e] == l.pen) {
              a0 = std::min(a0, e);
              a1 = std::max(a1, e);
            }
          eq = l.amb == a1 - a0;
        }
        if (eq && o.loci.size() < ol.size()) eq = o.list_bound <= ol[o.loci.size()].pen;
        if (eq && o.list_bound == o.yh) eq = ol.size() == o.loci.size();
        completeness[k] += eq;
        CHECK(eq, "%s read %zu (family %d): %zu loci below YH %d, the oracle has %zu (first: %d/%d vs %d/%d)",
              settings[k].name, i, t.family, o.loci.size(), o.yh, ol.size(), o.loci.empty() ? -1 : o.loci[0].pen,
              o.loci.empty() ? -1 : static_cast<int>(o.loci[0].end), ol.empty() ? -1 : ol[0].pen,
              ol.empty() ? -1 : static_cast<int>(ol[0].end));
        int tied = 0;
        for (const OLocus& l : ol) tied += l.pen == gmin;
        CHECK(o.n_best == tied, "%s read %zu: %d loci tied, oracle %d", settings[k].name, i, o.n_best, tied);
        const int p2 = tied > 1 ? gmin : (ol.size() > 1 ? ol[1].pen : -1);
        CHECK(o.p2 == p2, "%s read %zu: second %d, oracle %d", settings[k].name, i, o.p2, p2);
        if (k == 0) {
          n_ties += o.n_best > 1;
          n_moved += o.moved;
        }
        // Canonical: no gap moves one base left at the same penalty or less.
        for (int c = 0; c < P.r.n_cigar; ++c) {
          const uint32_t op = P.r.cigar[c] & 0xF;
          if (op != kOpI && op != kOpD) continue;
          Result s;
          if (!shift_left(P.r, c, s)) continue;
          const int ps = bio_rescore(views[k], models[k], br[x], s);
          if (ps < 0) continue;
          ++n_shift_checked;
          CHECK(ps > pen, "%s read %zu: gap %d moves left at penalty %d (reported %d)", settings[k].name, i, c, ps,
                pen);
        }
      }
      // Truth, default setting.
      if (!pair || true) {
        const Placed& P = res[0][i];
        const int f = t.family;
        if (P.ok) {
          ++n_fam[f];
          at_truth_bio[f] += at_truth(t, P.r);
          // The fewest-edits placement of the same read (default passes' core).
          Params pe;
          pe.k = KMAX;
          pe.budget = HOST_BUDGET_MAX;
          pe.mapq_limit = 8;
          std::unique_ptr<Workspace> w(new Workspace);
          Result re;
          process_read(ix, pe, t.codes.data(), L, 0x9e3779b97f4a7c15ULL * (i + 1), *w, re);
          at_truth_edit[f] += re.certified == 1 && at_truth(t, re);
        }
      }
    }
    if (!pair) return;
    // Pairs: a joint placement must be the oracle's best pair of loci.
    for (size_t k = 0; k < ns; ++k) {
      const Placed &A = res[k][i0], &B = res[k][i0 + 1];
      auto contig_of = [&](int64_t pos) { return G.ref.contig_of(pos); };
      const BioPair bp = bio_pair(settings[k].bp, A.o, L, B.o, L, 850, contig_of, 0x1234567ULL * (i0 + 1));
      if (!bp.ok) continue;
      if (k == 0) ++n_pairs_ok;
      const int g1 = global_min(V[0]), g2 = global_min(V[1]);
      const std::vector<OLocus> l1 = oracle_loci(V[0], L, bp.pen - g2 + 1), l2 = oracle_loci(V[1], L, bp.pen - g1 + 1);
      int best = kBioInf;
      for (const OLocus& a : l1)
        for (const OLocus& b : l2) {
          const BioLocus la{a.pen, static_cast<uint8_t>(a.strand), a.end, 0, 0, 0};
          const BioLocus lb{b.pen, static_cast<uint8_t>(b.strand), b.end, 0, 0, 0};
          best = std::min(best, a.pen + b.pen + bio_insert(settings[k].bp, la, L, lb, L, 850, contig_of));
        }
      CHECK(bp.pen == best, "%s pair %zu: joint placement %d, oracle's best pair %d", settings[k].name, i0, bp.pen,
            best);
      CHECK(A.o.loci[bp.x].pen + B.o.loci[bp.y].pen + bp.ins == bp.pen, "%s pair %zu: chosen loci do not add up",
            settings[k].name, i0);
    }
  });

  for (size_t k = 0; k < ns; ++k)
    std::fprintf(stderr, "%s: %d of %zu reads ML-certified, loci complete on %d\n", settings[k].name,
                 certified[k].load(), nr, completeness[k].load());
  std::fprintf(stderr, "default: %d ties, %d placed elsewhere than the fewest edits; %d gap shifts checked; %d pairs "
                       "placed jointly\n",
               n_ties.load(), n_moved.load(), n_shift_checked.load(), n_pairs_ok.load());
  static const char* const fam[10] = {"sampled",        "qualities decide", "at the horizon",
                                      "indel in a run", "repeat copies",    "pair mates",
                                      "band reach",     "wide band",        "list-bound pairs",
                                      "skipped cluster"};
  for (int f = 0; f < 10; ++f)
    if (n_fam[f])
      std::fprintf(stderr, "  %-16s certified %3d: at truth ML %3d, fewest edits %3d\n", fam[f], n_fam[f].load(),
                   at_truth_bio[f].load(), at_truth_edit[f].load());
  CHECK(n_fam[1] >= 20 && at_truth_bio[1] == n_fam[1] && at_truth_edit[1] == 0,
        "qualities-decide reads: ML should place all %d at the truth (%d) and fewest edits none (%d)", n_fam[1].load(),
        at_truth_bio[1].load(), at_truth_edit[1].load());
  CHECK(n_ties > 0, "no tie certified: the test does not exercise ties");
  CHECK(n_shift_checked > 20, "too few gaps checked for left-alignment (%d)", n_shift_checked.load());
  CHECK(n_pairs_ok > npairs / 2, "only %d of %d pairs placed jointly", n_pairs_ok.load(), npairs);
  CHECK(certified[1] > 0 && certified[1] < certified[0], "tight setting: %d certified vs %d (pruning not exercised)",
        certified[1].load(), certified[0].load());
  if (failures) {
    std::fprintf(stderr, "%d check(s) FAILED\n", failures.load());
    return 1;
  }
  std::fprintf(stderr, "all bio checks passed\n");
  return 0;
}
