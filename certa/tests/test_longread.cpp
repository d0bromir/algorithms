// Oracle tests of the long-read certificates (certa/longread.h): on a
// synthetic genome with planted structural variants, every proven local
// optimum and every split (jump) optimum or interval is checked against an
// exact brute-force DP over the whole genome, both strands.
//
//   test_longread [seed]
//
// The brute force for the jump objective: one pass over the read, for each
// genome position on each strand a Gotoh row, a piece starting from the best
// prefix value minus P (or 0 / -5 for the first piece), the prefix value
// updated after every read base. With P infinite it is the local optimum.
#include <algorithm>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "certa/core.h"
#include "certa/longread.h"
#include "index.h"

using namespace certa;

namespace {

int failures = 0;
#define CHECK(cond, ...)                                          \
  do {                                                            \
    if (!(cond)) {                                                \
      if (++failures <= 20) {                                     \
        std::fprintf(stderr, "FAIL %s:%d: ", __FILE__, __LINE__); \
        std::fprintf(stderr, __VA_ARGS__);                        \
        std::fprintf(stderr, "\n");                               \
      }                                                           \
    }                                                             \
  } while (0)

std::vector<uint8_t> random_seq(std::mt19937_64& g, size_t n) {
  std::vector<uint8_t> s(n);
  for (auto& c : s) c = static_cast<uint8_t>(g() % 4);
  return s;
}

std::vector<uint8_t> rc(const std::vector<uint8_t>& v) {
  std::vector<uint8_t> o(v.size());
  for (size_t i = 0; i < v.size(); ++i) o[i] = lr_comp(v[v.size() - 1 - i]);
  return o;
}

// Substitutions, short insertions and deletions at about `rate` per base.
std::vector<uint8_t> noisy(std::mt19937_64& g, const std::vector<uint8_t>& s, double rate) {
  std::vector<uint8_t> o;
  std::uniform_real_distribution<double> u(0, 1);
  for (size_t i = 0; i < s.size(); ++i) {
    const double x = u(g);
    if (x < rate / 3) o.push_back(static_cast<uint8_t>((s[i] + 1 + g() % 3) % 4));
    else if (x < 2 * rate / 3) continue;
    else if (x < rate) {
      o.push_back(s[i]);
      o.push_back(static_cast<uint8_t>(g() % 4));
    } else {
      o.push_back(s[i]);
    }
  }
  return o;
}

// The exact optimum of the jump objective over the whole sequence `ref`
// (both strands: the reverse-complemented reference read forward), by a
// full DP. P >= 1e8: the local optimum.
long brute_jump(const std::vector<uint8_t>& ref, const std::vector<uint8_t>& rd, long P) {
  const int L = static_cast<int>(rd.size());
  const int64_t N = static_cast<int64_t>(ref.size());
  const std::vector<uint8_t> t1 = rc(ref);
  const std::vector<uint8_t>* T[2] = {&ref, &t1};
  const int NEG = INT_MIN / 4;
  std::vector<int> H[2], I[2], Hn(N + 1), In(N + 1);
  for (int s = 0; s < 2; ++s) {
    H[s].assign(N + 1, NEG);
    I[s].assign(N + 1, NEG);
  }
  long M = LONG_MIN / 4, end_at_L = LONG_MIN / 4;
  for (int f = 0; f < L; ++f) {
    const long Gl = std::max(f > 0 ? -(long)kClip : 0L, P >= 100000000 ? LONG_MIN / 4 : M - P);
    const int G = static_cast<int>(std::max<long>(Gl, NEG));
    long F = LONG_MIN / 4;
    for (int s = 0; s < 2; ++s) {
      const std::vector<uint8_t>& t = *T[s];
      // column j (1-based): target base j - 1; H[j] at row f - 1, diagonal moves.
      int d = NEG;
      Hn[0] = NEG;
      In[0] = NEG;
      for (int64_t j = 1; j <= N; ++j) {
        const int m = std::max(H[s][j - 1], G) + (lr_same(rd[f], t[j - 1]) ? kMatch : -kMismatch);
        const int ins = f > 0 ? std::max(H[s][j] - kGapOpen - kGapExt, I[s][j] - kGapExt) : NEG;
        d = std::max(Hn[j - 1] - kGapOpen - kGapExt, d - kGapExt);
        const int h = std::max(m, std::max(ins, d));
        Hn[j] = h;
        In[j] = ins;
        F = std::max(F, (long)m);
      }
      std::swap(H[s], Hn);
      std::swap(I[s], In);
      Hn.assign(N + 1, NEG);
      In.assign(N + 1, NEG);
    }
    M = std::max(M, F);
    if (f == L - 1) end_at_L = F;
  }
  return std::max(M - kClip, end_at_L);
}

// The score of the planted pieces' own alignments: each piece (a read
// interval and its source window) aligned locally there, summed as a jump
// alignment. A real jump alignment, so a valid lb.
long known_jump(const std::vector<uint8_t>& ref, const std::vector<uint8_t>& rd,
                const std::vector<std::pair<std::pair<int, int>, std::pair<int64_t, int>>>& pieces, long P) {
  // pieces: ((read a, b), (window start, strand))
  long total = 0;
  int k = 0;
  const int L = static_cast<int>(rd.size());
  for (size_t pi = 0; pi < pieces.size(); ++pi) {
    const auto& pc = pieces[pi];
    // The read's ends: the first piece pays 5 if its alignment starts after
    // read base 0, the last if it ends before the read's end.
    const bool is_first = pi == 0, is_last = pi + 1 == pieces.size();
    const int a = pc.first.first, b = pc.first.second;
    const int64_t w0 = std::max<int64_t>(0, pc.second.first - 50);
    const int64_t w1 = std::min<int64_t>(static_cast<int64_t>(ref.size()), pc.second.first + (b - a) + 300);
    std::vector<uint8_t> win(ref.begin() + w0, ref.begin() + w1);
    std::vector<uint8_t> seg(rd.begin() + a, rd.begin() + b);
    if (pc.second.second) win = rc(win);
    // best local score of seg in win (inner piece ends cost nothing)
    const int n = static_cast<int>(seg.size()), W = static_cast<int>(win.size()), NEG = INT_MIN / 4;
    std::vector<int> Hp(W + 1, NEG), Ip(W + 1, NEG), Hc(W + 1), Ic(W + 1);
    long best = LONG_MIN / 4;
    for (int i = 0; i < n; ++i) {
      int d = NEG;
      Hc[0] = Ic[0] = NEG;
      for (int j = 1; j <= W; ++j) {
        const int start = is_first && a + i > 0 ? -kClip : 0;
        const int m = std::max(Hp[j - 1], start) + (lr_same(seg[i], win[j - 1]) ? kMatch : -kMismatch);
        const int ins = i > 0 ? std::max(Hp[j] - kGapOpen - kGapExt, Ip[j] - kGapExt) : NEG;
        d = std::max(Hc[j - 1] - kGapOpen - kGapExt, d - kGapExt);
        Hc[j] = std::max(m, std::max(ins, d));
        Ic[j] = ins;
        best = std::max(best, (long)m - (is_last && a + i + 1 < L ? kClip : 0));
      }
      std::swap(Hp, Hc);
      std::swap(Ip, Ic);
    }
    total += best;
    ++k;
  }
  return total - P * (k - 1);
}

}  // namespace

int main(int argc, char** argv) {
  const uint64_t seed = argc > 1 ? std::strtoull(argv[1], nullptr, 10) : 1;
  std::mt19937_64 g(seed);
  std::fprintf(stderr, "seed %llu\n", static_cast<unsigned long long>(seed));
  // A genome of two contigs, with a segment duplicated (two loci for some reads).
  std::vector<uint8_t> c1 = random_seq(g, 40000), c2 = random_seq(g, 25000);
  {
    const std::vector<uint8_t> seg(c1.begin() + 5000, c1.begin() + 7000);
    std::copy(seg.begin(), seg.end(), c2.begin() + 15000);
    c2[15500] = static_cast<uint8_t>((c2[15500] + 1) % 4);  // one difference between the copies
  }
  Reference ref;
  ref.seq.owned.assign(kContigPad, 4);
  const std::vector<uint8_t>* cs[2] = {&c1, &c2};
  for (int c = 0; c < 2; ++c) {
    ref.names.push_back("c" + std::to_string(c + 1));
    ref.offsets.push_back(ref.seq.owned.size());
    ref.lengths.push_back(cs[c]->size());
    ref.seq.owned.insert(ref.seq.owned.end(), cs[c]->begin(), cs[c]->end());
    ref.seq.owned.insert(ref.seq.owned.end(), kContigPad, 4);
  }
  const std::vector<uint8_t>& R = ref.seq.owned;
  const int64_t o1 = static_cast<int64_t>(ref.offsets[0]), o2 = static_cast<int64_t>(ref.offsets[1]);
  Index ix_ = Index::build(ref, 16, 1, 4);
  const IndexView ix = ix_.view(ref);

  // Reads: (sequence, planted pieces). Kinds: plain, deletion, novel
  // insertion, inversion, translocation, duplicated segment.
  struct Read {
    std::vector<uint8_t> seq;
    std::vector<std::pair<std::pair<int, int>, std::pair<int64_t, int>>> pieces;
    const char* kind;
  };
  std::vector<Read> reads;
  auto from = [&](int64_t a, int64_t n) { return std::vector<uint8_t>(R.begin() + a, R.begin() + a + n); };
  const double rate = 0.004;
  for (int n = 0; n < 6; ++n) {
    const int64_t a = o1 + 1000 + static_cast<int64_t>(g() % 30000);
    const int len = 700 + static_cast<int>(g() % 400);
    Read r;
    r.kind = "plain";
    r.seq = noisy(g, from(a, len), rate);
    r.pieces = {{{0, (int)r.seq.size()}, {a, 0}}};
    reads.push_back(r);
  }
  for (int n = 0; n < 6; ++n) {  // deletion of D reference bases
    const int64_t a = o1 + 1000 + static_cast<int64_t>(g() % 30000);
    const int D = n % 2 ? 60 : 400, h = 400;
    Read r;
    r.kind = "deletion";
    std::vector<uint8_t> p1 = noisy(g, from(a, h), rate), p2 = noisy(g, from(a + h + D, h), rate);
    r.seq = p1;
    r.seq.insert(r.seq.end(), p2.begin(), p2.end());
    r.pieces = {{{0, (int)p1.size()}, {a, 0}}, {{(int)p1.size(), (int)r.seq.size()}, {a + h + D, 0}}};
    reads.push_back(r);
  }
  for (int n = 0; n < 4; ++n) {  // novel insertion
    const int64_t a = o1 + 1000 + static_cast<int64_t>(g() % 30000);
    const int I = n % 2 ? 40 : 250, h = 400;
    Read r;
    r.kind = "insertion";
    std::vector<uint8_t> p1 = noisy(g, from(a, h), rate), ins = random_seq(g, I), p2 = noisy(g, from(a + h, h), rate);
    r.seq = p1;
    r.seq.insert(r.seq.end(), ins.begin(), ins.end());
    r.seq.insert(r.seq.end(), p2.begin(), p2.end());
    const int b2 = (int)(p1.size() + ins.size());
    r.pieces = {{{0, (int)p1.size()}, {a, 0}}, {{b2, (int)r.seq.size()}, {a + h, 0}}};
    reads.push_back(r);
  }
  for (int n = 0; n < 5; ++n) {  // inversion of V bases
    const int64_t a = o1 + 1000 + static_cast<int64_t>(g() % 30000);
    const int V = 300, h = 350;
    Read r;
    r.kind = "inversion";
    std::vector<uint8_t> p1 = noisy(g, from(a, h), rate), p2 = noisy(g, rc(from(a + h, V)), rate),
                         p3 = noisy(g, from(a + h + V, h), rate);
    r.seq = p1;
    r.seq.insert(r.seq.end(), p2.begin(), p2.end());
    r.seq.insert(r.seq.end(), p3.begin(), p3.end());
    const int b1 = (int)p1.size(), b2 = (int)(p1.size() + p2.size());
    r.pieces = {{{0, b1}, {a, 0}}, {{b1, b2}, {a + h, 1}}, {{b2, (int)r.seq.size()}, {a + h + V, 0}}};
    reads.push_back(r);
  }
  for (int n = 0; n < 5; ++n) {  // translocation between the contigs
    const int64_t a = o1 + 1000 + static_cast<int64_t>(g() % 30000), b = o2 + 1000 + static_cast<int64_t>(g() % 12000);
    Read r;
    r.kind = "translocation";
    std::vector<uint8_t> p1 = noisy(g, from(a, 450), rate), p2 = noisy(g, from(b, 450), rate);
    r.seq = p1;
    r.seq.insert(r.seq.end(), p2.begin(), p2.end());
    r.pieces = {{{0, (int)p1.size()}, {a, 0}}, {{(int)p1.size(), (int)r.seq.size()}, {b, 0}}};
    reads.push_back(r);
  }
  for (int n = 0; n < 3; ++n) {  // ragged ends: 6 random bases on each, so the optimum clips both
    const int64_t a = o1 + 1000 + static_cast<int64_t>(g() % 30000);
    Read r;
    r.kind = "ragged";
    const std::vector<uint8_t> head = random_seq(g, 6), body = from(a, 700), tail = random_seq(g, 6);
    r.seq = head;
    r.seq.insert(r.seq.end(), body.begin(), body.end());
    r.seq.insert(r.seq.end(), tail.begin(), tail.end());
    r.pieces = {{{0, (int)r.seq.size()}, {a - 6, 0}}};
    reads.push_back(r);
  }
  for (int n = 0; n < 4; ++n) {  // a 20-bp deletion, noise-free: cheaper as a gap than as a jump (P = 30)
    const int64_t a = o1 + 1000 + static_cast<int64_t>(g() % 30000);
    const int h = 400;
    Read r;
    r.kind = "small-del";
    r.seq = from(a, h);
    const std::vector<uint8_t> p2 = from(a + h + 20, h);
    r.seq.insert(r.seq.end(), p2.begin(), p2.end());
    r.pieces = {{{0, (int)r.seq.size()}, {a, 0}}};
    reads.push_back(r);
  }
  for (int n = 0; n < 4; ++n) {  // insertion of a diverged copy: good, but no exact part
    const int64_t a = o1 + 1000 + static_cast<int64_t>(g() % 30000), b = o2 + 1000 + static_cast<int64_t>(g() % 12000);
    const int h = 400, I = 260;
    Read r;
    r.kind = "diverged-ins";
    std::vector<uint8_t> p1 = noisy(g, from(a, h), rate), ins = from(b, I), p2 = noisy(g, from(a + h, h), rate);
    for (int k = 7; k < I; k += 13) ins[k] = static_cast<uint8_t>((ins[k] + 1 + g() % 3) % 4);  // no exact 16-mer
    r.seq = p1;
    r.seq.insert(r.seq.end(), ins.begin(), ins.end());
    r.seq.insert(r.seq.end(), p2.begin(), p2.end());
    const int b1 = (int)p1.size(), b2 = b1 + I;
    r.pieces = {{{0, b1}, {a, 0}}, {{b2, (int)r.seq.size()}, {a + h, 0}}};  // the copy is not in the known alignment
    reads.push_back(r);
  }
  for (int n = 0; n < 4; ++n) {  // in the duplicated segment (two loci)
    const int64_t a = o1 + 5000 + static_cast<int64_t>(g() % 1200);
    Read r;
    r.kind = "duplicate";
    r.seq = noisy(g, from(a, 700), rate);
    r.pieces = {{{0, (int)r.seq.size()}, {a, 0}}};
    reads.push_back(r);
  }
  // Half the reads reverse-complemented (their pieces then too).
  for (auto& r : reads)
    if (g() % 2) {
      const int L = static_cast<int>(r.seq.size());
      r.seq = rc(r.seq);
      for (auto& pc : r.pieces) {
        const int a = pc.first.first, b = pc.first.second;
        pc.first = {L - b, L - a};
        pc.second.second ^= 1;
      }
      std::reverse(r.pieces.begin(), r.pieces.end());
    }

  PartitionParams pp;
  pp.keep_occ = true;
  LrScratch sc;
  BandDp dp;
  const long P = 30;
  int single_checked = 0, single_proven = 0, jump_proven = 0, jump_interval = 0, jump_none = 0, weak_checked = 0,
      weak_proven = 0, tight_checked = 0;
  for (size_t i = 0; i < reads.size(); ++i) {
    const Read& rr = reads[i];
    const int L = static_cast<int>(rr.seq.size());
    const std::vector<uint8_t> rcs = rc(rr.seq);
    StrandParts so[2];
    partition_strand(ix, rr.seq.data(), L, pp, sc, so[0]);
    partition_strand(ix, rcs.data(), L, pp, sc, so[1]);
    const long U = std::max(lemma_l_span(L, so[0].parts), lemma_l_span(L, so[1].parts));
    // The local optimum: lb from the first planted piece alone (a real
    // alignment, clipped where the read goes on).
    const long opt1 = brute_jump(R, rr.seq, 100000000);
    {
      const auto& pc = rr.pieces[0];
      std::vector<std::pair<std::pair<int, int>, std::pair<int64_t, int>>> one = {pc};
      const long lb1 = known_jump(R, rr.seq, one, P);
      const uint8_t* rds[2] = {rr.seq.data(), rcs.data()};
      const LrCert c1 = certify_single(ix, rds, L, so, lb1, U, 4e9, dp);
      ++single_checked;
      CHECK(lb1 <= opt1, "read %zu (%s): the known piece scores %ld above the local optimum %ld", i, rr.kind, lb1, opt1);
      if (c1.status == 1) {
        ++single_proven;
        CHECK(c1.opt == opt1, "read %zu (%s): proven local optimum %ld, brute force %ld", i, rr.kind, c1.opt, opt1);
      }
      // From lb = OPT - 1 (the narrowest band that must still hold the
      // optimum): it must be found.
      if (opt1 - 1 > U) {
        const LrCert ct = certify_single(ix, rds, L, so, opt1 - 1, U, 4e9, dp);
        ++tight_checked;
        CHECK(ct.status == 1 && ct.opt == opt1, "read %zu (%s): from lb = OPT - 1, status %d, OPT %ld, brute force %ld", i,
              rr.kind, ct.status, ct.opt, opt1);
      }
    }
    // The jump optimum: lb from all planted pieces.
    const long optJ = brute_jump(R, rr.seq, P);
    const long lbJ = known_jump(R, rr.seq, rr.pieces, P);
    CHECK(lbJ <= optJ, "read %zu (%s): the planted jump alignment scores %ld above OPT_J %ld", i, rr.kind, lbJ, optJ);
    CHECK(opt1 <= optJ, "read %zu (%s): the local optimum %ld above OPT_J %ld", i, rr.kind, opt1, optJ);
    std::vector<ReadPart> fp[2];
    forward_parts(L, so, fp);
    const long UJ = lemma_l_jump(L, fp, P);
    const LrCert cj = certify_jump(ix, rr.seq.data(), L, so, lbJ, P, 4e9);
    // The same from the weakest admissible lb, U_J + 1 (<= lb_J <= OPT_J):
    // the search and the pruning must then find the optimum themselves.
    if (lbJ > UJ) {
      const LrCert cw = certify_jump(ix, rr.seq.data(), L, so, UJ + 1, P, 4e9);
      ++weak_checked;
      CHECK(cw.status == 1 || cw.status == 3, "read %zu (%s): weak lb not certified (status %d)", i, rr.kind, cw.status);
      if (cw.status == 1)
        CHECK(cw.opt == optJ, "read %zu (%s): from the weak lb, proven OPT_J %ld, brute force %ld", i, rr.kind, cw.opt,
              optJ);
      if (cw.status == 3)
        CHECK(cw.opt <= optJ && optJ <= cw.hi, "read %zu (%s): from the weak lb, OPT_J %ld outside [%ld, %ld]", i,
              rr.kind, optJ, cw.opt, cw.hi);
      weak_proven += cw.status == 1;
    }
    if (optJ - 1 > UJ) {
      const LrCert ct = certify_jump(ix, rr.seq.data(), L, so, optJ - 1, P, 4e9);
      ++tight_checked;
      CHECK((ct.status == 1 && ct.opt == optJ) || (ct.status == 3 && ct.opt <= optJ && optJ <= ct.hi),
            "read %zu (%s): from lb = OPT_J - 1, status %d [%ld, %ld], brute force %ld", i, rr.kind, ct.status, ct.opt,
            ct.hi, optJ);
    }
    if (cj.status == 1) {
      ++jump_proven;
      CHECK(cj.opt == optJ, "read %zu (%s): proven OPT_J %ld, brute force %ld (lb %ld, U_J %ld)", i, rr.kind, cj.opt,
            optJ, lbJ, UJ);
    } else if (cj.status == 3) {
      ++jump_interval;
      CHECK(cj.opt <= optJ && optJ <= cj.hi, "read %zu (%s): OPT_J %ld outside the interval [%ld, %ld]", i, rr.kind,
            optJ, cj.opt, cj.hi);
    } else {
      ++jump_none;
      CHECK(lbJ <= UJ || cj.status == 2, "read %zu (%s): fits (lb %ld > U_J %ld) but not attempted", i, rr.kind, lbJ,
            UJ);
    }
    std::fprintf(stderr, "  read %2zu %-13s L %4d  OPT %5ld  OPT_J %5ld  lb_J %5ld  U_J %5ld  -> %s [%ld, %ld], %d of %d "
                 "chains searched\n", i, rr.kind, L, opt1, optJ, lbJ, UJ,
                 cj.status == 1 ? "proven" : cj.status == 3 ? "interval" : cj.status == 2 ? "capped" : "no fit",
                 cj.opt, cj.hi, cj.searched, cj.chains);
  }
  std::fprintf(stderr, "local optimum: %d proven of %d (all equal to the brute force); split optimum: %d proven, %d "
               "intervals (all containing it), %d not attempted; from the weak lb U_J + 1: %d proven of %d; from lb "
               "= OPT - 1: %d checked\n",
               single_proven, single_checked, jump_proven, jump_interval, jump_none, weak_proven, weak_checked,
               tight_checked);
  if (failures) {
    std::fprintf(stderr, "%d check(s) FAILED\n", failures);
    return 1;
  }
  std::fprintf(stderr, "all long-read checks passed\n");
  return 0;
}
