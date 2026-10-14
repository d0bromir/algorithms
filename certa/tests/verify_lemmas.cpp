// Exhaustive checks of the finite lemmas behind CERTA's bounds (Aim 5,
// THEORY.md). Each lemma is a statement about alignment transcripts or short
// strings. Here it is checked on every instance up to a size, and each
// hypothesis is shown necessary by a counterexample just outside it.
//
//  Lemma L  parts >= 3 long: every part inside the span that does not occur
//           exactly costs >= 5, disjointly (fails for length 2).
//  Lemma E  a local alignment whose end-to-end completion has d >= 1 edits
//           scores <= L - d - (d >= 2 ? 5 : 4), and this is attained.
//  Lemma H  parts >= 4 long, one-mismatch probes: every unwitnessed part
//           costs >= 7, disjointly (fails for length 3).
//  Lemma S  parts >= 5 long: Lemma L's maximum over spans is the full span
//           (fails for length 4); the span candidates of span_bound match a
//           brute-force maximum over all spans.
//  Prop. A  the constructions give unwitnessed parts of cost e + 7.
//  Thm G    the greedy partition is optimal (random instances, against DP).
//  XQ       BWA-MEM's MAPQ formula (mapq_se) is non-increasing in SUB and
//           non-decreasing in the score, so XQ bounds it from below.
//  Lemma L-J  the one-pass bound of jump (split) alignments, lemma_l_jump,
//           equals the best decomposition into pieces found by brute force
//           (random two-strand part layouts), and the single-span bound
//           when jumps cost too much to use.
//  Thm L-J  the pruning values of certify_jump: the prefix values (a pass
//           over the read), the suffix values (the same pass over the
//           reversed read, on the mirrored parts), pot and the read window
//           of a chain, against their definitions by brute force (random
//           layouts).
//
// Transcripts: each read base is a match, a mismatch or inserted; between
// consecutive bases lies a deletion of 0, 1 or 2 reference bases. Loss
// (span minus score): mismatch 5, insertion run of k bases 6 + 2k, deletion
// of k bases 6 + k, clipped end of c bases c + 5. A part occurs exactly when
// all its bases are matches and no deletion lies between two of them.
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

#include "certa/interval.h"
#include "certa/longread.h"
#include "certa/pair.h"

namespace {

int failures = 0;
#define EXPECT(cond, ...)                                        \
  do {                                                           \
    if (!(cond)) {                                               \
      if (++failures <= 20) {                                    \
        std::fprintf(stderr, "FAIL %s:%d: ", __FILE__, __LINE__); \
        std::fprintf(stderr, __VA_ARGS__);                       \
        std::fprintf(stderr, "\n");                              \
      }                                                          \
    }                                                            \
  } while (0)

enum : int { kMatchB = 0, kMis = 1, kIns = 2 };
struct Part { int off, end; };

// All layouts of disjoint parts, each at least `minlen` long, in [0, L).
void layouts(int L, int minlen, int from, std::vector<Part>& cur, std::vector<std::vector<Part>>& out) {
  if (!cur.empty()) out.push_back(cur);
  for (int o = from; o + minlen <= L; ++o)
    for (int e = o + minlen; e <= L; ++e) {
      cur.push_back({o, e});
      layouts(L, minlen, e, cur, out);
      cur.pop_back();
    }
}

// Visits every transcript of n bases: states (match / mismatch / insert) and
// the deletions between consecutive bases (0..maxdel), restricted to the gap
// positions allowed by `delok`.
template <class F>
void transcripts(int n, int maxdel, const std::vector<char>& delok, F f) {
  std::vector<int> st(n, 0), del(std::max(0, n - 1), 0);
  for (;;) {
    f(st, del);
    int i = 0;  // odometer over the deletions, then the states
    for (; i < n - 1; ++i) {
      if (!delok[i]) continue;
      if (++del[i] <= maxdel) break;
      del[i] = 0;
    }
    if (i < n - 1) continue;
    int j = 0;
    for (; j < n; ++j) {
      if (++st[j] <= kIns) break;
      st[j] = 0;
    }
    if (j == n) return;
  }
}

// Loss of a transcript (no clipping).
int loss_of(const std::vector<int>& st, const std::vector<int>& del) {
  const int n = static_cast<int>(st.size());
  int loss = 0;
  for (int i = 0; i < n; ++i) {
    if (st[i] == kMis) loss += 5;
    if (st[i] == kIns && (i == 0 || st[i - 1] != kIns || del[i - 1])) {  // a new insertion run
      int k = 0;
      for (int j = i; j < n && st[j] == kIns && (j == i || !del[j - 1]); ++j) ++k;
      loss += 6 + 2 * k;
    }
  }
  for (size_t i = 0; i < del.size(); ++i)
    if (del[i]) loss += 6 + del[i];
  return loss;
}

bool exact(const Part& p, const std::vector<int>& st, const std::vector<int>& del) {
  for (int i = p.off; i < p.end; ++i)
    if (st[i] != kMatchB || (i + 1 < p.end && del[i])) return false;
  return true;
}

// Witnessed by a one-mismatch probe (Lemma H), read off the transcript: no
// deletion inside, and the part's bases are all matches except one mismatch,
// or except one inserted base at an end of the part (its other bases lie on
// one diagonal, so the ungapped placement there has <= 1 mismatch). Parts
// witnessed only through sequence coincidences count as unwitnessed, which
// only makes the check stricter.
bool witnessed1(const Part& p, const std::vector<int>& st, const std::vector<int>& del) {
  int mis = 0, ins = 0, insat = -1;
  for (int i = p.off; i < p.end; ++i) {
    if (i + 1 < p.end && del[i]) return false;
    if (st[i] == kMis) ++mis;
    if (st[i] == kIns) { ++ins; insat = i; }
  }
  if (mis + ins == 0) return true;
  if (ins == 0) return mis == 1;
  return ins == 1 && mis == 0 && (insat == p.off || insat == p.end - 1);
}

// Lemma L / H on all layouts of a read of length L with parts >= minlen:
// loss >= per * (#parts failing `ok`). Returns the number of violations
// (0 expected) and the first one in `example`.
template <class Ok>
long long part_lemma(int L, int minlen, int per, Ok ok, std::string* example, long long* checked,
                     bool densest = false) {
  std::vector<std::vector<Part>> lay;
  std::vector<Part> cur;
  if (densest) {  // only the layout of consecutive parts of length minlen
    for (int o = 0; o + minlen <= L; o += minlen) cur.push_back({o, o + minlen});
    lay.push_back(cur);
  } else {
    layouts(L, minlen, 0, cur, lay);
  }
  long long bad = 0;
  // Deletions only inside parts matter: elsewhere they add loss and touch
  // no part, so dropping them only makes the check stricter.
  for (const auto& parts : lay) {
    std::vector<char> delok(L - 1, 0);
    for (const Part& p : parts)
      for (int i = p.off; i + 1 < p.end; ++i) delok[i] = 1;
    transcripts(L, 1, delok, [&](const std::vector<int>& st, const std::vector<int>& del) {
      ++*checked;
      int n = 0;
      for (const Part& p : parts) n += !ok(p, st, del);
      const int loss = loss_of(st, del);
      if (loss < per * n) {
        if (++bad == 1) {
          char buf[256];
          int k = std::snprintf(buf, sizeof buf, "L=%d parts", L);
          for (const Part& p : parts) k += std::snprintf(buf + k, sizeof buf - k, " [%d,%d)", p.off, p.end);
          k += std::snprintf(buf + k, sizeof buf - k, " states ");
          for (int i = 0; i < L; ++i) k += std::snprintf(buf + k, sizeof buf - k, "%c%s", "=XI"[st[i]], i + 1 < L && del[i] ? "D" : "");
          std::snprintf(buf + k, sizeof buf - k, ": loss %d < %d x %d", loss, per, n);
          *example = buf;
        }
      }
    });
  }
  return bad;
}

// Lemma L's bound for one span: (y - x) - 5 max(0, #counted parts inside - extra) - clips.
int span_value(int L, const std::vector<Part>& p, uint32_t counted, int extra, int x, int y) {
  int k = 0;
  for (size_t i = 0; i < p.size(); ++i) k += ((counted >> i) & 1) && p[i].off >= x && p[i].end <= y;
  k = std::max(0, k - extra);
  return (y - x) - 5 * k - (x > 0 ? 5 : 0) - (y < L ? 5 : 0);
}

}  // namespace

int main() {
  // ---- Lemma L
  {
    std::string ex;
    long long n = 0;
    const long long bad = part_lemma(9, 3, 5, exact, &ex, &n);
    EXPECT(bad == 0, "Lemma L fails for parts >= 3: %s", ex.c_str());
    std::fprintf(stderr, "Lemma L, parts >= 3, L = 9: %lld transcript-layout pairs, %lld violations\n", n, bad);
    long long n2 = 0;
    std::string ex2;
    const long long bad2 = part_lemma(8, 2, 5, exact, &ex2, &n2);
    EXPECT(bad2 > 0, "Lemma L should fail for parts of length 2");
    std::fprintf(stderr, "Lemma L, parts >= 2, L = 8: %lld violations, e.g. %s\n", bad2, ex2.c_str());
  }
  // ---- Lemma H (one-mismatch probes, bound 7)
  {
    std::string ex;
    long long n = 0;
    long long bad = part_lemma(10, 4, 7, witnessed1, &ex, &n);
    bad += part_lemma(12, 4, 7, witnessed1, &ex, &n, true);  // three consecutive parts of 4
    EXPECT(bad == 0, "Lemma H fails for parts >= 4: %s", ex.c_str());
    std::fprintf(stderr, "Lemma H, parts >= 4 (all layouts, L = 10; three parts of 4, L = 12): %lld transcript-layout "
                 "pairs, %lld violations\n", n, bad);
    long long n2 = 0;
    std::string ex2;
    const long long bad2 = part_lemma(9, 3, 7, witnessed1, &ex2, &n2, true);
    EXPECT(bad2 > 0, "Lemma H should fail for parts of length 3");
    std::fprintf(stderr, "Lemma H, parts of 3, L = 9: %lld violations, e.g. %s\n", bad2, ex2.c_str());
  }
  // ---- Lemma E: clipped ends, interior transcript; the completion's edits
  // are at most mismatches + inserted + deleted + clipped bases (worst case:
  // every clipped base mismatches), and edit_bound decreases in d.
  {
    const int L = 8;
    long long n = 0, bad = 0;
    int tight[L + 2] = {};  // tight[d]: some transcript attains the bound
    for (int c1 = 0; c1 < L; ++c1)
      for (int c2 = 0; c1 + c2 < L; ++c2) {
        const int m = L - c1 - c2;
        std::vector<char> delok(std::max(0, m - 1), 1);
        transcripts(m, 2, delok, [&](const std::vector<int>& st, const std::vector<int>& del) {
          if (st[0] == kIns || st[m - 1] == kIns) return;  // an alignment starts and ends with an aligned base
          ++n;
          int matched = 0, mis = 0, ins = 0, dels = 0, gap = 0;
          for (int i = 0; i < m; ++i) {
            matched += st[i] == kMatchB;
            mis += st[i] == kMis;
            ins += st[i] == kIns;
            if (st[i] == kIns && (st[i - 1] != kIns || del[i - 1])) {
              int k = 0;
              for (int j = i; j < m && st[j] == kIns && (j == i || !del[j - 1]); ++j) ++k;
              gap += 6 + k;
            }
          }
          for (int x : del)
            if (x) { dels += x; gap += 6 + x; }
          const int score = matched - 4 * mis - gap - (c1 ? 5 : 0) - (c2 ? 5 : 0);
          const int d = mis + ins + dels + c1 + c2;
          if (d < 1) return;
          const int b = certa::edit_bound(L, d);
          if (score > b) ++bad;
          if (score == b && d <= L) tight[d] = 1;
        });
      }
    EXPECT(bad == 0, "Lemma E: %lld transcripts score above edit_bound", bad);
    int ntight = 0;
    for (int d = 1; d <= 4; ++d) ntight += tight[d];
    EXPECT(ntight == 4, "Lemma E should be attained for d = 1..4");
    std::fprintf(stderr, "Lemma E, L = %d: %lld clipped transcripts, %lld violations; attained for d = 1..4: %s\n", L, n,
                 bad, ntight == 4 ? "yes" : "no");
  }
  // ---- Lemma S: the full span maximizes Lemma L's bound for parts >= 5,
  // and span_bound (pair.h) equals the brute-force maximum over all spans.
  {
    long long n = 0, notfull = 0, mism = 0;
    for (int minlen = 4; minlen <= 5; ++minlen) {
      const int L = 24;
      std::vector<std::vector<Part>> lay;
      std::vector<Part> cur;
      layouts(L, minlen, 0, cur, lay);
      long long nf = 0;
      for (const auto& p : lay) {
        if (p.size() > 6) continue;
        for (uint32_t counted = 0; counted < (1u << p.size()); ++counted)
          for (int extra = 0; extra <= 2; ++extra) {
            int best = -1000;
            for (int x = 0; x < L; ++x)
              for (int y = x + 1; y <= L; ++y) best = std::max(best, span_value(L, p, counted, extra, x, y));
            nf += best > span_value(L, p, counted, extra, 0, L);
            if (minlen == 5) ++n;
            if (extra == 0) {
              int off[certa::kWinParts], end[certa::kWinParts];
              for (size_t i = 0; i < p.size(); ++i) { off[i] = p[i].off; end[i] = p[i].end; }
              mism += certa::span_bound(L, static_cast<int>(p.size()), off, end, counted) != best;
            }
          }
      }
      if (minlen == 5) notfull = nf;
      else EXPECT(nf > 0, "Lemma S should fail for parts of length 4");
      if (minlen == 4) std::fprintf(stderr, "Lemma S, parts >= 4: the full span is not the maximum in %lld cases\n", nf);
    }
    EXPECT(notfull == 0, "Lemma S fails for parts >= 5 in %lld cases", notfull);
    EXPECT(mism == 0, "span_bound differs from the brute-force maximum in %lld cases", mism);
    std::fprintf(stderr, "Lemma S, parts >= 5, L = 24: %lld (layout, counted set, extra) cases, full span always the "
                 "maximum: %s; span_bound = brute force on all layouts: %s\n", n, notfull ? "no" : "yes", mism ? "no" : "yes");
  }
  // ---- Proposition A: p = (AC)^m, the reference holds its left half, G^(e+1),
  // its right half, T-flanks; the alignment deletes the G's (cost e + 7). No
  // reference substring is within e edits of p (brute-force edit distance).
  {
    auto edit = [](const std::string& a, const std::string& b) {
      std::vector<int> prev(b.size() + 1), cur(b.size() + 1);
      for (size_t j = 0; j <= b.size(); ++j) prev[j] = static_cast<int>(j);
      for (size_t i = 1; i <= a.size(); ++i) {
        cur[0] = static_cast<int>(i);
        for (size_t j = 1; j <= b.size(); ++j)
          cur[j] = std::min({prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (a[i - 1] != b[j - 1])});
        std::swap(prev, cur);
      }
      return prev[b.size()];
    };
    for (int e = 0; e <= 3; ++e) {
      const int m = 2 * (e + 1);
      std::string p;
      for (int i = 0; i < m; ++i) p += "AC";
      const std::string ref = std::string(12, 'T') + p.substr(0, m) + std::string(e + 1, 'G') + p.substr(m) + std::string(12, 'T');
      int closest = 1 << 20;
      for (size_t a = 0; a < ref.size(); ++a)
        for (size_t b = a + 1; b <= ref.size(); ++b) closest = std::min(closest, edit(p, ref.substr(a, b - a)));
      // Ungapped placements (Hamming), for e = 0.
      int ham = 1 << 20;
      for (size_t a = 0; a + p.size() <= ref.size(); ++a) {
        int h = 0;
        for (size_t i = 0; i < p.size(); ++i) h += p[i] != ref[a + i];
        ham = std::min(ham, h);
      }
      EXPECT(closest > e, "Prop. A, e = %d: a substring within %d edits", e, closest);
      if (e == 0) EXPECT(ham >= 2, "Prop. A, Hamming: a placement with %d mismatches", ham);
      std::fprintf(stderr, "Prop. A, e = %d: part (AC)^%d costs %d (= e + 7); nearest substring %d edits; nearest ungapped "
                   "placement %d mismatches\n", e, m, 6 + e + 1, closest, ham);
    }
  }
  // ---- Theorem G: greedy vs DP on random instances.
  {
    std::mt19937_64 g(7);
    long long n = 0, bad = 0;
    for (int it = 0; it < 20000; ++it) {
      const int TL = 40 + static_cast<int>(g() % 260), L = 8 + static_cast<int>(g() % 24), tau = static_cast<int>(g() % 9);
      const int sigma = 2 + static_cast<int>(g() % 2);
      std::string T(TL, 'a'), r(L, 'a');
      for (char& c : T) c = static_cast<char>('a' + g() % sigma);
      if (g() % 2) {  // a read from the text, with an edit or two
        const int at = static_cast<int>(g() % std::max(1, TL - L));
        r = T.substr(at, std::min(L, TL - at));
        r.resize(L, 'a');
        for (int k = static_cast<int>(g() % 3); k > 0; --k) r[g() % L] = static_cast<char>('a' + g() % sigma);
      } else {
        for (char& c : r) c = static_cast<char>('a' + g() % sigma);
      }
      auto occ = [&](int i, int j) {
        int c = 0;
        const std::string s = r.substr(i, j - i);
        for (size_t a = 0; a + s.size() <= T.size(); ++a) c += T.compare(a, s.size(), s) == 0;
        return c;
      };
      std::vector<std::vector<char>> ok(L + 1, std::vector<char>(L + 1, 0));
      for (int i = 0; i < L; ++i)
        for (int j = i + 1; j <= L; ++j) ok[i][j] = occ(i, j) <= tau;
      int greedy = 0;
      for (int s = 0; s < L;) {
        int e = s + 1;
        while (e <= L && !ok[s][e]) ++e;
        if (e > L) break;
        ++greedy;
        s = e;
      }
      std::vector<int> best(L + 1, 0);
      for (int i = L - 1; i >= 0; --i) {
        best[i] = best[i + 1];
        for (int j = i + 1; j <= L; ++j)
          if (ok[i][j]) best[i] = std::max(best[i], 1 + best[j]);
      }
      ++n;
      bad += greedy != best[0];
    }
    EXPECT(bad == 0, "Theorem G: greedy below the optimum in %lld instances", bad);
    std::fprintf(stderr, "Theorem G: %lld random instances, greedy = DP optimum in all: %s\n", n, bad ? "no" : "yes");
  }
  // ---- Lemma L-J: lemma_l_jump against a brute force over all
  // decompositions into pieces (start, end, strand), each bounded by Lemma L.
  {
    std::mt19937_64 g(7);
    long long n = 0, bad = 0, bad1 = 0, n2 = 0, bad2 = 0;
    for (int it = 0; it < 30000; ++it) {
      const int L = 3 + static_cast<int>(g() % 22);
      std::vector<certa::ReadPart> part[2];
      for (int s = 0; s < 2; ++s) {  // disjoint parts >= 3 long, gaps allowed
        int x = static_cast<int>(g() % 3);
        for (;;) {
          const int len = 3 + static_cast<int>(g() % 4);
          if (x + len > L) break;
          if (g() % 4) part[s].push_back({x, x + len, 0});
          x += len + static_cast<int>(g() % 3);
        }
      }
      auto V = [&](int s, int a, int b) {  // Lemma L for one piece [a, b) on strand s
        long v = b - a;
        for (const auto& q : part[s]) v -= 5 * (q.off >= a && q.end <= b);
        return v;
      };
      for (long P : {0L, 3L, 5L, 9L, 10L, 30L, 1000000L}) {
        // f[b]: the best decomposition whose last piece ends at b; M[a]: max f[e], e <= a.
        std::vector<long> f(L + 1, LONG_MIN / 4), M(L + 1, LONG_MIN / 4);
        for (int b = 1; b <= L; ++b) {
          for (int a = 0; a < b; ++a) {
            const long pre = std::max(a > 0 ? -5L : 0L, M[a] - P);
            for (int s = 0; s < 2; ++s) f[b] = std::max(f[b], pre + V(s, a, b));
          }
          M[b] = std::max(M[b - 1], f[b]);
        }
        long brute = LONG_MIN;
        for (int b = 1; b <= L; ++b) brute = std::max(brute, f[b] - (b < L ? 5 : 0));
        ++n;
        bad += certa::lemma_l_jump(L, part, P) != brute;
      }
      bad1 += certa::lemma_l_jump(L, part, 1000000) !=
              std::max(certa::lemma_l_span(L, part[0]), certa::lemma_l_span(L, part[1]));
      // jump_bound in general: piece types with free parts, and a type that
      // must be used (the pruning bound of certify_jump).
      std::vector<char> fr[2];
      for (int s = 0; s < 2; ++s)
        for (size_t k = 0; k < part[s].size(); ++k) fr[s].push_back(static_cast<char>(g() % 2));
      const std::vector<certa::PieceType> types = {{0, nullptr}, {1, nullptr}, {0, &fr[0]}, {1, &fr[1]}};
      auto VT = [&](int t, int a, int b) {
        const int s = types[t].strand;
        long v = b - a;
        for (size_t k = 0; k < part[s].size(); ++k)
          v -= 5 * (part[s][k].off >= a && part[s][k].end <= b && !(types[t].free && (*types[t].free)[k]));
        return v;
      };
      const int req = static_cast<int>(g() % 5) - 1;
      const long P = static_cast<long>(g() % 12);
      // f[b][u]: the best decomposition whose last piece ends at b, u = req used.
      const long NEG = LONG_MIN / 4;
      std::vector<std::array<long, 2>> f(L + 1, {NEG, NEG}), M(L + 1, {NEG, NEG});
      for (int b = 1; b <= L; ++b) {
        for (int a = 0; a < b; ++a)
          for (int t = 0; t < 4; ++t)
            for (int u = 0; u < 2; ++u) {  // u: used before this piece
              long pre = u == 0 ? (a > 0 ? -5L : 0L) : NEG;
              if (M[a][u] > NEG) pre = std::max(pre, M[a][u] - P);
              if (pre <= NEG) continue;
              const int u2 = u | (t == req);
              f[b][u2] = std::max(f[b][u2], pre + VT(t, a, b));
            }
        for (int u = 0; u < 2; ++u) M[b][u] = std::max(M[b - 1][u], f[b][u]);
      }
      long brute = NEG;
      for (int b = 1; b <= L; ++b)
        for (int u = 0; u < 2; ++u)
          if (req < 0 || u == 1) brute = std::max(brute, f[b][u] - (b < L ? 5 : 0));
      ++n2;
      bad2 += certa::jump_bound(L, part, types, P, req) != brute;
    }
    EXPECT(bad == 0, "Lemma L-J: the one-pass bound differs from the brute force in %lld cases", bad);
    EXPECT(bad1 == 0, "Lemma L-J: without jumps, not the single-span bound in %lld cases", bad1);
    EXPECT(bad2 == 0, "jump_bound (free parts, a required type) differs from the brute force in %lld cases", bad2);
    std::fprintf(stderr, "Lemma L-J: %lld (layout, P) cases, one-pass bound = brute force: %s; = single-span bound "
                 "without jumps: %s; jump_bound with free parts and a required type = brute force on %lld: %s\n", n,
                 bad ? "no" : "yes", bad1 ? "no" : "yes", n2, bad2 ? "no" : "yes");
  }
  // ---- Theorem L-J: the pruning values (pre, suf, pot) against their
  // definitions over all decompositions, each piece bounded by Lemma L.
  {
    std::mt19937_64 g(11);
    long long n = 0, bad_pre = 0, bad_suf = 0, bad_pot = 0, bad_win = 0;
    const certa::IndexView ix{};  // no bands: the passes read no reference
    const long NEG = LONG_MIN / 4;
    for (int it = 0; it < 20000; ++it) {
      const int L = 3 + static_cast<int>(g() % 22);
      std::vector<certa::ReadPart> part[2], rpart[2];
      for (int s = 0; s < 2; ++s) {
        int x = static_cast<int>(g() % 3);
        for (;;) {
          const int len = 3 + static_cast<int>(g() % 4);
          if (x + len > L) break;
          if (g() % 4) part[s].push_back({x, x + len, 0});
          x += len + static_cast<int>(g() % 3);
        }
      }
      std::vector<char> fr[2], fc;  // the free parts of the any-types; of the chain (strand cs)
      for (int s = 0; s < 2; ++s)
        for (size_t k = 0; k < part[s].size(); ++k) fr[s].push_back(static_cast<char>(g() % 2));
      const int cs = static_cast<int>(g() % 2);
      for (size_t k = 0; k < part[cs].size(); ++k) fc.push_back(static_cast<char>(g() % 2));
      const long P = static_cast<long>(g() % 12);
      certa::PartTables tb[2], tbr[2];
      for (int s = 0; s < 2; ++s) {
        tb[s].build(L, part[s]);
        certa::mirror_parts(L, part[s], rpart[s]);
        tbr[s].build(L, rpart[s]);
      }
      const std::vector<uint8_t> rd(L, 0);
      std::vector<long> pre, mrev, suf(L + 1);
      double cells = 0;
      certa::jump_pass(ix, rd.data(), L, {}, tb, {{0, nullptr}, {1, nullptr}, {0, &fr[0]}, {1, &fr[1]}}, P, &pre,
                       &cells);
      certa::jump_pass(ix, rd.data(), L, {}, tbr, {{0, nullptr}, {1, nullptr}}, P, &mrev, &cells);
      for (int b = 0; b <= L; ++b) suf[b] = mrev[L - b];
      auto V = [&](int s, const std::vector<char>* fl, int a, int b) {  // Lemma L for a piece [a, b)
        long v = b - a;
        for (size_t k = 0; k < part[s].size(); ++k)
          v -= 5 * (part[s][k].off >= a && part[s][k].end <= b && !(fl && (*fl)[k]));
        return v;
      };
      // M[a]: the best decomposition with every piece ending at or before a
      // (pieces of all four types); S[b]: with every piece starting at or
      // after b (unwitnessed pieces only), the last one's clip included.
      std::vector<long> f(L + 1, NEG), M(L + 1, NEG), h(L + 1, NEG), S(L + 1, NEG);
      for (int b = 1; b <= L; ++b) {
        for (int a = 0; a < b; ++a) {
          const long st = std::max(a > 0 ? -5L : 0L, M[a] - P);
          for (int s = 0; s < 2; ++s) f[b] = std::max({f[b], st + V(s, nullptr, a, b), st + V(s, &fr[s], a, b)});
        }
        M[b] = std::max(M[b - 1], f[b]);
      }
      for (int a = L - 1; a >= 0; --a) {
        for (int b = a + 1; b <= L; ++b) {
          const long en = std::max(b < L ? -5L : 0L, S[b] - P);
          for (int s = 0; s < 2; ++s) h[a] = std::max(h[a], V(s, nullptr, a, b) + en);
        }
        S[a] = std::max(S[a + 1], h[a]);
      }
      // A piece of the chain contains one of its free parts.
      auto holds = [&](int a, int b) {
        for (size_t k = 0; k < part[cs].size(); ++k)
          if (fc[k] && part[cs][k].off >= a && part[cs][k].end <= b) return true;
        return false;
      };
      long pot = NEG;
      for (int a = 0; a < L; ++a)
        for (int b = a + 1; b <= L; ++b)
          if (holds(a, b))
            pot = std::max(pot, std::max(a > 0 ? -5L : 0L, M[a] - P) + V(cs, &fc, a, b) +
                                    std::max(b < L ? -5L : 0L, S[b] - P));
      ++n;
      bad_pre += pre != M;
      bad_suf += suf != S;
      bad_pot += certa::jump_pot(L, tb[cs], fc, pre, suf, P) != pot;
      // The window: the smallest [w0, w1) holding every [a, b) above lb.
      const long lb = pot > NEG ? pot - static_cast<long>(g() % 12) : 0;
      int a0 = L, b1 = -1;
      for (int a = 0; a < L; ++a)
        for (int b = a + 1; b <= L; ++b)
          if (holds(a, b) &&
              std::max(a > 0 ? -5L : 0L, M[a] - P) + V(cs, &fc, a, b) + std::max(b < L ? -5L : 0L, S[b] - P) > lb) {
            a0 = std::min(a0, a);
            b1 = std::max(b1, b);
          }
      const std::vector<char> fc_rev(fc.rbegin(), fc.rend());
      int w0 = -1, w1 = -1;
      certa::jump_window(L, tb[cs], tbr[cs], fc, fc_rev, pre, suf, P, lb, &w0, &w1);
      if (b1 < 0) bad_win += !(w0 == 0 && w1 == L);
      else bad_win += !(w0 == a0 && w1 == b1);
    }
    EXPECT(bad_pre == 0, "Theorem L-J: prefix values differ from their definition in %lld cases", bad_pre);
    EXPECT(bad_suf == 0, "Theorem L-J: suffix values (reversed pass) differ from their definition in %lld cases",
           bad_suf);
    EXPECT(bad_pot == 0, "Theorem L-J: pot differs from its definition in %lld cases", bad_pot);
    EXPECT(bad_win == 0, "Theorem L-J: the read window differs from its definition in %lld cases", bad_win);
    std::fprintf(stderr, "Theorem L-J: %lld random layouts, prefix values = definition: %s; suffix values (reversed "
                 "pass, mirrored parts) = definition: %s; pot = definition: %s; read window = definition: %s\n", n,
                 bad_pre ? "no" : "yes", bad_suf ? "no" : "yes", bad_pot ? "no" : "yes", bad_win ? "no" : "yes");
  }
  // ---- XQ: mapq_se is monotone (non-increasing in SUB, non-decreasing in the score).
  {
    long long n = 0, bad = 0;
    for (int L = 30; L <= 320; L += 5)
      for (int score = 0; score <= L; ++score)
        for (int sub = -5; sub <= L; ++sub) {
          ++n;
          const int q = certa::mapq_se(score, sub, L);
          bad += certa::mapq_se(score, sub + 1, L) > q;
          if (score < L) bad += certa::mapq_se(score + 1, sub, L) < q;
        }
    EXPECT(bad == 0, "mapq_se not monotone in %lld cases", bad);
    std::fprintf(stderr, "XQ: mapq_se monotone on %lld (L, score, SUB) triples: %s\n", n, bad ? "no" : "yes");
  }
  if (failures) {
    std::fprintf(stderr, "%d check(s) FAILED\n", failures);
    return 1;
  }
  std::fprintf(stderr, "all lemma checks passed\n");
  return 0;
}
