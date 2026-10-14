// Aim 6 (EXACT_ALIGNMENT_PROPOSAL.md §6.6): does the exact-part certificate
// reach long reads? For each read, the Lemma L ceiling U1 (t = 1) over a
// greedy partition into parts that occur <= tau times (Theorem G, counted by
// anchor extension as in Lemma R, both strands), compared with the score of
// the read's best known alignment (minimap2 --eqx PAF, rescored under CERTA's
// scores). A read "fits" when sigma > U1.
//
//   lr_ceiling INDEX.cidx READS.fa ALN.paf [--tau 16] [--cmax 4096]
//              [--window 40] [--step 8] [-t threads] [--verify]
//              [--certify [--max-cells 4e9]] [--jump 10,30,60] > out.tsv
//
// --certify completes the certificate of each fitting read (certify_read):
// OPT is then proven exactly, and the columns cert, opt, chains, searched,
// cells, us give it and its cost. With --verify, a full DP around the best
// known alignment's locus must not exceed the proven OPT.
//
// --jump P,...: the split (jump) objective of certa/longread.h at each jump
// penalty P: U_J (Lemma L-J), the best jump alignment assembled from
// minimap2's lines (lb_J: pieces of its primary, supplementary and secondary
// alignments, each also split at its gaps longer than P - 6, where a jump is
// cheaper; in read order, a later piece trimmed where it overlaps the earlier
// one), and fits_J (lb_J > U_J).
//
// The index must have s = 1. Unlike tier R (repeat.h), reads may be of any
// length, and a stretch where no part of <= tau occurrences starts is skipped
// instead of ending the partition: Lemma L needs disjoint enumerated parts,
// not a tiling.
#include <algorithm>
#include <atomic>
#include <chrono>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#include "certa/core.h"
#include "certa/longread.h"
#include "index.h"

using namespace certa;

namespace {

constexpr int kBig = INT_MAX / 2;

using Part = ReadPart;

struct Opts {
  int tau = 16, cmax = 4096, window = 40, step = 8, threads = 0;
  bool verify = false;
  bool certify = false;  // complete the certificate of fitting reads (banded DP)
  double max_cells = 4e9;  // per read; beyond it the read is left uncertified
  std::vector<long> jump;  // jump penalties P of the split objective (--jump)
};

struct Line {  // one PAF line (an alignment of a read interval [qs, qe))
  int qs, qe, strand;
  std::string cigar;
};

struct Aln {  // best known alignment of a read (minimap2)
  bool mapped = false;
  int strand = 0, qs = 0, qe = 0;
  int sigma = INT_MIN;      // the primary line, under CERTA's scores
  int sigma_all = INT_MIN;  // best over all lines
  int nlines = 0;
  std::string tname, cigar;  // the primary line's target and --eqx CIGAR
  int64_t tstart = 0;
  std::vector<Line> lines;  // every line (for the jump objective)
};

inline bool same(uint8_t a, uint8_t b) { return a < 4 && a == b; }

// CERTA's local score of a PAF line from its --eqx CIGAR: match +1,
// mismatch -4, gap of k bases -(6 + k), 5 per clipped end.
int score_paf(const char* cg, int qlen, int qs, int qe) {
  long s = 0;
  long n = 0;
  for (const char* p = cg; *p && *p != '\t' && *p != '\n'; ++p) {
    if (*p >= '0' && *p <= '9') { n = n * 10 + (*p - '0'); continue; }
    switch (*p) {
      case '=': s += n * kMatch; break;
      case 'X': s -= n * kMismatch; break;
      case 'I': case 'D': s -= kGapOpen + n * kGapExt; break;
      case 'M': return INT_MIN;  // needs --eqx
      default: break;
    }
    n = 0;
  }
  if (qs > 0) s -= kClip;
  if (qe < qlen) s -= kClip;
  return (int)s;
}

std::unordered_map<std::string, Aln> read_paf(const char* path) {
  std::unordered_map<std::string, Aln> m;
  std::ifstream f(path);
  if (!f) { std::fprintf(stderr, "cannot open %s\n", path); std::exit(1); }
  std::string line;
  while (std::getline(f, line)) {
    std::vector<const char*> col;
    col.push_back(line.c_str());
    for (size_t i = 0; i < line.size(); ++i)
      if (line[i] == '\t') col.push_back(line.c_str() + i + 1);
    if (col.size() < 12) continue;
    const std::string name(col[0], col[1] - 1);
    const int qlen = std::atoi(col[1]), qs = std::atoi(col[2]), qe = std::atoi(col[3]);
    const int strand = col[4][0] == '-';
    const char* cg = nullptr;
    char tp = 'P';
    for (size_t c = 12; c < col.size(); ++c) {
      if (!std::strncmp(col[c], "cg:Z:", 5)) cg = col[c] + 5;
      if (!std::strncmp(col[c], "tp:A:", 5)) tp = col[c][5];
    }
    if (!cg) continue;
    const int s = score_paf(cg, qlen, qs, qe);
    if (s == INT_MIN) { std::fprintf(stderr, "PAF needs --eqx CIGARs\n"); std::exit(1); }
    Aln& a = m[name];
    a.sigma_all = std::max(a.sigma_all, s);
    ++a.nlines;
    {
      const char* ce = cg;
      while (*ce && *ce != '\t') ++ce;
      a.lines.push_back({qs, qe, strand, std::string(cg, ce)});
    }
    if (tp == 'P' && !a.mapped) {  // minimap2 writes the primary first
      a.mapped = true;
      a.strand = strand;
      a.qs = qs;
      a.qe = qe;
      a.sigma = s;
      a.tname.assign(col[5], col[6] - 1);
      a.tstart = std::atoll(col[7]);
      const char* ce = cg;
      while (*ce && *ce != '\t') ++ce;
      a.cigar.assign(cg, ce);
    }
  }
  return m;
}

// Walks the primary CIGAR over the reference: every '=' base must match and
// every 'X' base differ (else the PAF's score is not to be trusted). Also
// returns the best sub-alignment of the path (from the start of a '=' run to
// the end of one), itself an alignment, so also <= OPT.
struct PathCheck {
  long bad = 0;
  long trim = LONG_MIN;
};
PathCheck check_path(const IndexView& ix, int64_t tpos, const uint8_t* rd, int L, int r0, const std::string& cg) {
  PathCheck pc;
  long cur = LONG_MIN / 4, n = 0;
  int r = r0;
  int64_t t = tpos;
  for (char ch : cg) {
    if (ch >= '0' && ch <= '9') { n = n * 10 + (ch - '0'); continue; }
    if (ch == '=' || ch == 'X') {
      for (long k = 0; k < n; ++k) pc.bad += (ch == '=') != same(rd[r + k], ref_at(ix, t + k));
    }
    if (ch == '=') {
      cur = std::max(cur + n * kMatch, n * kMatch - (r > 0 ? kClip : 0));
      r += n;
      t += n;
      pc.trim = std::max(pc.trim, cur - (r < L ? kClip : 0));
    } else if (ch == 'X') {
      cur -= n * kMismatch;
      r += n;
      t += n;
    } else if (ch == 'I') {
      cur -= kGapOpen + n * kGapExt;
      r += n;
    } else if (ch == 'D') {
      cur -= kGapOpen + n * kGapExt;
      t += n;
    }
    n = 0;
  }
  return pc;
}

// A line's path as its aligned bases, in path order: the forward read
// position of each, the path's score after it (gaps included), and its own
// score.
struct PathBases {
  std::vector<int> pos, own;
  std::vector<long> cum;
  bool rev = false;  // positions decrease along the path (reverse strand)
  std::vector<std::pair<int, long>> gaps;  // (index of the aligned base before it, length)
};
void path_bases(const Line& l, int L, PathBases& pb) {
  pb.pos.clear();
  pb.own.clear();
  pb.cum.clear();
  pb.gaps.clear();
  pb.rev = l.strand != 0;
  long cum = 0, n = 0;
  int r = l.strand ? L - l.qe : l.qs;  // read position in the aligned orientation
  for (char ch : l.cigar) {
    if (ch >= '0' && ch <= '9') { n = n * 10 + (ch - '0'); continue; }
    if (ch == '=' || ch == 'X') {
      const int d = ch == '=' ? kMatch : -kMismatch;
      for (long k = 0; k < n; ++k, ++r) {
        cum += d;
        pb.pos.push_back(l.strand ? L - 1 - r : r);
        pb.own.push_back(d);
        pb.cum.push_back(cum);
      }
    } else if (ch == 'I') {
      cum -= kGapOpen + n * kGapExt;
      r += (int)n;
      pb.gaps.push_back({(int)pb.pos.size() - 1, n});
    } else if (ch == 'D') {
      cum -= kGapOpen + n * kGapExt;
      pb.gaps.push_back({(int)pb.pos.size() - 1, n});
    }
    n = 0;
  }
}
// The path restricted to forward read positions [u, v): from its first
// aligned base there to its last, gaps between included. A real alignment
// (no clip penalty: it is a piece). LONG_MIN if it has no base there.
long piece_score(const PathBases& pb, int u, int v) {
  const size_t n = pb.pos.size();
  size_t k1, k2;
  if (!pb.rev) {
    k1 = std::lower_bound(pb.pos.begin(), pb.pos.end(), u) - pb.pos.begin();
    const size_t e = std::lower_bound(pb.pos.begin(), pb.pos.end(), v) - pb.pos.begin();
    if (k1 >= e) return LONG_MIN;
    k2 = e - 1;
  } else {
    k1 = std::lower_bound(pb.pos.begin(), pb.pos.end(), v - 1, std::greater<int>()) - pb.pos.begin();
    const size_t e = std::lower_bound(pb.pos.begin(), pb.pos.end(), u - 1, std::greater<int>()) - pb.pos.begin();
    if (k1 >= e || k1 >= n) return LONG_MIN;
    k2 = e - 1;
  }
  return pb.cum[k2] - (pb.cum[k1] - pb.own[k1]);
}
// lb_J: the best jump alignment made of pieces of minimap2's lines: each
// line whole, and its segments between gaps longer than P - 6 (a jump costs
// less than such a gap). Pieces in read order; a later piece is trimmed to
// start where the earlier one ends.
long jump_lb(const std::vector<Line>& lines, int L, long P, std::vector<PathBases>& pbs) {
  const long NEG = LONG_MIN / 4;
  const int n = (int)lines.size();
  pbs.resize(n);
  struct Cand {
    int line, lo, hi;  // forward read interval [lo, hi)
  };
  std::vector<Cand> cand;
  for (int i = 0; i < n; ++i) {
    PathBases& pb = pbs[i];
    path_bases(lines[i], L, pb);
    if (pb.pos.empty()) continue;
    auto add = [&](int k1, int k2) {  // aligned bases k1..k2 (path order)
      const int a = pb.rev ? pb.pos[k2] : pb.pos[k1], b = (pb.rev ? pb.pos[k1] : pb.pos[k2]) + 1;
      cand.push_back({i, a, b});
    };
    add(0, (int)pb.pos.size() - 1);
    int k1 = 0;
    for (const auto& gp : pb.gaps)
      if (gp.second > P - kGapOpen && gp.first >= k1 && gp.first + 1 < (int)pb.pos.size()) {
        add(k1, gp.first);
        k1 = gp.first + 1;
      }
    if (k1 > 0) add(k1, (int)pb.pos.size() - 1);
  }
  std::sort(cand.begin(), cand.end(), [](const Cand& x, const Cand& y) { return x.lo < y.lo; });
  const int m = (int)cand.size();
  std::vector<long> be(m, NEG);  // the best decomposition whose last piece is (part of) cand j
  long best = NEG;
  for (int j = 0; j < m; ++j) {
    const Cand& c = cand[j];
    const long whole = piece_score(pbs[c.line], c.lo, c.hi);
    if (whole > NEG) be[j] = whole - (c.lo > 0 ? kClip : 0);
    for (int i = 0; i < j; ++i) {
      if (be[i] <= NEG || cand[i].hi >= c.hi) continue;
      const long w = piece_score(pbs[c.line], std::max(c.lo, cand[i].hi), c.hi);
      if (w > NEG) be[j] = std::max(be[j], be[i] - P + w);
    }
    if (be[j] > NEG) best = std::max(best, be[j] - (c.hi < L ? kClip : 0));
  }
  return best;
}

struct Scratch {
  std::vector<uint64_t> lo;
  std::vector<uint32_t> cnt;
  std::vector<int> nz;
  std::vector<int64_t> elig;
  std::vector<int> ext;
};

struct StrandOut {
  std::vector<Part> parts;
  long skipped = 0;   // read bases where no part could start
  long occ_sum = 0;   // enumerated occurrences
  int absent = 0;
  int longest_free = 0;  // longest stretch covered by no part
  // --certify: every exact occurrence of a rare part, as (diagonal = the
  // reference position of read base 0, part index).
  std::vector<std::pair<int64_t, int>> occ;
};

// Greedy partition of one strand (Theorem G at tau, Lemma R).
void partition(const IndexView& ix, const uint8_t* rd, int L, const Opts& o, Scratch& sc, StrandOut& out) {
  const int q = ix.q;
  out = StrandOut();
  if (L < q) { out.skipped = L; out.longest_free = L; return; }
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
        for (int k = c - 1; k >= s && ok; --k) ok = same(rd[k], ref_at(ix, h - (c - k)));
        if (ok) sc.elig.push_back(h);
      }
      if ((int)sc.elig.size() <= o.tau) {
        banc = c + q;
      } else {
        for (int64_t h : sc.elig) {
          int R = 0;
          while (c + q + R < L && same(rd[c + q + R], ref_at(ix, h + q + R))) ++R;
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
          if (o.certify) out.occ.push_back({sc.elig[k] - c, part});
        }
    }
    out.parts.push_back({s, b, occ});
    out.occ_sum += occ;
    out.absent += !rare;
    s = b;
  }
  if (s < L) out.skipped += L - s;  // a tail shorter than q
  int prev = 0;
  for (const Part& p : out.parts) {
    out.longest_free = std::max(out.longest_free, p.off - prev);
    prev = p.end;
  }
  out.longest_free = std::max(out.longest_free, L - prev);
}

long u1_brute(int L, const std::vector<Part>& p) {
  std::vector<long> xs{0}, ys{L};
  for (const Part& q : p) { xs.push_back(q.off + 1); ys.push_back(q.end - 1); }
  long best = LONG_MIN;
  for (long x : xs)
    for (long y : ys) {
      if (y <= x) continue;
      long in = 0;
      for (const Part& q : p) in += q.off >= x && q.end <= y;
      best = std::max(best, (y - x) - 5 * in - (x > 0 ? kClip : 0) - (y < L ? kClip : 0));
    }
  return best;
}

// The best local score (BWA-MEM's scores, -5 per clipped end) of rd over
// alignments whose cells all lie on diagonals [dlo, dhi] (reference
// position minus read position). Gotoh, one row of the band at a time.
struct Dp {
  std::vector<int> H, I, Hn, In;
};
long band_best(const IndexView& ix, const uint8_t* rd, int L, int64_t dlo, int64_t dhi, Dp& dp, double* cells) {
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
      const int m = std::max(dp.H[b], start) + (same(c, ref_at(ix, j0 + b)) ? kMatch : -kMismatch);
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

// The rest of the certificate of a fitting read (lb = the best known
// alignment's score > U1): every alignment above lb has gap total <= pad =
// L - lb - 7 (Lemma W) and, having an exact part (it beats U1), lies in a
// chain of part occurrences on that strand (split where diagonals are more
// than pad apart), within its diagonals +- pad (Lemma C), scoring at most
// the chain's Lemma L bound with its parts free. So OPT = max(lb, the band
// maxima of the chains whose bound exceeds the best so far).
struct Cert {
  int status = 0;  // 0: not attempted (no fit); 1: OPT proven; 2: over the cell cap
  long opt = 0;
  int chains = 0, searched = 0;
  double cells = 0, us = 0;
};
Cert certify_read(const IndexView& ix, const uint8_t* const* rd, int L, StrandOut* so, long lb, long U, const Opts& o,
                  Dp& dp) {
  Cert c;
  if (lb <= U) return c;
  const auto t0 = std::chrono::steady_clock::now();
  const long pad = std::max(0L, (long)L - lb - 7);
  long best = lb;
  std::vector<char> free;
  for (int st = 0; st < 2; ++st) {
    auto& occ = so[st].occ;
    std::sort(occ.begin(), occ.end());
    const std::vector<Part>& P = so[st].parts;
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
      if (c.cells + (double)L * (double)(dhi - dlo + 1) > o.max_cells) {
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

// --verify: recount each rare part from a different anchor (its first
// q-mer) and check that U1 equals the brute-force maximum.
int verify_strand(const IndexView& ix, const uint8_t* rd, int L, const StrandOut& so) {
  int bad = 0;
  for (const Part& p : so.parts) {
    uint64_t key, lo = 0, hi = 0;
    if (!kmer_key(rd + p.off, ix.q, &key)) { bad += p.occ != 0; continue; }
    lookup(ix, key, &lo, &hi);
    if (hi - lo > 200000) continue;  // too many to recount
    uint32_t n = 0;
    for (uint64_t e = lo; e < hi; ++e) {
      const int64_t h = ix.pos[e];
      bool ok = true;
      for (int k = ix.q; ok && k < p.end - p.off; ++k) ok = same(rd[p.off + k], ref_at(ix, h + k));
      n += ok;
    }
    if (n != p.occ) {
      ++bad;
      std::fprintf(stderr, "verify: part [%d, %d) counted %u, recount %u\n", p.off, p.end, p.occ, n);
    }
  }
  if (lemma_l_span(L, so.parts) != u1_brute(L, so.parts)) {
    ++bad;
    std::fprintf(stderr, "verify: U1 %ld != brute force %ld\n", lemma_l_span(L, so.parts), u1_brute(L, so.parts));
  }
  return bad;
}

}  // namespace

int main(int argc, char** argv) {
  Opts o;
  std::vector<const char*> pos;
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    if (a == "--tau" && i + 1 < argc) o.tau = std::atoi(argv[++i]);
    else if (a == "--cmax" && i + 1 < argc) o.cmax = std::atoi(argv[++i]);
    else if (a == "--window" && i + 1 < argc) o.window = std::atoi(argv[++i]);
    else if (a == "--step" && i + 1 < argc) o.step = std::atoi(argv[++i]);
    else if (a == "-t" && i + 1 < argc) o.threads = std::atoi(argv[++i]);
    else if (a == "--verify") o.verify = true;
    else if (a == "--certify") o.certify = true;
    else if (a == "--max-cells" && i + 1 < argc) o.max_cells = std::atof(argv[++i]);
    else if (a == "--jump" && i + 1 < argc) {
      for (const char* q = argv[++i]; *q;) {
        o.jump.push_back(std::strtol(q, const_cast<char**>(&q), 10));
        if (*q == ',') ++q;
        else if (*q) break;
      }
    }
    else pos.push_back(argv[i]);
  }
  if (pos.size() != 3) {
    std::fprintf(stderr, "usage: lr_ceiling INDEX.cidx READS.fa ALN.paf [--tau 16] [--cmax 4096] "
                         "[--window 40] [--step 8] [-t threads] [--verify] [--certify [--max-cells 4e9]] [--jump P,...]\n");
    return 2;
  }
  Reference ref;
  Index index;
  load_index(pos[0], ref, index);
  const IndexView ix = index.view(ref);
  if (ix.s != 1) { std::fprintf(stderr, "the index must have s = 1 (has s = %d)\n", ix.s); return 2; }
  const auto alns = read_paf(pos[2]);
  std::unordered_map<std::string, int64_t> contig;
  for (size_t c = 0; c < ref.names.size(); ++c) contig[ref.names[c]] = (int64_t)ref.offsets[c];

  std::vector<std::string> names;
  std::vector<std::vector<uint8_t>> reads;
  {
    std::ifstream f(pos[1]);
    if (!f) { std::fprintf(stderr, "cannot open %s\n", pos[1]); return 2; }
    std::string line;
    while (std::getline(f, line)) {
      if (line.empty()) continue;
      if (line[0] == '>') {
        names.push_back(line.substr(1, line.find_first_of(" \t") - 1));
        reads.emplace_back();
      } else if (!reads.empty()) {
        for (char ch : line) reads.back().push_back(encode_base(ch));
      }
    }
  }
  const int nt = o.threads > 0 ? o.threads : (int)std::thread::hardware_concurrency();
  std::vector<std::string> rows(reads.size());
  std::atomic<size_t> next{0};
  std::atomic<long> fits{0}, verify_bad{0}, path_bad{0}, fits_trim{0};
  std::atomic<long> certified{0}, capped{0}, improved{0}, local_checked{0}, local_equal{0}, local_above{0};
  std::vector<std::atomic<long>> fits_j(o.jump.size()), fits_j_new(o.jump.size());
  auto work = [&]() {
    Scratch sc;
    StrandOut so[2];
    Dp dp;
    std::vector<PathBases> pbs;
    std::vector<Part> jparts[2];
    std::vector<uint8_t> rc;
    for (size_t r; (r = next++) < reads.size();) {
      const std::vector<uint8_t>& fw = reads[r];
      const int L = (int)fw.size();
      rc.resize(L);
      for (int i = 0; i < L; ++i) rc[i] = fw[L - 1 - i] < 4 ? 3 - fw[L - 1 - i] : 4;
      partition(ix, fw.data(), L, o, sc, so[0]);
      partition(ix, rc.data(), L, o, sc, so[1]);
      const long U0 = lemma_l_span(L, so[0].parts), U1 = lemma_l_span(L, so[1].parts), U = std::max(U0, U1);
      if (o.verify) verify_bad += verify_strand(ix, fw.data(), L, so[0]) + verify_strand(ix, rc.data(), L, so[1]);
      auto it = alns.find(names[r]);
      const Aln a = it == alns.end() ? Aln() : it->second;
      const bool fit = a.mapped && a.sigma > U;
      fits += fit;
      long trim = a.mapped ? a.sigma : 0;
      if (a.mapped) {
        auto ct = contig.find(a.tname);
        if (ct == contig.end()) {
          ++path_bad;
        } else {
          const uint8_t* rs = a.strand ? rc.data() : fw.data();
          const PathCheck pc = check_path(ix, ct->second + a.tstart, rs, L, a.strand ? L - a.qe : a.qs, a.cigar);
          path_bad += pc.bad > 0;
          trim = std::max(trim, pc.trim);
        }
      }
      fits_trim += a.mapped && trim > U;
      Cert ct;
      if (o.certify && a.mapped) {
        const long lb = std::max<long>(a.sigma, trim);  // the best known real alignment
        const uint8_t* rds[2] = {fw.data(), rc.data()};
        ct = certify_read(ix, rds, L, so, lb, U, o, dp);
        certified += ct.status == 1;
        capped += ct.status == 2;
        improved += ct.status == 1 && ct.opt > lb;
        auto cg = contig.find(a.tname);
        if (o.verify && ct.status == 1 && cg != contig.end()) {
          // Every alignment near the best known one, on its strand, by a
          // full DP: none may beat the proven OPT.
          const int64_t t0 = cg->second + a.tstart;
          double cells = 0;
          const long near = band_best(ix, rds[a.strand], L, t0 - 1000 - L, t0 + L + 1000, dp, &cells);
          ++local_checked;
          local_equal += near == ct.opt;
          if (near > ct.opt) {
            ++verify_bad;
            ++local_above;
            std::fprintf(stderr, "verify: %s: a full DP near the known locus scores %ld, above the proven OPT %ld\n",
                         names[r].c_str(), near, ct.opt);
          }
        }
      }
      std::string jcols;
      if (!o.jump.empty()) {
        // Strand 1's parts in forward read coordinates.
        jparts[0] = so[0].parts;
        jparts[1].clear();
        for (auto q = so[1].parts.rbegin(); q != so[1].parts.rend(); ++q)
          jparts[1].push_back({L - q->end, L - q->off, q->occ});
        for (size_t k = 0; k < o.jump.size(); ++k) {
          const long P = o.jump[k];
          const long UJ = lemma_l_jump(L, jparts, P);
          long lbJ = a.mapped ? std::max<long>(a.sigma, trim) : LONG_MIN / 4;
          if (!a.lines.empty()) lbJ = std::max(lbJ, jump_lb(a.lines, L, P, pbs));
          const bool fj = lbJ > UJ;
          fits_j[k] += fj;
          fits_j_new[k] += fj && !fit;
          char jb[96];
          std::snprintf(jb, sizeof jb, "\t%ld\t%ld\t%d", UJ, lbJ > LONG_MIN / 8 ? lbJ : 0L, (int)fj);
          jcols += jb;
        }
      }
      const StrandOut& sa = so[a.mapped ? a.strand : 0];
      char buf[640];
      std::snprintf(buf, sizeof buf,
                    "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%ld\t%ld\t%ld\t%zu\t%zu\t%d\t%d\t%ld\t%ld\t%d\t%ld\t%d\t%ld"
                    "\t%d\t%ld\t%d\t%d\t%.0f\t%.0f",
                    names[r].c_str(), L, a.mapped, a.strand, a.qs, a.qe, a.mapped ? a.sigma : 0,
                    a.mapped ? a.sigma_all : 0, U, U0, U1, so[0].parts.size(), so[1].parts.size(), so[0].absent,
                    so[1].absent, sa.skipped, so[0].occ_sum + so[1].occ_sum, sa.longest_free,
                    a.mapped ? (long)a.sigma - U : 0, (int)fit, trim, ct.status, ct.opt, ct.chains, ct.searched,
                    ct.cells, ct.us);
      rows[r] = std::string(buf) + jcols + "\n";
    }
  };
  std::vector<std::thread> pool;
  for (int t = 0; t < nt; ++t) pool.emplace_back(work);
  for (auto& t : pool) t.join();
  std::printf("read\tlen\tmapped\tstrand\tqs\tqe\tsigma\tsigma_all\tU1\tU1_fwd\tU1_rev\tparts_fwd\tparts_rev"
              "\tabsent_fwd\tabsent_rev\tskipped\tocc_sum\tlongest_free\tmargin\tfits\tsigma_trim"
              "\tcert\topt\tchains\tsearched\tcells\tus");
  for (long P : o.jump) std::printf("\tUJ_%ld\tlbJ_%ld\tfitsJ_%ld", P, P, P);
  std::printf("\n");
  for (const auto& row : rows) std::fputs(row.c_str(), stdout);
  std::fprintf(stderr, "reads %zu, fit %ld (%.2f %%), with the trimmed path %ld; CIGARs failing the base check %ld; "
               "tau %d, cmax %d%s\n", reads.size(), fits.load(),
               reads.empty() ? 0.0 : 100.0 * fits.load() / reads.size(), fits_trim.load(), path_bad.load(), o.tau, o.cmax,
               o.verify ? (verify_bad == 0 ? ", verify: all parts and U1 recounted equal"
                                           : ", verify: MISMATCHES (see above)")
                        : "");
  for (size_t k = 0; k < o.jump.size(); ++k)
    std::fprintf(stderr, "jump P = %ld: fit %ld (%.2f %%), of which %ld do not fit the single-piece objective\n",
                 o.jump[k], fits_j[k].load(), reads.empty() ? 0.0 : 100.0 * fits_j[k].load() / reads.size(),
                 fits_j_new[k].load());
  if (o.certify)
    std::fprintf(stderr, "certify: OPT proven for %ld reads (%.2f %%), above the best known alignment for %ld; over the "
                 "cell cap %ld%s\n", certified.load(), reads.empty() ? 0.0 : 100.0 * certified.load() / reads.size(),
                 improved.load(), capped.load(),
                 o.verify ? (" ; full DP near the known locus: " + std::to_string(local_checked.load()) +
                             " reads, equal to OPT for " + std::to_string(local_equal.load()) + ", above it for " +
                             std::to_string(local_above.load()))
                                .c_str()
                          : "");
  return (o.verify && verify_bad) || path_bad ? 1 : 0;
}
