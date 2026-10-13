// Aim 6 (EXACT_ALIGNMENT_PROPOSAL.md §6.6): does the exact-part certificate
// reach long reads? For each read, the Lemma L ceiling U1 (t = 1) over a
// greedy partition into parts that occur <= tau times (Theorem G, counted by
// anchor extension as in Lemma R, both strands), compared with the score of
// the read's best known alignment (minimap2 --eqx PAF, rescored under CERTA's
// scores). A read "fits" when sigma > U1.
//
//   lr_ceiling INDEX.cidx READS.fa ALN.paf [--tau 16] [--cmax 4096]
//              [--window 40] [--step 8] [-t threads] [--verify] > out.tsv
//
// The index must have s = 1. Unlike tier R (repeat.h), reads may be of any
// length, and a stretch where no part of <= tau occurrences starts is skipped
// instead of ending the partition: Lemma L needs disjoint enumerated parts,
// not a tiling.
#include <algorithm>
#include <atomic>
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
#include "index.h"

using namespace certa;

namespace {

constexpr int kBig = INT_MAX / 2;

struct Part {
  int off, end;
  uint32_t occ;  // exact occurrences (0: absent)
};

struct Opts {
  int tau = 16, cmax = 4096, window = 40, step = 8, threads = 0;
  bool verify = false;
};

struct Aln {  // best known alignment of a read (minimap2)
  bool mapped = false;
  int strand = 0, qs = 0, qe = 0;
  int sigma = INT_MIN;      // the primary line, under CERTA's scores
  int sigma_all = INT_MIN;  // best over all lines
  int nlines = 0;
  std::string tname, cigar;  // the primary line's target and --eqx CIGAR
  int64_t tstart = 0;
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
      if (sc.ext.empty()) occ = (uint32_t)sc.elig.size();
      else
        for (int R : sc.ext) occ += c + q + R >= b;
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

// Lemma L, t = 1: every alignment with no exact occurrence of a part scores
// <= max over spans [x, y) of (y - x) - 5 (parts inside) - 5 (clipped ends).
// Span candidates x in {0, off + 1}, y in {L, end - 1} (Lemma S); O(P).
long u1(int L, const std::vector<Part>& p) {
  const int P = (int)p.size();
  const long c5 = kMismatch + kMatch;  // 5: the least loss of an inexact part
  long G = 0;                          // i = -1: x = 0
  long best = LONG_MIN;
  for (int j = 0; j <= P; ++j) {
    const long y = j == P ? L : p[j].end - 1;
    const long f = y - c5 * j - (y < L ? kClip : 0);
    best = std::max(best, f + G);
    if (j < P) {
      const long x = p[j].off + 1;
      G = std::max(G, -x + c5 * (j + 1) - kClip);
    }
  }
  return best;
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
  if (u1(L, so.parts) != u1_brute(L, so.parts)) {
    ++bad;
    std::fprintf(stderr, "verify: U1 %ld != brute force %ld\n", u1(L, so.parts), u1_brute(L, so.parts));
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
    else pos.push_back(argv[i]);
  }
  if (pos.size() != 3) {
    std::fprintf(stderr, "usage: lr_ceiling INDEX.cidx READS.fa ALN.paf [--tau 16] [--cmax 4096] "
                         "[--window 40] [--step 8] [-t threads] [--verify]\n");
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
  auto work = [&]() {
    Scratch sc;
    StrandOut so[2];
    std::vector<uint8_t> rc;
    for (size_t r; (r = next++) < reads.size();) {
      const std::vector<uint8_t>& fw = reads[r];
      const int L = (int)fw.size();
      rc.resize(L);
      for (int i = 0; i < L; ++i) rc[i] = fw[L - 1 - i] < 4 ? 3 - fw[L - 1 - i] : 4;
      partition(ix, fw.data(), L, o, sc, so[0]);
      partition(ix, rc.data(), L, o, sc, so[1]);
      const long U0 = u1(L, so[0].parts), U1 = u1(L, so[1].parts), U = std::max(U0, U1);
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
      const StrandOut& sa = so[a.mapped ? a.strand : 0];
      char buf[512];
      std::snprintf(buf, sizeof buf,
                    "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%ld\t%ld\t%ld\t%zu\t%zu\t%d\t%d\t%ld\t%ld\t%d\t%ld\t%d\t%ld\n",
                    names[r].c_str(), L, a.mapped, a.strand, a.qs, a.qe, a.mapped ? a.sigma : 0,
                    a.mapped ? a.sigma_all : 0, U, U0, U1, so[0].parts.size(), so[1].parts.size(), so[0].absent,
                    so[1].absent, sa.skipped, so[0].occ_sum + so[1].occ_sum, sa.longest_free,
                    a.mapped ? (long)a.sigma - U : 0, (int)fit, trim);
      rows[r] = buf;
    }
  };
  std::vector<std::thread> pool;
  for (int t = 0; t < nt; ++t) pool.emplace_back(work);
  for (auto& t : pool) t.join();
  std::printf("read\tlen\tmapped\tstrand\tqs\tqe\tsigma\tsigma_all\tU1\tU1_fwd\tU1_rev\tparts_fwd\tparts_rev"
              "\tabsent_fwd\tabsent_rev\tskipped\tocc_sum\tlongest_free\tmargin\tfits\tsigma_trim\n");
  for (const auto& row : rows) std::fputs(row.c_str(), stdout);
  std::fprintf(stderr, "reads %zu, fit %ld (%.2f %%), with the trimmed path %ld; CIGARs failing the base check %ld; "
               "tau %d, cmax %d%s\n", reads.size(), fits.load(),
               reads.empty() ? 0.0 : 100.0 * fits.load() / reads.size(), fits_trim.load(), path_bad.load(), o.tau, o.cmax,
               o.verify ? (verify_bad == 0 ? ", verify: all parts and U1 recounted equal"
                                           : ", verify: MISMATCHES (see above)")
                        : "");
  return (o.verify && verify_bad) || path_bad ? 1 : 0;
}
