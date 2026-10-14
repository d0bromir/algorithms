// Aim 6 (EXACT_ALIGNMENT_PROPOSAL.md §6.6): does the exact-part certificate
// reach long reads? For each read, the Lemma L ceiling U1 (t = 1) over a
// greedy partition into parts that occur <= tau times (Theorem G, counted by
// anchor extension as in Lemma R, both strands), compared with the score of
// the read's best known alignment (minimap2 --eqx PAF, rescored under CERTA's
// scores). A read "fits" when sigma > U1.
//
//   lr_ceiling INDEX.cidx READS.fa ALN.paf [--tau 16] [--cmax 4096]
//              [--window 40] [--step 8] [-t threads] [--verify]
//              [--certify [--max-cells 4e9]] [--jump 10,30,60]
//              [--jump-certify P] > out.tsv
//
// --certify completes the certificate of each fitting read (certify_single):
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
// --jump-certify P: certify the split optimum OPT_J at penalty P (also in
// --jump) from lb_J (certify_jump in certa/longread.h): columns jcert (1
// proven, 3 an interval, 2 over the cell cap, 0 not attempted), jlo, jhi,
// jchains, jsearched, jcells, jus. With --jump-certify-new, only for the
// reads that fit the split objective but not the local one.
//
// The index must have s = 1. The parts, bounds and certificates are in
// certa/longread.h; this file reads the alignments and writes the table.
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

using Part = ReadPart;

struct Opts {
  int tau = 16, cmax = 4096, window = 40, step = 8, threads = 0;
  bool verify = false;
  bool certify = false;  // complete the certificate of fitting reads (banded DP)
  double max_cells = 4e9;  // per read; beyond it the read is left uncertified
  std::vector<long> jump;  // jump penalties P of the split objective (--jump)
  long jump_certify = -1;  // certify OPT_J at this P (--jump-certify)
  bool jump_certify_new = false;  // only for reads that fit J but not the local objective
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

inline bool same(uint8_t a, uint8_t b) { return lr_same(a, b); }

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
int verify_strand(const IndexView& ix, const uint8_t* rd, int L, const StrandParts& so) {
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
    else if (a == "--jump-certify" && i + 1 < argc) o.jump_certify = std::atol(argv[++i]);
    else if (a == "--jump-certify-new") o.jump_certify_new = true;
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
                         "[--window 40] [--step 8] [-t threads] [--verify] [--certify [--max-cells 4e9]] [--jump P,...] [--jump-certify P [--jump-certify-new]]\n");
    return 2;
  }
  Reference ref;
  Index index;
  load_index(pos[0], ref, index);
  const IndexView ix = index.view(ref);
  if (ix.s != 1) { std::fprintf(stderr, "the index must have s = 1 (has s = %d)\n", ix.s); return 2; }
  if (o.jump_certify >= 0 && std::find(o.jump.begin(), o.jump.end(), o.jump_certify) == o.jump.end())
    o.jump.push_back(o.jump_certify);
  PartitionParams pp;
  pp.tau = o.tau;
  pp.cmax = o.cmax;
  pp.window = o.window;
  pp.step = o.step;
  pp.keep_occ = o.certify || o.jump_certify >= 0;
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
  std::atomic<long> jproven{0}, jinterval{0}, jcapped{0}, jproven_new{0};
  auto work = [&]() {
    LrScratch sc;
    StrandParts so[2];
    BandDp dp;
    std::vector<PathBases> pbs;
    std::vector<Part> jparts[2];
    std::vector<uint8_t> rc;
    for (size_t r; (r = next++) < reads.size();) {
      const std::vector<uint8_t>& fw = reads[r];
      const int L = (int)fw.size();
      rc.resize(L);
      for (int i = 0; i < L; ++i) rc[i] = fw[L - 1 - i] < 4 ? 3 - fw[L - 1 - i] : 4;
      partition_strand(ix, fw.data(), L, pp, sc, so[0]);
      partition_strand(ix, rc.data(), L, pp, sc, so[1]);
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
      LrCert ct;
      if (o.certify && a.mapped) {
        const long lb = std::max<long>(a.sigma, trim);  // the best known real alignment
        const uint8_t* rds[2] = {fw.data(), rc.data()};
        ct = certify_single(ix, rds, L, so, lb, U, o.max_cells, dp);
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
      LrCert jc;
      if (!o.jump.empty()) {
        forward_parts(L, so, jparts);
        for (size_t k = 0; k < o.jump.size(); ++k) {
          const long P = o.jump[k];
          const long UJ = lemma_l_jump(L, jparts, P);
          long lbJ = a.mapped ? std::max<long>(a.sigma, trim) : LONG_MIN / 4;
          if (!a.lines.empty()) lbJ = std::max(lbJ, jump_lb(a.lines, L, P, pbs));
          const bool fj = lbJ > UJ;
          fits_j[k] += fj;
          fits_j_new[k] += fj && !fit;
          if (P == o.jump_certify && fj && !(o.jump_certify_new && fit)) {
            jc = certify_jump(ix, fw.data(), L, so, lbJ, P, o.max_cells);
            jproven += jc.status == 1;
            jinterval += jc.status == 3;
            jcapped += jc.status == 2;
            jproven_new += jc.status == 1 && !fit;
          }
          char jb[96];
          std::snprintf(jb, sizeof jb, "\t%ld\t%ld\t%d", UJ, lbJ > LONG_MIN / 8 ? lbJ : 0L, (int)fj);
          jcols += jb;
        }
      }
      const StrandParts& sa = so[a.mapped ? a.strand : 0];
      char buf[640];
      std::snprintf(buf, sizeof buf,
                    "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%ld\t%ld\t%ld\t%zu\t%zu\t%d\t%d\t%ld\t%ld\t%d\t%ld\t%d\t%ld"
                    "\t%d\t%ld\t%d\t%d\t%.0f\t%.0f",
                    names[r].c_str(), L, a.mapped, a.strand, a.qs, a.qe, a.mapped ? a.sigma : 0,
                    a.mapped ? a.sigma_all : 0, U, U0, U1, so[0].parts.size(), so[1].parts.size(), so[0].absent,
                    so[1].absent, sa.skipped, so[0].occ_sum + so[1].occ_sum, sa.longest_free,
                    a.mapped ? (long)a.sigma - U : 0, (int)fit, trim, ct.status, ct.opt, ct.chains, ct.searched,
                    ct.cells, ct.us);
      if (o.jump_certify >= 0) {
        char jb[160];
        std::snprintf(jb, sizeof jb, "\t%d\t%ld\t%ld\t%d\t%d\t%.0f\t%.0f", jc.status, jc.opt, jc.hi, jc.chains,
                      jc.searched, jc.cells, jc.us);
        jcols += jb;
      }
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
  if (o.jump_certify >= 0) std::printf("\tjcert\tjlo\tjhi\tjchains\tjsearched\tjcells\tjus");
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
  if (o.jump_certify >= 0)
    std::fprintf(stderr, "jump-certify P = %ld: OPT_J proven for %ld reads (%.2f %%; %ld of them do not fit the "
                 "single-piece objective), an interval for %ld, over the cell cap %ld\n", o.jump_certify,
                 jproven.load(), reads.empty() ? 0.0 : 100.0 * jproven.load() / reads.size(), jproven_new.load(),
                 jinterval.load(), jcapped.load());
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
