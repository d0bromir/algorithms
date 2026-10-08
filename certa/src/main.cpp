// certa: prototype of the certified short-read fast path.
//
//   certa index ref.fa -o ref.cidx [-q 22] [-s 8] [-t threads]
//   certa map ref.cidx reads.fq[.gz] [-k 4] [--budget 256] [-t threads]
//             [--gpu [--device 0]] [-o out.sam] [-u uncertified.fq]
//             [--stats stats.json] [--batch N] [--io-threads N]
//
// Certified reads are written to SAM. All other reads are written unchanged
// to the uncertified FASTQ, for a full aligner (minibwa, BWA-MEM2, ...).
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <map>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "certa/pair.h"
#include "index.h"
#include "mapper.h"
#include "seqio.h"

using namespace certa;

namespace {

using Clock = std::chrono::steady_clock;
double secs(Clock::time_point a) {
  return std::chrono::duration<double>(Clock::now() - a).count();
}

const char* const kReasonNames[kNumReasons] = {
    "certified", "bad_length", "radius_negative", "not_found",
    "cluster_too_wide", "cross_contig"};

struct Stats {
  uint64_t reads = 0, bases = 0, certified = 0, s0 = 0, s1 = 0, sr = 0, sl = 0, s2 = 0, ties = 0;
  uint64_t reason[kNumReasons] = {};
  uint64_t by_radius[KMAX + 2] = {};  // certified reads per radius R
  uint64_t by_d1[KMAX + 1] = {};
  uint64_t unc_by_used[PMAX + 1] = {};  // uncertified reads per parts searched |S|
  uint64_t pr = 0;                      // reads certified by the pair certificate (PR)
  uint64_t pairs = 0, pairs_cert = 0, pairs_proper = 0;
  uint64_t cands = 0;                   // candidate diagonals enumerated (last pass per read)
  void add(const Stats& o) {
    reads += o.reads; bases += o.bases; certified += o.certified;
    s0 += o.s0; s1 += o.s1; sr += o.sr; sl += o.sl; s2 += o.s2; ties += o.ties;
    for (int i = 0; i < kNumReasons; ++i) reason[i] += o.reason[i];
    for (int i = 0; i < KMAX + 2; ++i) by_radius[i] += o.by_radius[i];
    for (int i = 0; i < KMAX + 1; ++i) by_d1[i] += o.by_d1[i];
    for (int i = 0; i <= PMAX; ++i) unc_by_used[i] += o.unc_by_used[i];
    cands += o.cands;
    pr += o.pr; pairs += o.pairs; pairs_cert += o.pairs_cert; pairs_proper += o.pairs_proper;
  }
};

// Provisional MAPQ from the certified distance gap (to be calibrated).
// Provisional MAPQ from the distance gap (to be calibrated). Certified tiers:
// the second-best distance is exact (or > R). S2: only candidates were
// searched, so the scale is halved and capped at 40.
int g_s2_limit = 0;
// BWA-MEM's single-end MAPQ, 30 (1 - S2/S1) ln(L), with alignment scores
// approximated from edit distances (each edit costs match + mismatch = 5).
int bwa_mapq(int L, int d1, int d2) {
  const double s1 = L * kMatch - d1 * (kMatch + kMismatch);
  const double s2 = L * kMatch - d2 * (kMatch + kMismatch);
  if (s1 <= 0 || s2 >= s1) return 0;
  const int q = static_cast<int>(30.0 * (1.0 - s2 / s1) * std::log(static_cast<double>(L)) + 0.499);
  return std::max(0, std::min(60, q));
}

int mapq_of(const Result& r, int L) {
  if (r.n_best > 1) return 0;
  if (r.tier == kTierS2) {
    // Second best: exact among the evaluated candidates, else at least D + 1.
    int q = bwa_mapq(L, r.d1, r.d2 >= 0 ? r.d2 : g_s2_limit + 1);
    // Clusters beyond the S2_TOP best-supported ones were not evaluated, so
    // an unseen tie is possible: keep MAPQ below the usual caller cut-off.
    if (r.n_clusters > S2_TOP) q = std::min(q, 10);
    return q;
  }
  // Certified: the second locus is exact within R, or beyond R the best other
  // candidate found (d2x), as bwa-mem uses its suboptimal hit. MAPQ follows
  // bwa-mem's default single-end formula (mem_approx_mapq_se with
  // mapQ_coef_len = 50) on affine scores, with sub = 19 (min_seed_len x a)
  // when no second locus was found.
  if (r.score > 0) {
    int sub = r.sub_score != kNoSub ? r.sub_score : 19;
    if (sub < 19) sub = 19;
    if (sub >= r.score) return 0;
    const double l = L, identity = 1.0 - (l * kMatch - r.score) / (kMatch + kMismatch) / l;
    double t = l < 50 ? 1.0 : std::log(50.0) / std::log(l);
    t *= identity * identity;
    int q = static_cast<int>(6.02 * (r.score - sub) / kMatch * t * t + 0.499);
    // bwa subtracts 4.343 ln(sub_n + 1) for near-optimal alternative hits;
    // lower bound: one when the second locus is within one edit of the best.
    const int second = r.d2 >= 0 ? r.d2 : r.d2x;
    if (second >= 0 && second <= r.d1 + 1) q -= static_cast<int>(4.343 * std::log(2.0) + 0.499);
    return std::max(0, std::min(60, q));
  }
  const int second = r.d2 >= 0 ? r.d2 : r.d2x;
  return second >= 0 ? bwa_mapq(L, r.d1, second) : 60;
}

char comp(char c) {
  switch (c) {
    case 'A': case 'a': return 'T';
    case 'C': case 'c': return 'G';
    case 'G': case 'g': return 'C';
    case 'T': case 't': return 'A';
    default: return 'N';
  }
}

int64_t ref_span(const Result& r) {
  int64_t span = 0;
  for (int i = 0; i < r.n_cigar; ++i) {
    const uint32_t op = r.cigar[i] & 0xF;
    if (op == kOpM || op == kOpD) span += r.cigar[i] >> 4;
  }
  return span;
}

// Host-side check: a certified alignment must lie within one contig.
void check_contig(Result& r, const Reference& ref) {
  if (!r.certified) return;
  const int c = ref.contig_of(r.ref_pos);
  if (c < 0 || ref.contig_of(r.ref_pos + ref_span(r) - 1) != c) {
    r.certified = 0;
    r.reason = kCrossContig;
  }
}

void fastq_record(const FastqRecord& rec, std::string& fq) {
  fq += '@'; fq += rec.name; fq += '\n';
  fq += rec.seq; fq += "\n+\n";
  fq += rec.qual; fq += '\n';
}

// Mate fields of a SAM line (single-end: none).
struct MateInfo {
  int flag = 0;
  std::string rnext = "*";
  int64_t pnext = 0, tlen = 0;
};

void sam_line(const Result& r, const FastqRecord& rec, const Reference& ref, int mapq,
              const MateInfo& mi, std::string& sam);

// Appends a SAM line (certified) or a FASTQ record (fallback). Returns the
// final reason code after the host-side contig check.
int emit(const Result& r0, const FastqRecord& rec, const Reference& ref,
         std::string& sam, std::string& fq) {
  Result r = r0;
  check_contig(r, ref);
  if (!r.certified) {
    fastq_record(rec, fq);
    return r.reason;
  }
  sam_line(r, rec, ref, mapq_of(r, static_cast<int>(rec.seq.size())), MateInfo(), sam);
  return kOk;
}

void sam_line(const Result& r, const FastqRecord& rec, const Reference& ref, int mapq,
              const MateInfo& mi, std::string& sam) {
  int c = ref.contig_of(r.ref_pos);
  char buf[96];
  sam += rec.name;
  std::snprintf(buf, sizeof buf, "\t%d\t", mi.flag | (r.strand ? 16 : 0));
  sam += buf;
  sam += ref.names[c];
  std::snprintf(buf, sizeof buf, "\t%lld\t%d\t",
                static_cast<long long>(r.ref_pos - ref.offsets[c] + 1), mapq);
  sam += buf;
  for (int i = 0; i < r.n_cigar; ++i) {
    std::snprintf(buf, sizeof buf, "%u%c", r.cigar[i] >> 4, "MIDNS"[r.cigar[i] & 0xF]);
    sam += buf;
  }
  std::snprintf(buf, sizeof buf, "\t%s\t%lld\t%lld\t", mi.rnext.c_str(),
                static_cast<long long>(mi.pnext), static_cast<long long>(mi.tlen));
  sam += buf;
  if (r.strand) {
    for (size_t i = rec.seq.size(); i-- > 0;) sam += comp(rec.seq[i]);
    sam += '\t';
    sam.append(rec.qual.rbegin(), rec.qual.rend());
  } else {
    sam += rec.seq; sam += '\t'; sam += rec.qual;
  }
  static const char* const kTierNames[] = {"S0", "S1", "SR", "S2", "SL", "PR"};
  // NM: edits in the reported (affine, possibly clipped) alignment.
  // XE: the certified minimum end-to-end edit distance over the reference
  // (-1 for SL). SL: AS (clip penalties included) is the certified maximum
  // local score over the reference, and XF the threshold it exceeds.
  char tags[160];
  int n = std::snprintf(tags, sizeof tags, "\tNM:i:%d\tAS:i:%d\tXE:i:%d\tXT:Z:%s\tXR:i:%d\tXD:i:%d\tXB:i:%d",
                        r.nm, r.score, r.d1, kTierNames[r.tier], r.radius, r.d2, r.n_best);
  if (r.tier == kTierSL) n += std::snprintf(tags + n, sizeof tags - n, "\tXF:i:%d", r.floor);
  // PR: XP is the certified maximum pair score, XF its threshold, XS a bound
  // on every other proper pair.
  if (r.tier == kTierPR) std::snprintf(tags + n, sizeof tags - n, "\tXF:i:%d\tXS:i:%d", r.floor, r.sub_score);
  sam += tags;
  sam += '\n';
}

// Proper pair (as certified): one contig, opposite strands, and the reverse
// mate's first aligned base 0..max_dist after the forward mate's.
bool is_proper(const Result& a, const Result& b, const Reference& ref, int max_dist) {
  if (a.strand == b.strand || ref.contig_of(a.ref_pos) != ref.contig_of(b.ref_pos)) return false;
  const Result& f = a.strand ? b : a;
  const Result& r = a.strand ? a : b;
  return r.ref_pos >= f.ref_pos && r.ref_pos - f.ref_pos <= max_dist;
}

// BWA-MEM's pairing: q_pe from the pair score margin; a mate's MAPQ is raised
// towards it (by at most 40 over its single-end MAPQ).
int mapq_pair(const Result& r, const Result& mate, const Result& se, int L) {
  if (r.tier != kTierPR) return mapq_of(r, L);
  int q_se = se.certified == 1 ? mapq_of(se, L) : 0;
  int q_pe = 0;
  if (r.n_best == 1) {
    const int margin = r.score + mate.score - r.sub_score;
    q_pe = std::min(60, static_cast<int>(6.02 * margin / kMatch + 0.499));
  }
  if (q_pe > q_se) q_se = std::min(q_pe, q_se + 40);
  return q_se;
}

// Both mates as SAM lines if both are certified, else to the fallback.
// Writes each mate's final reason to reason[0..1]. With `single` (not null),
// a pair with one certified mate sends only the other mate, single-end, there.
void emit_pair(const Result* r0, const Result* se, const FastqRecord* rec, const Reference& ref,
               int max_dist, std::string& sam, std::string& fq1, std::string& fq2, std::string* single,
               int* reason) {
  Result r[2] = {r0[0], r0[1]};
  for (int x = 0; x < 2; ++x) {
    if (r[x].certified == 2) { r[x].certified = 0; r[x].reason = kNotFound; }
    check_contig(r[x], ref);
  }
  if (!r[0].certified && !r[1].certified) {
    fastq_record(rec[0], fq1);
    fastq_record(rec[1], fq2);
    for (int x = 0; x < 2; ++x) reason[x] = r[x].reason;
    return;
  }
  if ((!r[0].certified || !r[1].certified) && single) {
    // One mate keeps its single-end certificate; the other is mapped
    // single-end by the fallback, so both records are written unpaired.
    const int k = r[0].certified ? 0 : 1;
    sam_line(r[k], rec[k], ref, mapq_of(r[k], static_cast<int>(rec[k].seq.size())), MateInfo(), sam);
    fastq_record(rec[1 - k], *single);
    reason[k] = kOk;
    reason[1 - k] = r[1 - k].reason;
    return;
  }
  if (!r[0].certified || !r[1].certified) {
    // One mate keeps its single-end certificate. Both mates still go to the
    // fallback (so it can pair them), named "name~k" where k is the certified
    // mate, whose fallback record is dropped downstream
    // (bench/drop_certified_mates.awk). Mate fields are unknown here.
    const int k = r[0].certified ? 0 : 1;
    MateInfo mi;
    mi.flag = 1 | (k == 0 ? 64 : 128);
    sam_line(r[k], rec[k], ref, mapq_of(r[k], static_cast<int>(rec[k].seq.size())), mi, sam);
    FastqRecord m0 = rec[0], m1 = rec[1];
    m0.name += k == 0 ? "~1" : "~2";
    m1.name = m0.name;
    fastq_record(m0, fq1);
    fastq_record(m1, fq2);
    reason[k] = kOk;
    reason[1 - k] = r[1 - k].reason;
    return;
  }
  const bool proper = is_proper(r[0], r[1], ref, max_dist);
  const int c0 = ref.contig_of(r[0].ref_pos), c1 = ref.contig_of(r[1].ref_pos);
  for (int x = 0; x < 2; ++x) {
    const Result& a = r[x];
    const Result& b = r[1 - x];
    MateInfo mi;
    mi.flag = 1 | (proper ? 2 : 0) | (b.strand ? 32 : 0) | (x == 0 ? 64 : 128);
    const int cb = x == 0 ? c1 : c0;
    mi.rnext = c0 == c1 ? "=" : ref.names[cb];
    mi.pnext = b.ref_pos - static_cast<int64_t>(ref.offsets[cb]) + 1;
    if (c0 == c1) {
      const int64_t lo = std::min(a.ref_pos, b.ref_pos);
      const int64_t hi = std::max(a.ref_pos + ref_span(a), b.ref_pos + ref_span(b));
      const bool left = a.ref_pos < b.ref_pos || (a.ref_pos == b.ref_pos && x == 0);
      mi.tlen = left ? hi - lo : -(hi - lo);
    }
    sam_line(a, rec[x], ref, mapq_pair(a, b, se[x], static_cast<int>(rec[x].seq.size())), mi, sam);
    reason[x] = kOk;
  }
}

struct Args {
  std::vector<std::string> pos;
  std::map<std::string, std::string> opt;
  bool has(const std::string& k) const { return opt.count(k) > 0; }
  std::string get(const std::string& k, const std::string& d) const {
    auto it = opt.find(k);
    return it == opt.end() ? d : it->second;
  }
  int geti(const std::string& k, int d) const {
    return has(k) ? std::atoi(opt.at(k).c_str()) : d;
  }
};

Args parse(int argc, char** argv, int start) {
  static const char* flags[] = {"--gpu"};
  Args a;
  for (int i = start; i < argc; ++i) {
    std::string s = argv[i];
    if (s.size() > 1 && s[0] == '-') {
      bool flag = std::find(std::begin(flags), std::end(flags), s) != std::end(flags);
      if (flag) { a.opt[s] = "1"; continue; }
      if (i + 1 >= argc) throw std::runtime_error("missing value for " + s);
      a.opt[s] = argv[++i];
    } else {
      a.pos.push_back(s);
    }
  }
  return a;
}

int default_threads() {
  unsigned n = std::thread::hardware_concurrency();
  return n ? static_cast<int>(n) : 1;
}

int cmd_index(const Args& a) {
  if (a.pos.size() != 1 || !a.has("-o"))
    throw std::runtime_error("usage: certa index ref.fa -o ref.cidx [-q 22] [-s 8] [-t N]");
  auto t0 = Clock::now();
  Reference ref = Reference::load_fasta(a.pos[0]);
  std::fprintf(stderr, "[index] %zu contigs, %llu bases (with padding) loaded in %.1f s\n",
               ref.names.size(), static_cast<unsigned long long>(ref.seq.size()), secs(t0));
  auto t1 = Clock::now();
  Index ix = Index::build(ref, a.geti("-q", 22), a.geti("-s", 8), a.geti("-t", default_threads()));
  std::fprintf(stderr, "[index] q=%d s=%d dir_bits=%d, %llu entries, built in %.1f s\n",
               ix.q, ix.s, ix.dir_bits, static_cast<unsigned long long>(ix.pos.size()), secs(t1));
  auto t2 = Clock::now();
  save_index(a.get("-o", ""), ref, ix);
  std::fprintf(stderr, "[index] written to %s in %.1f s\n", a.get("-o", "").c_str(), secs(t2));
  return 0;
}

int cmd_map(const Args& a) {
  if (a.pos.size() != 2 && a.pos.size() != 3)
    throw std::runtime_error(
        "usage: certa map ref.cidx reads.fq[.gz] [mates.fq[.gz]] [-k 5] [--k1 2] [--budget 16] [--budget2 256] "
        "[--budget3 0] [--support3 2] [--local 0] [--mapq-limit 8] [-t N] [--gpu] [--device 0] [-o out.sam] "
        "[-u uncertified.fq] [-U uncertified_mates.fq] [--fallback-mates pair|single -S single_mates.fq] [--max-dist 850] [--stats s.json] [--batch N] [--io-threads N]");
  Params p;
  p.k = a.geti("-k", 5);
  p.budget = a.geti("--budget", a.geti("--cap", 16));  // --cap: the v0.1 name
  if (p.k < 0 || p.k > KMAX) throw std::runtime_error("-k must be in [0, 5]");
  if (p.budget < 1 || p.budget > GPU_BUDGET_MAX) throw std::runtime_error("--budget must be in [1, 256]");
  // Escalation: reads that fail because the budget limited how many parts
  // they could search are re-run with this budget (0 disables).
  p.mapq_limit = a.geti("--mapq-limit", 8);  // MAPQ second-best estimate (0 = off)
  // S2 is off by default: on synthetic truth most S2 placements with MAPQ >= 20
  // were wrong, and on HG002 S2 reads were enriched at false-positive calls.
  p.s2_limit = a.geti("--s2", 0);
  if (p.s2_limit < 0 || p.s2_limit > S2_MAX) throw std::runtime_error("--s2 must be in [0, 10]");
  g_s2_limit = p.s2_limit;
  // Pass 1 runs every read with the small budget and without S2 (uniform,
  // fast GPU warps). Pass 2 re-runs only the reads pass 1 could not certify,
  // as a compacted batch with --budget2 and S2, so the expensive work stays
  // among hard reads. Each re-run read gets exactly its single-pass result.
  Params p2 = p;
  p2.budget = a.geti("--budget2", 256);
  if (p2.budget != 0 && (p2.budget < p.budget || p2.budget > GPU_BUDGET_MAX))
    throw std::runtime_error("--budget2 must be 0 or in [budget, 256]");
  if (p2.budget == 0) p2.budget = p.budget;
  p.s2_limit = 0;
  // Pass 1 also uses a smaller radius (--k1): wide bands for R up to k are
  // only needed by the few reads pass 1 cannot certify.
  p.k = std::min(p2.k, a.geti("--k1", 2));
  if (p.k < 0) throw std::runtime_error("--k1 must be >= 0");
  const bool second_pass = p2.budget > p.budget || p2.s2_limit > 0 || p2.k > p.k;
  // Pass 3 (CPU, opt-in): larger budgets for reads whose parts are all
  // repetitive. It certifies ~0.9% more reads but costs more CPU than the
  // fallback mapper spends on them, so it is off by default; use
  // --budget3 4096 for maximum certification.
  Params p3 = p2;
  p3.budget = a.geti("--budget3", 0);
  if (p3.budget != 0 && (p3.budget < p2.budget || p3.budget > HOST_BUDGET_MAX))
    throw std::runtime_error("--budget3 must be 0 or in [budget2, 4096]");
  if (p3.budget == 0) p3.budget = p2.budget;
  // Pass 3 verifies only clusters hit by >= t distinct parts (radius |S| - t),
  // which skips the flood of single-part hits from repetitive parts.
  p3.min_support = a.geti("--support3", 2);
  if (p3.min_support < 1 || p3.min_support > PMAX) throw std::runtime_error("--support3 must be in [1, 16]");
  // Pass L (CPU, opt-in): certified local alignment (tier SL) for reads still
  // uncertified. --local t verifies chains hit by >= t parts; 2 is the
  // cheap setting (t = 1 costs ~10x more on HG002 for +0.16 points).
  Params pl = p2;
  pl.budget = std::max(p2.budget, p3.budget);
  pl.local = a.geti("--local", 0);
  pl.local_only = 1;
  if (pl.local < 0 || pl.local > PMAX) throw std::runtime_error("--local must be in [0, 16]");
  // Paired-end: a second FASTQ with the mates, in the same order.
  const bool paired = a.pos.size() == 3;
  PairParams pp;
  pp.se = p2;
  pp.se.budget = std::max(p2.budget, p3.budget);
  pp.se.min_support = 1;
  pp.max_dist = a.geti("--max-dist", 850);
  if (pp.max_dist < 0) throw std::runtime_error("--max-dist must be >= 0");
  const bool gpu = a.has("--gpu");
  const int threads = a.geti("-t", default_threads());
  // Threads that parse and encode the next batch while the current one maps.
  const int io_threads = a.geti("--io-threads", std::max(1, std::min(16, threads / 4)));
  const size_t batch_size = static_cast<size_t>(a.geti("--batch", gpu ? 1000000 : 200000));

  auto t_all = Clock::now();
  Reference ref;
  Index ix;
  load_index(a.pos[0], ref, ix);
  const double t_load = secs(t_all);
  std::fprintf(stderr, "[map] index %s: q=%d s=%d, %llu entries, loaded in %.1f s\n",
               a.pos[0].c_str(), ix.q, ix.s, static_cast<unsigned long long>(ix.pos.size()), t_load);
  std::fprintf(stderr,
               "[map] pass 1: k=%d budget=%d; pass 2: k=%d budget=%d S2<=%d; parts are >= %d bases\n",
               p.k, p.budget, p2.k, p2.budget, p2.s2_limit, min_read_length(ix.q, ix.s));
  const IndexView view = ix.view(ref);

  // Open input and outputs first, so a failure cannot leave the upload
  // thread running. Uncompressed regular files are memory-mapped; gzip and
  // pipes are streamed.
  std::unique_ptr<MappedFastq> mapped, mapped2;
  std::unique_ptr<LineReader> stream, stream2;
  if (MappedFastq::usable(a.pos[1])) mapped.reset(new MappedFastq(a.pos[1]));
  else stream.reset(new LineReader(a.pos[1]));
  if (paired) {
    if (MappedFastq::usable(a.pos[2])) mapped2.reset(new MappedFastq(a.pos[2]));
    else stream2.reset(new LineReader(a.pos[2]));
  }
  std::FILE* sam = a.has("-o") ? std::fopen(a.get("-o", "").c_str(), "wb") : stdout;
  if (!sam) throw std::runtime_error("cannot open SAM output");
  std::FILE* fq = a.has("-u") ? std::fopen(a.get("-u", "").c_str(), "wb") : nullptr;
  if (a.has("-u") && !fq) throw std::runtime_error("cannot open uncertified FASTQ output");
  std::FILE* fq2 = a.has("-U") ? std::fopen(a.get("-U", "").c_str(), "wb") : nullptr;
  if (a.has("-U") && !fq2) throw std::runtime_error("cannot open uncertified mates output");
  // --fallback-mates single: the uncertified mate of a half-certified pair goes
  // single-end to -S (default "pair": both mates to -u/-U, renamed name~k).
  const std::string fb_mode = a.get("--fallback-mates", "pair");
  if (fb_mode != "pair" && fb_mode != "single") throw std::runtime_error("--fallback-mates must be pair or single");
  std::FILE* fqs = nullptr;
  if (fb_mode == "single") {
    if (!a.has("-S")) throw std::runtime_error("--fallback-mates single needs -S single_mates.fq");
    fqs = std::fopen(a.get("-S", "").c_str(), "wb");
    if (!fqs) throw std::runtime_error("cannot open single-mate fallback output");
  }
  if (paired && fq && !fq2) throw std::runtime_error("paired input: -u needs -U for the mates");
  std::fprintf(sam, "@HD\tVN:1.6\tSO:unsorted\n");
  for (size_t c = 0; c < ref.names.size(); ++c)
    std::fprintf(sam, "@SQ\tSN:%s\tLN:%llu\n", ref.names[c].c_str(),
                 static_cast<unsigned long long>(ref.lengths[c]));
  std::fprintf(sam, "@PG\tID:certa\tPN:certa\tVN:0.5\tCL:certa map -k %d --k1 %d --budget %d --budget2 %d --s2 %d%s\n",
               p2.k, p.k, p.budget, p2.budget, p2.s2_limit, gpu ? " --gpu" : "");

  // The GPU upload runs in the background, overlapped with reading the
  // first batch; the first map call waits for it.
  std::unique_ptr<GpuMapper> gm;
  double t_upload = 0;
  std::exception_ptr upload_error;
  std::thread uploader;
  if (gpu) {
    uploader = std::thread([&] {
      try {
        auto t = Clock::now();
        gm.reset(new GpuMapper(ref, ix, a.geti("--device", 0)));
        t_upload = secs(t);
      } catch (...) {
        upload_error = std::current_exception();
      }
    });
  }
  auto wait_for_gpu = [&] {
    if (!uploader.joinable()) return;
    uploader.join();
    if (upload_error) std::rethrow_exception(upload_error);
    std::fprintf(stderr, "[map] GPU %s, index uploaded in %.1f s (overlapped with input)\n",
                 gm->device_name().c_str(), t_upload);
  };

  // Double buffering: a reader thread parses and encodes batch n+1 while
  // batch n is mapped and written, so input parsing overlaps with mapping.
  std::vector<FastqRecord> recs_buf[2];
  ReadBatch batch_buf[2];
  std::exception_ptr read_error;
  double t_io = 0, t_wait = 0, t_map = 0, t_out = 0;
  uint64_t escalated = 0;  // reads re-run with --budget2
  uint64_t third = 0;      // reads re-run on the CPU with --budget3
  uint64_t local = 0;      // reads re-run on the CPU for tier SL
  uint64_t pair_tried = 0; // pairs given to the pair certificate
  std::vector<Result> se_res;  // single-end results (paired mode, for MAPQ)
  auto load = [&](int slot) {
    try {
      auto t = Clock::now();
      if (!paired) {
        if (mapped) mapped->next_batch(recs_buf[slot], batch_size, io_threads);
        else read_fastq_batch(*stream, recs_buf[slot], batch_size, io_threads);
      } else {
        // Mates interleaved: read 2i is mate 1, read 2i + 1 its mate 2.
        std::vector<FastqRecord> m1, m2;
        const size_t half = std::max<size_t>(1, batch_size / 2);
        if (mapped) mapped->next_batch(m1, half, io_threads);
        else read_fastq_batch(*stream, m1, half, io_threads);
        if (mapped2) mapped2->next_batch(m2, half, io_threads);
        else read_fastq_batch(*stream2, m2, half, io_threads);
        if (m1.size() != m2.size()) throw std::runtime_error("mate files have different numbers of reads");
        std::vector<FastqRecord>& out = recs_buf[slot];
        out.resize(2 * m1.size());
        for (size_t i = 0; i < m1.size(); ++i) {
          if (m1[i].name != m2[i].name)
            throw std::runtime_error("mate names differ: " + m1[i].name + " / " + m2[i].name);
          out[2 * i] = std::move(m1[i]);
          out[2 * i + 1] = std::move(m2[i]);
        }
      }
      encode_batch(recs_buf[slot], batch_buf[slot], io_threads);
      t_io += secs(t);
    } catch (...) {
      recs_buf[slot].clear();
      read_error = std::current_exception();
    }
  };
  std::vector<Result> res;
  Stats total;
  load(0);
  wait_for_gpu();
  for (int cur = 0; !recs_buf[cur].empty(); cur ^= 1) {
    std::thread prefetch(load, cur ^ 1);
    // Join the reader even if mapping throws, so the error is reported
    // instead of std::terminate on a joinable thread.
    struct Joiner {
      std::thread& t;
      ~Joiner() {
        if (t.joinable()) t.join();
      }
    } joiner{prefetch};
    const std::vector<FastqRecord>& recs = recs_buf[cur];
    const ReadBatch& batch = batch_buf[cur];

    auto t = Clock::now();
    if (gm) gm->map(p, batch, res);
    else map_cpu(view, p, batch, res, threads);
    // Re-runs the selected reads as a compacted batch with other parameters.
    auto rerun = [&](const std::vector<size_t>& redo, const Params& pp, bool on_gpu) {
      if (redo.empty()) return;
      ReadBatch sub;
      sub.offs.resize(redo.size());
      sub.lens.resize(redo.size());
      sub.hashes.resize(redo.size());
      for (size_t x = 0; x < redo.size(); ++x) {
        const size_t i = redo[x];
        sub.offs[x] = sub.codes.size();
        sub.lens[x] = batch.lens[i];
        sub.hashes[x] = batch.hashes[i];
        sub.codes.insert(sub.codes.end(), batch.codes.begin() + batch.offs[i],
                         batch.codes.begin() + batch.offs[i] + batch.lens[i]);
      }
      std::vector<Result> res2;
      if (on_gpu) gm->map(pp, sub, res2);
      else map_cpu(view, pp, sub, res2, threads);
      for (size_t x = 0; x < redo.size(); ++x) res[redo[x]] = res2[x];
    };
    if (second_pass) {
      // Pass 2: every read pass 1 could not certify (except too-long/short
      // reads). Certified reads keep pass-1 results.
      std::vector<size_t> redo;
      for (size_t i = 0; i < res.size(); ++i)
        if (res[i].certified == 0 && res[i].reason != kBadLength) redo.push_back(i);
      rerun(redo, p2, gm != nullptr);
      escalated += redo.size();
    }
    if (p3.budget > p2.budget) {
      // Pass 3, always on the CPU (identical in CPU and GPU runs): reads that
      // failed because the budget limited how many parts they could search.
      std::vector<size_t> redo;
      for (size_t i = 0; i < res.size(); ++i) {
        const Result& r = res[i];
        const int want = r.parts < p3.k + p3.min_support ? r.parts : p3.k + p3.min_support;
        if (r.certified == 0 && (r.reason == kRadiusNegative || r.reason == kNotFound) && r.used < want)
          redo.push_back(i);
      }
      rerun(redo, p3, false);
      third += redo.size();
    }
    if (pl.local > 0) {
      // Pass L, on the CPU: certified local alignment (tier SL).
      std::vector<size_t> redo;
      for (size_t i = 0; i < res.size(); ++i)
        if (res[i].certified == 0 && (res[i].reason == kNotFound || res[i].reason == kClusterTooWide))
          redo.push_back(i);
      rerun(redo, pl, false);
      local += redo.size();
    }
    if (paired) {
      // Pair pass, on the CPU: pairs with one certified mate and an
      // uncertified or tied mate.
      se_res = res;
      std::vector<size_t> todo;
      for (size_t i = 0; 2 * i + 1 < res.size(); ++i) {
        const Result& x = res[2 * i];
        const Result& y = res[2 * i + 1];
        // Pairs with both mates uncertified are skipped: on HG002 under 1%
        // of them get a pair certificate, at a fifth of the pass's cost.
        const bool one = x.certified == 1 || y.certified == 1;
        if (one && (x.certified != 1 || y.certified != 1 || x.n_best > 1 || y.n_best > 1)) todo.push_back(i);
      }
      pair_tried += todo.size();
      std::atomic<size_t> next{0};
      auto worker = [&] {
        std::unique_ptr<Workspace> w1(new Workspace), w2(new Workspace);
        PairScratch sc;
        for (size_t z; (z = next.fetch_add(1)) < todo.size();) {
          const size_t i = todo[z];
          Result r1 = res[2 * i], r2 = res[2 * i + 1];
          const PairOut po = certify_pair(view, pp, batch.codes.data() + batch.offs[2 * i], batch.lens[2 * i],
                                          batch.codes.data() + batch.offs[2 * i + 1], batch.lens[2 * i + 1],
                                          *w1, *w2, sc, r1, r2);
          if (po.certified) { res[2 * i] = r1; res[2 * i + 1] = r2; }
        }
      };
      std::vector<std::thread> pw;
      for (int w = 1; w < threads; ++w) pw.emplace_back(worker);
      worker();
      for (auto& th : pw) th.join();
    }
    t_map += secs(t);

    // Format in parallel over contiguous slices, then write in input order.
    t = Clock::now();
    const int nt = std::max(1, std::min<int>(threads, static_cast<int>(recs.size() / 1024) + 1));
    std::vector<std::string> sam_out(nt), fq_out(nt), fq2_out(nt), fqs_out(nt);
    // Paired fallback: pair-aligned cut points, so the two mate files can be
    // written in alternating small chunks (a consumer reading both pipes in
    // lockstep, record by record, never waits on one while we block on the other).
    std::vector<std::vector<std::pair<size_t, size_t>>> cuts(nt);
    std::vector<Stats> st(nt);
    std::vector<std::thread> pool;
    for (int w = 0; w < nt; ++w) {
      pool.emplace_back([&, w] {
        size_t a0 = recs.size() * w / nt, b0 = recs.size() * (w + 1) / nt;
        if (paired) { a0 &= ~size_t(1); b0 = w == nt - 1 ? recs.size() : (b0 & ~size_t(1)); }
        int pair_reason[2] = {0, 0};
        for (size_t i = a0; i < b0; ++i) {
          int reason;
          if (!paired) {
            reason = emit(res[i], recs[i], ref, sam_out[w], fq_out[w]);
          } else {
            if (i % 2 == 0) {
              emit_pair(&res[i], &se_res[i], &recs[i], ref, pp.max_dist, sam_out[w], fq_out[w], fq2_out[w],
                        fqs ? &fqs_out[w] : nullptr,
                        pair_reason);
              Stats& s = st[w];
              ++s.pairs;
              if (pair_reason[0] == kOk && pair_reason[1] == kOk) {
                ++s.pairs_cert;
                s.pairs_proper += is_proper(res[i], res[i + 1], ref, pp.max_dist);
              }
              if (s.pairs % 64 == 0) cuts[w].push_back({fq_out[w].size(), fq2_out[w].size()});
            }
            reason = pair_reason[i % 2];
          }
          Stats& s = st[w];
          ++s.reads;
          s.cands += res[i].ncand;
          s.bases += recs[i].seq.size();
          ++s.reason[reason];
          if (reason != kOk) ++s.unc_by_used[res[i].used];
          if (reason == kOk && res[i].certified == 2) {
            ++s.s2;  // placed without certificate
          } else if (reason == kOk) {
            ++s.certified;
            const uint8_t tier = res[i].tier;
            ++(tier == kTierS0 ? s.s0 : tier == kTierS1 ? s.s1 : tier == kTierSR ? s.sr :
               tier == kTierSL ? s.sl : s.pr);
            s.ties += res[i].n_best > 1;
            if (res[i].radius >= 0 && res[i].radius <= KMAX + 1) ++s.by_radius[res[i].radius];
            if (res[i].d1 >= 0) ++s.by_d1[res[i].d1];
          }
        }
      });
    }
    for (auto& th : pool) th.join();
    for (int w = 0; w < nt; ++w) {
      std::fwrite(sam_out[w].data(), 1, sam_out[w].size(), sam);
      if (fq && fq2) {
        size_t a1 = 0, a2 = 0;
        cuts[w].push_back({fq_out[w].size(), fq2_out[w].size()});
        for (const auto& c : cuts[w]) {
          if (c.first == a1 && c.second == a2) continue;
          std::fwrite(fq_out[w].data() + a1, 1, c.first - a1, fq);
          std::fflush(fq);
          std::fwrite(fq2_out[w].data() + a2, 1, c.second - a2, fq2);
          std::fflush(fq2);
          a1 = c.first;
          a2 = c.second;
        }
      } else if (fq) {
        std::fwrite(fq_out[w].data(), 1, fq_out[w].size(), fq);
      }
      if (fqs) std::fwrite(fqs_out[w].data(), 1, fqs_out[w].size(), fqs);
      total.add(st[w]);
    }
    t_out += secs(t);

    t = Clock::now();
    prefetch.join();
    t_wait += secs(t);
  }
  if (read_error) std::rethrow_exception(read_error);
  if (sam != stdout) std::fclose(sam);
  if (fq) std::fclose(fq);
  if (fq2) std::fclose(fq2);
  if (fqs) std::fclose(fqs);
  const double t_total = secs(t_all);

  auto pct = [&](uint64_t x) { return total.reads ? 100.0 * x / total.reads : 0.0; };
  std::fprintf(stderr,
               "[map] %llu reads: certified %.2f%% (S0 exact %.2f%%, S1 <=R edits %.2f%%, "
               "SR repeat %.2f%%, SL local %.2f%%), ties %.2f%%\n",
               static_cast<unsigned long long>(total.reads), pct(total.certified),
               pct(total.s0), pct(total.s1), pct(total.sr), pct(total.sl), pct(total.ties));
  for (int i = 1; i < kNumReasons; ++i)
    std::fprintf(stderr, "[map]   uncertified/%s: %.2f%%\n", kReasonNames[i], pct(total.reason[i]));
  std::fprintf(stderr, "[map] placed without certificate (S2, <= %d edits): %.2f%%; to fallback: %.2f%%\n",
               p2.s2_limit, pct(total.s2), pct(total.reads - total.certified - total.s2));
  std::fprintf(stderr, "[map] second pass (budget %d, S2 <= %d): %.2f%% of reads; third pass (CPU, budget %d): %.2f%%\n",
               p2.budget, p2.s2_limit, pct(escalated), p3.budget, pct(third));
  std::fprintf(stderr, "[map] local pass (CPU, support %d): %.2f%% of reads\n", pl.local, pct(local));
  if (paired) {
    const double np = total.pairs ? 100.0 / total.pairs : 0.0;
    std::fprintf(stderr,
                 "[map] pairs %llu: both mates certified %.2f%% (proper %.2f%%); PR pair certificate "
                 "%.2f%% of reads; pair pass tried %.2f%% of pairs\n",
                 static_cast<unsigned long long>(total.pairs), np * total.pairs_cert, np * total.pairs_proper,
                 pct(total.pr), np * pair_tried);
  }
  std::fprintf(stderr, "[map] candidates enumerated per read (last pass): %.1f\n",
               total.reads ? double(total.cands) / total.reads : 0.0);
  std::fprintf(stderr, "[map] uncertified by parts searched |S|:");
  for (int i = 0; i <= PMAX; ++i) std::fprintf(stderr, " %d:%.2f%%", i, pct(total.unc_by_used[i]));
  std::fprintf(stderr, "\n");
  std::fprintf(stderr,
               "[map] time: load %.2f s, upload %.2f s, read+encode %.2f s (overlapped; "
               "waited %.2f s), map %.2f s (%.0f reads/s), output %.2f s, total %.2f s "
               "(%.0f reads/s end to end)\n",
               t_load, t_upload, t_io, t_wait, t_map, t_map > 0 ? total.reads / t_map : 0.0,
               t_out, t_total, t_total > 0 ? total.reads / t_total : 0.0);

  if (a.has("--stats")) {
    std::FILE* js = std::fopen(a.get("--stats", "").c_str(), "wb");
    if (!js) throw std::runtime_error("cannot open stats output");
    std::fprintf(js, "{\n  \"mode\": \"%s\",\n  \"device\": \"%s\",\n  \"threads\": %d,\n",
                 gpu ? "gpu" : "cpu", gm ? gm->device_name().c_str() : "", threads);
    std::fprintf(js, "  \"k\": %d, \"k1\": %d, \"budget\": %d, \"budget2\": %d, \"second_pass\": %llu, \"q\": %d, \"s\": %d,\n",
                 p2.k, p.k, p.budget, p2.budget, static_cast<unsigned long long>(escalated), ix.q, ix.s);
    std::fprintf(js, "  \"reads\": %llu, \"bases\": %llu, \"certified\": %llu, \"s0\": %llu, \"s1\": %llu, \"sr\": %llu, \"sl\": %llu, \"s2\": %llu, \"s2_limit\": %d, \"ties\": %llu,\n",
                 static_cast<unsigned long long>(total.reads), static_cast<unsigned long long>(total.bases),
                 static_cast<unsigned long long>(total.certified), static_cast<unsigned long long>(total.s0),
                 static_cast<unsigned long long>(total.s1), static_cast<unsigned long long>(total.sr),
                 static_cast<unsigned long long>(total.sl), static_cast<unsigned long long>(total.s2), p2.s2_limit,
                 static_cast<unsigned long long>(total.ties));
    std::fprintf(js, "  \"uncertified\": {");
    for (int i = 1; i < kNumReasons; ++i)
      std::fprintf(js, "%s\"%s\": %llu", i > 1 ? ", " : "", kReasonNames[i],
                   static_cast<unsigned long long>(total.reason[i]));
    std::fprintf(js, "},\n  \"certified_by_d1\": [");
    for (int i = 0; i <= KMAX; ++i)
      std::fprintf(js, "%s%llu", i ? ", " : "", static_cast<unsigned long long>(total.by_d1[i]));
    std::fprintf(js, "],\n  \"certified_by_radius\": [");
    for (int i = 0; i <= KMAX; ++i)
      std::fprintf(js, "%s%llu", i ? ", " : "", static_cast<unsigned long long>(total.by_radius[i]));
    std::fprintf(js, "],\n  \"seconds\": {\"load\": %.3f, \"upload\": %.3f, \"read_encode\": %.3f, "
                     "\"read_wait\": %.3f, \"map\": %.3f, \"output\": %.3f, \"total\": %.3f}\n}\n",
                 t_load, t_upload, t_io, t_wait, t_map, t_out, t_total);
    std::fclose(js);
  }
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc < 2) {
      std::fprintf(stderr,
                   "certa 0.7 - certified short-read fast path (prototype)\n"
                   "  certa index ref.fa -o ref.cidx [-q 22] [-s 8] [-t N]\n"
                   "  certa map ref.cidx reads.fq[.gz] [-k 5] [--k1 2] [--budget 16] [--budget2 256]\n"
                   "            [--budget3 0|4096] [--support3 2] [--local 0|2] [-t N] [--gpu] [--device 0]\n"
                   "            [-o out.sam] [-u uncertified.fq] [--stats s.json] [--batch N] [--io-threads N]\n"
                   "GPU support compiled in: %s\n",
                   GpuMapper::compiled_in() ? "yes" : "no");
      return 1;
    }
    std::string cmd = argv[1];
    Args a = parse(argc, argv, 2);
    if (cmd == "index") return cmd_index(a);
    if (cmd == "map") return cmd_map(a);
    throw std::runtime_error("unknown command: " + cmd);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "certa: error: %s\n", e.what());
    return 1;
  }
}
