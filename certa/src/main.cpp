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
  uint64_t reads = 0, bases = 0, certified = 0, s0 = 0, s1 = 0, sr = 0, s2 = 0, ties = 0;
  uint64_t reason[kNumReasons] = {};
  uint64_t by_radius[KMAX + 2] = {};  // certified reads per radius R
  uint64_t by_d1[KMAX + 1] = {};
  void add(const Stats& o) {
    reads += o.reads; bases += o.bases; certified += o.certified;
    s0 += o.s0; s1 += o.s1; sr += o.sr; s2 += o.s2; ties += o.ties;
    for (int i = 0; i < kNumReasons; ++i) reason[i] += o.reason[i];
    for (int i = 0; i < KMAX + 2; ++i) by_radius[i] += o.by_radius[i];
    for (int i = 0; i < KMAX + 1; ++i) by_d1[i] += o.by_d1[i];
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
  // candidate found (d2x), as bwa-mem uses its suboptimal hit; no second: 60.
  // bwa-mem compares affine scores, so use them when available.
  if (r.sub_score != kNoSub && r.score > 0) {
    if (r.sub_score >= r.score) return 0;
    const int q = static_cast<int>(30.0 * (1.0 - static_cast<double>(r.sub_score) / r.score) *
                                       std::log(static_cast<double>(L)) + 0.499);
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

// Appends a SAM line (certified) or a FASTQ record (fallback). Returns the
// final reason code after the host-side contig check.
int emit(const Result& r0, const FastqRecord& rec, const Reference& ref,
         std::string& sam, std::string& fq) {
  Result r = r0;
  if (r.certified) {
    int64_t span = 0;
    for (int i = 0; i < r.n_cigar; ++i) {
      uint32_t op = r.cigar[i] & 0xF;
      if (op == kOpM || op == kOpD) span += r.cigar[i] >> 4;
    }
    int c = ref.contig_of(r.ref_pos);
    if (c < 0 || ref.contig_of(r.ref_pos + span - 1) != c) {
      r.certified = 0;
      r.reason = kCrossContig;
    }
  }
  if (!r.certified) {
    fq += '@'; fq += rec.name; fq += '\n';
    fq += rec.seq; fq += "\n+\n";
    fq += rec.qual; fq += '\n';
    return r.reason;
  }
  int c = ref.contig_of(r.ref_pos);
  char buf[64];
  sam += rec.name;
  std::snprintf(buf, sizeof buf, "\t%d\t", r.strand ? 16 : 0);
  sam += buf;
  sam += ref.names[c];
  std::snprintf(buf, sizeof buf, "\t%lld\t%d\t",
                static_cast<long long>(r.ref_pos - ref.offsets[c] + 1),
                mapq_of(r, static_cast<int>(rec.seq.size())));
  sam += buf;
  for (int i = 0; i < r.n_cigar; ++i) {
    std::snprintf(buf, sizeof buf, "%u%c", r.cigar[i] >> 4, "MIDNS"[r.cigar[i] & 0xF]);
    sam += buf;
  }
  sam += "\t*\t0\t0\t";
  if (r.strand) {
    for (size_t i = rec.seq.size(); i-- > 0;) sam += comp(rec.seq[i]);
    sam += '\t';
    sam.append(rec.qual.rbegin(), rec.qual.rend());
  } else {
    sam += rec.seq; sam += '\t'; sam += rec.qual;
  }
  static const char* const kTierNames[] = {"S0", "S1", "SR", "S2"};
  // NM: edits in the reported (affine, possibly clipped) alignment.
  // XE: the certified minimum end-to-end edit distance over the reference.
  char tags[128];
  std::snprintf(tags, sizeof tags, "\tNM:i:%d\tAS:i:%d\tXE:i:%d\tXT:Z:%s\tXR:i:%d\tXD:i:%d\tXB:i:%d\n",
                r.nm, r.score, r.d1, kTierNames[r.tier], r.radius, r.d2, r.n_best);
  sam += tags;
  return kOk;
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
  if (a.pos.size() != 2)
    throw std::runtime_error(
        "usage: certa map ref.cidx reads.fq[.gz] [-k 4] [--budget 256] [-t N] [--gpu] "
        "[--device 0] [-o out.sam] [-u uncertified.fq] [--stats s.json] [--batch N] [--io-threads N]");
  Params p;
  p.k = a.geti("-k", 5);
  p.budget = a.geti("--budget", a.geti("--cap", 16));  // --cap: the v0.1 name
  if (p.k < 0 || p.k > KMAX) throw std::runtime_error("-k must be in [0, 5]");
  if (p.budget < 1 || p.budget > BUDGET_MAX) throw std::runtime_error("--budget must be in [1, 256]");
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
  if (p2.budget != 0 && (p2.budget < p.budget || p2.budget > BUDGET_MAX))
    throw std::runtime_error("--budget2 must be 0 or in [budget, 256]");
  if (p2.budget == 0) p2.budget = p.budget;
  p.s2_limit = 0;
  // Pass 1 also uses a smaller radius (--k1): wide bands for R up to k are
  // only needed by the few reads pass 1 cannot certify.
  p.k = std::min(p2.k, a.geti("--k1", 2));
  if (p.k < 0) throw std::runtime_error("--k1 must be >= 0");
  const bool second_pass = p2.budget > p.budget || p2.s2_limit > 0 || p2.k > p.k;
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
  std::unique_ptr<MappedFastq> mapped;
  std::unique_ptr<LineReader> stream;
  if (MappedFastq::usable(a.pos[1])) mapped.reset(new MappedFastq(a.pos[1]));
  else stream.reset(new LineReader(a.pos[1]));
  std::FILE* sam = a.has("-o") ? std::fopen(a.get("-o", "").c_str(), "wb") : stdout;
  if (!sam) throw std::runtime_error("cannot open SAM output");
  std::FILE* fq = a.has("-u") ? std::fopen(a.get("-u", "").c_str(), "wb") : nullptr;
  if (a.has("-u") && !fq) throw std::runtime_error("cannot open uncertified FASTQ output");
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
  auto load = [&](int slot) {
    try {
      auto t = Clock::now();
      if (mapped) mapped->next_batch(recs_buf[slot], batch_size, io_threads);
      else read_fastq_batch(*stream, recs_buf[slot], batch_size, io_threads);
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
    if (second_pass) {
      // Pass 2 on a compacted batch: every read pass 1 could not certify
      // (except too-long/short reads). Certified reads keep pass-1 results.
      std::vector<size_t> redo;
      for (size_t i = 0; i < res.size(); ++i) {
        const Result& r = res[i];
        if (r.certified == 0 && r.reason != kBadLength) redo.push_back(i);
      }
      if (!redo.empty()) {
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
        if (gm) gm->map(p2, sub, res2);
        else map_cpu(view, p2, sub, res2, threads);
        for (size_t x = 0; x < redo.size(); ++x) res[redo[x]] = res2[x];
        escalated += redo.size();
      }
    }
    t_map += secs(t);

    // Format in parallel over contiguous slices, then write in input order.
    t = Clock::now();
    const int nt = std::max(1, std::min<int>(threads, static_cast<int>(recs.size() / 1024) + 1));
    std::vector<std::string> sam_out(nt), fq_out(nt);
    std::vector<Stats> st(nt);
    std::vector<std::thread> pool;
    for (int w = 0; w < nt; ++w) {
      pool.emplace_back([&, w] {
        size_t a0 = recs.size() * w / nt, b0 = recs.size() * (w + 1) / nt;
        for (size_t i = a0; i < b0; ++i) {
          int reason = emit(res[i], recs[i], ref, sam_out[w], fq_out[w]);
          Stats& s = st[w];
          ++s.reads;
          s.bases += recs[i].seq.size();
          ++s.reason[reason];
          if (reason == kOk && res[i].certified == 2) {
            ++s.s2;  // placed without certificate
          } else if (reason == kOk) {
            ++s.certified;
            ++(res[i].tier == kTierS0 ? s.s0 : res[i].tier == kTierS1 ? s.s1 : s.sr);
            s.ties += res[i].n_best > 1;
            if (res[i].radius >= 0) ++s.by_radius[res[i].radius];
            ++s.by_d1[res[i].d1];
          }
        }
      });
    }
    for (auto& th : pool) th.join();
    for (int w = 0; w < nt; ++w) {
      std::fwrite(sam_out[w].data(), 1, sam_out[w].size(), sam);
      if (fq) std::fwrite(fq_out[w].data(), 1, fq_out[w].size(), fq);
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
  const double t_total = secs(t_all);

  auto pct = [&](uint64_t x) { return total.reads ? 100.0 * x / total.reads : 0.0; };
  std::fprintf(stderr,
               "[map] %llu reads: certified %.2f%% (S0 exact %.2f%%, S1 <=R edits %.2f%%, "
               "SR repeat %.2f%%), ties %.2f%%\n",
               static_cast<unsigned long long>(total.reads), pct(total.certified),
               pct(total.s0), pct(total.s1), pct(total.sr), pct(total.ties));
  for (int i = 1; i < kNumReasons; ++i)
    std::fprintf(stderr, "[map]   uncertified/%s: %.2f%%\n", kReasonNames[i], pct(total.reason[i]));
  std::fprintf(stderr, "[map] placed without certificate (S2, <= %d edits): %.2f%%; to fallback: %.2f%%\n",
               p2.s2_limit, pct(total.s2), pct(total.reads - total.certified - total.s2));
  std::fprintf(stderr, "[map] second pass (budget %d, S2 <= %d): %.2f%% of reads\n", p2.budget,
               p2.s2_limit, pct(escalated));
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
    std::fprintf(js, "  \"reads\": %llu, \"bases\": %llu, \"certified\": %llu, \"s0\": %llu, \"s1\": %llu, \"sr\": %llu, \"s2\": %llu, \"s2_limit\": %d, \"ties\": %llu,\n",
                 static_cast<unsigned long long>(total.reads), static_cast<unsigned long long>(total.bases),
                 static_cast<unsigned long long>(total.certified), static_cast<unsigned long long>(total.s0),
                 static_cast<unsigned long long>(total.s1), static_cast<unsigned long long>(total.sr),
                 static_cast<unsigned long long>(total.s2), p2.s2_limit,
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
                   "certa 0.1 - certified short-read fast path (prototype)\n"
                   "  certa index ref.fa -o ref.cidx [-q 22] [-s 8] [-t N]\n"
                   "  certa map ref.cidx reads.fq[.gz] [-k 4] [--budget 256] [-t N] [--gpu] [--device 0]\n"
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
