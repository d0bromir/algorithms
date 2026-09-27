// certa: prototype of the certified short-read fast path.
//
//   certa index ref.fa -o ref.cidx [-q 22] [-s 8] [-t threads]
//   certa map ref.cidx reads.fq[.gz] [-k 2] [--cap 32] [-t threads]
//             [--gpu [--device 0]] [-o out.sam] [-u uncertified.fq]
//             [--stats stats.json] [--batch N]
//
// Certified reads are written to SAM. All other reads are written unchanged
// to the uncertified FASTQ, for a full aligner (minibwa, BWA-MEM2, ...).
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
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
  uint64_t reads = 0, bases = 0, certified = 0, s0 = 0, s1 = 0, ties = 0;
  uint64_t reason[kNumReasons] = {};
  uint64_t by_radius[KMAX + 2] = {};  // certified reads per radius R
  uint64_t by_d1[KMAX + 1] = {};
  void add(const Stats& o) {
    reads += o.reads; bases += o.bases; certified += o.certified;
    s0 += o.s0; s1 += o.s1; ties += o.ties;
    for (int i = 0; i < kNumReasons; ++i) reason[i] += o.reason[i];
    for (int i = 0; i < KMAX + 2; ++i) by_radius[i] += o.by_radius[i];
    for (int i = 0; i < KMAX + 1; ++i) by_d1[i] += o.by_d1[i];
  }
};

// Provisional MAPQ from the certified distance gap (to be calibrated).
int mapq_of(const Result& r) {
  if (r.n_best > 1) return 0;
  int second = r.d2 >= 0 ? r.d2 : r.radius + 1;  // lower bound when unseen
  return std::min(60, 20 * (second - r.d1));
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
                static_cast<long long>(r.ref_pos - ref.offsets[c] + 1), mapq_of(r));
  sam += buf;
  for (int i = 0; i < r.n_cigar; ++i) {
    std::snprintf(buf, sizeof buf, "%u%c", r.cigar[i] >> 4, "MID"[r.cigar[i] & 0xF]);
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
  std::snprintf(buf, sizeof buf, "\tNM:i:%d\tXT:Z:S%d\tXR:i:%d\tXD:i:%d\tXB:i:%d\n",
                r.d1, r.tier, r.radius, r.d2, r.n_best);
  sam += buf;
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
               ix.q, ix.s, ix.dir_bits, static_cast<unsigned long long>(ix.keys.size()), secs(t1));
  auto t2 = Clock::now();
  save_index(a.get("-o", ""), ref, ix);
  std::fprintf(stderr, "[index] written to %s in %.1f s\n", a.get("-o", "").c_str(), secs(t2));
  return 0;
}

int cmd_map(const Args& a) {
  if (a.pos.size() != 2)
    throw std::runtime_error(
        "usage: certa map ref.cidx reads.fq[.gz] [-k 2] [--cap 32] [-t N] [--gpu] "
        "[--device 0] [-o out.sam] [-u uncertified.fq] [--stats s.json] [--batch N]");
  Params p;
  p.k = a.geti("-k", 2);
  p.cap = a.geti("--cap", 32);
  if (p.k < 0 || p.k > KMAX) throw std::runtime_error("-k must be in [0, 5]");
  if (p.cap < 1 || p.cap > CAP_MAX) throw std::runtime_error("--cap must be in [1, 32]");
  const bool gpu = a.has("--gpu");
  const int threads = a.geti("-t", default_threads());
  const size_t batch_size = static_cast<size_t>(a.geti("--batch", gpu ? 1000000 : 200000));

  auto t_all = Clock::now();
  Reference ref;
  Index ix;
  load_index(a.pos[0], ref, ix);
  const double t_load = secs(t_all);
  std::fprintf(stderr, "[map] index %s: q=%d s=%d, %llu entries, loaded in %.1f s\n",
               a.pos[0].c_str(), ix.q, ix.s, static_cast<unsigned long long>(ix.keys.size()), t_load);
  std::fprintf(stderr, "[map] k=%d cap=%d: reads shorter than %d bases cannot be certified\n",
               p.k, p.cap, min_read_length(p.k, ix.q, ix.s));
  const IndexView view = ix.view(ref);

  std::unique_ptr<GpuMapper> gm;
  double t_upload = 0;
  if (gpu) {
    auto t = Clock::now();
    gm.reset(new GpuMapper(ref, ix, a.geti("--device", 0)));
    t_upload = secs(t);
    std::fprintf(stderr, "[map] GPU %s, index uploaded in %.1f s\n", gm->device_name().c_str(), t_upload);
  }

  std::FILE* sam = a.has("-o") ? std::fopen(a.get("-o", "").c_str(), "wb") : stdout;
  if (!sam) throw std::runtime_error("cannot open SAM output");
  std::FILE* fq = a.has("-u") ? std::fopen(a.get("-u", "").c_str(), "wb") : nullptr;
  if (a.has("-u") && !fq) throw std::runtime_error("cannot open uncertified FASTQ output");
  std::fprintf(sam, "@HD\tVN:1.6\tSO:unsorted\n");
  for (size_t c = 0; c < ref.names.size(); ++c)
    std::fprintf(sam, "@SQ\tSN:%s\tLN:%llu\n", ref.names[c].c_str(),
                 static_cast<unsigned long long>(ref.lengths[c]));
  std::fprintf(sam, "@PG\tID:certa\tPN:certa\tVN:0.1\tCL:certa map -k %d --cap %d%s\n",
               p.k, p.cap, gpu ? " --gpu" : "");

  LineReader in(a.pos[1]);
  std::vector<FastqRecord> recs;
  ReadBatch batch;
  std::vector<Result> res;
  Stats total;
  double t_io = 0, t_map = 0, t_out = 0;
  for (;;) {
    auto t = Clock::now();
    if (read_fastq_batch(in, recs, batch_size) == 0) break;
    encode_batch(recs, batch);
    t_io += secs(t);

    t = Clock::now();
    if (gm) gm->map(p, batch, res);
    else map_cpu(view, p, batch, res, threads);
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
          if (reason == kOk) {
            ++s.certified;
            ++(res[i].tier == 0 ? s.s0 : s.s1);
            s.ties += res[i].n_best > 1;
            ++s.by_radius[res[i].radius];
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
  }
  if (sam != stdout) std::fclose(sam);
  if (fq) std::fclose(fq);
  const double t_total = secs(t_all);

  auto pct = [&](uint64_t x) { return total.reads ? 100.0 * x / total.reads : 0.0; };
  std::fprintf(stderr,
               "[map] %llu reads: certified %.2f%% (S0 exact %.2f%%, S1 <=k %.2f%%), "
               "ties %.2f%%\n",
               static_cast<unsigned long long>(total.reads), pct(total.certified),
               pct(total.s0), pct(total.s1), pct(total.ties));
  for (int i = 1; i < kNumReasons; ++i)
    std::fprintf(stderr, "[map]   uncertified/%s: %.2f%%\n", kReasonNames[i], pct(total.reason[i]));
  std::fprintf(stderr,
               "[map] time: load %.2f s, upload %.2f s, read+encode %.2f s, map %.2f s "
               "(%.0f reads/s), output %.2f s, total %.2f s\n",
               t_load, t_upload, t_io, t_map, t_map > 0 ? total.reads / t_map : 0.0, t_out, t_total);

  if (a.has("--stats")) {
    std::FILE* js = std::fopen(a.get("--stats", "").c_str(), "wb");
    if (!js) throw std::runtime_error("cannot open stats output");
    std::fprintf(js, "{\n  \"mode\": \"%s\",\n  \"device\": \"%s\",\n  \"threads\": %d,\n",
                 gpu ? "gpu" : "cpu", gm ? gm->device_name().c_str() : "", threads);
    std::fprintf(js, "  \"k\": %d, \"cap\": %d, \"q\": %d, \"s\": %d,\n", p.k, p.cap, ix.q, ix.s);
    std::fprintf(js, "  \"reads\": %llu, \"bases\": %llu, \"certified\": %llu, \"s0\": %llu, \"s1\": %llu, \"ties\": %llu,\n",
                 static_cast<unsigned long long>(total.reads), static_cast<unsigned long long>(total.bases),
                 static_cast<unsigned long long>(total.certified), static_cast<unsigned long long>(total.s0),
                 static_cast<unsigned long long>(total.s1), static_cast<unsigned long long>(total.ties));
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
                     "\"map\": %.3f, \"output\": %.3f, \"total\": %.3f}\n}\n",
                 t_load, t_upload, t_io, t_map, t_out, t_total);
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
                   "  certa map ref.cidx reads.fq[.gz] [-k 2] [--cap 32] [-t N] [--gpu] [--device 0]\n"
                   "            [-o out.sam] [-u uncertified.fq] [--stats s.json] [--batch N]\n"
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
