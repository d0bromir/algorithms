// certa: prototype of the certified short-read fast path.
//
//   certa index ref.fa -o ref.cidx [-q 22] [-s 8] [-t threads]
//   certa map ref.cidx reads.fq[.gz] [mates.fq[.gz]] [options]   (cmd_map.cpp)
//   certa repeat ref.cidx reads.fq [--tau 16] [--tsv out.tsv]
//   certa check-index ref.cidx ref.fa
//   certa check ref.cidx aln.sam|- [--sample N] [--failures f.txt]
//
// Run without arguments for the options of each.
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "certa/repeat.h"
#include "check.h"
#include "cli.h"
#include "index.h"
#include "mapper.h"
#include "seqio.h"

namespace certa {

// Tier R (repeats) on a FASTQ of reads, e.g. those left uncertified by
// `certa map -u`. Prints a summary; --tsv writes one line per read.
int cmd_repeat(const Args& a) {
  if (a.pos.size() != 2)
    throw std::runtime_error(
        "usage: certa repeat ref.cidx reads.fq[.gz] [-t N] [--tau 16] [--cmax 4096] [--samples 8] [--chains 0] "
        "[--fallback 0] [--tsv out.tsv]");
  RepeatParams rp;
  rp.tau = a.geti("--tau", 16);
  rp.cmax = a.geti("--cmax", 4096);
  rp.samples = a.geti("--samples", 8);
  rp.lb_chains = a.geti("--chains", 0);
  rp.fallback_samples = a.geti("--fallback", 0);
  if (rp.tau < 0 || rp.cmax < 1 || rp.samples < 1 || rp.lb_chains < 0 || rp.fallback_samples < 0)
    throw std::runtime_error("bad tier R parameters");
  const int threads = a.geti("-t", default_threads());
  Reference ref;
  Index ix;
  load_index(a.pos[0], ref, ix);
  if (ix.s != 1) throw std::runtime_error("tier R needs an index with s = 1");
  const IndexView view = ix.view(ref);
  LineReader in(a.pos[1]);
  std::vector<FastqRecord> recs, chunk;
  while (read_fastq_batch(in, chunk, 1 << 20, threads) > 0)
    for (auto& r : chunk) recs.push_back(std::move(r));
  ReadBatch batch;
  encode_batch(recs, batch, threads);
  std::vector<RepeatOut> out(recs.size());
  auto t0 = Clock::now();
  std::atomic<size_t> next{0};
  auto worker = [&] {
    std::unique_ptr<Workspace> ws(new Workspace);
    RepeatScratch sc;
    for (size_t i; (i = next.fetch_add(1)) < recs.size();)
      out[i] = repeat_tier(view, rp, batch.codes.data() + batch.offs[i], batch.lens[i], *ws, sc);
  };
  std::vector<std::thread> pool;
  for (int t = 1; t < threads; ++t) pool.emplace_back(worker);
  worker();
  for (auto& th : pool) th.join();
  const double t_run = secs(t0);
  // Summary. "Ambiguous" means two loci provably score >= lo while OPT <= hi,
  // so OPT - SUB <= hi - lo; at <= 5 that keeps BWA-MEM's MAPQ below 20.
  uint64_t n = recs.size(), claim = 0, exact = 0, tie = 0, amb = 0, nolo = 0, w5 = 0, w25 = 0;
  for (const auto& o : out) {
    if (!o.claim) continue;
    ++claim;
    if (o.lo == kNoSub) { ++nolo; continue; }
    if (o.exact) { ++exact; tie += o.copies >= 2; continue; }
    const int w = o.hi - o.lo;
    w5 += w <= 5;
    w25 += w <= 25;
    amb += o.copies >= 2 && w <= 5;
  }
  auto pct = [&](uint64_t x) { return n ? 100.0 * x / n : 0.0; };
  std::fprintf(stderr,
               "[repeat] %llu reads (tau %d, cmax %d): exact %.2f%% (proven tie %.2f%%); not exact: "
               "interval <= 5 %.2f%%, <= 25 %.2f%%, proven ambiguous %.2f%%, no alignment %.2f%%\n",
               static_cast<unsigned long long>(n), rp.tau, rp.cmax, pct(exact), pct(tie), pct(w5), pct(w25),
               pct(amb), pct(nolo));
  std::fprintf(stderr, "[repeat] time %.2f s on %d threads (%.3f ms CPU-equivalent per read)\n", t_run, threads,
               n ? t_run * threads * 1000.0 / n : 0.0);
  if (a.has("--tsv")) {
    std::FILE* f = std::fopen(a.get("--tsv", "").c_str(), "wb");
    if (!f) throw std::runtime_error("cannot open --tsv output");
    std::fprintf(f, "read\texact\tlo\thi\tcopies\tstrand\tcontig\tpos\tparts\tabsent\toccurrences\n");
    for (size_t i = 0; i < recs.size(); ++i) {
      const RepeatOut& o = out[i];
      std::string contig = "*";
      long long pos = 0;
      if (o.lo != kNoSub) {
        const int c = ref.contig_of(o.aln.ref_pos);
        if (c >= 0) { contig = ref.names[c]; pos = o.aln.ref_pos - static_cast<long long>(ref.offsets[c]) + 1; }
      }
      std::fprintf(f, "%s\t%d\t%d\t%d\t%d\t%d\t%s\t%lld\t%d\t%d\t%d\n", recs[i].name.c_str(), o.exact ? 1 : 0,
                   o.lo == kNoSub ? -1 : o.lo, o.hi, o.copies, o.lo == kNoSub ? -1 : o.aln.strand, contig.c_str(),
                   pos, o.nparts, o.nabsent, o.nocc);
    }
    std::fclose(f);
  }
  return 0;
}

// Validates an index file against its reference and the FASTA, once per
// reference build (the trust base of certa check, Theorem V).
int cmd_check_index(const Args& a) {
  if (a.pos.size() != 2) throw std::runtime_error("usage: certa check-index ref.cidx ref.fa [-t N]");
  auto t0 = Clock::now();
  Reference ref;
  Index ix;
  load_index(a.pos[0], ref, ix);
  std::string why;
  const bool ok = check::validate_index(ref, ix, a.pos[1], a.geti("-t", default_threads()), &why);
  std::fprintf(stderr, "[check-index] %s: q=%d s=%d, %llu entries, %zu contigs: %s (%.1f s)\n", a.pos[0].c_str(),
               ix.q, ix.s, static_cast<unsigned long long>(ix.pos.size()), ref.names.size(),
               ok ? "valid" : ("INVALID: " + why).c_str(), secs(t0));
  return ok ? 0 : 1;
}

// Checks the claims of certa map --intervals --witness records (check.h):
// every record, or a uniform sample of --sample N.
int cmd_check(const Args& a) {
  if (a.pos.size() != 2)
    throw std::runtime_error("usage: certa check ref.cidx aln.sam|- [-t N] [--sample N] [--seed 1] [--failures f.txt]");
  const int threads = a.geti("-t", default_threads());
  const int64_t sample = a.geti("--sample", 0);
  Reference ref;
  Index ix;
  load_index(a.pos[0], ref, ix);
  const check::Checker ck(ref, ix);
  std::ifstream file;
  std::istream* in = &std::cin;
  if (a.pos[1] != "-") {
    file.open(a.pos[1]);
    if (!file) throw std::runtime_error("cannot open " + a.pos[1]);
    in = &file;
  }
  std::FILE* fails = a.has("--failures") ? std::fopen(a.get("--failures", "").c_str(), "wb") : nullptr;
  uint64_t verdict[check::kNumVerdicts] = {}, records = 0, xu = 0, xh = 0, pair_xq = 0, pair_need = 0,
           xs = 0, occ = 0, windows = 0, cells = 0, lines = 0, groups = 0;
  double t_check = 0;
  const std::clock_t c0 = std::clock();
  // A group: consecutive records with one QNAME (a pair's two records are
  // checked together: a PR mate's XQ rests on both). Kept whole by batches
  // and by --sample, which samples groups.
  struct Group {
    std::vector<std::string> line;
    std::vector<uint64_t> number;  // record number (1-based, headers excluded)
  };
  auto run = [&](std::vector<Group>& batch) {
    std::vector<std::vector<check::Verdict>> v(batch.size());
    std::vector<std::vector<check::Work>> w(batch.size());
    auto t = Clock::now();
    std::atomic<size_t> next{0};
    auto worker = [&] {
      for (size_t g; (g = next.fetch_add(1)) < batch.size();) {
        const size_t n = batch[g].line.size();
        v[g].assign(n, check::kMalformed);
        w[g].assign(n, check::Work());
        std::vector<check::Claim> c(n);
        std::vector<char> ok(n);
        for (size_t i = 0; i < n; ++i) ok[i] = ck.parse_sam(batch[g].line[i], &c[i]);
        if (n == 2 && ok[0] && ok[1]) {
          ck.check_pair(c[0], c[1], &w[g][0], &w[g][1], &v[g][0], &v[g][1]);
        } else {
          for (size_t i = 0; i < n; ++i)
            if (ok[i]) v[g][i] = ck.check(c[i], &w[g][i]);
        }
      }
    };
    std::vector<std::thread> pool;
    for (int k = 1; k < threads; ++k) pool.emplace_back(worker);
    worker();
    for (auto& th : pool) th.join();
    t_check += secs(t);
    for (size_t g = 0; g < batch.size(); ++g)
      for (size_t i = 0; i < batch[g].line.size(); ++i) {
        const check::Work& wi = w[g][i];
        ++records;
        ++verdict[v[g][i]];
        xu += wi.xu_checked;
        xh += wi.xh_checked;
        pair_xq += wi.pair_xq;
        pair_need += wi.needs_pair;
        xs += wi.xs_checked;
        occ += wi.occurrences;
        windows += wi.windows;
        cells += wi.cells;
        if (v[g][i] != check::kAccepted && fails)
          std::fprintf(fails, "%s\t%llu\t%s\n", check::kVerdictNames[v[g][i]],
                       static_cast<unsigned long long>(batch[g].number[i]), batch[g].line[i].c_str());
      }
    batch.clear();
  };
  std::vector<Group> batch;
  Group cur;
  std::string cur_name;
  std::mt19937_64 rng(static_cast<uint64_t>(a.geti("--seed", 1)));
  auto flush = [&]() {
    if (cur.line.empty()) return;
    ++groups;
    if (sample > 0) {  // reservoir sample of groups
      if (static_cast<int64_t>(batch.size()) < sample) {
        batch.push_back(std::move(cur));
      } else {
        const uint64_t j = rng() % groups;
        if (j < static_cast<uint64_t>(sample)) batch[j] = std::move(cur);
      }
    } else {
      batch.push_back(std::move(cur));
      if (batch.size() == (1u << 17)) run(batch);
    }
    cur = Group();
  };
  std::string line;
  while (std::getline(*in, line)) {
    if (line.empty() || line[0] == '@') continue;
    ++lines;
    const std::string name = line.substr(0, line.find('\t'));
    if (name != cur_name) {
      flush();
      cur_name = name;
    }
    cur.line.push_back(line);
    cur.number.push_back(lines);
  }
  flush();
  if (!batch.empty()) run(batch);
  if (fails) std::fclose(fails);
  const double cpu = static_cast<double>(std::clock() - c0) / CLOCKS_PER_SEC;
  const double n = records ? static_cast<double>(records) : 1.0;
  std::fprintf(stderr, "[check] %llu of %llu records checked%s: accepted %.4f%%\n",
               static_cast<unsigned long long>(records), static_cast<unsigned long long>(lines),
               sample > 0 ? " (uniform sample)" : "", 100.0 * verdict[check::kAccepted] / n);
  for (int i = 1; i < check::kNumVerdicts; ++i)
    std::fprintf(stderr, "[check]   %s: %llu\n", check::kVerdictNames[i], static_cast<unsigned long long>(verdict[i]));
  std::fprintf(stderr,
               "[check] claims proven from witnesses: XU < L %llu, XH %llu; PR mates with XQ > 0 %llu, of which "
               "%llu rest on the pair certificate and %llu have XS proven\n",
               static_cast<unsigned long long>(xu), static_cast<unsigned long long>(xh),
               static_cast<unsigned long long>(pair_xq), static_cast<unsigned long long>(pair_need),
               static_cast<unsigned long long>(xs));
  std::fprintf(stderr,
               "[check] per record: %.1f probe occurrences, %.2f DP windows, %.0f DP cells; %.1f us CPU (%.1f s wall "
               "on %d threads)\n",
               occ / n, windows / n, cells / n, 1e6 * cpu / n, t_check, threads);
  if (sample > 0 && verdict[check::kAccepted] == records)
    std::fprintf(stderr,
                 "[check] audit: if more than %.3g%% of all records were wrong, a sample of %llu would have shown one "
                 "with 95%% probability\n",
                 100.0 * (1.0 - std::pow(0.05, 1.0 / n)), static_cast<unsigned long long>(records));
  return verdict[check::kAccepted] == records ? 0 : 2;
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

}  // namespace certa

int main(int argc, char** argv) {
  using namespace certa;
  try {
    if (argc < 2) {
      std::fprintf(stderr,
                   "certa 0.7 - certified short-read fast path (prototype)\n"
                   "  certa index ref.fa -o ref.cidx [-q 22] [-s 8] [-t N]\n"
                   "  certa map ref.cidx reads.fq[.gz] [-k 5] [--k1 2] [--budget 16] [--budget2 256]\n"
                   "            [--budget3 0|4096] [--support3 2] [--local 0|2] [--full-cert] [--adaptive-parts 0|1|2]\n"
                   "            [--index2 second.cidx] [-t N] [--gpu] [--device 0]\n"
                   "            [-o out.sam] [-u uncertified.fq] [--stats s.json] [--batch N] [--io-threads N]\n"
                   "            [--intervals [--witness] [--iv-mapq point|cons]]\n"
                   "  certa repeat ref.cidx reads.fq [--tau 16] [--tsv out.tsv]\n"
                   "  certa check-index ref.cidx ref.fa       (once per reference build)\n"
                   "  certa check ref.cidx aln.sam|- [--sample N] [--failures f.txt]\n"
                   "GPU support compiled in: %s\n",
                   GpuMapper::compiled_in() ? "yes" : "no");
      return 1;
    }
    std::string cmd = argv[1];
    Args a = Args::parse(argc, argv, 2);
    if (cmd == "index") return cmd_index(a);
    if (cmd == "map") return cmd_map(a);
    if (cmd == "repeat") return cmd_repeat(a);
    if (cmd == "check") return cmd_check(a);
    if (cmd == "check-index") return cmd_check_index(a);
    throw std::runtime_error("unknown command: " + cmd);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "certa: error: %s\n", e.what());
    return 1;
  }
}
