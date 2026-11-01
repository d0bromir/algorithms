// bio_probe: why `certa map --bio` certifies a read or sends it to the
// fallback, and what it costs. Runs bio_place (certa/bio.h) with certa map's
// --bio settings on single reads and prints one line per read:
//   bio_probe INDEX.cidx reads.fq [N [THREADS]] > probe.tsv
// Columns: name, L, parts P, the edit certificate (certified, tier, d1, d2,
// reason) and its alignment's --bio penalty (-1 if none or not end to end),
// |S| per strand, H, G, the horizon if every part without an N were in S
// (h_all), the hits of the rarest part left out of S (next, per strand),
// yh, ok, p1, p2, MAPQ, reason, wide, candidates, clusters, bands aligned,
// DP cells, then microseconds in the edit certificate alone and in
// bio_place (which runs the edit certificate too).
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "certa/bio.h"
#include "index.h"

using namespace certa;

namespace {

struct Read {
  std::string name, seq, qual;
};

std::vector<Read> read_fastq(const char* path, long n) {
  std::FILE* f = std::fopen(path, "r");
  if (!f) throw std::runtime_error(std::string("cannot open ") + path);
  std::vector<Read> out;
  std::vector<char> b(1 << 16);
  auto line = [&](std::string& s) {
    if (!std::fgets(b.data(), static_cast<int>(b.size()), f)) return false;
    s = b.data();
    while (!s.empty() && (s.back() == '\n' || s.back() == '\r')) s.pop_back();
    return true;
  };
  std::string h, plus;
  Read r;
  while ((n < 0 || static_cast<long>(out.size()) < n) && line(h) && line(r.seq) && line(plus) && line(r.qual)) {
    r.name = h.substr(1, h.find_first_of(" \t") == std::string::npos ? std::string::npos : h.find_first_of(" \t") - 1);
    out.push_back(r);
  }
  std::fclose(f);
  return out;
}

uint64_t name_hash(const std::string& s) {  // as certa map hashes read names
  uint64_t h = 1469598103934665603ULL;
  for (char c : s) { h ^= static_cast<uint8_t>(c); h *= 1099511628211ULL; }
  return h;
}

// The rarest part's hits among those not in the first m of strand st.
long next_count(const Workspace& ws, int st, int P, int m) {
  return m < P ? static_cast<long>(ws.count[st][ws.order[st][m]]) : -1;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 3 || argc > 5) {
    std::fprintf(stderr, "usage: bio_probe INDEX.cidx reads.fq [N [THREADS]] > probe.tsv\n");
    return 2;
  }
  const long N = argc > 3 ? std::atol(argv[3]) : -1;
  const int T = argc > 4 ? std::atoi(argv[4]) : 1;
  Reference ref;
  Index index;
  load_index(argv[1], ref, index);
  const IndexView ix = index.view(ref);
  const std::vector<Read> reads = read_fastq(argv[2], N);
  // certa map --bio's pass (map_config.cpp): radius KMAX, budget 4096.
  Params p{};
  p.k = KMAX;
  p.budget = HOST_BUDGET_MAX;
  p.min_support = 1;
  const BioModel model;
  std::vector<std::string> lines(reads.size());
  auto work = [&](int t) {
    std::unique_ptr<Workspace> ws(new Workspace);
    std::unique_ptr<BioScratch> sc(new BioScratch);
    std::vector<uint8_t> codes(LMAX);
    BioRead br;
    using Clock = std::chrono::steady_clock;
    for (size_t i = t; i < reads.size(); i += T) {
      const Read& rd = reads[i];
      const int L = static_cast<int>(rd.seq.size());
      if (L > LMAX || static_cast<int>(rd.qual.size()) != L) continue;
      for (int j = 0; j < L; ++j) codes[j] = encode_base(rd.seq[j]);
      const uint64_t h = name_hash(rd.name);
      Result e;
      auto t0 = Clock::now();
      process_read(ix, p, codes.data(), L, h, *ws, e);
      auto t1 = Clock::now();
      int p_edit = -1;
      if (e.certified == 1 && e.n_cigar > 0 && e.ref_pos >= 0) {
        bio_read(model, codes.data(), rd.qual.c_str(), L, br);
        p_edit = bio_rescore(ix, model, br, e);
      }
      Result r;
      BioOut o;
      auto t2 = Clock::now();
      bio_place(ix, model, p, codes.data(), rd.qual.c_str(), L, h, *ws, *sc, r, o);
      auto t3 = Clock::now();
      // The horizon with every part that has no N.
      const int P = e.parts;
      int h_all = -1;
      long nx[2] = {-1, -1};
      if (P > 0) {
        bio_read(model, codes.data(), rd.qual.c_str(), L, br);
        int hc = kBioInf;
        for (int st = 0; st < 2; ++st) {
          uint32_t req = 0;
          for (int j = 0; j < P; ++j)
            if (ws->count[st][j] < static_cast<uint32_t>(BIG)) req |= 1u << j;
          hc = std::min(hc, bio_cover(model, br, st, P, req));
          nx[st] = next_count(*ws, st, P, o.n_s[st]);
        }
        h_all = std::min(hc, model.go_min() + model.params().gap_cap * model.ge());
      }
      char b[1024];
      std::snprintf(b, sizeof b,
                    "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%ld\t%ld\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d"
                    "\t%lld\t%lld\t%lld\n",
                    rd.name.c_str(), L, P, e.certified, e.tier, e.d1, e.d2, e.reason, p_edit, o.n_s[0], o.n_s[1],
                    o.horizon, o.reach, h_all, nx[0], nx[1], o.yh, o.ok ? 1 : 0, o.p1, o.p2, o.mapq, r.reason,
                    o.wide ? 1 : 0, o.n_cand, o.n_clusters, o.n_bands, static_cast<long long>(o.cells),
                    static_cast<long long>(std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count()),
                    static_cast<long long>(std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count()));
      lines[i] = b;
    }
  };
  std::vector<std::thread> th;
  for (int t = 0; t < T; ++t) th.emplace_back(work, t);
  for (auto& x : th) x.join();
  std::printf("name\tL\tP\tecert\tetier\td1\td2\tereason\tp_edit\tns0\tns1\tH\tG\th_all\tnext0\tnext1\tyh\tok\tp1\tp2\tmapq"
              "\treason\twide\tncand\tncl\tnbands\tcells\tt_edit\tt_bio\n");
  for (const std::string& s : lines) std::fputs(s.c_str(), stdout);
  return 0;
}
