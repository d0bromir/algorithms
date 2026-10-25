// payoff: what part-selection techniques taken from the slow mappers would
// certify on the reads CERTA leaves uncertified (bench/headroom.py measured
// what they could reach). For each read it certifies the minimum end-to-end
// edit distance over the reference (both strands) with a pigeonhole
// certificate, its parts chosen one of four ways:
//   fixed  CERTA's parts (L / q equal parts, each looked up by its first q-mer);
//   rare   disjoint q-mers anywhere in the read, fewest hits in all (DP);
//   long   disjoint exact parts of any length q..q + 2E (longer parts are
//          rarer: repeats), counted from a q-mer inside them;
//   err1   disjoint q-mers, each searched with up to one edit;
//   err1r  the same, but only the rarest disjoint q-mers are expanded (cheap).
// The parts whose hits fit the budget (B per strand) form S. Every alignment
// with at most R = |S| - 1 edits (err1, err1r: R = 2 |S| - 1) leaves one part
// error-free (err1, err1r: with at most one edit), so its path lies within R + 1
// diagonals of a hit of that part, where a banded edit distance finds it. The
// read is certified when the best distance found on either strand is <= R,
// R being the smaller of the two strands' radii (capped by --k).
//   payoff INDEX.cidx READS.fq [-t N] [--budget 4096] [--k 15] [--cap 20000]
//          [--ext 40] [--modes fixed,rare,long,err1,err1r] [--truth tool.sam]
//          [--out per_read.tsv] [--of N]
// --truth: a full-sensitivity mapper's SAM (Yara) to check every certified
// distance against (its NM is the fewest edits within its error rate).
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cstdint>
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

constexpr uint64_t INF = UINT64_MAX / 4;
const char* const kModes[] = {"fixed", "rare", "long", "err1", "err1r"};
constexpr int kNumModes = 5;

struct Opts {
  int threads = 0, budget = 4096, k = 15, ext = 40;
  uint64_t cap = 20000;
  bool mode[kNumModes] = {true, true, true, true, true};
  std::string truth, out;
  long of = 0;
};

uint8_t code_of(char c) {
  switch (c) {
    case 'A': case 'a': return 0;
    case 'C': case 'c': return 1;
    case 'G': case 'g': return 2;
    case 'T': case 't': return 3;
    default: return 4;
  }
}

struct Read {
  std::string name;
  std::vector<uint8_t> seq;
};

std::vector<Read> read_fastq(const char* path) {
  std::ifstream f(path);
  if (!f) { std::fprintf(stderr, "cannot open %s\n", path); std::exit(1); }
  std::vector<Read> out;
  std::string h, s, p, q;
  while (std::getline(f, h) && std::getline(f, s) && std::getline(f, p) && std::getline(f, q)) {
    Read r;
    r.name = h.substr(1, h.find_first_of(" \t") == std::string::npos ? std::string::npos : h.find_first_of(" \t") - 1);
    for (char c : s) r.seq.push_back(code_of(c));
    out.push_back(std::move(r));
  }
  return out;
}

// name -> NM of the primary alignment, -1 when unmapped.
std::unordered_map<std::string, int> read_truth(const std::string& path) {
  std::unordered_map<std::string, int> m;
  std::ifstream f(path);
  if (!f) { std::fprintf(stderr, "cannot open %s\n", path.c_str()); std::exit(1); }
  std::string line;
  while (std::getline(f, line)) {
    if (line.empty() || line[0] == '@') continue;
    std::vector<std::string> t;
    size_t a = 0;
    for (size_t i = 0; i <= line.size(); ++i)
      if (i == line.size() || line[i] == '\t') { t.push_back(line.substr(a, i - a)); a = i + 1; }
    const int flag = std::atoi(t[1].c_str());
    if (flag & 0x900) continue;
    int nm = -1;
    if (!(flag & 4))
      for (size_t i = 11; i < t.size(); ++i)
        if (t[i].compare(0, 5, "NM:i:") == 0) nm = std::atoi(t[i].c_str() + 5);
    m[t[0]] = nm;
  }
  return m;
}

// One strand of a read: its q-mers' index ranges.
struct Strand {
  const uint8_t* s = nullptr;
  int L = 0, n = 0;  // n = L - q + 1 q-mer starts
  std::vector<uint64_t> lo, hi, c;
};

void prepare(const IndexView& ix, Strand& st) {
  const int q = ix.q;
  st.n = st.L >= q ? st.L - q + 1 : 0;
  st.lo.assign(st.n, 0);
  st.hi.assign(st.n, 0);
  st.c.assign(st.n, INF);
  for (int a = 0; a < st.n; ++a) {
    uint64_t key;
    if (!kmer_key(st.s + a, q, &key)) continue;
    lookup(ix, key, &st.lo[a], &st.hi[a]);
    st.c[a] = st.hi[a] - st.lo[a];
  }
}

// The most disjoint parts [a, a + len), q <= len <= maxlen, whose costs sum
// to <= B: f[k][i] = least total cost of k parts within [0, i).
std::vector<std::pair<int, int>> best_parts(int L, int q, int maxlen, const std::vector<std::vector<uint64_t>>& cost,
                                            uint64_t B) {
  const int K = L / q;
  std::vector<std::vector<uint64_t>> f(K + 1, std::vector<uint64_t>(L + 1, INF));
  std::vector<std::vector<int>> from(K + 1, std::vector<int>(L + 1, -1));
  for (int i = 0; i <= L; ++i) f[0][i] = 0;
  for (int k = 1; k <= K; ++k)
    for (int i = 1; i <= L; ++i) {
      f[k][i] = f[k][i - 1];
      for (int len = q; len <= maxlen && len <= i; ++len) {
        const int a = i - len;
        if (f[k - 1][a] >= INF || a >= (int)cost.size() || len - q >= (int)cost[a].size()) continue;
        const uint64_t c = cost[a][len - q];
        if (c >= INF) continue;
        if (f[k - 1][a] + c < f[k][i]) { f[k][i] = f[k - 1][a] + c; from[k][i] = a; }
      }
    }
  int m = 0;
  for (int k = K; k >= 1; --k)
    if (f[k][L] <= B) { m = k; break; }
  std::vector<std::pair<int, int>> parts;
  for (int i = L, k = m; k > 0;) {
    if (from[k][i] < 0) { --i; continue; }  // f[k][i] = f[k][i - 1]
    const int a = from[k][i];
    parts.push_back({a, i});
    i = a;
    --k;
  }
  return parts;
}

// Distinct q-mer keys within one edit of the q-mer at s: substitutions; one
// read base missing in the reference (the q - 1 others, then any base); one
// reference base inserted (at 0..q-1, the result cut to q bases).
void neighbours(const uint8_t* s, int q, std::vector<uint64_t>& keys) {
  keys.clear();
  uint8_t w[64];
  auto add = [&](const uint8_t* v) {
    uint64_t key;
    if (kmer_key(v, q, &key)) keys.push_back(key);
  };
  add(s);
  for (int i = 0; i < q; ++i)
    for (uint8_t x = 0; x < 4; ++x) {
      if (x == s[i]) continue;
      std::memcpy(w, s, q);
      w[i] = x;
      add(w);
    }
  for (int i = 0; i < q; ++i)
    for (uint8_t x = 0; x < 4; ++x) {
      int n = 0;
      for (int j = 0; j < q; ++j) if (j != i) w[n++] = s[j];
      w[n] = x;
      add(w);
    }
  for (int i = 0; i < q; ++i)
    for (uint8_t x = 0; x < 4; ++x) {
      int n = 0;
      for (int j = 0; j < i; ++j) w[n++] = s[j];
      w[n++] = x;
      for (int j = i; n < q; ++j) w[n++] = s[j];
      add(w);
    }
  std::sort(keys.begin(), keys.end());
  keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
}

struct ModeStats {
  long certified = 0, by_d1[64] = {}, truth_agree = 0, truth_fewer = 0, truth_more = 0, truth_unmapped = 0;
  double secs = 0;
  uint64_t verified = 0;  // band checks
  std::vector<std::string> bad;  // reads whose certified distance a truth SAM contradicts
};

struct Worker {
  const IndexView& ix;
  const Opts& o;
  std::vector<int> prev, cur;
  std::vector<int64_t> diags;
  std::vector<uint64_t> keys;
  explicit Worker(const IndexView& v, const Opts& op) : ix(v), o(op) {}

  // Best distance (<= R) over the hits' bands, else R + 1. Each hit at
  // diagonal g puts any alignment with <= R edits through it within
  // [g - R - 1, g + R + 1] (one more for err1's indel hits).
  int verify(const Strand& st, std::vector<int64_t>& d, int R, uint64_t* checks) {
    std::sort(d.begin(), d.end());
    d.erase(std::unique(d.begin(), d.end()), d.end());
    const int B = 2 * R + 5;
    prev.resize(B);
    cur.resize(B);
    int best = R + 1;
    int64_t last = INT64_MIN / 2;
    for (int64_t g : d) {
      if (g - last <= 1) continue;  // the last band already holds every such path
      last = g;
      ++*checks;
      const int v = band_distance_scalar(ix, st.s, st.L, g - R - 2, B, R, prev.data(), cur.data());
      if (v < best) best = v;
      if (best == 0) break;
    }
    return best;
  }

  // One read, one mode: certified distance, or -1.
  int run(int mode, Strand* sd, uint64_t* checks) {
    const int q = ix.q;
    const uint64_t B = (uint64_t)o.budget;
    int m[2] = {0, 0};
    std::vector<std::pair<int, int>> parts[2];
    std::vector<std::vector<uint64_t>> nb_lo[2], nb_hi[2];  // err1: per part, its neighbours' ranges
    for (int s = 0; s < 2; ++s) {
      Strand& st = sd[s];
      if (st.n == 0) return -1;
      if (mode == 0) {  // fixed: CERTA's parts, rarest first while they fit
        const int P = part_count(st.L, q, ix.s);
        std::vector<std::pair<uint64_t, int>> by;
        for (int j = 0; j < P; ++j) {
          int off, len;
          part_geometry(st.L, P, j, &off, &len);
          by.push_back({st.c[off], j});
        }
        std::stable_sort(by.begin(), by.end(), [](auto& x, auto& y) { return x.first < y.first; });
        uint64_t used = 0;
        for (auto& [c, j] : by) {
          if (c >= INF || used + c > B) break;
          used += c;
          int off, len;
          part_geometry(st.L, P, j, &off, &len);
          parts[s].push_back({off, off + q});
        }
      } else if (mode == 1) {  // rare: q-mers anywhere
        std::vector<std::vector<uint64_t>> cost(st.n, std::vector<uint64_t>(1));
        for (int a = 0; a < st.n; ++a) cost[a][0] = st.c[a];
        parts[s] = best_parts(st.L, q, q, cost, B);
      } else if (mode == 2) {  // long: exact parts of length q..q + 2E, counted through their rarest q-mer
        const int E = o.ext, W = E + 1;
        // Per anchor q-mer j with <= cap hits: S[j][x][y] = hits whose exact match extends >= x bases
        // left and >= y bases right of the q-mer.
        std::vector<std::vector<uint32_t>> S(st.n);
        for (int j = 0; j < st.n; ++j) {
          if (st.c[j] > o.cap) continue;
          std::vector<uint32_t> H((size_t)W * W, 0);
          for (uint64_t e = st.lo[j]; e < st.hi[j]; ++e) {
            const int64_t p = ix.pos[e];
            int el = 0, er = 0;
            while (el < E && j - 1 - el >= 0 && sub_cost(st.s[j - 1 - el], ref_at(ix, p - 1 - el)) == 0) ++el;
            while (er < E && j + q + er < st.L && sub_cost(st.s[j + q + er], ref_at(ix, p + q + er)) == 0) ++er;
            ++H[(size_t)el * W + er];
          }
          for (int x = E; x >= 0; --x)
            for (int y = E; y >= 0; --y) {
              uint32_t v = H[(size_t)x * W + y];
              if (x < E) v += H[(size_t)(x + 1) * W + y];
              if (y < E) v += H[(size_t)x * W + y + 1];
              if (x < E && y < E) v -= H[(size_t)(x + 1) * W + y + 1];
              H[(size_t)x * W + y] = v;
            }
          S[j] = std::move(H);
        }
        const int maxlen = q + 2 * E;
        std::vector<std::vector<uint64_t>> cost(st.n, std::vector<uint64_t>(maxlen - q + 1, INF));
        for (int a = 0; a < st.n; ++a)
          for (int len = q; len <= maxlen && a + len <= st.L; ++len) {
            const int b = a + len;
            int jb = -1;  // the anchor with the fewest hits that the part can be counted through
            for (int j = std::max(a, b - q - E); j <= std::min(b - q, a + E); ++j)
              if (!S[j].empty() && (jb < 0 || st.c[j] < st.c[jb])) jb = j;
            if (jb >= 0) cost[a][len - q] = S[jb][(size_t)(jb - a) * W + (b - jb - q)];
          }
        parts[s] = best_parts(st.L, q, maxlen, cost, B);
        // Each part's hits are its anchor's hits whose exact match spans the part.
        anchors[s].clear();
        for (auto& [a, b] : parts[s]) {
          int jb = -1;
          for (int j = std::max(a, b - q - E); j <= std::min(b - q, a + E); ++j)
            if (!S[j].empty() && (jb < 0 || st.c[j] < st.c[jb])) jb = j;
          anchors[s].push_back(jb);
        }
      } else {  // err1: q-mers with their one-edit neighbourhoods; err1r: only the rarest disjoint q-mers'
        std::vector<std::vector<uint64_t>> cost(st.n, std::vector<uint64_t>(1, INF));
        std::vector<std::vector<uint64_t>> lo(st.n), hi(st.n);
        auto expand = [&](int a) {
          neighbours(st.s + a, q, keys);
          uint64_t tot = 0;
          for (uint64_t key : keys) {
            uint64_t l, h;
            lookup(ix, key, &l, &h);
            if (h > l) { lo[a].push_back(l); hi[a].push_back(h); tot += h - l; }
          }
          cost[a][0] = tot;
        };
        if (mode == 3) {
          for (int a = 0; a < st.n; ++a)
            if (st.c[a] < INF) expand(a);
          parts[s] = best_parts(st.L, q, q, cost, B);
        } else {
          // The most disjoint q-mers with the fewest exact hits, then the
          // cheapest of their neighbourhoods while they fit the budget.
          std::vector<std::vector<uint64_t>> exact(st.n, std::vector<uint64_t>(1));
          for (int a = 0; a < st.n; ++a) exact[a][0] = st.c[a];
          std::vector<std::pair<uint64_t, int>> by;
          for (auto& [a, b] : best_parts(st.L, q, q, exact, INF - 1)) { expand(a); by.push_back({cost[a][0], a}); }
          std::sort(by.begin(), by.end());
          uint64_t used = 0;
          for (auto& [c, a] : by) {
            if (used + c > B) break;
            used += c;
            parts[s].push_back({a, a + q});
          }
        }
        nb_lo[s].clear();
        nb_hi[s].clear();
        for (auto& [a, b] : parts[s]) { nb_lo[s].push_back(lo[a]); nb_hi[s].push_back(hi[a]); }
      }
      m[s] = (int)parts[s].size();
    }
    const int ms = std::min(m[0], m[1]);
    int R = mode >= 3 ? 2 * ms - 1 : ms - 1;
    R = std::min(R, o.k);
    if (R < 0) return -1;
    int best = R + 1;
    for (int s = 0; s < 2; ++s) {
      Strand& st = sd[s];
      diags.clear();
      for (size_t x = 0; x < parts[s].size(); ++x) {
        const int a = parts[s][x].first, b = parts[s][x].second;
        if (mode >= 3) {
          for (size_t y = 0; y < nb_lo[s][x].size(); ++y)
            for (uint64_t e = nb_lo[s][x][y]; e < nb_hi[s][x][y]; ++e) diags.push_back((int64_t)ix.pos[e] - a);
        } else if (mode == 2) {
          const int j = anchors[s][x], E = o.ext;
          for (uint64_t e = st.lo[j]; e < st.hi[j]; ++e) {
            const int64_t p = ix.pos[e];
            int el = 0, er = 0;
            while (el < E && j - 1 - el >= 0 && sub_cost(st.s[j - 1 - el], ref_at(ix, p - 1 - el)) == 0) ++el;
            while (er < E && j + q + er < st.L && sub_cost(st.s[j + q + er], ref_at(ix, p + q + er)) == 0) ++er;
            if (el >= j - a && er >= b - j - q) diags.push_back(p - j);
          }
        } else {
          for (uint64_t e = st.lo[a]; e < st.hi[a]; ++e) diags.push_back((int64_t)ix.pos[e] - a);
        }
      }
      best = std::min(best, verify(st, diags, R, checks));
      if (best == 0) break;
    }
    return best <= R ? best : -1;
  }

  std::vector<int> anchors[2];  // long: each part's anchor q-mer
};

}  // namespace

int main(int argc, char** argv) {
  Opts o;
  std::vector<const char*> pos;
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    auto next = [&]() { if (i + 1 >= argc) { std::fprintf(stderr, "%s needs a value\n", a.c_str()); std::exit(2); } return std::string(argv[++i]); };
    if (a == "-t") o.threads = std::atoi(next().c_str());
    else if (a == "--budget") o.budget = std::atoi(next().c_str());
    else if (a == "--k") o.k = std::atoi(next().c_str());
    else if (a == "--cap") o.cap = std::strtoull(next().c_str(), nullptr, 10);
    else if (a == "--ext") o.ext = std::atoi(next().c_str());
    else if (a == "--truth") o.truth = next();
    else if (a == "--out") o.out = next();
    else if (a == "--of") o.of = std::atol(next().c_str());
    else if (a == "--modes") {
      const std::string v = "," + next() + ",";
      for (int k = 0; k < kNumModes; ++k) o.mode[k] = v.find("," + std::string(kModes[k]) + ",") != std::string::npos;
    } else pos.push_back(argv[i]);
  }
  if (pos.size() != 2) {
    std::fprintf(stderr, "usage: payoff INDEX.cidx READS.fq [-t N] [--budget 4096] [--k 15] [--cap 20000] [--ext 40] "
                         "[--modes fixed,rare,long,err1,err1r] [--truth tool.sam] [--out per_read.tsv] [--of N]\n");
    return 2;
  }
  Reference ref;
  Index index;
  load_index(pos[0], ref, index);
  const IndexView ix = index.view(ref);
  if (ix.s != 1) { std::fprintf(stderr, "the index must have s = 1\n"); return 2; }
  const std::vector<Read> reads = read_fastq(pos[1]);
  std::unordered_map<std::string, int> truth;
  if (!o.truth.empty()) truth = read_truth(o.truth);
  const int nt = o.threads > 0 ? o.threads : (int)std::thread::hardware_concurrency();
  std::fprintf(stderr, "payoff: %zu reads, q = %d, budget %d, k %d, cap %llu, ext %d, %d threads\n", reads.size(), ix.q,
               o.budget, o.k, (unsigned long long)o.cap, o.ext, nt);
  std::vector<std::vector<ModeStats>> stats(nt, std::vector<ModeStats>(kNumModes));
  std::vector<std::array<int8_t, kNumModes>> dist(reads.size());  // per read and mode: certified distance or -1
  std::vector<std::array<float, kNumModes>> usec(reads.size());   // and the time it took (microseconds)
  std::atomic<size_t> next{0};
  std::vector<std::thread> th;
  for (int w = 0; w < nt; ++w)
    th.emplace_back([&, w]() {
      Worker wk(ix, o);
      std::vector<uint8_t> rc;
      for (size_t i; (i = next.fetch_add(1)) < reads.size();) {
        const Read& r = reads[i];
        const int L = (int)r.seq.size();
        rc.resize(L);
        for (int x = 0; x < L; ++x) rc[x] = r.seq[L - 1 - x] < 4 ? (uint8_t)(3 - r.seq[L - 1 - x]) : (uint8_t)4;
        Strand sd[2];
        sd[0].s = r.seq.data();
        sd[1].s = rc.data();
        for (int s = 0; s < 2; ++s) { sd[s].L = L; prepare(ix, sd[s]); }
        auto t_it = truth.find(r.name);
        dist[i].fill(-1);
        usec[i].fill(0);
        for (int mode = 0; mode < kNumModes; ++mode) {
          if (!o.mode[mode]) continue;
          ModeStats& ms = stats[w][mode];
          const auto t0 = std::chrono::steady_clock::now();
          const int d = wk.run(mode, sd, &ms.verified);
          const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
          ms.secs += sec;
          dist[i][mode] = (int8_t)d;
          usec[i][mode] = (float)(1e6 * sec);
          if (d < 0) continue;
          ++ms.certified;
          ++ms.by_d1[std::min(d, 63)];
          if (t_it == truth.end()) continue;
          const int nm = t_it->second;
          if (nm < 0) ++ms.truth_unmapped;
          else if (nm == d) ++ms.truth_agree;
          else if (nm < d) { ++ms.truth_fewer; if (ms.bad.size() < 5) ms.bad.push_back(r.name); }
          else ++ms.truth_more;
        }
      }
    });
  for (auto& t : th) t.join();
  if (!o.out.empty()) {  // per read and mode: certified distance (-1: not certified or not run), microseconds
    std::FILE* f = std::fopen(o.out.c_str(), "w");
    if (!f) { std::fprintf(stderr, "cannot write %s\n", o.out.c_str()); return 1; }
    std::fprintf(f, "read");
    for (int mode = 0; mode < kNumModes; ++mode) std::fprintf(f, "\t%s\t%s_us", kModes[mode], kModes[mode]);
    std::fprintf(f, "\n");
    for (size_t i = 0; i < reads.size(); ++i) {
      std::fprintf(f, "%s", reads[i].name.c_str());
      for (int mode = 0; mode < kNumModes; ++mode) std::fprintf(f, "\t%d\t%.0f", dist[i][mode], usec[i][mode]);
      std::fprintf(f, "\n");
    }
    std::fclose(f);
  }
  const double of = o.of > 0 ? (double)o.of : (double)reads.size();
  for (int mode = 0; mode < kNumModes; ++mode) {
    if (!o.mode[mode]) continue;
    ModeStats a;
    for (int w = 0; w < nt; ++w) {
      const ModeStats& s = stats[w][mode];
      a.certified += s.certified;
      for (int d = 0; d < 64; ++d) a.by_d1[d] += s.by_d1[d];
      a.truth_agree += s.truth_agree; a.truth_fewer += s.truth_fewer; a.truth_more += s.truth_more;
      a.truth_unmapped += s.truth_unmapped;
      a.secs += s.secs;
      a.verified += s.verified;
      a.bad.insert(a.bad.end(), s.bad.begin(), s.bad.end());
    }
    std::printf("%-6s certified %ld of %zu (%.2f %% of %.0f reads); CPU %.1f s (%.2f ms per read), %.1f band checks per read\n",
                kModes[mode], a.certified, reads.size(), 100.0 * a.certified / of, of, a.secs,
                1e3 * a.secs / std::max<size_t>(1, reads.size()), (double)a.verified / std::max<size_t>(1, reads.size()));
    std::printf("       distance:");
    for (int d = 0; d < 64; ++d) if (a.by_d1[d]) std::printf(" %d: %ld", d, a.by_d1[d]);
    std::printf("\n");
    if (!truth.empty())
      std::printf("       vs truth: equal %ld, truth fewer edits %ld%s, truth more %ld, truth unmapped %ld\n", a.truth_agree,
                  a.truth_fewer, a.truth_fewer ? " (A CERTIFICATE IS WRONG)" : "", a.truth_more, a.truth_unmapped);
    for (const auto& b : a.bad) std::printf("       contradicted: %s\n", b.c_str());
  }
  return 0;
}
