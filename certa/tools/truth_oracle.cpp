// truth_oracle: the best local score of a read at given loci of the reference
// (bench/truth_oracle.py; bench/VALIDATION_RESULTS.md, "How far from the
// limit"). Shares no code with the mapper.
//
//   truth_oracle REF.fa R1.fq R2.fq JOBS.bin OUT.bin [-t THREADS] [-w WINDOW]
//
// JOBS.bin: 16-byte records sorted by pair: int32 pair (0-based index in the
// FASTQ files), int32 pos (0-based position on the contig of the read's first
// base as it aligns, clips included), int16 contig (index in REF.fa.fai),
// int8 mate (0 or 1), int8 strand (1: the read's reverse complement aligns),
// int32 unused. OUT.bin: one 8-byte record per job, in order: int16 score,
// int16 read bases aligned (clips excluded), int32 0-based contig position of
// the first aligned base; score -32768 when the window is empty.
//
// The score is CERTA's objective (BWA-MEM's defaults): +1 a match, -4 a
// mismatch, a gap of k bases 6 + k, 5 per clipped read end, maximized over
// every alignment that lies in [pos - W, pos + L + W) by a full DP of the
// window, with no band. Ties go to the first end in row-major order.
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

namespace {

constexpr int kMatch = 1, kMismatch = 4, kGapOpen = 6, kGapExt = 1, kClip = 5;
constexpr int NEG = -(1 << 20);

struct Job {
  int32_t pair, pos;
  int16_t contig;
  int8_t mate, strand;
  int32_t unused;
};
static_assert(sizeof(Job) == 16, "job layout");

struct Out {
  int16_t score, aligned;
  int32_t start;
};
static_assert(sizeof(Out) == 8, "out layout");

struct Contig {
  std::string name;
  int64_t len, off, line_bases, line_width;
  std::vector<uint8_t> seq;  // 0..3, 4 = N
};

uint8_t code(char c) {
  switch (c) {
    case 'A': case 'a': return 0;
    case 'C': case 'c': return 1;
    case 'G': case 'g': return 2;
    case 'T': case 't': return 3;
    default: return 4;
  }
}

std::vector<Contig> load_reference(const std::string& fa) {
  std::vector<Contig> cs;
  std::ifstream fai(fa + ".fai");
  if (!fai) { std::cerr << "cannot open " << fa << ".fai\n"; std::exit(1); }
  Contig c;
  while (fai >> c.name >> c.len >> c.off >> c.line_bases >> c.line_width) cs.push_back(c);
  FILE* f = std::fopen(fa.c_str(), "rb");
  if (!f) { std::cerr << "cannot open " << fa << "\n"; std::exit(1); }
  std::vector<char> buf;
  for (Contig& x : cs) {
    const int64_t lines = (x.len + x.line_bases - 1) / x.line_bases;
    const int64_t bytes = x.len + lines * (x.line_width - x.line_bases);
    buf.resize(bytes);
#if defined(_WIN32)
    _fseeki64(f, x.off, SEEK_SET);
#else
    fseeko(f, x.off, SEEK_SET);
#endif
    const size_t got = std::fread(buf.data(), 1, bytes, f);
    x.seq.resize(x.len);
    int64_t n = 0;
    for (size_t i = 0; i < got && n < x.len; ++i)
      if (buf[i] != '\n' && buf[i] != '\r') x.seq[n++] = code(buf[i]);
    if (n != x.len) { std::cerr << "short contig " << x.name << "\n"; std::exit(1); }
  }
  std::fclose(f);
  return cs;
}

// The reads of the pairs in `want` (sorted, unique), coded 0..4.
void load_reads(const std::string& fq, const std::vector<int32_t>& want, std::vector<std::vector<uint8_t>>& out) {
  std::ifstream f(fq);
  if (!f) { std::cerr << "cannot open " << fq << "\n"; std::exit(1); }
  out.assign(want.size(), {});
  std::string h, s, p, q;
  int64_t i = 0;
  size_t w = 0;
  while (w < want.size() && std::getline(f, h) && std::getline(f, s) && std::getline(f, p) && std::getline(f, q)) {
    if (i == want[w]) {
      if (!s.empty() && s.back() == '\r') s.pop_back();
      out[w].resize(s.size());
      for (size_t k = 0; k < s.size(); ++k) out[w][k] = code(s[k]);
      ++w;
    }
    ++i;
  }
  if (w != want.size()) { std::cerr << fq << ": pair " << want[w] << " not found\n"; std::exit(1); }
}

struct Scratch {
  std::vector<int> h[2], e[2], si[2], sc_[2];
  std::vector<int> es[2], ec[2];
};

// The full DP of read rd[0, L) against ref[0, W). Returns the best score, the
// number of read bases aligned and the window column of the first aligned
// base.
Out align(const uint8_t* rd, int L, const uint8_t* ref, int W, std::vector<int>* H, std::vector<int>* E,
          std::vector<int>* SI, std::vector<int>* SC, std::vector<int>* ESI, std::vector<int>* ESC) {
  // H: best alignment ending at (i, c) in any state; E: ending with an
  // insertion (read base i, no reference base). SI/SC: the read base and the
  // column where that alignment starts (ESI/ESC for E).
  for (int x = 0; x < 2; ++x) {
    H[x].assign(W, NEG); E[x].assign(W, NEG);
    SI[x].assign(W, 0); SC[x].assign(W, 0); ESI[x].assign(W, 0); ESC[x].assign(W, 0);
  }
  int best = NEG, bi = 0, bsi = 0, bsc = 0;
  int cur = 1, prv = 0;
  for (int i = 1; i <= L; ++i) {
    const uint8_t r = rd[i - 1];
    const int start = i == 1 ? 0 : -kClip;
    const int clip_end = i < L ? kClip : 0;
    int f = NEG, fsi = 0, fsc = 0, hl = NEG, hlsi = 0, hlsc = 0;
    int* hp = H[prv].data(); int* hc = H[cur].data();
    int* ep = E[prv].data(); int* ec = E[cur].data();
    int* sip = SI[prv].data(); int* sic = SI[cur].data();
    int* scp = SC[prv].data(); int* scc = SC[cur].data();
    int* esip = ESI[prv].data(); int* esic = ESI[cur].data();
    int* escp = ESC[prv].data(); int* escc = ESC[cur].data();
    for (int c = 0; c < W; ++c) {
      int diag = c > 0 ? hp[c - 1] : NEG, dsi = c > 0 ? sip[c - 1] : 0, dsc = c > 0 ? scp[c - 1] : 0;
      if (start >= diag) { diag = start; dsi = i; dsc = c; }
      const bool mis = r > 3 || ref[c] > 3 || r != ref[c];
      const int m = diag + (mis ? -kMismatch : kMatch);
      int e = NEG, esi = 0, esc = 0;
      if (i > 1) {
        const int open = hp[c] - kGapOpen - kGapExt, ext = ep[c] - kGapExt;
        if (open >= ext) { e = open; esi = sip[c]; esc = scp[c]; }
        else { e = ext; esi = esip[c]; esc = escp[c]; }
      }
      int fo = NEG, fosi = 0, fosc = 0;
      if (c > 0) {
        const int open = hl - kGapOpen - kGapExt, ext = f - kGapExt;
        if (open >= ext) { fo = open; fosi = hlsi; fosc = hlsc; }
        else { fo = ext; fosi = fsi; fosc = fsc; }
      }
      f = fo; fsi = fosi; fsc = fosc;
      int h = m, hsi = dsi, hsc = dsc;
      bool is_m = true;
      if (e > h) { h = e; hsi = esi; hsc = esc; is_m = false; }
      if (f > h) { h = f; hsi = fsi; hsc = fsc; is_m = false; }
      hc[c] = h; sic[c] = hsi; scc[c] = hsc;
      ec[c] = e; esic[c] = esi; escc[c] = esc;
      hl = h; hlsi = hsi; hlsc = hsc;
      if (is_m) {
        const int end = h - clip_end;
        if (end > best) { best = end; bi = i; bsi = hsi; bsc = hsc; }
      }
    }
    std::swap(cur, prv);
  }
  Out o;
  if (best < NEG / 2) { o.score = INT16_MIN; o.aligned = 0; o.start = 0; return o; }
  o.score = (int16_t)best;
  o.aligned = (int16_t)(bi - bsi + 1);
  o.start = bsc;
  return o;
}

}  // namespace

int main(int argc, char** argv) {
  int threads = 8, W = 150;
  std::vector<std::string> pos;
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    if (a == "-t" && i + 1 < argc) threads = std::atoi(argv[++i]);
    else if (a == "-w" && i + 1 < argc) W = std::atoi(argv[++i]);
    else pos.push_back(a);
  }
  if (pos.size() != 5) {
    std::cerr << "usage: truth_oracle REF.fa R1.fq R2.fq JOBS.bin OUT.bin [-t THREADS] [-w WINDOW]\n";
    return 1;
  }
  std::vector<Job> jobs;
  {
    std::ifstream f(pos[3], std::ios::binary | std::ios::ate);
    if (!f) { std::cerr << "cannot open " << pos[3] << "\n"; return 1; }
    const std::streamsize n = f.tellg();
    f.seekg(0);
    jobs.resize(n / sizeof(Job));
    f.read(reinterpret_cast<char*>(jobs.data()), jobs.size() * sizeof(Job));
  }
  std::vector<int32_t> want;
  for (const Job& j : jobs)
    if (want.empty() || want.back() != j.pair) {
      if (!want.empty() && j.pair < want.back()) { std::cerr << "jobs are not sorted by pair\n"; return 1; }
      want.push_back(j.pair);
    }
  std::cerr << "[oracle] " << jobs.size() << " jobs, " << want.size() << " pairs\n";
  std::vector<Contig> ref;
  std::vector<std::vector<uint8_t>> reads[2];
  {
    std::thread t0([&] { ref = load_reference(pos[0]); });
    std::thread t1([&] { load_reads(pos[1], want, reads[0]); });
    std::thread t2([&] { load_reads(pos[2], want, reads[1]); });
    t0.join(); t1.join(); t2.join();
  }
  std::cerr << "[oracle] reference and reads loaded\n";
  std::vector<Out> out(jobs.size());
  std::vector<std::thread> pool;
  for (int t = 0; t < threads; ++t)
    pool.emplace_back([&, t] {
      std::vector<int> H[2], E[2], SI[2], SC[2], ESI[2], ESC[2];
      std::vector<uint8_t> rd;
      for (size_t x = t; x < jobs.size(); x += threads) {
        const Job& j = jobs[x];
        const size_t w = std::lower_bound(want.begin(), want.end(), j.pair) - want.begin();
        const std::vector<uint8_t>& fw = reads[j.mate][w];
        const int L = (int)fw.size();
        Out o;
        o.score = INT16_MIN; o.aligned = 0; o.start = 0;
        if (j.contig >= 0 && j.contig < (int)ref.size() && L > 0) {
          const Contig& c = ref[j.contig];
          rd.resize(L);
          for (int i = 0; i < L; ++i) {
            const uint8_t b = fw[j.strand ? L - 1 - i : i];
            rd[i] = j.strand ? (b < 4 ? (uint8_t)(3 - b) : (uint8_t)4) : b;
          }
          const int64_t lo = std::max<int64_t>(0, (int64_t)j.pos - W);
          const int64_t hi = std::min<int64_t>(c.len, (int64_t)j.pos + L + W);
          if (hi > lo) {
            o = align(rd.data(), L, c.seq.data() + lo, (int)(hi - lo), H, E, SI, SC, ESI, ESC);
            if (o.score != INT16_MIN) o.start += (int32_t)lo;
          }
        }
        out[x] = o;
      }
    });
  for (std::thread& t : pool) t.join();
  std::ofstream f(pos[4], std::ios::binary);
  f.write(reinterpret_cast<const char*>(out.data()), out.size() * sizeof(Out));
  std::cerr << "[oracle] done\n";
  return 0;
}
