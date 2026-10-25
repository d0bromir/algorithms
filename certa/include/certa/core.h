// CERTA certified short-read fast path: per-read core.
//
// This header is compiled both by the host C++ compiler (CPU mode, any
// architecture) and by nvcc (GPU mode). Every function is CERTA_HD, so the
// CPU and GPU back-ends run the *same* code and must produce bit-identical
// results.
//
// Certificate (see RESEARCH_ROADMAP.md, section 3.3), adaptive version:
//   A read of length L is split into P disjoint parts, as many as fit
//   (each part needs length >= q + s - 1). Under unit-cost edit distance a
//   locus with d edits touches at most d parts. So if every occurrence of
//   each part in a set S is enumerated, every locus with d < |S| edits has
//   an untouched part in S, occurs exactly there, and is found:
//   the certified radius is R = |S| - 1, whichever parts S contains.
//   The index stores the q-mer at every s-th reference position; a part is
//   queried at all s shifts, so each exact occurrence contains a sampled
//   position. Hit counts are known from the index ranges before anything is
//   enumerated, so S is chosen as the *rarest* parts that fit a per-strand
//   budget, up to k + 1 of them. The same |S| is used on both strands, and
//   every candidate cluster is verified by exact banded edit distance.
//   An 'N' never matches (it costs 1, even against 'N'), so a part containing
//   'N' is never chosen; it is an edit for every locus, which the argument
//   already allows.
//   Certified repeat (tier SR): when no part fits the budget, hits of the
//   rarest part are sampled. If >= 2 distinct exact full-read matches are
//   found, the optimum is 0 edits and the read is a proven multi-mapper
//   (MAPQ 0) even though not every copy was enumerated.
//   Certified local alignment (tier SL, host only): see certify_local.
#pragma once
#include <climits>
#include <cstdint>
#include <cstring>
#ifndef __CUDA_ARCH__
#include <algorithm>
#endif

#ifdef __CUDACC__
#define CERTA_HD __host__ __device__
#else
#define CERTA_HD
#endif

namespace certa {

constexpr int KMAX = 5;                         // max certified radius (k)
// Max parts per read. A q = 22 index gives 6 parts for 150 bp; finer
// indexes (q = 15: 10 parts) lift the certifiable ceiling (Aim 1).
constexpr int PMAX = 16;
// Max hits enumerated per strand. The GPU keeps its workspace in per-thread
// stack memory, so its budget is small; on the host (CPU passes, including the
// CPU third pass of GPU runs) the workspace is heap-allocated.
constexpr int GPU_BUDGET_MAX = 256;
constexpr int HOST_BUDGET_MAX = 4096;
#ifdef __CUDA_ARCH__
constexpr int BUDGET_MAX = GPU_BUDGET_MAX;
#else
constexpr int BUDGET_MAX = HOST_BUDGET_MAX;
#endif
constexpr int SMAX = 16;                        // max index sampling step
constexpr int LMAX = 320;                       // max read length handled
// Max band width (diagonals). The host allows wider bands, for tiers SL and
// PR. The end-to-end certificate (S0/S1/S2) always uses BAND_E2E, so CPU and
// GPU runs give byte-identical results.
constexpr int BAND_E2E = 48;
#ifdef __CUDA_ARCH__
constexpr int BMAX = BAND_E2E;
#else
constexpr int BMAX = 96;
#endif
constexpr int MAX_CAND = BUDGET_MAX;            // candidates per strand
constexpr int MAX_CLUST = 2 * MAX_CAND;         // clusters over both strands
constexpr int MAX_CIGAR = 32;                   // S2_MAX edits, clips, affine gaps
constexpr int REP_SAMPLE = 64;                  // hits examined for tier SR
constexpr int REP_KEEP = 8;                     // exact copies kept for tier SR
constexpr int BIG = 1 << 20;

// S0/S1/SR are certified. S2 is a heuristic placement, reported with
// Result::certified = 2 and a conservative MAPQ.
// SL: certified local (clipped) alignment; PR: certified proper pair (both host only).
// Interval output (host only, certa/interval.h): R exact by tier R, IV an
// interval [lo, hi] with lo < hi, MR a mate placed by the pair window search.
enum Tier : uint8_t { kTierS0 = 0, kTierS1 = 1, kTierSR = 2, kTierS2 = 3, kTierSL = 4, kTierPR = 5,
                      kTierR = 6, kTierIV = 7, kTierMR = 8 };
constexpr int S2_MAX = 10;  // largest --s2 edit limit (band must fit BMAX)
constexpr int S2_TOP = 8;   // S2 evaluates only this many best-supported clusters
static_assert(MAX_CIGAR >= 2 * S2_MAX + 1, "CIGAR buffer too small for S2");
static_assert(4 * S2_MAX + 1 <= BAND_E2E, "S2 band does not fit BAND_E2E");
// Adaptive parts (pass V, host only; certify_vparts): largest certified
// radius, and the most one-edit neighbours of a q-mer (q <= 32).
constexpr int VMAX = 12;
constexpr int NBMAX = 1 + 11 * 32;
static_assert(MAX_CIGAR >= 2 * VMAX + 3, "CIGAR buffer too small for pass V");

enum Reason : uint8_t {
  kOk = 0,
  kBadLength = 1,        // read shorter than one part needs, or longer than LMAX
  kRadiusNegative = 2,   // every part exceeds the budget and no SR proof
  kNotFound = 3,         // no locus within the certified radius
  kClusterTooWide = 4,   // candidate cluster wider than BAND_E2E (tandem repeat)
  kCrossContig = 5,      // set by the host: alignment spans a contig boundary
  kNumReasons = 6
};

// BAM CIGAR op codes.
constexpr uint32_t kOpM = 0, kOpI = 1, kOpD = 2, kOpS = 4;

// Reported alignments use BWA-MEM's default affine scoring (the certificate
// itself is about unit-cost edit distance and picks the locus).
constexpr int kMatch = 1, kMismatch = 4, kGapOpen = 6, kGapExt = 1, kClip = 5;
constexpr int16_t kNoSub = -32768;  // Result::sub_score when no other locus was found

// Lemma E. A local alignment of a read of length L whose end-to-end
// completion (clipped ends aligned base by base) has d >= 1 edits scores at
// most edit_bound(L, d). Its loss L - score is a sum over events, each
// costing more than the edits it adds: a mismatch 5 for 1 edit, a clipped end
// c + 5 for <= c, a gap of k 6 + k (deletion) or 6 + 2k (insertion) for k.
// One edit: >= 5. Two or more: one event costs >= d + 5 (a clip or gap),
// several cost >= d + 8.
CERTA_HD inline int edit_bound(int L, int d) {
  static_assert(kMatch == 1 && kMismatch == 4 && kGapOpen == 6 && kGapExt == 1 && kClip == 5,
                "Lemma E is stated for BWA-MEM's default scores");
  return d <= 0 ? L * kMatch : L * kMatch - d - (d >= 2 ? 5 : 4);
}

struct IndexView {
  const uint8_t* ref;     // concatenated reference, codes 0..3, 4 = N
  uint64_t ref_len;
  const uint64_t* keys;   // sorted q-mer keys of sampled positions, or null:
  const uint16_t* keys16; // low 16 bits only, when 2q - dir_bits <= 16 (the
                          // bucket directory already fixes the high bits)
  const uint32_t* pos;    // positions (multiples of s), parallel to keys
  const uint32_t* dir;    // bucket directory, size 2^dir_bits + 1 (n < 2^32)
  uint64_t n;
  int q, s, dir_bits;
};

struct Params {
  int k;       // maximum certified radius sought (at most k + 1 parts used)
  int budget;  // max hits enumerated per strand
  int reverse_order = 0;  // testing only: verify clusters least-supported first
  // Tier S2 (not certified): when no locus lies within the certified radius,
  // align the enumerated candidates with this larger edit limit (0 = off).
  int s2_limit = 0;
  // For MAPQ only: when the second-best locus lies beyond the certified
  // radius, estimate it among the other candidates up to this many edits.
  int mapq_limit = 0;
  // q-gram-lemma filter: only clusters hit by at least this many distinct
  // enumerated parts are verified; the radius becomes |S| - min_support.
  int min_support = 1;
  // Tier SL (host only): when no locus lies within R, certify the best local
  // (clipped) alignment; chains need hits of this many parts to be verified
  // (0 = off). local_only skips the end-to-end path (re-runs of reads that
  // already failed it).
  int local = 0;
  int local_only = 0;
  // Pass V (host only, certify_vparts): 1 = the rarest disjoint q-mers as
  // parts; 2 = also those q-mers searched with up to one edit. Needs s = 1.
  int vparts = 0;
};

struct Result {
  int64_t ref_pos;              // 0-based start in concatenated reference
  uint32_t cigar[MAX_CIGAR];    // BAM encoding: (len << 4) | op
  uint8_t n_cigar;
  uint8_t certified;            // 1 certified, 2 emitted without certificate (S2)
  uint8_t tier;                 // Tier: S0 exact, S1 1..R edits, SR repeat
  uint8_t strand;               // 0 forward, 1 reverse complement
  uint8_t reason;               // Reason when not certified
  int8_t radius;                // certified radius R (-1 for SR / uncertified)
  int8_t d1;                    // best (certified) edit distance
  uint8_t nm;                   // edits in the reported alignment (SAM NM)
  int16_t score;                // affine score of the reported alignment
  int16_t sub_score;            // affine score of the best other locus found (MAPQ), or kNoSub
  int8_t d2;                    // second-best distance within R, -1 if none
  int8_t d2x;                   // MAPQ only: best other candidate beyond R, -1 if none
  uint16_t n_best;              // loci tied at d1 (a lower bound for SR)
  uint16_t n_clusters;          // loci verified (diagnostic)
  uint8_t parts;                // parts the read was split into
  uint8_t used;                 // parts enumerated per strand (|S|)
  int16_t floor;                // SL: certified threshold (score > floor proves optimality)
  uint32_t ncand;               // candidate diagonals enumerated, both strands (diagnostic)
  // Certified tiers: an upper bound on SUB, the best local score at a locus
  // not overlapping the reported one (interval output, Aim 3), or kNoSub.
  int16_t sub_hi;
  // Witness (Theorem V): the parts enumerated per strand (bit j = part j)
  // and the support t of the certificate that used them.
  uint16_t probes[2];
  uint8_t support;
};

struct Cluster {
  int64_t lo;       // lowest diagonal (reference start) in band
  uint8_t width;    // number of diagonals in band
  uint8_t strand;
  int8_t dist;      // exact distance, or radius + 1 if it cannot beat the best two
  uint16_t beg, end;  // member diagonals: ws.cand[strand][beg, end)
  uint16_t pmask;     // distinct parts with a hit in the cluster (bit j = part j)
};

// Candidates store diagonal * 16 + part index (PMAX <= 16), so sorting orders
// by diagonal and the part survives the sort.
static_assert(PMAX <= 16, "part index must fit in 4 bits");
CERTA_HD inline int64_t diag_of(int64_t c) { return c >> 4; }  // floor(c / 16)
CERTA_HD inline int part_of(int64_t c) { return (int)(c & 15); }
CERTA_HD inline int popcount16(uint16_t v) {
  int n = 0;
  for (; v; v &= (uint16_t)(v - 1)) ++n;
  return n;
}

struct Workspace {
  uint8_t seq[2][LMAX];
  uint64_t rlo[2][PMAX][SMAX], rhi[2][PMAX][SMAX];  // index ranges per part/shift
  uint32_t count[2][PMAX];                          // hits per part (BIG if 'N')
  uint8_t order[2][PMAX];                           // parts, rarest first
  int64_t cand[2][MAX_CAND];
  int ncand[2];
  Cluster clusters[MAX_CLUST];
  uint16_t corder[MAX_CLUST];  // evaluation order of clusters
  int count_buf[MAX_CAND + 2];  // counting sort of clusters by support
  int row_a[BMAX], row_b[BMAX];
  int aff[4][BMAX];  // affine DP rows: H prev/cur, E prev/cur
  uint8_t tb[(LMAX + 1) * BMAX];
#ifndef __CUDA_ARCH__
  int16_t lub[MAX_CLUST];     // SL: score bound of each chain
  int16_t lscore[MAX_CLUST];  // SL: best score in each chain's bands, or kNoSub
  // SL/PR: Lemma L bound per (strand, part mask); kNoSub = not computed.
  // Initialised once (ub_magic), then only the entries a call touched are
  // reset at the start of the next call. ub_magic starts at 0: heap memory
  // that last held a Workspace still carries the magic while other
  // allocations may have overwritten its touched list since.
  uint64_t ub_magic = 0;
  int16_t ubmemo[2][1 << PMAX];
  uint32_t ub_touched[MAX_CLUST];
  int ub_ntouched;
  // Pass V: the index range of the q-mer at each read offset (lo > hi: it
  // holds an 'N'), the DP over disjoint q-mers, the chosen parts, and the
  // index ranges of their one-edit neighbourhoods.
  uint64_t vlo[2][LMAX], vhi[2][LMAX];
  uint64_t vf[2][PMAX + 1][LMAX + 1];
  int16_t vfrom[2][PMAX + 1][LMAX + 1];
  int16_t voff[2][PMAX];
  uint64_t nblo[2][PMAX][NBMAX], nbhi[2][PMAX][NBMAX];
  uint16_t nbn[2][PMAX];
  uint64_t nbcost[2][PMAX];
#endif
};

CERTA_HD inline uint8_t ref_at(const IndexView& ix, int64_t j) {
  return (j < 0 || (uint64_t)j >= ix.ref_len) ? 4 : ix.ref[j];
}

CERTA_HD inline int sub_cost(uint8_t a, uint8_t b) {
  return (a == b && a < 4) ? 0 : 1;
}

CERTA_HD inline bool kmer_key(const uint8_t* s, int q, uint64_t* key) {
  uint64_t v = 0;
  for (int i = 0; i < q; ++i) {
    uint8_t c = s[i];
    if (c > 3) return false;
    v = (v << 2) | c;
  }
  *key = v;
  return true;
}

// Half-open range [*lo, *hi) of the entries whose key equals `key` within
// its bucket [l, r) = [dir[b], dir[b + 1]), b = key >> (2q - dir_bits).
CERTA_HD inline void lookup_bucket(const IndexView& ix, uint64_t key, uint64_t l, uint64_t r, uint64_t* lo,
                                   uint64_t* hi) {
  const int shift = 2 * ix.q - ix.dir_bits;
  const uint64_t r0 = r;
  if (ix.keys16) {  // within a bucket, keys are ordered by their low bits
    const uint16_t k = (uint16_t)(key & ((1ull << shift) - 1));
#ifndef __CUDA_ARCH__
    // A small bucket (most are; they fit a cache line or two): a linear scan
    // finds the same first key >= k and first key > k without the binary
    // search's unpredictable branches.
    if (r - l <= 64) {
      uint64_t a = l;
      while (a < r && ix.keys16[a] < k) ++a;
      *lo = a;
      while (a < r && ix.keys16[a] == k) ++a;
      *hi = a;
      return;
    }
#endif
    while (l < r) {
      uint64_t m = l + (r - l) / 2;
      if (ix.keys16[m] < k) l = m + 1; else r = m;
    }
    *lo = l;
    r = r0;
    while (l < r) {
      uint64_t m = l + (r - l) / 2;
      if (ix.keys16[m] <= k) l = m + 1; else r = m;
    }
    *hi = l;
    return;
  }
  while (l < r) {
    uint64_t m = l + (r - l) / 2;
    if (ix.keys[m] < key) l = m + 1; else r = m;
  }
  *lo = l;
  r = r0;
  while (l < r) {
    uint64_t m = l + (r - l) / 2;
    if (ix.keys[m] <= key) l = m + 1; else r = m;
  }
  *hi = l;
}

// Half-open range [*lo, *hi) of index entries whose key equals `key`.
CERTA_HD inline void lookup(const IndexView& ix, uint64_t key, uint64_t* lo,
                            uint64_t* hi) {
  const uint64_t b = key >> (2 * ix.q - ix.dir_bits);
  lookup_bucket(ix, key, ix.dir[b], ix.dir[b + 1], lo, hi);
}

// Part j of P parts over a read of length L: offset and length.
CERTA_HD inline void part_geometry(int L, int P, int j, int* off, int* len) {
  int base = L / P, rem = L % P;
  *off = j * base + (j < rem ? j : rem);
  *len = base + (j < rem ? 1 : 0);
}

// Parts a read of length L is split into: as many as fit, each >= q + s - 1.
CERTA_HD inline int part_count(int L, int q, int s) {
  int p = L / (q + s - 1);
  return p < PMAX ? p : PMAX;
}

// Shortest read with at least one part; radius R needs R + 1 parts.
CERTA_HD inline int min_read_length(int q, int s) { return q + s - 1; }

// Shell sort (Ciura gaps): O(n^1.3) instead of insertion sort's O(n^2), which
// stalled whole GPU warps on repetitive reads with hundreds of candidates.
CERTA_HD inline void sort_i64(int64_t* a, int n) {
  const int gaps[] = {132, 57, 23, 10, 4, 1};
  for (int g : gaps) {
    for (int i = g; i < n; ++i) {
      int64_t v = a[i];
      int j = i;
      while (j >= g && a[j - g] > v) { a[j] = a[j - g]; j -= g; }
      a[j] = v;
    }
  }
}

// Look up every part of one strand (index ranges only, nothing enumerated)
// and order the parts rarest first. Parts containing 'N' get count BIG.
CERTA_HD inline void count_parts(const IndexView& ix, const uint8_t* seq, int L,
                                 int P, Workspace& ws, int st) {
  for (int j = 0; j < P; ++j) {
    int off, len;
    part_geometry(L, P, j, &off, &len);
    bool has_n = false;
    for (int i = off; i < off + len; ++i) has_n |= (seq[i] > 3);
    uint64_t total = 0;
    // The part's s q-mers, rolled (each is kmer_key(seq + off + t)): q + s - 1
    // steps instead of s q. key holds the first q - 1 bases before shift 0.
    const uint64_t mask = ix.q >= 32 ? ~0ull : (1ull << (2 * ix.q)) - 1;
    uint64_t key = 0;
    if (!has_n)
      for (int i = 0; i < ix.q - 1; ++i) key = (key << 2) | seq[off + i];
    for (int t = 0; t < ix.s; ++t) {
      if (has_n) { ws.rlo[st][j][t] = ws.rhi[st][j][t] = 0; continue; }
      key = ((key << 2) | seq[off + t + ix.q - 1]) & mask;
      lookup(ix, key, &ws.rlo[st][j][t], &ws.rhi[st][j][t]);
      total += ws.rhi[st][j][t] - ws.rlo[st][j][t];
    }
    ws.count[st][j] = has_n ? BIG : (total > (uint64_t)BIG - 1 ? BIG - 1 : (uint32_t)total);
    // Insertion into the rarest-first order (stable: ties keep part order).
    int x = j;
    while (x > 0 && ws.count[st][ws.order[st][x - 1]] > ws.count[st][j]) {
      ws.order[st][x] = ws.order[st][x - 1];
      --x;
    }
    ws.order[st][x] = (uint8_t)j;
  }
}

// Largest m <= limit such that the m rarest parts together fit the budget.
CERTA_HD inline int parts_that_fit(const Workspace& ws, int st, int P, int limit,
                                   int budget) {
  uint32_t used = 0;
  int m = 0;
  while (m < P && m < limit) {
    uint32_t c = ws.count[st][ws.order[st][m]];
    if (c >= (uint32_t)BIG || used + c > (uint32_t)budget) break;
    used += c;
    ++m;
  }
  return m;
}

// Enumerate every hit of the m rarest parts as candidate diagonals.
CERTA_HD inline void enumerate_parts(const IndexView& ix, int L, int P, int m,
                                     Workspace& ws, int st) {
  int n = 0;
  for (int x = 0; x < m; ++x) {
    int j = ws.order[st][x], off, len;
    part_geometry(L, P, j, &off, &len);
    for (int t = 0; t < ix.s; ++t)
      for (uint64_t e = ws.rlo[st][j][t]; e < ws.rhi[st][j][t]; ++e)
        ws.cand[st][n++] = ((int64_t)ix.pos[e] - (off + t)) * 16 + j;  // diagonal, part
  }
  ws.ncand[st] = n;
}

// Semi-global unit-cost edit distance of the whole read against the
// reference, restricted to diagonals [lo, lo + B). Returns `limit + 1` as
// soon as every cell exceeds `limit` (exact whenever the result <= limit).
CERTA_HD inline int band_distance_scalar(const IndexView& ix, const uint8_t* rd, int L,
                                         int64_t lo, int B, int limit, int* prev,
                                         int* cur) {
  for (int b = 0; b < B; ++b) prev[b] = 0;
  for (int i = 1; i <= L; ++i) {
    const uint8_t c = rd[i - 1];
    int row_min = BIG;
    for (int b = 0; b < B; ++b) {
      int v = prev[b] + sub_cost(c, ref_at(ix, (int64_t)i - 1 + lo + b));
      if (b + 1 < B && prev[b + 1] + 1 < v) v = prev[b + 1] + 1;  // insertion
      if (b > 0 && cur[b - 1] + 1 < v) v = cur[b - 1] + 1;        // deletion
      cur[b] = v;
      if (v < row_min) row_min = v;
    }
    if (row_min > limit) return limit + 1;
    int* t = prev; prev = cur; cur = t;
  }
  int best = BIG;
  for (int b = 0; b < B; ++b) if (prev[b] < best) best = prev[b];
  return best;
}

CERTA_HD inline int hamming_scalar(const IndexView& ix, const uint8_t* rd, int L,
                                   int64_t start, int limit) {
  int d = 0;
  for (int i = 0; i < L && d <= limit; ++i)
    d += sub_cost(rd[i], ref_at(ix, start + i));
  return d;
}

#if defined(__GNUC__) && !defined(__CUDACC__) && !defined(CERTA_NO_VEC) && (defined(__clang__) || __GNUC__ >= 12)
#define CERTA_EDIT_VEC 1
// band_distance_scalar, 8 cells at a time (16-bit lanes; NEON or SSE), with
// the same result for every input, early stop included:
// - the match and insertion terms depend on the previous row only;
// - the deletion step cur[b] = min(v[b], cur[b-1] + 1) unrolls to
//   cur[b] = b + min over b' <= b of (v[b'] - b'), an inclusive prefix
//   minimum: three in-register shifts per vector and a carry across vectors
//   (the minimum form of band_row's identity, lean/Certa/BandRow.lean);
// - the row's minimum is min_b v[b], since cur[b] >= that minimum and the
//   cell that attains it keeps it.
// Distances stay below L + B, far under the 16-bit sentinel.
inline int band_distance_vec16(const IndexView& ix, const uint8_t* rd, int L, int64_t lo, int B, int limit) {
  typedef int16_t v8 __attribute__((vector_size(16)));
  typedef uint8_t u8x8 __attribute__((vector_size(8)));
  typedef int16_t v8u __attribute__((vector_size(16), aligned(2), may_alias));
  typedef uint8_t u8x8u __attribute__((vector_size(8), aligned(1), may_alias));
  constexpr int W = 8;
  constexpr int16_t HI = 16000;
  auto vmin = [](v8 a, v8 b) { return a < b ? a : b; };
  const v8 hi = {HI, HI, HI, HI, HI, HI, HI, HI}, one = {1, 1, 1, 1, 1, 1, 1, 1};
  const v8 lane = {0, 1, 2, 3, 4, 5, 6, 7};
  const int Bp = (B + W - 1) / W * W;
  // Rows: lanes [0, B) hold distances; [B, Bp + W) the sentinel, so no
  // insertion enters from beyond the band and padding lanes never win.
  int16_t rows[2][BMAX + 2 * W];
  int16_t* prev = rows[0];
  int16_t* cur = rows[1];
  for (int b = 0; b < Bp + W; ++b) { prev[b] = b < B ? 0 : HI; cur[b] = HI; }
  // The band's reference window, once: row i reads it from i - 1 (4 outside the reference).
  uint8_t refw[LMAX + BMAX + 2 * W];
  const int nw = L - 1 + Bp;
  for (int k = 0; k < nw; ++k) refw[k] = ref_at(ix, lo + k);
  for (int k = nw; k < nw + W; ++k) refw[k] = 4;
  for (int i = 1; i <= L; ++i) {
    const uint8_t c = rd[i - 1];
    const int16_t cc = c < 4 ? (int16_t)c : (int16_t)-1;  // an N matches nothing
    const v8 cv = {cc, cc, cc, cc, cc, cc, cc, cc};
    v8 carry = hi, rmin = hi;
    for (int b = 0; b < Bp; b += W) {
      const int16_t b16 = (int16_t)b;
    const v8 bv = lane + v8{b16, b16, b16, b16, b16, b16, b16, b16};
      const v8 p = *(const v8u*)(prev + b), pn = *(const v8u*)(prev + b + 1);
      const v8 r = __builtin_convertvector(*(const u8x8u*)(refw + (i - 1) + b), v8);
      v8 v = p + one + (r == cv);  // a match adds 0 (the comparison is -1), anything else 1
      v = vmin(v, pn + one);
      v8 x = v - bv;
      x = vmin(x, __builtin_shufflevector(x, hi, 8, 0, 1, 2, 3, 4, 5, 6));
      x = vmin(x, __builtin_shufflevector(x, hi, 8, 9, 0, 1, 2, 3, 4, 5));
      x = vmin(x, __builtin_shufflevector(x, hi, 8, 9, 10, 11, 0, 1, 2, 3));
      x = vmin(x, carry);
      carry = __builtin_shufflevector(x, x, 7, 7, 7, 7, 7, 7, 7, 7);
      v8 out = x + bv;
      if (b + W > B) {
        const int16_t bb = (int16_t)B;
        const v8 in = bv < v8{bb, bb, bb, bb, bb, bb, bb, bb};
        out = in ? out : hi;
        v = in ? v : hi;
      }
      *(v8u*)(cur + b) = out;
      rmin = vmin(rmin, v);
    }
    int row_min = HI;
    for (int l = 0; l < W; ++l) row_min = rmin[l] < row_min ? rmin[l] : row_min;
    if (row_min > limit) return limit + 1;
    int16_t* t = prev; prev = cur; cur = t;
  }
  int best = BIG;
  for (int b = 0; b < B; ++b) if (prev[b] < best) best = prev[b];
  return best;
}

// The same DP on values capped at cap = limit + 1, 16 cells at a time (8-bit
// lanes), for small limits and a band whose reference window lies inside the
// reference. The result is unchanged:
// - every recurrence term is non-decreasing in its argument and adds >= 0,
//   so min(D, cap) satisfies the recurrence with every term capped, and
//   cells outside the band (infinite) count as cap;
// - band_distance_scalar returns min(final row minimum, limit + 1), since
//   row minima never decrease (each cell is >= some cell of the row above),
//   and so does this.
// It stops once every cell of a row is at the cap.
inline int band_distance_vec8(const IndexView& ix, const uint8_t* rd, int L, int64_t lo, int B, int limit) {
  typedef int8_t v16 __attribute__((vector_size(16)));
  typedef int8_t v16u __attribute__((vector_size(16), aligned(1), may_alias));
  typedef uint64_t u64x2 __attribute__((vector_size(16)));
  constexpr int W = 16;
  auto vmin = [](v16 a, v16 b) { return a < b ? a : b; };
  auto splat = [](int x) {
    const int8_t s = (int8_t)x;
    return v16{s, s, s, s, s, s, s, s, s, s, s, s, s, s, s, s};
  };
  const int cap = limit + 1;
  const v16 capv = splat(cap), one = splat(1);
  const v16 lane = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  const int Bp = (B + W - 1) / W * W;
  int8_t rows[2][BMAX + 2 * W];
  int8_t* prev = rows[0];
  int8_t* cur = rows[1];
  for (int b = 0; b < Bp + W; ++b) { prev[b] = (int8_t)(b < B ? 0 : cap); cur[b] = (int8_t)cap; }
  const uint8_t* rf = ix.ref + lo;  // row i compares rd[i - 1] with rf[i - 1 + b]
  for (int i = 1; i <= L; ++i) {
    const uint8_t c = rd[i - 1];
    const v16 cv = splat(c < 4 ? c : -1);  // an N matches nothing
    v16 carry = capv, rmin = capv;
    for (int b = 0; b < Bp; b += W) {
      const v16 bv = lane + splat(b);
      const v16 p = *(const v16u*)(prev + b), pn = *(const v16u*)(prev + b + 1);
      const v16 r = *(const v16u*)(rf + (i - 1) + b);
      v16 v = p + one + (r == cv);  // a match adds 0 (the comparison is -1), anything else 1
      v = vmin(v, pn + one);
      v16 x = v - bv;
      x = vmin(x, __builtin_shufflevector(x, capv, 16, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14));
      x = vmin(x, __builtin_shufflevector(x, capv, 16, 17, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13));
      x = vmin(x, __builtin_shufflevector(x, capv, 16, 17, 18, 19, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11));
      x = vmin(x, __builtin_shufflevector(x, capv, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7));
      x = vmin(x, carry);
      carry = __builtin_shufflevector(x, x, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15);
      v16 out = vmin(x + bv, capv);
      if (b + W > B) out = bv < splat(B) ? out : capv;
      *(v16u*)(cur + b) = out;
      rmin = vmin(rmin, out);
    }
    const u64x2 below = (u64x2)(rmin < capv);
    if ((below[0] | below[1]) == 0) return cap;  // every cell at the cap: the result is limit + 1
    int8_t* t = prev; prev = cur; cur = t;
  }
  int best = cap;
  for (int b = 0; b < B; ++b) if (prev[b] < best) best = prev[b];
  return best;
}

// Dispatch: the 8-bit form when its window fits the reference and the cap
// fits its lanes (every caller's limit is <= S2_MAX), else the 16-bit form.
inline int band_distance_vec(const IndexView& ix, const uint8_t* rd, int L, int64_t lo, int B, int limit) {
  constexpr int W = 16;
  const int Bp = (B + W - 1) / W * W;
  if (limit >= 0 && limit <= 60 && lo >= 0 && (uint64_t)(lo + L - 1 + Bp + W) <= ix.ref_len)
    return band_distance_vec8(ix, rd, L, lo, B, limit);
  return band_distance_vec16(ix, rd, L, lo, B, limit);
}

// hamming_scalar, 16 positions at a time when the window lies inside the
// reference, with the same result: the exact count while it stays <= limit,
// else limit + 1.
inline int hamming_vec(const IndexView& ix, const uint8_t* rd, int L, int64_t start, int limit) {
  if (start < 0 || (uint64_t)(start + L) > ix.ref_len) return hamming_scalar(ix, rd, L, start, limit);
  typedef uint8_t u8x16 __attribute__((vector_size(16)));
  typedef uint8_t u8x16u __attribute__((vector_size(16), aligned(1), may_alias));
  const uint8_t* rf = ix.ref + start;
  const u8x16 four = {4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};
  int d = 0, i = 0;
  for (; i + 16 <= L; i += 16) {
    const u8x16 a = *(const u8x16u*)(rd + i), b = *(const u8x16u*)(rf + i);
    const auto ok = (a == b) & (a < four);  // -1 where the bases match
    int m = 0;
    for (int l = 0; l < 16; ++l) m += ok[l] ? 0 : 1;
    d += m;
    if (d > limit) return limit + 1;
  }
  for (; i < L && d <= limit; ++i) d += sub_cost(rd[i], rf[i]);
  return d <= limit ? d : limit + 1;
}
#endif

CERTA_HD inline int band_distance(const IndexView& ix, const uint8_t* rd, int L,
                                  int64_t lo, int B, int limit, int* prev, int* cur) {
#ifdef CERTA_EDIT_VEC
  (void)prev; (void)cur;
  return band_distance_vec(ix, rd, L, lo, B, limit);
#else
  return band_distance_scalar(ix, rd, L, lo, B, limit, prev, cur);
#endif
}

CERTA_HD inline int hamming(const IndexView& ix, const uint8_t* rd, int L,
                            int64_t start, int limit) {
#ifdef CERTA_EDIT_VEC
  return hamming_vec(ix, rd, L, start, limit);
#else
  return hamming_scalar(ix, rd, L, start, limit);
#endif
}

CERTA_HD inline void push_op(Result& r, uint32_t op) {
  if (r.n_cigar > 0 && (r.cigar[r.n_cigar - 1] & 0xF) == op)
    r.cigar[r.n_cigar - 1] += 1u << 4;
  else
    r.cigar[r.n_cigar++] = (1u << 4) | op;
}

// Exact alignment of the chosen band with traceback. Tie preference:
// substitution/match > insertion > deletion, and the smallest end diagonal.
CERTA_HD inline void band_traceback(const IndexView& ix, const uint8_t* rd,
                                    int L, int64_t lo, int B, Workspace& ws,
                                    Result& r) {
  int* prev = ws.row_a;
  int* cur = ws.row_b;
  for (int b = 0; b < B; ++b) prev[b] = 0;
  for (int i = 1; i <= L; ++i) {
    const uint8_t c = rd[i - 1];
    for (int b = 0; b < B; ++b) {
      int v = prev[b] + sub_cost(c, ref_at(ix, (int64_t)i - 1 + lo + b));
      uint8_t dir = 0;
      if (b + 1 < B && prev[b + 1] + 1 < v) { v = prev[b + 1] + 1; dir = 1; }
      if (b > 0 && cur[b - 1] + 1 < v) { v = cur[b - 1] + 1; dir = 2; }
      cur[b] = v;
      ws.tb[i * B + b] = dir;
    }
    int* t = prev; prev = cur; cur = t;
  }
  int b = 0;
  for (int x = 1; x < B; ++x) if (prev[x] < prev[b]) b = x;
  // Walk back, collecting ops in reverse.
  uint32_t rev[LMAX + BMAX];
  int nrev = 0, i = L;
  while (i > 0) {
    uint8_t dir = ws.tb[i * B + b];
    if (dir == 0) { rev[nrev++] = kOpM; --i; }
    else if (dir == 1) { rev[nrev++] = kOpI; --i; ++b; }
    else { rev[nrev++] = kOpD; --b; }
  }
  r.ref_pos = lo + b;  // row 0: reference consumed up to diagonal lo + b
  r.n_cigar = 0;
  for (int x = nrev - 1; x >= 0; --x) push_op(r, rev[x]);
}

// Groups the sorted candidate diagonals of both strands into clusters whose
// bands (+-pad) overlap. Returns the number of clusters, or -1 when a cluster
// is wider than max_width and skip_wide is false (with skip_wide it is
// dropped). Only host passes (pass V) widen max_width beyond BAND_E2E.
CERTA_HD inline int build_clusters(Workspace& ws, int pad, bool skip_wide, int max_width = BAND_E2E) {
  int nclust = 0;
  for (int st = 0; st < 2; ++st) {
    const int64_t* c = ws.cand[st];
    const int n = ws.ncand[st];
    int i = 0;
    while (i < n) {
      int64_t dmin = diag_of(c[i]), dmax = dmin;
      uint16_t pmask = (uint16_t)(1u << part_of(c[i]));
      int j = i + 1;
      while (j < n && diag_of(c[j]) - dmax <= 2 * pad) {
        dmax = diag_of(c[j]);
        pmask |= (uint16_t)(1u << part_of(c[j]));
        ++j;
      }
      const int64_t width = dmax - dmin + 2 * pad + 1;
      if (width > max_width) {
        if (!skip_wide) return -1;
        i = j;
        continue;
      }
      Cluster& cl = ws.clusters[nclust++];
      cl.lo = dmin - pad;
      cl.width = (uint8_t)width;
      cl.strand = (uint8_t)st;
      cl.beg = (uint16_t)i;
      cl.end = (uint16_t)j;
      cl.pmask = pmask;
      i = j;
    }
  }
  return nclust;
}

// ws.corder = cluster indices ordered by member count, descending; ties keep
// index order (stable). Insertion sort for small n, counting sort for the
// large candidate sets of host passes; both give the same order.
CERTA_HD inline void order_by_support(Workspace& ws, int n) {
  auto support = [&](int x) { return ws.clusters[x].end - ws.clusters[x].beg; };
  if (n <= 64) {
    for (int x = 0; x < n; ++x) ws.corder[x] = (uint16_t)x;
    for (int x = 1; x < n; ++x) {
      uint16_t v = ws.corder[x];
      int sv = support(v), y = x - 1;
      while (y >= 0 && support(ws.corder[y]) < sv) {
        ws.corder[y + 1] = ws.corder[y];
        --y;
      }
      ws.corder[y + 1] = v;
    }
    return;
  }
  // Counting sort on support in [1, MAX_CAND]: start offsets for each count,
  // highest count first, then place clusters in index order (stable).
  int maxs = 0;
  for (int x = 0; x < n; ++x) if (support(x) > maxs) maxs = support(x);
  int* start = ws.count_buf;  // size MAX_CAND + 2
  for (int s = 0; s <= maxs + 1; ++s) start[s] = 0;
  for (int x = 0; x < n; ++x) ++start[support(x)];
  int at = 0;
  for (int s = maxs; s >= 0; --s) {
    const int c = start[s];
    start[s] = at;
    at += c;
  }
  for (int x = 0; x < n; ++x) ws.corder[start[support(x)]++] = (uint16_t)x;
}

// Verify clusters, best-supported first, tracking the two smallest distances
// b1 <= b2 (initially cap + 1). Only a distance below b2 can change (b1, b2),
// so each cluster is evaluated with limit b2 - 1: b1 is exact, and when
// b1 < b2 the second-best distance b2 (<= cap) is exact too. Once two loci tie
// at b1 the rest are skipped, so the tie count becomes a lower bound (MAPQ 0
// either way). A Hamming check at the member diagonals gives an upper bound
// first, so exact matches never reach the DP. Cluster.dist is exact or cap + 1.
CERTA_HD inline void evaluate_clusters(const IndexView& ix, int L, int nclust, int cap,
                                       int reverse, int max_eval, int min_support, Workspace& ws,
                                       int* pb1, int* pb2, int lb_parts = 0, int lb_per = 0) {
  order_by_support(ws, nclust);
  if (reverse)  // results must not depend on the order (tested)
    for (int x = 0, y = nclust - 1; x < y; ++x, --y) {
      uint16_t t = ws.corder[x];
      ws.corder[x] = ws.corder[y];
      ws.corder[y] = t;
    }
  int b1 = cap + 1, b2 = cap + 1;
#ifndef __CUDA_ARCH__
  // Each cluster's first reference access is a random read of a 3 GB
  // reference; start those of the clusters a few places ahead (results do
  // not change).
  // Only for the many clusters of host passes (pass 3, pass V).
  constexpr int kAhead = 6;
  const bool pf = nclust > 16;
  auto prefetch_window = [&](int y) {
    const int64_t lo = ws.clusters[ws.corder[y]].lo;
    if (lo < 0 || (uint64_t)lo + L + BMAX > ix.ref_len) return;
    for (int b = 0; b < L + 64; b += 64) __builtin_prefetch(ix.ref + lo + b);
  };
  if (pf)
    for (int y = 0; y < kAhead && y < nclust; ++y) prefetch_window(y);
#endif
  for (int x = 0; x < nclust; ++x) {
#ifndef __CUDA_ARCH__
    if (pf && x + kAhead < nclust) prefetch_window(x + kAhead);
#endif
    Cluster& cl = ws.clusters[ws.corder[x]];
    cl.dist = (int8_t)(cap + 1);  // "cannot improve (b1, b2)", or not evaluated
    const int limit = b2 - 1;
    if (limit < 0 || x >= max_eval) continue;
    if (popcount16(cl.pmask) < min_support) continue;  // q-gram lemma: no locus within R here
    if (lb_per * (lb_parts - popcount16(cl.pmask)) > limit) continue;  // pass V: none below b2 here
    const uint8_t* rd = ws.seq[cl.strand];
    int u = limit + 1;  // Hamming upper bound, only tracked below limit + 1
    for (int e = cl.beg; e < cl.end && u > 0; ++e) {
      int h = hamming(ix, rd, L, diag_of(ws.cand[cl.strand][e]), u - 1);
      if (h < u) u = h;
    }
    int d;
    if (u == 0) {
      d = 0;
    } else if (u == 1) {
      // The band's distance is <= 1 (or, with limit 0, only whether it is 0
      // matters), and it is 0 only on a diagonal without a mismatch (no
      // edits: an ungapped path). So a Hamming check per diagonal gives what
      // the DP below would (tested in test_certa).
      d = 1;
      for (int b = 0; b < cl.width && d > 0; ++b)
        if (hamming(ix, rd, L, cl.lo + b, 0) == 0) d = 0;
      if (d > limit) continue;  // distance >= b2: irrelevant
    } else {
      const int lim = u <= limit ? u : limit;
      d = band_distance(ix, rd, L, cl.lo, cl.width, lim, ws.row_a, ws.row_b);
      if (d > lim) continue;  // distance >= b2: irrelevant
    }
    cl.dist = (int8_t)d;
    if (d < b1) { b2 = b1; b1 = d; } else if (d < b2) { b2 = d; }
  }
  *pb1 = b1;
  *pb2 = b2;
}

// The DP fill of affine_align (below): Gotoh recurrences on the same band
// indexing as band_distance, cell (i, b) aligning read base i-1 to reference
// i-1+lo+b. Writes ws.tb[i * B + b] for every cell and returns the best end
// score (NEG if there is none), with its cell in (*pbi, *pbb), the first in
// row-major order.
// tb bits: 0-1 H source (0 match/mismatch, 1 insertion E, 2 deletion F),
// 2 the diagonal step starts the alignment (read prefix clipped or none),
// 3 E extends E, 4 F extends F.
CERTA_HD inline int affine_fill_scalar(const IndexView& ix, const uint8_t* rd, int L, int64_t lo, int B,
                                       Workspace& ws, int64_t smin, int64_t smax, int* pbi, int* pbb) {
  const int NEG = -(1 << 20);
  int* hp = ws.aff[0];
  int* hc = ws.aff[1];
  int* ep = ws.aff[2];
  int* ec = ws.aff[3];
  for (int b = 0; b < B; ++b) { hp[b] = NEG; ep[b] = NEG; }
  int best = NEG, best_i = 0, best_b = 0;
  // The band's reference window, read directly when it lies inside the
  // reference (ref_at's bounds check otherwise).
  const bool inside = lo >= 0 && lo + L - 1 + B <= (int64_t)ix.ref_len;
  for (int i = 1; i <= L; ++i) {
    const uint8_t c = rd[i - 1];
    const int start = i == 1 ? 0 : -kClip;
    const int clip_end = i < L ? kClip : 0;
    // Cells b in [blo, bhi] put their first aligned base in [smin, smax].
    const int64_t base = (int64_t)i - 1 + lo;
    const int64_t blo64 = smin == INT64_MIN ? 0 : smin - base, bhi64 = smax == INT64_MAX ? B : smax - base;
    const int blo = blo64 < 0 ? 0 : (blo64 > B ? B : (int)blo64);
    const int bhi = bhi64 < -1 ? -1 : (bhi64 >= B ? B - 1 : (int)bhi64);
    const uint8_t* rw = inside ? ix.ref + base : nullptr;
    int f = NEG;
    for (int b = 0; b < B; ++b) {
      uint8_t bits = 0;
      int diag = hp[b];
      if (start >= diag && b >= blo && b <= bhi) { diag = start; bits |= 4; }
      const uint8_t t = rw ? rw[b] : ref_at(ix, base + b);
      const int m = diag + (sub_cost(c, t) ? -kMismatch : kMatch);
      int e = NEG;
      if (b + 1 < B && i > 1) {  // no alignment may start with a gap
        const int open = hp[b + 1] - kGapOpen - kGapExt, ext = ep[b + 1] - kGapExt;
        e = open >= ext ? open : ext;
        if (ext > open) bits |= 8;
      }
      int fo = NEG;
      if (b > 0) {
        const int open = hc[b - 1] - kGapOpen - kGapExt, ext = f - kGapExt;
        fo = open >= ext ? open : ext;
        if (ext > open) bits |= 16;
      }
      f = fo;
      int h = m;
      if (e > h) { h = e; bits = (uint8_t)((bits & ~3) | 1); }
      if (f > h) { h = f; bits = (uint8_t)((bits & ~3) | 2); }
      hc[b] = h;
      ec[b] = e;
      ws.tb[i * B + b] = bits;
      // An alignment may end after a match/mismatch (a select, not a branch
      // on the cell's state; the first best cell still wins).
      const int end = (bits & 3) == 0 ? h - clip_end : NEG;
      if (end > best) { best = end; best_i = i; best_b = b; }
    }
    int* t = hp; hp = hc; hc = t;
    t = ep; ep = ec; ec = t;
  }
  *pbi = best_i;
  *pbb = best_b;
  return best;
}

#if defined(__GNUC__) && !defined(__CUDACC__) && !defined(CERTA_NO_VEC) && (defined(__clang__) || __GNUC__ >= 12)
#define CERTA_AFFINE_VEC 1
// The same fill, 4 cells at a time (32-bit lanes; NEON or SSE), for a band
// whose padded reference window lies inside the reference. Every cell gets
// the scalar loop's values and traceback bits, sentinels included:
// - E and the match term depend on the previous row only.
// - F, the deletion state, satisfies F[b] = max(hc[b-1] - 7, F[b-1] - 1)
//   with hc = max(h0, F) and h0 = max(m, E); opening from an F cell is
//   dominated by extending it, so F[b] = max(h0[b-1] - 7, F[b-1] - 1), and
//   unrolled with F[0] = NEG: F[b] = max(NEG + 6, max over b' < b of
//   h0[b'] + b') - b - 6, an exclusive prefix maximum (two in-register shifts
//   per vector and a carry across vectors; as band_row in longread.h).
// - "F extends F" at b is F[b-1] - 1 > hc[b-1] - 7, from the previous lane.
// The best end cell is taken lane by lane, in order, so ties resolve as in
// the scalar loop.
inline int affine_fill_vec(const IndexView& ix, const uint8_t* rd, int L, int64_t lo, int B, Workspace& ws,
                           int64_t smin, int64_t smax, int* pbi, int* pbb) {
  static_assert(kMatch == 1 && kGapExt == 1, "affine_fill_vec assumes +1 matches and gap extension 1");
  static_assert(__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__, "affine_fill_vec packs the low byte of each lane");
  typedef int32_t v4 __attribute__((vector_size(16)));
  typedef uint8_t u8x16 __attribute__((vector_size(16)));
  typedef uint64_t u64x2 __attribute__((vector_size(16)));
  typedef uint32_t u32x4 __attribute__((vector_size(16)));
  // Unaligned loads and stores as plain accesses (memcpy at a variable
  // offset becomes a checked library call under _FORTIFY_SOURCE).
  typedef int32_t v4u __attribute__((vector_size(16), aligned(4), may_alias));
  typedef uint32_t u32u __attribute__((aligned(1), may_alias));
  constexpr int W = 4;
  const int NEG = -(1 << 20);
  auto vmax = [](v4 a, v4 b) { return a > b ? a : b; };
  auto splat = [](int x) { return v4{x, x, x, x}; };
  const v4 vneg = splat(NEG), lane = {0, 1, 2, 3};
  const int Bp = (B + W - 1) / W * W;
  int32_t hb[2][BMAX + 2 * W], eb[2][BMAX + 2 * W];
  int32_t* hp = hb[0];
  int32_t* hc = hb[1];
  int32_t* ep = eb[0];
  int32_t* ec = eb[1];
  for (int b = 0; b < Bp + W; ++b) hp[b] = hc[b] = ep[b] = ec[b] = NEG;
  // The band's reference window as 32-bit lanes, once: row i reads it from i - 1.
  int32_t refw[LMAX + BMAX + W];
  for (int k = 0; k < L - 1 + Bp; ++k) refw[k] = ix.ref[lo + k];
  int best = NEG, best_i = 0, best_b = 0;
  for (int i = 1; i <= L; ++i) {
    const v4 c = splat(rd[i - 1]);
    const v4 start = splat(i == 1 ? 0 : -kClip);
    const int clip_end = i < L ? kClip : 0;
    const int64_t base = (int64_t)i - 1 + lo;
    const int64_t blo64 = smin == INT64_MIN ? 0 : smin - base, bhi64 = smax == INT64_MAX ? B : smax - base;
    const int blo = blo64 < 0 ? 0 : (blo64 > B ? B : (int)blo64);
    const int bhi = bhi64 < -1 ? -1 : (bhi64 >= B ? B - 1 : (int)bhi64);
    const int32_t* rw = refw + (i - 1);
    v4 carry = splat(NEG + 6);  // running max of h0[b'] + b'
    v4 fprev = vneg, hprev = vneg;  // the previous vector's F and H (lane 3 used)
    for (int b = 0; b < Bp; b += W) {
      const v4 bv = lane + splat(b);
      const v4 hpv = *(const v4u*)(hp + b), hp1 = *(const v4u*)(hp + b + 1), ep1 = *(const v4u*)(ep + b + 1);
      const v4 st = (start >= hpv) & (bv >= splat(blo)) & (bv <= splat(bhi));
      const v4 diag = st ? start : hpv;
      const v4 t = *(const v4u*)(rw + b);
      const v4 mis = (t != c) | (c >= splat(4));
      const v4 m = diag + (mis ? splat(-kMismatch) : splat(kMatch));
      v4 e = vneg, ebit = splat(0);
      if (i > 1) {
        const v4 open = hp1 - splat(kGapOpen + kGapExt), ext = ep1 - splat(kGapExt);
        const v4 valid = bv + splat(1) < splat(B);
        e = valid ? (open >= ext ? open : ext) : vneg;
        ebit = (ext > open) & valid & splat(8);
      }
      const v4 esel = e > m;
      const v4 h0 = esel ? e : m;
      v4 src = esel & splat(1);
      // F from the exclusive prefix maximum of h0[b] + b.
      const v4 ninf = splat(INT32_MIN);
      const v4 g = h0 + bv;
      v4 inc = vmax(g, __builtin_shufflevector(g, ninf, 4, 0, 1, 2));
      inc = vmax(inc, __builtin_shufflevector(inc, ninf, 4, 5, 0, 1));
      inc = vmax(inc, carry);
      const v4 excl = __builtin_shufflevector(carry, inc, 3, 4, 5, 6);
      carry = __builtin_shufflevector(inc, inc, 3, 3, 3, 3);
      const v4 f = bv > splat(0) ? excl - bv - splat(kGapOpen + kGapExt - 1) : vneg;
      const v4 fsel = f > h0;
      const v4 h = fsel ? f : h0;
      src = fsel ? splat(2) : src;
      const v4 fm1 = __builtin_shufflevector(fprev, f, 3, 4, 5, 6);
      const v4 hm1 = __builtin_shufflevector(hprev, h, 3, 4, 5, 6);
      const v4 fbit = (fm1 + splat(kGapOpen + kGapExt - 1) > hm1) & (bv > splat(0)) & splat(16);
      fprev = f;
      hprev = h;
      *(v4u*)(hc + b) = h;
      *(v4u*)(ec + b) = e;
      const v4 bits = (st & splat(4)) | ebit | fbit | src;
      // The low byte of each lane (little-endian), packed into the first four
      // and stored whole: past the row's end it lands on cells of row i + 1,
      // written later, or, after row L, past the cells the traceback reads
      // (tb holds (LMAX + 1) * BMAX bytes and B < BMAX whenever B is not a
      // multiple of 4).
      const u8x16 bytes = (u8x16)bits;
      const u8x16 b8 = __builtin_shufflevector(bytes, bytes, 0, 4, 8, 12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
      *(u32u*)(ws.tb + (size_t)i * B + b) = ((u32x4)b8)[0];
      // An alignment may end after a match/mismatch.
      const v4 end = (src == splat(0)) & (bv < splat(B)) ? h - splat(clip_end) : vneg;
      const u64x2 over = (u64x2)(end > splat(best));
      if (over[0] | over[1])
        for (int l = 0; l < W; ++l)
          if (end[l] > best) { best = end[l]; best_i = i; best_b = b + l; }
    }
    int32_t* tmp = hp; hp = hc; hc = tmp;
    tmp = ep; ep = ec; ec = tmp;
  }
  *pbi = best_i;
  *pbb = best_b;
  return best;
}

// affine_fill_vec with 16-bit lanes, 8 cells at a time, and the same values
// relative to "minus infinity": NEG is -16000 here instead of -(1 << 20).
// Every value the fill compares is either a real score (above -9 LMAX - 5: a
// cell loses at most 9 points per row) or NEG plus a bounded offset (below
// NEG + LMAX, and above NEG - 9 LMAX - BMAX, so no lane wraps), and the same operations give the same offsets from either
// NEG; so every comparison, traceback bit and best cell is the same. A best
// derived from NEG (no allowed start) is mapped back to the 32-bit NEG plus
// its offset, as the 32-bit fill returns it.
inline int affine_fill_vec16(const IndexView& ix, const uint8_t* rd, int L, int64_t lo, int B, Workspace& ws,
                             int64_t smin, int64_t smax, int* pbi, int* pbb) {
  static_assert(kMatch == 1 && kGapExt == 1, "affine_fill_vec16 assumes +1 matches and gap extension 1");
  static_assert(__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__, "affine_fill_vec16 packs the low byte of each lane");
  static_assert(-16000 - 9 * LMAX - BMAX > INT16_MIN && LMAX * 9 + 5 < 16000 / 2, "16-bit affine fill");
  typedef int16_t v8 __attribute__((vector_size(16)));
  typedef uint8_t u8x16 __attribute__((vector_size(16)));
  typedef uint64_t u64x2 __attribute__((vector_size(16)));
  typedef int16_t v8u __attribute__((vector_size(16), aligned(2), may_alias));
  typedef uint8_t u8x8u __attribute__((vector_size(8), aligned(1), may_alias));
  typedef uint64_t u64u __attribute__((aligned(1), may_alias));
  constexpr int W = 8;
  const int NEG = -16000;
  auto vmax = [](v8 a, v8 b) { return a > b ? a : b; };
  auto splat = [](int x) {
    const int16_t s = (int16_t)x;
    return v8{s, s, s, s, s, s, s, s};
  };
  const v8 vneg = splat(NEG), lane = {0, 1, 2, 3, 4, 5, 6, 7};
  const v8 ninf = splat(INT16_MIN);
  const int Bp = (B + W - 1) / W * W;
  int16_t hb[2][BMAX + 2 * W], eb[2][BMAX + 2 * W];
  int16_t* hp = hb[0];
  int16_t* hc = hb[1];
  int16_t* ep = eb[0];
  int16_t* ec = eb[1];
  for (int b = 0; b < Bp + W; ++b) hp[b] = hc[b] = ep[b] = ec[b] = (int16_t)NEG;
  const uint8_t* refw = ix.ref + lo;  // row i reads it from i - 1 (inside the reference: see the caller)
  int best = NEG, best_i = 0, best_b = 0;
  for (int i = 1; i <= L; ++i) {
    const v8 c = splat(rd[i - 1]);
    const v8 start = splat(i == 1 ? 0 : -kClip);
    const int clip_end = i < L ? kClip : 0;
    const int64_t base = (int64_t)i - 1 + lo;
    const int64_t blo64 = smin == INT64_MIN ? 0 : smin - base, bhi64 = smax == INT64_MAX ? B : smax - base;
    const int blo = blo64 < 0 ? 0 : (blo64 > B ? B : (int)blo64);
    const int bhi = bhi64 < -1 ? -1 : (bhi64 >= B ? B - 1 : (int)bhi64);
    const uint8_t* rw = refw + (i - 1);
    v8 carry = splat(NEG + 6);  // running max of h0[b'] + b'
    v8 fprev = vneg, hprev = vneg;  // the previous vector's F and H (lane 7 used)
    for (int b = 0; b < Bp; b += W) {
      const v8 bv = lane + splat(b);
      const v8 hpv = *(const v8u*)(hp + b), hp1 = *(const v8u*)(hp + b + 1), ep1 = *(const v8u*)(ep + b + 1);
      const v8 st = (start >= hpv) & (bv >= splat(blo)) & (bv <= splat(bhi));
      const v8 diag = st ? start : hpv;
      const v8 t = __builtin_convertvector(*(const u8x8u*)(rw + b), v8);
      const v8 mis = (t != c) | (c >= splat(4));
      const v8 m = diag + (mis ? splat(-kMismatch) : splat(kMatch));
      v8 e = vneg, ebit = splat(0);
      if (i > 1) {
        const v8 open = hp1 - splat(kGapOpen + kGapExt), ext = ep1 - splat(kGapExt);
        const v8 valid = bv + splat(1) < splat(B);
        e = valid ? (open >= ext ? open : ext) : vneg;
        ebit = (ext > open) & valid & splat(8);
      }
      const v8 esel = e > m;
      const v8 h0 = esel ? e : m;
      v8 src = esel & splat(1);
      // F from the exclusive prefix maximum of h0[b] + b.
      const v8 g = h0 + bv;
      v8 inc = vmax(g, __builtin_shufflevector(g, ninf, 8, 0, 1, 2, 3, 4, 5, 6));
      inc = vmax(inc, __builtin_shufflevector(inc, ninf, 8, 9, 0, 1, 2, 3, 4, 5));
      inc = vmax(inc, __builtin_shufflevector(inc, ninf, 8, 9, 10, 11, 0, 1, 2, 3));
      inc = vmax(inc, carry);
      const v8 excl = __builtin_shufflevector(carry, inc, 7, 8, 9, 10, 11, 12, 13, 14);
      carry = __builtin_shufflevector(inc, inc, 7, 7, 7, 7, 7, 7, 7, 7);
      const v8 f = bv > splat(0) ? excl - bv - splat(kGapOpen + kGapExt - 1) : vneg;
      const v8 fsel = f > h0;
      const v8 h = fsel ? f : h0;
      src = fsel ? splat(2) : src;
      const v8 fm1 = __builtin_shufflevector(fprev, f, 7, 8, 9, 10, 11, 12, 13, 14);
      const v8 hm1 = __builtin_shufflevector(hprev, h, 7, 8, 9, 10, 11, 12, 13, 14);
      const v8 fbit = (fm1 + splat(kGapOpen + kGapExt - 1) > hm1) & (bv > splat(0)) & splat(16);
      fprev = f;
      hprev = h;
      *(v8u*)(hc + b) = h;
      *(v8u*)(ec + b) = e;
      const v8 bits = (st & splat(4)) | ebit | fbit | src;
      // The low byte of each lane, packed into 8 bytes and stored whole (past
      // the row's end it lands on cells of row i + 1, written later, or after
      // row L inside tb's (LMAX + 1) * BMAX bytes, past what traceback reads).
      const u8x16 bytes = (u8x16)bits;
      const u8x16 b8 = __builtin_shufflevector(bytes, bytes, 0, 2, 4, 6, 8, 10, 12, 14, 0, 0, 0, 0, 0, 0, 0, 0);
      *(u64u*)(ws.tb + (size_t)i * B + b) = ((u64x2)b8)[0];
      // An alignment may end after a match/mismatch.
      const v8 end = (src == splat(0)) & (bv < splat(B)) ? h - splat(clip_end) : vneg;
      const u64x2 over = (u64x2)(end > splat(best));
      if (over[0] | over[1])
        for (int l = 0; l < W; ++l)
          if (end[l] > best) { best = end[l]; best_i = i; best_b = b + l; }
    }
    int16_t* tmp = hp; hp = hc; hc = tmp;
    tmp = ep; ep = ec; ec = tmp;
  }
  *pbi = best_i;
  *pbb = best_b;
  return best < NEG / 2 ? best - NEG - (1 << 20) : best;
}
#endif

// Affine-optimal alignment of the read within diagonals [lo, lo + B) with
// BWA-MEM's default scores (+1 / -4 / gap -6 - 1 per base / clip -5 per end),
// reference ends free. Writes CIGAR (with soft clips), ref_pos, nm and score
// into `r`; returns false (leaving `r` unchanged) if the CIGAR would not fit.
// Optionally only alignments whose first aligned reference base lies in
// [smin, smax] count (false if there is none).
CERTA_HD inline bool affine_align(const IndexView& ix, const uint8_t* rd, int L,
                                  int64_t lo, int B, Workspace& ws, Result& r,
                                  bool score_only = false, int64_t smin = INT64_MIN,
                                  int64_t smax = INT64_MAX) {
  const int NEG = -(1 << 20);
  int best_i = 0, best_b = 0, best;
#ifdef CERTA_AFFINE_VEC
  if (lo >= 0 && lo + L - 1 + (B + 7) / 8 * 8 <= (int64_t)ix.ref_len)
    best = affine_fill_vec16(ix, rd, L, lo, B, ws, smin, smax, &best_i, &best_b);
  else if (lo >= 0 && lo + L - 1 + (B + 3) / 4 * 4 <= (int64_t)ix.ref_len)
    best = affine_fill_vec(ix, rd, L, lo, B, ws, smin, smax, &best_i, &best_b);
  else
#endif
    best = affine_fill_scalar(ix, rd, L, lo, B, ws, smin, smax, &best_i, &best_b);
  if (best < NEG / 2) return false;  // no allowed start
  if (score_only) {
    r.score = (int16_t)best;
    return true;
  }
  // Traceback from (best_i, best_b) in state H.
  uint32_t rev[MAX_CIGAR];
  int nrev = 0, nm = 0, i = best_i, b = best_b, state = 0;  // 0 H, 1 E, 2 F
  uint32_t run_op = 99, run_len = 0;
  auto emit = [&](uint32_t op) -> bool {
    if (op == run_op) { ++run_len; return true; }
    if (run_len) {
      if (nrev == MAX_CIGAR) return false;
      rev[nrev++] = (run_len << 4) | run_op;
    }
    run_op = op;
    run_len = 1;
    return true;
  };
  int first_i = 0, first_b = 0;
  for (;;) {
    const uint8_t bits = ws.tb[i * B + b];
    if (state == 0) state = bits & 3;
    if (state == 0) {  // match or mismatch
      if (!emit(kOpM)) return false;
      nm += sub_cost(rd[i - 1], ref_at(ix, (int64_t)i - 1 + lo + b));
      if (bits & 4) { first_i = i; first_b = b; break; }
      --i;
    } else if (state == 1) {  // insertion: read base without reference
      if (!emit(kOpI)) return false;
      ++nm;
      state = (bits & 8) ? 1 : 0;
      --i;
      ++b;
    } else {  // deletion: reference base without read
      if (!emit(kOpD)) return false;
      ++nm;
      state = (bits & 16) ? 2 : 0;
      --b;
    }
  }
  if (run_len) {
    if (nrev == MAX_CIGAR) return false;
    rev[nrev++] = (run_len << 4) | run_op;
  }
  const int lead = first_i - 1, trail = L - best_i;
  if (nrev + (lead > 0) + (trail > 0) > MAX_CIGAR) return false;
  r.n_cigar = 0;
  if (lead > 0) r.cigar[r.n_cigar++] = ((uint32_t)lead << 4) | kOpS;
  for (int x = nrev - 1; x >= 0; --x) r.cigar[r.n_cigar++] = rev[x];
  if (trail > 0) r.cigar[r.n_cigar++] = ((uint32_t)trail << 4) | kOpS;
  r.ref_pos = (int64_t)first_i - 1 + lo + first_b;
  r.nm = (uint8_t)nm;
  r.score = (int16_t)best;
  return true;
}

// Chooses among the clusters at distance d1 (deterministically, seeded by the
// read name), aligns the read there and returns how many clusters tie.
// *width receives the chosen cluster's band width.
CERTA_HD inline int align_best(const IndexView& ix, int L, int nclust, int d1,
                               uint64_t name_hash, Workspace& ws, Result& r, int* width) {
  int n_best = 0;
  for (int x = 0; x < nclust; ++x) n_best += ws.clusters[x].dist == d1;
  int pick = (int)(name_hash % (uint64_t)n_best), chosen = -1;
  for (int x = 0; x < nclust; ++x)
    if (ws.clusters[x].dist == d1 && pick-- == 0) { chosen = x; break; }
  const Cluster& cl = ws.clusters[chosen];
  *width = cl.width;
  const uint8_t* rd = ws.seq[cl.strand];
  // Prefer an ungapped alignment when it is optimal.
  int64_t ungapped = -1;
  for (int b = 0; b < cl.width && ungapped < 0; ++b)
    if (hamming(ix, rd, L, cl.lo + b, d1) == d1) ungapped = cl.lo + b;
  if (ungapped >= 0) {
    r.ref_pos = ungapped;
    r.n_cigar = 1;
    r.cigar[0] = ((uint32_t)L << 4) | kOpM;
  } else {
    band_traceback(ix, rd, L, cl.lo, cl.width, ws, r);
  }
  r.nm = (uint8_t)d1;
  r.score = (int16_t)(L * kMatch - d1 * (kMatch + kMismatch));  // exact for ungapped
  // Report the locus with BWA-MEM's affine scoring (clipping, gap placement),
  // so callers see the same alignment representation as from bwa-mem. The
  // locus and d1 (the certified quantity) do not change.
  // With d1 = 1 the ungapped alignment found above (on the first diagonal with
  // one mismatch) is already what affine_align reports: d1 is the band's
  // minimum edit distance, so it scores L - 5, while a gapped alignment scores
  // <= L - 7, a clipped one <= L - 6 and an ungapped one with more mismatches
  // less. The fill's first best end cell is then in row L on that diagonal,
  // and its traceback is the ungapped path (tested in test_certa).
  if (d1 > 0 && !(d1 == 1 && ungapped >= 0)) {
    Result tmp = r;
    if (affine_align(ix, rd, L, cl.lo, cl.width, ws, tmp)) {
      for (int x = 0; x < tmp.n_cigar; ++x) r.cigar[x] = tmp.cigar[x];
      r.n_cigar = tmp.n_cigar;
      r.ref_pos = tmp.ref_pos;
      r.nm = tmp.nm;
      r.score = tmp.score;
    }
  }
  r.strand = cl.strand;
  return n_best;
}

// MAPQ only: affine score of the read aligned in cluster `c`'s band.
CERTA_HD inline int16_t cluster_score(const IndexView& ix, int L, const Cluster& c, Workspace& ws) {
  Result tmp;
  tmp.n_cigar = 0;
  if (!affine_align(ix, ws.seq[c.strand], L, c.lo, c.width, ws, tmp)) return kNoSub;
  return tmp.score;
}

// MAPQ only (not part of the certificate): the best distance, up to `limit`
// edits, among candidate clusters other than the reported locus, or -1; its
// affine score goes to *sub. Evaluates the S2_TOP best-supported clusters;
// overwrites ws.clusters.
CERTA_HD inline int other_candidates(const IndexView& ix, int L, int limit, const Result& r,
                                     Workspace& ws, int16_t* sub) {
  const int M = limit < S2_MAX ? limit : S2_MAX;
  const int nc = build_clusters(ws, 2 * M, true);
  order_by_support(ws, nc);
  int best = M + 1, evaluated = 0, best_c = -1;
  for (int x = 0; x < nc && evaluated < S2_TOP && best > 0; ++x) {
    const Cluster& cl = ws.clusters[ws.corder[x]];
    if (cl.strand == r.strand && r.ref_pos >= cl.lo - L && r.ref_pos <= cl.lo + cl.width + L)
      continue;  // the reported locus itself
    ++evaluated;
    const int d = band_distance(ix, ws.seq[cl.strand], L, cl.lo, cl.width, best - 1, ws.row_a, ws.row_b);
    if (d < best) { best = d; best_c = ws.corder[x]; }
  }
  if (best > M) return -1;
  *sub = cluster_score(ix, L, ws.clusters[best_c], ws);
  return best;
}

// Tier SR: sample hits of the rarest part on each strand; >= 2 distinct exact
// full-read matches prove d1 = 0 and a multi-mapping read.
CERTA_HD inline bool certify_repeat(const IndexView& ix, int L, int P,
                                    uint64_t name_hash, Workspace& ws,
                                    Result& r) {
  int64_t diag[REP_KEEP];
  uint8_t strand[REP_KEEP];
  int found = 0;
  for (int st = 0; st < 2; ++st) {
    const int j = ws.order[st][0];
    if (ws.count[st][j] >= (uint32_t)BIG) continue;  // every part contains 'N'
    int off, len, examined = 0;
    part_geometry(L, P, j, &off, &len);
    for (int t = 0; t < ix.s; ++t) {
      for (uint64_t e = ws.rlo[st][j][t]; e < ws.rhi[st][j][t]; ++e) {
        if (examined++ >= REP_SAMPLE || found >= REP_KEEP) break;
        int64_t d = (int64_t)ix.pos[e] - (off + t);
        if (hamming(ix, ws.seq[st], L, d, 0) != 0) continue;
        bool dup = false;
        for (int x = 0; x < found; ++x) dup |= diag[x] == d && strand[x] == st;
        if (!dup) { diag[found] = d; strand[found] = (uint8_t)st; ++found; }
      }
    }
  }
  if (found < 2) return false;
  const int pick = (int)(name_hash % (uint64_t)found);
  r.ref_pos = diag[pick];
  r.n_cigar = 1;
  r.cigar[0] = ((uint32_t)L << 4) | kOpM;
  r.certified = 1;
  r.tier = kTierSR;
  r.strand = strand[pick];
  r.d1 = 0;
  r.nm = 0;
  r.score = (int16_t)(L * kMatch);
  r.d2 = 0;
  r.n_best = (uint16_t)found;
  r.sub_hi = (int16_t)(L * kMatch);  // a second exact copy: SUB = OPT = L
  return true;
}

#ifndef __CUDA_ARCH__
// ---------------------------------------------------------------------------
// Tier SL: certified local alignment (host only).
//
// Objective: the score affine_align maximizes, i.e. BWA-MEM's defaults with
// soft clipping: +1 match, -4 mismatch (N included), gap -(6 + k), -5 per
// clipped read end. For an alignment A of the read span [x, y), its loss is
// (y - x) - score(A) - clip penalties >= 0: a match costs nothing, a
// mismatch 5, an insertion of k bases 6 + 2k (its bases also lose their +1),
// a deletion of k bases 6 + k.
//
// Lemma L. Every part of S lying inside [x, y) that does not occur exactly
// in A adds >= 5 to the loss, disjointly. (Such a part contains a mismatch,
// costing 5, or part of a gap. A deletion lies between two read bases, so it
// touches at most one part, for >= 7. An insertion touching j >= 2 parts
// covers the j - 2 between them and a base of each end part, so k >= 2 +
// (j - 2) * minimum part length, and 6 + 2k >= 5j.)
//
// So an alignment in which at most `extra` parts of S occur exactly scores
// at most local_bound(.., free_mask = 0, extra). An alignment whose exactly
// occurring parts all lie in free_mask scores at most local_bound(.., mask, 0).
//
// Theorem L (certified local optimum). Let floor = local_bound(0, t - 1), the
// most an alignment with < t exact parts of S can score, and
// pad = L - floor - 7. Group the sorted candidate diagonals into chains
// (consecutive diagonals <= pad apart), and let ub(C) =
// local_bound(pmask(C), 0). Then every alignment A scoring > floor:
//  - has a gap total G <= pad (its loss is >= 6 + G and < L - floor);
//  - has >= t exactly occurring parts of S, all enumerated, whose diagonals
//    lie within G of each other, so in one chain C; A stays within G of
//    each of them, so inside the band [d - pad, d + pad] of such a diagonal d;
//  - scores <= ub(C), since every part of S outside pmask(C) is touched.
// The chains with ub > floor are evaluated in any order, keeping the two best
// scores s1 >= s2 (initially floor) and skipping chains with ub <= s2. If
// s1 > floor, s1 is the maximum score over the whole reference: an unseen
// better alignment would have an evaluated chain (and be found) or a skipped
// one (and score <= s2 <= s1). The reported alignment attains s1.
constexpr uint64_t kUbMagic = 0x6365727461554221ull;

// Start a fresh set of Lemma L memo entries (see Workspace::ubmemo).
inline void ubmemo_begin(Workspace& ws) {
  if (ws.ub_magic != kUbMagic) {
    for (int st = 0; st < 2; ++st)
      for (int x = 0; x < (1 << PMAX); ++x) ws.ubmemo[st][x] = kNoSub;
    ws.ub_magic = kUbMagic;
    ws.ub_ntouched = 0;
  }
  for (int i = 0; i < ws.ub_ntouched; ++i) ws.ubmemo[ws.ub_touched[i] >> 16][ws.ub_touched[i] & 0xFFFF] = kNoSub;
  ws.ub_ntouched = 0;
}

inline int local_bound(const Workspace& ws, int st, int L, int P, int m, uint16_t free_mask,
                       int extra);

// The Lemma L bound of a chain's part mask, memoised for this call.
inline int ub_of_mask(Workspace& ws, int st, int L, int P, int m, uint16_t mask) {
  int16_t& u = ws.ubmemo[st][mask];
  if (u == kNoSub) {
    u = (int16_t)local_bound(ws, st, L, P, m, mask, 0);
    if (ws.ub_ntouched < MAX_CLUST) ws.ub_touched[ws.ub_ntouched++] = ((uint32_t)st << 16) | mask;
    else ws.ub_magic = 0;  // too many to track: full reset next time
  }
  return u;
}

inline int local_bound(const Workspace& ws, int st, int L, int P, int m, uint16_t free_mask,
                       int extra) {
  static_assert(kMatch == 1 && kMismatch == 4 && kGapOpen == 6 && kGapExt == 1 && kClip == 5,
                "Lemma L is stated for BWA-MEM's default scores");
  int off[PMAX], end[PMAX], xs[PMAX + 1], ys[PMAX + 1];
  bool counted[PMAX];
  int nx = 0, ny = 0;
  xs[nx++] = 0;
  ys[ny++] = L;
  for (int x = 0; x < m; ++x) {
    const int j = ws.order[st][x];
    int len;
    part_geometry(L, P, j, &off[x], &len);
    end[x] = off[x] + len;
    counted[x] = !((free_mask >> j) & 1);
    xs[nx++] = off[x] + 1;  // spans starting just inside a part exclude it
    ys[ny++] = end[x] - 1;
  }
  int best = -BIG;
  for (int a = 0; a < nx; ++a)
    for (int b = 0; b < ny; ++b) {
      const int x0 = xs[a], y0 = ys[b];
      if (y0 <= x0) continue;
      int n = 0;
      for (int x = 0; x < m; ++x) n += counted[x] && off[x] >= x0 && end[x] <= y0;
      n = n > extra ? n - extra : 0;
      const int v = (y0 - x0) * kMatch - (kMatch + kMismatch) * n - (x0 > 0 ? kClip : 0) -
                    (y0 < L ? kClip : 0);
      if (v > best) best = v;
    }
  return best;
}

// Chains: the sorted candidate diagonals of each strand, split where two
// consecutive ones are more than pad apart. Stored in ws.clusters with
// lo = first diagonal (width unused).
inline int build_chains(Workspace& ws, int pad) {
  int n = 0;
  for (int st = 0; st < 2; ++st) {
    const int64_t* c = ws.cand[st];
    const int nc = ws.ncand[st];
    int i = 0;
    while (i < nc) {
      int64_t dmax = diag_of(c[i]);
      uint16_t pmask = (uint16_t)(1u << part_of(c[i]));
      int j = i + 1;
      while (j < nc && diag_of(c[j]) - dmax <= pad) {
        dmax = diag_of(c[j]);
        pmask |= (uint16_t)(1u << part_of(c[j]));
        ++j;
      }
      Cluster& cl = ws.clusters[n++];
      cl.lo = diag_of(c[i]);
      cl.width = 0;
      cl.strand = (uint8_t)st;
      cl.beg = (uint16_t)i;
      cl.end = (uint16_t)j;
      cl.pmask = pmask;
      i = j;
    }
  }
  return n;
}

// Score-only twin of affine_align: the same objective, band and optional
// start range, returning the best score (or kNoSub when no start is allowed).
// No traceback; the reference window is copied once; 16-bit lanes (scores
// lie in [-5 L, L]; "minus infinity" is -16000 and stays above -32768 after
// a read's worth of gap penalties), so the row loops vectorize 8 cells per
// 128-bit register. Deletions: extending (-ext) always beats reopening
// (-open - ext), so F(b) = max over k < b of H'(k) - open - ext (b - k), with
// H' = max(M, E): a running prefix maximum of H'(k) + ext k, one register per
// row. Tested against affine_align on random bands.
inline void band_row16(int B, int blo, int bhi, int16_t start, int16_t eopen, int16_t NEG,
                       const int16_t* __restrict hp, const int16_t* __restrict ep,
                       const int16_t* __restrict sub, int16_t* __restrict m, int16_t* __restrict hq,
                       int16_t* __restrict en) {
  for (int b = 0; b < B; ++b) {
    const int16_t st = (b >= blo && b <= bhi) ? start : NEG;
    const int16_t d = hp[b] > st ? hp[b] : st;
    const int16_t mv = (int16_t)(d + sub[b]);
    m[b] = mv;
    const int16_t o = (int16_t)(hp[b + 1] - eopen), x = (int16_t)(ep[b + 1] - kGapExt);
    const int16_t e = o > x ? o : x;
    en[b] = e;
    hq[b] = mv > e ? mv : e;
  }
}

// The deletion pass of band_score's row: hn[b] = max(hq[b], F(b)), F(b) =
// max over b' < b of hq[b'] + b' (and NEG), minus kGapOpen + b (kGapExt = 1).
// The scalar reference:
inline void band_del16_scalar(int B, const int16_t* hq, int16_t* hn, int16_t NEG) {
  int pm = NEG;
  hn[0] = hq[0];
  for (int b = 1; b < B; ++b) {
    const int g = hq[b - 1] + kGapExt * (b - 1);
    pm = pm > g ? pm : g;
    const int f = pm - kGapOpen - kGapExt * b;
    hn[b] = (int16_t)(hq[b] > f ? hq[b] : f);
  }
}

#if defined(__GNUC__) && !defined(__CUDACC__) && !defined(CERTA_NO_VEC) && (defined(__clang__) || __GNUC__ >= 12)
// The same, 8 cells at a time (16-bit lanes; NEON or SSE). Taking b' <= b
// adds hq[b] - kGapOpen, which never beats hq[b], so F is an inclusive prefix
// maximum: three in-register shifts per vector and a carry across vectors
// (as band_row in longread.h; lean/Certa/BandRow.lean proves the identity).
// The values stay in int16: hq >= NEG = -16000, b < BMAX.
inline void band_del16(int B, const int16_t* hq, int16_t* hn, int16_t NEG) {
  static_assert(kGapExt == 1, "band_del16 assumes kGapExt = 1");
  typedef int16_t v8 __attribute__((vector_size(16)));
  auto vmax = [](v8 a, v8 b) { return a > b ? a : b; };
  const v8 lane = {0, 1, 2, 3, 4, 5, 6, 7};
  v8 carry = {NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG};
  int b = 0;
  for (; b + 8 <= B; b += 8) {
    v8 h;
    std::memcpy(&h, hq + b, sizeof h);
    const v8 idx = lane + (int16_t)b;
    v8 x = h + idx;
    x = vmax(x, __builtin_shufflevector(x, x, 0, 0, 1, 2, 3, 4, 5, 6));
    x = vmax(x, __builtin_shufflevector(x, x, 0, 1, 0, 1, 2, 3, 4, 5));
    x = vmax(x, __builtin_shufflevector(x, x, 0, 1, 2, 3, 0, 1, 2, 3));
    x = vmax(x, carry);
    const v8 f = x - idx - (int16_t)kGapOpen;
    const v8 out = vmax(h, f);
    std::memcpy(hn + b, &out, sizeof out);
    carry = __builtin_shufflevector(x, x, 7, 7, 7, 7, 7, 7, 7, 7);
  }
  int pm = carry[0];
  for (; b < B; ++b) {
    const int g = hq[b] + b;
    pm = pm > g ? pm : g;
    const int f = pm - kGapOpen - b;
    hn[b] = (int16_t)(hq[b] > f ? hq[b] : f);
  }
}
#else
inline void band_del16(int B, const int16_t* hq, int16_t* hn, int16_t NEG) { band_del16_scalar(B, hq, hn, NEG); }
#endif

inline int band_score_scalar(const IndexView& ix, const uint8_t* rd, int L, int64_t lo, int B,
                             int64_t smin = INT64_MIN, int64_t smax = INT64_MAX) {
  static_assert(LMAX * (kMismatch + kGapOpen + kGapExt) < 16000, "16-bit band scores");
  const int16_t NEG = -16000;
  alignas(16) int16_t buf[4][BMAX + 8], m[BMAX + 8], sub[BMAX + 8], hq[BMAX + 8];
  int16_t *hp = buf[0], *hn = buf[1], *ep = buf[2], *en = buf[3];
  uint8_t refw[LMAX + BMAX];
  for (int j = 0; j < L - 1 + B; ++j) refw[j] = ref_at(ix, lo + j);
  for (int b = 0; b <= B; ++b) { hp[b] = NEG; ep[b] = NEG; }
  int best = NEG;
  for (int i = 1; i <= L; ++i) {
    const uint8_t c = rd[i - 1];
    const int16_t start = i == 1 ? 0 : -kClip;
    const int64_t base = (int64_t)i - 1 + lo;
    const int64_t blo64 = smin == INT64_MIN ? 0 : smin - base, bhi64 = smax == INT64_MAX ? B : smax - base;
    const int blo = blo64 < 0 ? 0 : (blo64 > B ? B : (int)blo64);
    const int bhi = bhi64 < -1 ? -1 : (bhi64 >= B ? B - 1 : (int)bhi64);
    const uint8_t* rw = refw + (i - 1);
    const int16_t eopen = i > 1 ? kGapOpen + kGapExt : 8000;  // no alignment starts with a gap
    const uint8_t cc = c < 4 ? c : 0xFF;                      // N never matches
    for (int b = 0; b < B; ++b) sub[b] = rw[b] == cc ? kMatch : -kMismatch;
    band_row16(B, blo, bhi, start, eopen, NEG, hp, ep, sub, m, hq, en);
    band_del16(B, hq, hn, NEG);
    int16_t rb = NEG;
    for (int b = 0; b < B; ++b) rb = rb > m[b] ? rb : m[b];
    const int clip = i < L ? kClip : 0;
    if (rb - clip > best) best = rb - clip;
    hn[B] = NEG;
    en[B] = NEG;
    int16_t* t = hp; hp = hn; hn = t;
    t = ep; ep = en; en = t;
  }
  return best < NEG / 2 ? kNoSub : best;
}

#if defined(__GNUC__) && !defined(__CUDACC__) && !defined(CERTA_NO_VEC) && (defined(__clang__) || __GNUC__ >= 12)
#define CERTA_BAND_SCORE_VEC 1
// band_score_scalar with each row in one pass, 8 cells at a time: the same
// int16 values in every cell (substitution, start, insertion, the deletion
// prefix maximum of band_del16, the row's best), so the same result. Lanes
// past the band compute values nothing reads: the insertion term of the
// band's last cell takes cell B, reset to minus infinity every row, the
// deletion maximum runs left to right, and the row's best masks them.
// With `theta`, only scores above theta need to be exact: the DP stops after
// row i once no alignment can still exceed it. One continuing from a cell
// gains at most a point per remaining base, so it ends at most H + (L - i);
// one starting later starts at -kClip, so it ends at most L - i - kClip.
// The result is then the exact score when that is above theta, and some
// value at or below theta (the best so far) otherwise.
inline int band_score_vec(const IndexView& ix, const uint8_t* rd, int L, int64_t lo, int B,
                          int64_t smin, int64_t smax, int theta) {
  static_assert(kGapExt == 1, "band_score_vec assumes kGapExt = 1");
  typedef int16_t v8 __attribute__((vector_size(16)));
  typedef int16_t v8u __attribute__((vector_size(16), aligned(2), may_alias));
  constexpr int W = 8;
  const int16_t NEG = -16000;
  auto vmax = [](v8 a, v8 b) { return a > b ? a : b; };
  auto splat = [](int x) {
    const int16_t s = (int16_t)x;
    return v8{s, s, s, s, s, s, s, s};
  };
  const v8 lane = {0, 1, 2, 3, 4, 5, 6, 7}, vneg = splat(NEG);
  const int Bp = (B + W - 1) / W * W;
  int16_t buf[4][BMAX + 2 * W];
  int16_t *hp = buf[0], *hn = buf[1], *ep = buf[2], *en = buf[3];
  for (int b = 0; b < Bp + W; ++b) hp[b] = hn[b] = ep[b] = en[b] = NEG;
  // The window's bases as 16-bit lanes, so no row widens them.
  int16_t refw[LMAX + BMAX + 2 * W];
  for (int j = 0; j < L - 1 + B; ++j) refw[j] = ref_at(ix, lo + j);
  for (int j = L - 1 + B; j < L - 1 + Bp + W; ++j) refw[j] = 4;
  int best = NEG;
  for (int i = 1; i <= L; ++i) {
    const uint8_t c = rd[i - 1];
    const int16_t start = i == 1 ? 0 : -kClip;
    const int64_t base = (int64_t)i - 1 + lo;
    const int64_t blo64 = smin == INT64_MIN ? 0 : smin - base, bhi64 = smax == INT64_MAX ? B : smax - base;
    const int blo = blo64 < 0 ? 0 : (blo64 > B ? B : (int)blo64);
    const int bhi = bhi64 < -1 ? -1 : (bhi64 >= B ? B - 1 : (int)bhi64);
    const int16_t* rw = refw + (i - 1);
    const v8 eopen = splat(i > 1 ? kGapOpen + kGapExt : 8000);  // no alignment starts with a gap
    const v8 cv = splat(c < 4 ? c : -1);                         // N never matches
    const v8 st0 = splat(start), vlo = splat(blo), vhi = splat(bhi), vB = splat(B);
    // A row that excludes no start gives every lane st0, lanes past the band
    // included: their values reach no cell inside it (see above).
    const bool open = blo == 0 && bhi >= B - 1;
    v8 carry = vneg, rbv = vneg, hbv = vneg;
    v8 bv = lane;
    for (int b = 0; b < Bp; b += W, bv += splat(W)) {
      const v8 p = *(const v8u*)(hp + b), pn = *(const v8u*)(hp + b + 1), epn = *(const v8u*)(ep + b + 1);
      const v8 r = *(const v8u*)(rw + b);
      const v8 sub = r == cv ? splat(kMatch) : splat(-kMismatch);
      const v8 st = open ? st0 : ((bv >= vlo) & (bv <= vhi) ? st0 : vneg);
      const v8 m = vmax(p, st) + sub;
      const v8 e = vmax(pn - eopen, epn - splat(kGapExt));
      const v8 h0 = vmax(m, e);
      v8 x = h0 + bv;
      x = vmax(x, __builtin_shufflevector(x, x, 0, 0, 1, 2, 3, 4, 5, 6));
      x = vmax(x, __builtin_shufflevector(x, x, 0, 1, 0, 1, 2, 3, 4, 5));
      x = vmax(x, __builtin_shufflevector(x, x, 0, 1, 2, 3, 0, 1, 2, 3));
      x = vmax(x, carry);
      carry = __builtin_shufflevector(x, x, 7, 7, 7, 7, 7, 7, 7, 7);
      const v8 h = vmax(h0, x - bv - splat(kGapOpen));
      *(v8u*)(hn + b) = h;
      *(v8u*)(en + b) = e;
      if (b + W <= B) {  // every lane inside the band
        rbv = vmax(rbv, m);
        hbv = vmax(hbv, h);
      } else {
        const v8 in = bv < vB;
        rbv = vmax(rbv, in ? m : vneg);
        hbv = vmax(hbv, in ? h : vneg);
      }
    }
    int16_t rb = NEG, hb = NEG;
    for (int l = 0; l < W; ++l) {
      rb = rb > rbv[l] ? rb : rbv[l];
      hb = hb > hbv[l] ? hb : hbv[l];
    }
    const int clip = i < L ? kClip : 0;
    if (rb - clip > best) best = rb - clip;
    if (hb + (L - i) <= theta && L - i - kClip <= theta) break;  // nothing can still exceed theta
    hn[B] = NEG;
    en[B] = NEG;
    int16_t* t = hp; hp = hn; hn = t;
    t = ep; ep = en; en = t;
  }
  return best < NEG / 2 ? kNoSub : best;
}
#endif

// theta (optional): only a score above theta must be exact; at or below it,
// any value <= theta may come back (band_score_vec stops early).
inline int band_score(const IndexView& ix, const uint8_t* rd, int L, int64_t lo, int B,
                      int64_t smin = INT64_MIN, int64_t smax = INT64_MAX, int theta = INT_MIN) {
#ifdef CERTA_BAND_SCORE_VEC
  return band_score_vec(ix, rd, L, lo, B, smin, smax, theta);
#else
  (void)theta;
  return band_score_scalar(ix, rd, L, lo, B, smin, smax);
#endif
}

// Best local score in a chain: its member diagonals are covered by bands
// [d - pad, d + pad], merged into bands of at most BMAX diagonals. Returns
// the score (or kNoSub) and the band that attains it. Optionally only
// alignments whose first aligned reference base lies in [smin, smax] count.
inline int chain_best(const IndexView& ix, int L, const Cluster& cl, int pad, Workspace& ws,
                      int64_t* best_lo, int* best_B, int64_t smin = INT64_MIN,
                      int64_t smax = INT64_MAX) {
  const int64_t* c = ws.cand[cl.strand];
  int best = kNoSub;
  int a = cl.beg;
  while (a < cl.end) {
    const int64_t d0 = diag_of(c[a]);
    int e = a + 1;
    while (e < cl.end && diag_of(c[e]) - d0 <= BMAX - 2 * pad - 1) ++e;
    const int64_t lo = d0 - pad;
    const int B = (int)(diag_of(c[e - 1]) - d0) + 2 * pad + 1;
    const int sc = band_score(ix, ws.seq[cl.strand], L, lo, B, smin, smax);
    if (sc > best) {
      best = sc;
      *best_lo = lo;
      *best_B = B;
    }
    a = e;
  }
  return best;
}

inline bool certify_local(const IndexView& ix, const Params& p, int L, int P, int m,
                          uint64_t name_hash, Workspace& ws, Result& r) {
  const int t = p.local;
  if (m < t) return false;
  int floor = -BIG;
  for (int st = 0; st < 2; ++st) {
    const int b = local_bound(ws, st, L, P, m, 0, t - 1);
    if (b > floor) floor = b;
  }
  if (floor >= L * kMatch) return false;
  int pad = L * kMatch - floor - 1 - kGapOpen;  // gap total G <= pad (kGapExt = 1)
  if (pad < 0) pad = 0;
  if (2 * pad + 1 > BMAX) return false;
  const int n = build_chains(ws, pad);
  ubmemo_begin(ws);
  int ne = 0;
  for (int x = 0; x < n; ++x) {
    const Cluster& cl = ws.clusters[x];
    const int ub = ub_of_mask(ws, cl.strand, L, P, m, cl.pmask);
    ws.lub[x] = (int16_t)ub;
    ws.lscore[x] = kNoSub;
    if (ub > floor) ws.corder[ne++] = (uint16_t)x;
  }
  // Most promising first (fewer evaluations); the result does not depend on it.
  std::sort(ws.corder, ws.corder + ne, [&](uint16_t a, uint16_t b) {
    return ws.lub[a] != ws.lub[b] ? ws.lub[a] > ws.lub[b] : a < b;
  });
  if (p.reverse_order) std::reverse(ws.corder, ws.corder + ne);
  int s1 = floor, s2 = floor, n_best = 0;
  for (int y = 0; y < ne; ++y) {
    const int x = ws.corder[y];
    if (ws.lub[x] <= s2) continue;
    int64_t lo;
    int B;
    const int s = chain_best(ix, L, ws.clusters[x], pad, ws, &lo, &B);
    ws.lscore[x] = (int16_t)s;
    if (s > s1) { s2 = s1; s1 = s; n_best = 1; }
    else if (s == s1 && s1 > floor) { ++n_best; s2 = s1; }
    else if (s > s2) { s2 = s; }
  }
  if (s1 <= floor) return false;
  // Report one of the chains attaining s1 (seeded by the read name).
  int pick = (int)(name_hash % (uint64_t)n_best), chosen = -1;
  for (int x = 0; x < n; ++x)
    if (ws.lscore[x] == s1 && pick-- == 0) { chosen = x; break; }
  int64_t lo = 0;
  int B = 0;
  chain_best(ix, L, ws.clusters[chosen], pad, ws, &lo, &B);
  Result tmp = r;
  if (!affine_align(ix, ws.seq[ws.clusters[chosen].strand], L, lo, B, ws, tmp) || tmp.score != s1)
    return false;  // CIGAR does not fit
  r = tmp;
  // SUB bound: every alignment outside the chosen chain scores <= s2 (an
  // evaluated chain, a skipped one with ub <= s2, or no chain: <= floor <=
  // s2). Inside it, one not overlapping [p, p + span) starts before p or at
  // or after p + span; both start ranges are searched in its bands.
  if (n_best > 1) {
    r.sub_hi = (int16_t)s1;
  } else {
    int64_t span = 0;
    for (int c = 0; c < r.n_cigar; ++c)
      if ((r.cigar[c] & 0xF) == kOpM || (r.cigar[c] & 0xF) == kOpD) span += r.cigar[c] >> 4;
    int64_t xlo;
    int xB;
    const Cluster& ch = ws.clusters[chosen];
    const int left = chain_best(ix, L, ch, pad, ws, &xlo, &xB, INT64_MIN, r.ref_pos - 1);
    const int right = chain_best(ix, L, ch, pad, ws, &xlo, &xB, r.ref_pos + span, INT64_MAX);
    int sh = s2;
    if (left > sh) sh = left;
    if (right > sh) sh = right;
    r.sub_hi = (int16_t)sh;
  }
  r.certified = 1;
  r.tier = kTierSL;
  r.support = (uint8_t)t;
  r.strand = ws.clusters[chosen].strand;
  r.d1 = -1;
  r.d2 = -1;
  r.d2x = -1;
  r.n_best = (uint16_t)n_best;
  r.n_clusters = (uint16_t)n;
  r.floor = (int16_t)floor;
  r.reason = kOk;
  // MAPQ only: the second-best score among evaluated chains, else the best of
  // the S2_TOP best-supported other chains (as bwa-mem uses its suboptimal hit).
  if (n_best > 1 || s2 > floor) {
    r.sub_score = (int16_t)s2;
  } else {
    for (int y = 0; y < n; ++y) ws.corder[y] = (uint16_t)y;
    std::stable_sort(ws.corder, ws.corder + n, [&](uint16_t a, uint16_t b) {
      return ws.clusters[a].end - ws.clusters[a].beg > ws.clusters[b].end - ws.clusters[b].beg;
    });
    int sub = kNoSub;
    for (int y = 0, done = 0; y < n && done < S2_TOP; ++y) {
      const int x = ws.corder[y];
      if (x == chosen) continue;
      ++done;
      const int s = ws.lscore[x] != kNoSub ? ws.lscore[x] : chain_best(ix, L, ws.clusters[x], pad, ws, &lo, &B);
      if (s > sub) sub = s;
    }
    r.sub_score = (int16_t)sub;
  }
  return true;
}
#endif

// The certified end-to-end result once clusters were evaluated with b1 <= R:
// the best locus (ties broken by name_hash), the second best within R, the
// SUB bound and, for MAPQ, the second locus's score.
CERTA_HD inline void finish_certified(const IndexView& ix, const Params& p, int L, int R, int nclust, int b1,
                                      int b2, uint64_t name_hash, Workspace& ws, Result& r) {
  int width;
  const int n_best = align_best(ix, L, nclust, b1, name_hash, ws, r, &width);
  r.certified = 1;
  r.tier = b1 == 0 ? kTierS0 : kTierS1;
  r.d1 = (int8_t)b1;
  r.d2 = (int8_t)(n_best > 1 ? b1 : (b2 <= R ? b2 : -1));
  r.n_best = (uint16_t)n_best;
  // SUB bound. Every other locus has end-to-end distance >= d* (the
  // exact second best, or R + 1 beyond the radius), so by Lemma E an
  // alignment there that scores above edit_bound(L, R + 1) completes
  // within R edits and lies in another cluster. Inside the reported band
  // (w diagonals), an alignment scoring s shares >= s + score - L read
  // bases with the reported one (each spans >= its score), and sharing
  // >= 2w + 5 forces their reference intervals to overlap (the paths
  // stay within w - 1 diagonals of each other); so a non-overlapping one
  // there scores at most L - score + 2w + 6 (and at most score, the band
  // optimum).
  {
    const int dstar = n_best > 1 ? b1 : (r.d2 >= 0 ? r.d2 : R + 1);
    const int other = edit_bound(L, dstar);
    int band = L * kMatch - r.score + 2 * width + 6;
    if (band > r.score) band = r.score;
    r.sub_hi = (int16_t)(other > band ? other : band);
  }
  if (n_best == 1 && r.d2 >= 0) {
    // The certified second-best locus: its affine score, for MAPQ.
    for (int x = 0; x < nclust; ++x)
      if (ws.clusters[x].dist == r.d2) { r.sub_score = cluster_score(ix, L, ws.clusters[x], ws); break; }
  } else if (n_best == 1 && p.mapq_limit > R) {
    r.d2x = (int8_t)other_candidates(ix, L, p.mapq_limit, r, ws, &r.sub_score);
  }
}

#ifndef __CUDA_ARCH__
// ---------------------------------------------------------------------------
// Pass V: adaptive parts (host only; THEORY.md §4b, Theorem A).
//
// The fixed partition looks each part up by one q-mer at a fixed offset, so
// one repetitive q-mer removes a part from S. Here the parts are chosen per
// read and strand, any disjoint q-mers of the read (s = 1: every reference
// position is indexed, so any q-mer can be looked up):
//  (A1) exact parts: the most disjoint q-mers whose hits fit the budget,
//       fewest hits in all among those (a DP over read offsets). Every
//       alignment with at most |S| - 1 edits leaves one of them untouched
//       (an edit touches at most one part), so R = |S| - 1, as in the fixed
//       certificate;
//  (A2) one-edit parts: the most disjoint q-mers with the fewest hits (no
//       budget), their one-edit neighbourhoods looked up, then the cheapest
//       neighbourhoods while they fit the budget. An alignment with at most
//       2|S| - 1 edits has a part with at most one edit; the reference q-mer
//       at that part's aligned start is then in its neighbourhood (Lemma N),
//       within one diagonal of the alignment there, so R = 2|S| - 1.
// Hits are clustered with pad R + 1 (an alignment with <= R edits stays
// within R diagonals of an untouched part's hit, R + 1 of a one-edit part's)
// into bands of up to BMAX diagonals, and verified exactly as in
// process_read. R is capped at min(p.k, VMAX).
constexpr uint64_t kVInf = UINT64_MAX / 4;

// One strand: the index range of the q-mer at every offset.
// lookup() of n keys, the same ranges: every key's bucket directory entry is
// prefetched, then every bucket, then each is searched, so the random reads
// of all n overlap instead of following one another.
inline void lookup_many(const IndexView& ix, const uint64_t* keys, int n, uint64_t* lo, uint64_t* hi) {
  const int shift = 2 * ix.q - ix.dir_bits;
  for (int i = 0; i < n; ++i) __builtin_prefetch(ix.dir + (keys[i] >> shift));
  for (int i = 0; i < n; ++i) {
    const uint64_t l = ix.dir[keys[i] >> shift];
    if (ix.keys16) __builtin_prefetch(ix.keys16 + l);
    else __builtin_prefetch(ix.keys + l);
  }
  for (int i = 0; i < n; ++i) lookup(ix, keys[i], lo + i, hi + i);
}

inline void vparts_ranges(const IndexView& ix, int L, Workspace& ws, int st) {
  const int n = L - ix.q + 1;
  uint64_t keys[LMAX], lo[LMAX], hi[LMAX];
  int at[LMAX], nk = 0;
  for (int a = 0; a < n; ++a) {
    if (kmer_key(ws.seq[st] + a, ix.q, &keys[nk])) at[nk++] = a;
    else { ws.vlo[st][a] = 1; ws.vhi[st][a] = 0; }  // an 'N': never a part
  }
  lookup_many(ix, keys, nk, lo, hi);
  for (int k = 0; k < nk; ++k) { ws.vlo[st][at[k]] = lo[k]; ws.vhi[st][at[k]] = hi[k]; }
}

// DP over disjoint q-mers of strand st with the given costs (kVInf: not a
// part): vf[st][k][i] is the least total cost of k of them within [0, i).
// Returns the most parts whose least total cost is <= budget.
inline int vparts_dp(int L, int q, const uint64_t* cost, uint64_t budget, Workspace& ws, int st) {
  const int K = L / q < PMAX ? L / q : PMAX;
  for (int i = 0; i <= L; ++i) { ws.vf[st][0][i] = 0; ws.vfrom[st][0][i] = -1; }
  for (int k = 1; k <= K; ++k)
    for (int i = 0; i <= L; ++i) {
      uint64_t best = i > 0 ? ws.vf[st][k][i - 1] : kVInf;
      int16_t from = -1;
      if (i >= q && ws.vf[st][k - 1][i - q] < kVInf && cost[i - q] < kVInf) {
        const uint64_t v = ws.vf[st][k - 1][i - q] + cost[i - q];
        if (v < best) { best = v; from = (int16_t)(i - q); }
      }
      ws.vf[st][k][i] = best;
      ws.vfrom[st][k][i] = from;
    }
  int m = 0;
  for (int k = K; k >= 1 && m == 0; --k)
    if (ws.vf[st][k][L] <= budget) m = k;
  return m;
}

// The m parts of least total cost found by vparts_dp, left to right, into voff.
inline void vparts_pick(int L, int m, Workspace& ws, int st) {
  for (int i = L, k = m; k > 0;) {
    const int a = ws.vfrom[st][k][i];
    if (a < 0) { --i; continue; }
    ws.voff[st][--k] = (int16_t)a;
    i = a;
  }
}

// Lemma N: the q-mers within one edit of the read q-mer s, as they start a
// reference stretch aligned to s: s itself; one substitution; one read base
// missing from the reference (the other q - 1, then any base); one reference
// base inserted at 0..q-1 (cut to q bases). Sorted, distinct.
inline int vparts_neighbours(const uint8_t* s, int q, uint64_t* keys) {
  uint8_t w[32];
  int n = 0;
  auto add = [&](const uint8_t* v) {
    uint64_t key;
    if (kmer_key(v, q, &key)) keys[n++] = key;
  };
  add(s);
  for (int i = 0; i < q; ++i)
    for (uint8_t x = 0; x < 4; ++x) {
      if (x == s[i]) continue;
      for (int j = 0; j < q; ++j) w[j] = s[j];
      w[i] = x;
      add(w);
    }
  for (int i = 0; i < q; ++i)
    for (uint8_t x = 0; x < 4; ++x) {
      int k = 0;
      for (int j = 0; j < q; ++j) if (j != i) w[k++] = s[j];
      w[k] = x;
      add(w);
    }
  for (int i = 0; i < q; ++i)
    for (uint8_t x = 0; x < 4; ++x) {
      int k = 0;
      for (int j = 0; j < i; ++j) w[k++] = s[j];
      w[k++] = x;
      for (int j = i; k < q; ++j) w[k++] = s[j];
      add(w);
    }
  std::sort(keys, keys + n);
  return (int)(std::unique(keys, keys + n) - keys);
}

// Clusters ws.cand (hits of the chosen parts, both strands) with pad R + 1
// and verifies them; certifies r when the best distance is <= R.
// lb_per: an alignment in a cluster hit by s of the m parts has at least
// lb_per * (m - s) edits (Theorem A), so a cluster whose bound already
// exceeds the current limit is skipped; (b1, b2) and every Cluster.dist come
// out as if it had been verified.
inline bool vparts_verify(const IndexView& ix, const Params& p, int L, int R, int m, int lb_per,
                          uint64_t name_hash, Workspace& ws, Result& r) {
  sort_i64(ws.cand[0], ws.ncand[0]);
  sort_i64(ws.cand[1], ws.ncand[1]);
  r.ncand = (uint32_t)(ws.ncand[0] + ws.ncand[1]);
  const int nclust = build_clusters(ws, R + 1, false, BMAX);
  if (nclust < 0) { r.reason = kClusterTooWide; return false; }
  r.n_clusters = (uint16_t)nclust;
  int b1 = R + 1, b2 = R + 1;
  evaluate_clusters(ix, L, nclust, R, p.reverse_order, nclust, 1, ws, &b1, &b2, m, lb_per);
  if (b1 > R) { r.reason = kNotFound; return false; }
  r.reason = kOk;
  finish_certified(ix, p, L, R, nclust, b1, b2, name_hash, ws, r);
  return true;
}

// Pass V for one read (ws.seq holds both strands): (A1), then (A2) if
// p.vparts >= 2. No witness: r.support = 0, r.probes = 0. r.used and
// r.radius describe the last step tried.
inline void certify_vparts(const IndexView& ix, const Params& p, int L, uint64_t name_hash, Workspace& ws,
                           Result& r) {
  const int q = ix.q;
  r.reason = kRadiusNegative;
  if (ix.s != 1 || L < q || q > 32) return;
  const int cap = p.k < VMAX ? p.k : VMAX;
  const uint64_t budget = (uint64_t)p.budget;
  uint64_t cost[2][LMAX];
  int mst[2];
  // (A1) the rarest disjoint q-mers that fit the budget.
  for (int st = 0; st < 2; ++st) {
    vparts_ranges(ix, L, ws, st);
    for (int a = 0; a + q <= L; ++a)
      cost[st][a] = ws.vlo[st][a] <= ws.vhi[st][a] ? ws.vhi[st][a] - ws.vlo[st][a] : kVInf;
    mst[st] = vparts_dp(L, q, cost[st], budget, ws, st);
  }
  int m = mst[0] < mst[1] ? mst[0] : mst[1];
  int R = m - 1 < cap ? m - 1 : cap;
  r.used = (uint8_t)m;
  r.radius = (int8_t)R;
  if (R >= 0) {
    for (int st = 0; st < 2; ++st) {
      vparts_pick(L, m, ws, st);  // the m of least total hits, within the budget
      int n = 0;
      for (int x = 0; x < m; ++x) {
        const int a = ws.voff[st][x];
        for (uint64_t e = ws.vlo[st][a]; e < ws.vhi[st][a]; ++e)
          ws.cand[st][n++] = ((int64_t)ix.pos[e] - a) * 16 + x;
      }
      ws.ncand[st] = n;
    }
    if (vparts_verify(ix, p, L, R, m, 1, name_hash, ws, r)) return;
  }
  if (p.vparts < 2) return;
  // (A2) one-edit parts: the most disjoint q-mers with the fewest hits, then
  // the cheapest one-edit neighbourhoods while they fit the budget (ties:
  // left to right).
  uint64_t keys[NBMAX], klo[NBMAX], khi[NBMAX];
  int order[2][PMAX];
  for (int st = 0; st < 2; ++st) {
    const int K = vparts_dp(L, q, cost[st], kVInf - 1, ws, st);
    vparts_pick(L, K, ws, st);
    for (int x = 0; x < K; ++x) {
      const int nk = vparts_neighbours(ws.seq[st] + ws.voff[st][x], q, keys);
      lookup_many(ix, keys, nk, klo, khi);
      uint64_t tot = 0;
      int nb = 0;
      for (int y = 0; y < nk; ++y)
        if (khi[y] > klo[y]) { ws.nblo[st][x][nb] = klo[y]; ws.nbhi[st][x][nb] = khi[y]; ++nb; tot += khi[y] - klo[y]; }
      ws.nbn[st][x] = (uint16_t)nb;
      ws.nbcost[st][x] = tot;
      order[st][x] = x;
      for (int y = x; y > 0 && ws.nbcost[st][order[st][y]] < ws.nbcost[st][order[st][y - 1]]; --y) {
        const int t = order[st][y];
        order[st][y] = order[st][y - 1];
        order[st][y - 1] = t;
      }
    }
    uint64_t used = 0;
    int mm = 0;
    while (mm < K && used + ws.nbcost[st][order[st][mm]] <= budget) used += ws.nbcost[st][order[st][mm++]];
    mst[st] = mm;
  }
  m = mst[0] < mst[1] ? mst[0] : mst[1];
  R = 2 * m - 1 < cap ? 2 * m - 1 : cap;
  r.used = (uint8_t)m;
  r.radius = (int8_t)R;
  if (R < 0) { r.reason = kRadiusNegative; return; }
  for (int st = 0; st < 2; ++st) {
    int n = 0;
    for (int y = 0; y < m; ++y) {  // part y: the y-th cheapest neighbourhood
      const int x = order[st][y], a = ws.voff[st][x];
      for (int z = 0; z < ws.nbn[st][x]; ++z)
        for (uint64_t e = ws.nblo[st][x][z]; e < ws.nbhi[st][x][z]; ++e)
          ws.cand[st][n++] = ((int64_t)ix.pos[e] - a) * 16 + y;
    }
    ws.ncand[st] = n;
  }
  vparts_verify(ix, p, L, R, m, 2, name_hash, ws, r);
}
#endif

// Map one read. `codes` holds L bases coded 0..3 (4 = N). `name_hash`
// breaks ties between equally good loci deterministically.
CERTA_HD inline void process_read(const IndexView& ix, const Params& p,
                                  const uint8_t* codes, int L,
                                  uint64_t name_hash, Workspace& ws,
                                  Result& r) {
  r.ref_pos = -1; r.n_cigar = 0; r.certified = 0; r.tier = kTierS0; r.strand = 0;
  r.reason = kOk; r.radius = -1; r.d1 = -1; r.d2 = -1; r.n_best = 0;
  r.n_clusters = 0; r.parts = 0; r.used = 0; r.nm = 0; r.score = 0; r.d2x = -1;
  r.sub_score = kNoSub; r.floor = 0; r.ncand = 0; r.sub_hi = kNoSub;
  r.probes[0] = r.probes[1] = 0; r.support = 0;

  const int P = L <= LMAX ? part_count(L, ix.q, ix.s) : 0;
  if (P < 1) {
    r.reason = kBadLength;
    return;
  }
  for (int i = 0; i < L; ++i) {
    ws.seq[0][i] = codes[i];
    uint8_t c = codes[L - 1 - i];
    ws.seq[1][i] = c < 4 ? (uint8_t)(3 - c) : (uint8_t)4;
  }
#ifndef __CUDA_ARCH__
  if (p.vparts > 0) {  // pass V: parts chosen per read (host only)
    r.parts = (uint8_t)P;
    certify_vparts(ix, p, L, name_hash, ws, r);
    return;
  }
#endif
  // Same number of rarest parts on both strands, within budget and k + 1.
  count_parts(ix, ws.seq[0], L, P, ws, 0);
  count_parts(ix, ws.seq[1], L, P, ws, 1);
  // With min_support t, |S| parts certify radius |S| - t (q-gram lemma).
  const int t = p.min_support < 1 ? 1 : p.min_support;
  // The end-to-end certificate needs at most k + t parts; a local-only pass
  // (tier SL) uses as many as fit the budget, which lowers its floor.
  const int lim = p.local_only ? PMAX : p.k + t;
  const int m0 = parts_that_fit(ws, 0, P, lim, p.budget);
  const int m1 = parts_that_fit(ws, 1, P, lim, p.budget);
  const int m = m0 < m1 ? m0 : m1;
  r.parts = (uint8_t)P;
  r.used = (uint8_t)m;
  for (int st = 0; st < 2; ++st)
    for (int x = 0; x < m; ++x) r.probes[st] |= (uint16_t)(1u << ws.order[st][x]);
  r.support = (uint8_t)t;
  const int R = m - t;
  r.radius = (int8_t)R;
  if (R < 0) {
    if (!certify_repeat(ix, L, P, name_hash, ws, r)) r.reason = kRadiusNegative;
    return;
  }
  enumerate_parts(ix, L, P, m, ws, 0);
  enumerate_parts(ix, L, P, m, ws, 1);
  r.ncand = (uint32_t)(ws.ncand[0] + ws.ncand[1]);

  sort_i64(ws.cand[0], ws.ncand[0]);
  sort_i64(ws.cand[1], ws.ncand[1]);

  // Certified path: every locus within R edits contains an enumerated part.
  const int nclust = p.local_only ? 0 : build_clusters(ws, 2 * R, false);  // indels shift <= R
  int b1 = R + 1, b2 = R + 1;
  if (nclust >= 0 && !p.local_only) {
    r.n_clusters = (uint16_t)nclust;
    evaluate_clusters(ix, L, nclust, R, p.reverse_order, nclust, p.min_support, ws, &b1, &b2);  // all: certificate
    if (b1 <= R) {
      finish_certified(ix, p, L, R, nclust, b1, b2, name_hash, ws, r);
      return;
    }
  }

#ifndef __CUDA_ARCH__
  if (p.local > 0 && certify_local(ix, p, L, P, m, name_hash, ws, r)) return;
#endif

  // Tier S2 (no certificate): best alignment among the enumerated candidates
  // with up to D edits. A better locus with R+1..D edits could exist where no
  // enumerated part matches exactly, so MAPQ is set conservatively by the host.
  const int D = p.s2_limit < S2_MAX ? p.s2_limit : S2_MAX;
  if (D > R) {
    const int nc2 = build_clusters(ws, 2 * D, true);  // skips over-wide clusters
    int c1 = D + 1, c2 = D + 1;
    // Heuristic: only the S2_TOP best-supported clusters, in a fixed order
    // (S2 makes no exactness claim, so the order-independence test does not
    // apply). The host caps MAPQ when clusters were left out.
    evaluate_clusters(ix, L, nc2, D, 0, S2_TOP, 1, ws, &c1, &c2);
    if (c1 <= D) {
      int width;
      const int n_best = align_best(ix, L, nc2, c1, name_hash, ws, r, &width);
      r.n_clusters = (uint16_t)nc2;
      r.certified = 2;
      r.tier = kTierS2;
      r.d1 = (int8_t)c1;
      r.d2 = (int8_t)(n_best > 1 ? c1 : (c2 <= D ? c2 : -1));
      r.n_best = (uint16_t)n_best;
      return;
    }
  }
  r.reason = nclust < 0 ? kClusterTooWide : kNotFound;
}

}  // namespace certa
