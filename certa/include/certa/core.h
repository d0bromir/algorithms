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
constexpr int PMAX = 8;                         // max parts per read
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
// Max band width (diagonals). The host allows wider bands, for tier SL.
#ifdef __CUDA_ARCH__
constexpr int BMAX = 48;
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
enum Tier : uint8_t { kTierS0 = 0, kTierS1 = 1, kTierSR = 2, kTierS2 = 3, kTierSL = 4, kTierPR = 5 };
constexpr int S2_MAX = 10;  // largest --s2 edit limit (band must fit BMAX)
constexpr int S2_TOP = 8;   // S2 evaluates only this many best-supported clusters
static_assert(MAX_CIGAR >= 2 * S2_MAX + 1, "CIGAR buffer too small for S2");
static_assert(4 * S2_MAX + 1 <= BMAX, "S2 band does not fit BMAX");

enum Reason : uint8_t {
  kOk = 0,
  kBadLength = 1,        // read shorter than one part needs, or longer than LMAX
  kRadiusNegative = 2,   // every part exceeds the budget and no SR proof
  kNotFound = 3,         // no locus within the certified radius
  kClusterTooWide = 4,   // candidate cluster wider than BMAX (tandem repeat)
  kCrossContig = 5,      // set by the host: alignment spans a contig boundary
  kNumReasons = 6
};

// BAM CIGAR op codes.
constexpr uint32_t kOpM = 0, kOpI = 1, kOpD = 2, kOpS = 4;

// Reported alignments use BWA-MEM's default affine scoring (the certificate
// itself is about unit-cost edit distance and picks the locus).
constexpr int kMatch = 1, kMismatch = 4, kGapOpen = 6, kGapExt = 1, kClip = 5;
constexpr int16_t kNoSub = -32768;  // Result::sub_score when no other locus was found

struct IndexView {
  const uint8_t* ref;     // concatenated reference, codes 0..3, 4 = N
  uint64_t ref_len;
  const uint64_t* keys;   // sorted q-mer keys of sampled positions, or null:
  const uint16_t* keys16; // low 16 bits only, when 2q - dir_bits <= 16 (the
                          // bucket directory already fixes the high bits)
  const uint32_t* pos;    // positions (multiples of s), parallel to keys
  const uint64_t* dir;    // bucket directory, size 2^dir_bits + 1
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
};

struct Cluster {
  int64_t lo;       // lowest diagonal (reference start) in band
  uint8_t width;    // number of diagonals in band
  uint8_t strand;
  int8_t dist;      // exact distance, or radius + 1 if it cannot beat the best two
  uint16_t beg, end;  // member diagonals: ws.cand[strand][beg, end)
  uint8_t pmask;      // distinct parts with a hit in the cluster (bit j = part j)
};

// Candidates store diagonal * 8 + part index (PMAX <= 8), so sorting orders
// by diagonal and the part survives the sort.
static_assert(PMAX <= 8, "part index must fit in 3 bits");
CERTA_HD inline int64_t diag_of(int64_t c) { return c >> 3; }  // floor(c / 8)
CERTA_HD inline int part_of(int64_t c) { return (int)(c & 7); }
CERTA_HD inline int popcount8(uint8_t v) {
  int n = 0;
  for (; v; v &= (uint8_t)(v - 1)) ++n;
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
  int16_t ubmemo[2][256];     // SL: bound per (strand, part mask)
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

// Half-open range [*lo, *hi) of index entries whose key equals `key`.
CERTA_HD inline void lookup(const IndexView& ix, uint64_t key, uint64_t* lo,
                            uint64_t* hi) {
  const int shift = 2 * ix.q - ix.dir_bits;
  uint64_t b = key >> shift;
  uint64_t l = ix.dir[b], r = ix.dir[b + 1];
  if (ix.keys16) {  // within a bucket, keys are ordered by their low bits
    const uint16_t k = (uint16_t)(key & ((1ull << shift) - 1));
    while (l < r) {
      uint64_t m = l + (r - l) / 2;
      if (ix.keys16[m] < k) l = m + 1; else r = m;
    }
    *lo = l;
    r = ix.dir[b + 1];
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
  r = ix.dir[b + 1];
  while (l < r) {
    uint64_t m = l + (r - l) / 2;
    if (ix.keys[m] <= key) l = m + 1; else r = m;
  }
  *hi = l;
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
    for (int t = 0; t < ix.s; ++t) {
      if (has_n) { ws.rlo[st][j][t] = ws.rhi[st][j][t] = 0; continue; }
      uint64_t key = 0;
      kmer_key(seq + off + t, ix.q, &key);
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
        ws.cand[st][n++] = ((int64_t)ix.pos[e] - (off + t)) * 8 + j;  // diagonal, part
  }
  ws.ncand[st] = n;
}

// Semi-global unit-cost edit distance of the whole read against the
// reference, restricted to diagonals [lo, lo + B). Returns `limit + 1` as
// soon as every cell exceeds `limit` (exact whenever the result <= limit).
CERTA_HD inline int band_distance(const IndexView& ix, const uint8_t* rd, int L,
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

CERTA_HD inline int hamming(const IndexView& ix, const uint8_t* rd, int L,
                            int64_t start, int limit) {
  int d = 0;
  for (int i = 0; i < L && d <= limit; ++i)
    d += sub_cost(rd[i], ref_at(ix, start + i));
  return d;
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
// is wider than BMAX and skip_wide is false (with skip_wide it is dropped).
CERTA_HD inline int build_clusters(Workspace& ws, int pad, bool skip_wide) {
  int nclust = 0;
  for (int st = 0; st < 2; ++st) {
    const int64_t* c = ws.cand[st];
    const int n = ws.ncand[st];
    int i = 0;
    while (i < n) {
      int64_t dmin = diag_of(c[i]), dmax = dmin;
      uint8_t pmask = (uint8_t)(1u << part_of(c[i]));
      int j = i + 1;
      while (j < n && diag_of(c[j]) - dmax <= 2 * pad) {
        dmax = diag_of(c[j]);
        pmask |= (uint8_t)(1u << part_of(c[j]));
        ++j;
      }
      const int64_t width = dmax - dmin + 2 * pad + 1;
      if (width > BMAX) {
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
                                       int* pb1, int* pb2) {
  order_by_support(ws, nclust);
  if (reverse)  // results must not depend on the order (tested)
    for (int x = 0, y = nclust - 1; x < y; ++x, --y) {
      uint16_t t = ws.corder[x];
      ws.corder[x] = ws.corder[y];
      ws.corder[y] = t;
    }
  int b1 = cap + 1, b2 = cap + 1;
  for (int x = 0; x < nclust; ++x) {
    Cluster& cl = ws.clusters[ws.corder[x]];
    cl.dist = (int8_t)(cap + 1);  // "cannot improve (b1, b2)", or not evaluated
    const int limit = b2 - 1;
    if (limit < 0 || x >= max_eval) continue;
    if (popcount8(cl.pmask) < min_support) continue;  // q-gram lemma: no locus within R here
    const uint8_t* rd = ws.seq[cl.strand];
    int u = limit + 1;  // Hamming upper bound, only tracked below limit + 1
    for (int e = cl.beg; e < cl.end && u > 0; ++e) {
      int h = hamming(ix, rd, L, diag_of(ws.cand[cl.strand][e]), u - 1);
      if (h < u) u = h;
    }
    int d;
    if (u == 0) {
      d = 0;
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

// Affine-optimal alignment of the read within diagonals [lo, lo + B) with
// BWA-MEM's default scores (+1 / -4 / gap -6 - 1 per base / clip -5 per end),
// reference ends free. Gotoh recurrences on the same band indexing as
// band_distance: cell (i, b) aligns read base i-1 to reference i-1+lo+b.
// Writes CIGAR (with soft clips), ref_pos, nm and score into `r`; returns
// false (leaving `r` unchanged) if the CIGAR would not fit. Optionally only
// alignments whose first aligned reference base lies in [smin, smax] count
// (false if there is none).
CERTA_HD inline bool affine_align(const IndexView& ix, const uint8_t* rd, int L,
                                  int64_t lo, int B, Workspace& ws, Result& r,
                                  bool score_only = false, int64_t smin = INT64_MIN,
                                  int64_t smax = INT64_MAX) {
  const int NEG = -(1 << 20);
  int* hp = ws.aff[0];
  int* hc = ws.aff[1];
  int* ep = ws.aff[2];
  int* ec = ws.aff[3];
  // tb bits: 0-1 H source (0 match/mismatch, 1 insertion E, 2 deletion F),
  // 2 the diagonal step starts the alignment (read prefix clipped or none),
  // 3 E extends E, 4 F extends F.
  for (int b = 0; b < B; ++b) { hp[b] = NEG; ep[b] = NEG; }
  int best = NEG, best_i = 0, best_b = 0;
  for (int i = 1; i <= L; ++i) {
    const uint8_t c = rd[i - 1];
    const int start = i == 1 ? 0 : -kClip;
    int f = NEG;
    for (int b = 0; b < B; ++b) {
      uint8_t bits = 0;
      int diag = hp[b];
      const int64_t pos = (int64_t)i - 1 + lo + b;
      if (start >= diag && pos >= smin && pos <= smax) { diag = start; bits |= 4; }
      const int m = diag + (sub_cost(c, ref_at(ix, pos)) ? -kMismatch : kMatch);
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
      if ((bits & 3) == 0) {  // an alignment may end after a match/mismatch
        const int end = h - (i < L ? kClip : 0);
        if (end > best) { best = end; best_i = i; best_b = b; }
      }
    }
    int* t = hp; hp = hc; hc = t;
    t = ep; ep = ec; ec = t;
  }
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
CERTA_HD inline int align_best(const IndexView& ix, int L, int nclust, int d1,
                               uint64_t name_hash, Workspace& ws, Result& r) {
  int n_best = 0;
  for (int x = 0; x < nclust; ++x) n_best += ws.clusters[x].dist == d1;
  int pick = (int)(name_hash % (uint64_t)n_best), chosen = -1;
  for (int x = 0; x < nclust; ++x)
    if (ws.clusters[x].dist == d1 && pick-- == 0) { chosen = x; break; }
  const Cluster& cl = ws.clusters[chosen];
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
  if (d1 > 0) {
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
inline int local_bound(const Workspace& ws, int st, int L, int P, int m, uint8_t free_mask,
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
      uint8_t pmask = (uint8_t)(1u << part_of(c[i]));
      int j = i + 1;
      while (j < nc && diag_of(c[j]) - dmax <= pad) {
        dmax = diag_of(c[j]);
        pmask |= (uint8_t)(1u << part_of(c[j]));
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

inline int band_score(const IndexView& ix, const uint8_t* rd, int L, int64_t lo, int B,
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
    int pm = NEG;
    hn[0] = hq[0];
    for (int b = 1; b < B; ++b) {
      const int g = hq[b - 1] + kGapExt * (b - 1);
      pm = pm > g ? pm : g;
      const int f = pm - kGapOpen - kGapExt * b;
      hn[b] = (int16_t)(hq[b] > f ? hq[b] : f);
    }
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

// Best local score in a chain: its member diagonals are covered by bands
// [d - pad, d + pad], merged into bands of at most BMAX diagonals. Returns
// the score (or kNoSub) and the band that attains it.
inline int chain_best(const IndexView& ix, int L, const Cluster& cl, int pad, Workspace& ws,
                      int64_t* best_lo, int* best_B) {
  const int64_t* c = ws.cand[cl.strand];
  int best = kNoSub;
  int a = cl.beg;
  while (a < cl.end) {
    const int64_t d0 = diag_of(c[a]);
    int e = a + 1;
    while (e < cl.end && diag_of(c[e]) - d0 <= BMAX - 2 * pad - 1) ++e;
    const int64_t lo = d0 - pad;
    const int B = (int)(diag_of(c[e - 1]) - d0) + 2 * pad + 1;
    const int sc = band_score(ix, ws.seq[cl.strand], L, lo, B);
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
  for (int st = 0; st < 2; ++st)
    for (int x = 0; x < 256; ++x) ws.ubmemo[st][x] = kNoSub;
  int ne = 0;
  for (int x = 0; x < n; ++x) {
    const Cluster& cl = ws.clusters[x];
    int16_t& ub = ws.ubmemo[cl.strand][cl.pmask];
    if (ub == kNoSub) ub = (int16_t)local_bound(ws, cl.strand, L, P, m, cl.pmask, 0);
    ws.lub[x] = ub;
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
  r.certified = 1;
  r.tier = kTierSL;
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

// Map one read. `codes` holds L bases coded 0..3 (4 = N). `name_hash`
// breaks ties between equally good loci deterministically.
CERTA_HD inline void process_read(const IndexView& ix, const Params& p,
                                  const uint8_t* codes, int L,
                                  uint64_t name_hash, Workspace& ws,
                                  Result& r) {
  r.ref_pos = -1; r.n_cigar = 0; r.certified = 0; r.tier = kTierS0; r.strand = 0;
  r.reason = kOk; r.radius = -1; r.d1 = -1; r.d2 = -1; r.n_best = 0;
  r.n_clusters = 0; r.parts = 0; r.used = 0; r.nm = 0; r.score = 0; r.d2x = -1;
  r.sub_score = kNoSub; r.floor = 0;

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
  // Same number of rarest parts on both strands, within budget and k + 1.
  count_parts(ix, ws.seq[0], L, P, ws, 0);
  count_parts(ix, ws.seq[1], L, P, ws, 1);
  // With min_support t, |S| parts certify radius |S| - t (q-gram lemma).
  const int t = p.min_support < 1 ? 1 : p.min_support;
  const int m0 = parts_that_fit(ws, 0, P, p.k + t, p.budget);
  const int m1 = parts_that_fit(ws, 1, P, p.k + t, p.budget);
  const int m = m0 < m1 ? m0 : m1;
  r.parts = (uint8_t)P;
  r.used = (uint8_t)m;
  const int R = m - t;
  r.radius = (int8_t)R;
  if (R < 0) {
    if (!certify_repeat(ix, L, P, name_hash, ws, r)) r.reason = kRadiusNegative;
    return;
  }
  enumerate_parts(ix, L, P, m, ws, 0);
  enumerate_parts(ix, L, P, m, ws, 1);

  sort_i64(ws.cand[0], ws.ncand[0]);
  sort_i64(ws.cand[1], ws.ncand[1]);

  // Certified path: every locus within R edits contains an enumerated part.
  const int nclust = p.local_only ? 0 : build_clusters(ws, 2 * R, false);  // indels shift <= R
  int b1 = R + 1, b2 = R + 1;
  if (nclust >= 0 && !p.local_only) {
    r.n_clusters = (uint16_t)nclust;
    evaluate_clusters(ix, L, nclust, R, p.reverse_order, nclust, p.min_support, ws, &b1, &b2);  // all: certificate
    if (b1 <= R) {
      const int n_best = align_best(ix, L, nclust, b1, name_hash, ws, r);
      r.certified = 1;
      r.tier = b1 == 0 ? kTierS0 : kTierS1;
      r.d1 = (int8_t)b1;
      r.d2 = (int8_t)(n_best > 1 ? b1 : (b2 <= R ? b2 : -1));
      r.n_best = (uint16_t)n_best;
      if (n_best == 1 && r.d2 >= 0) {
        // The certified second-best locus: its affine score, for MAPQ.
        for (int x = 0; x < nclust; ++x)
          if (ws.clusters[x].dist == r.d2) { r.sub_score = cluster_score(ix, L, ws.clusters[x], ws); break; }
      } else if (n_best == 1 && p.mapq_limit > R) {
        r.d2x = (int8_t)other_candidates(ix, L, p.mapq_limit, r, ws, &r.sub_score);
      }
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
      const int n_best = align_best(ix, L, nc2, c1, name_hash, ws, r);
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
