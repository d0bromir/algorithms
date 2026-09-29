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
#pragma once
#include <cstdint>

#ifdef __CUDACC__
#define CERTA_HD __host__ __device__
#else
#define CERTA_HD
#endif

namespace certa {

constexpr int KMAX = 5;                         // max certified radius (k)
constexpr int PMAX = 8;                         // max parts per read
constexpr int BUDGET_MAX = 256;                 // max hits enumerated per strand
constexpr int SMAX = 16;                        // max index sampling step
constexpr int LMAX = 320;                       // max read length handled
constexpr int BMAX = 48;                        // max band width (diagonals)
constexpr int MAX_CAND = BUDGET_MAX;            // candidates per strand
constexpr int MAX_CLUST = 2 * MAX_CAND;         // clusters over both strands
constexpr int MAX_CIGAR = 2 * KMAX + 2;
constexpr int REP_SAMPLE = 64;                  // hits examined for tier SR
constexpr int REP_KEEP = 8;                     // exact copies kept for tier SR
constexpr int BIG = 1 << 20;

enum Tier : uint8_t { kTierS0 = 0, kTierS1 = 1, kTierSR = 2 };

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
constexpr uint32_t kOpM = 0, kOpI = 1, kOpD = 2;

struct IndexView {
  const uint8_t* ref;     // concatenated reference, codes 0..3, 4 = N
  uint64_t ref_len;
  const uint64_t* keys;   // sorted q-mer keys of sampled positions
  const uint32_t* pos;    // positions (multiples of s), parallel to keys
  const uint64_t* dir;    // bucket directory, size 2^dir_bits + 1
  uint64_t n;
  int q, s, dir_bits;
};

struct Params {
  int k;       // maximum certified radius sought (at most k + 1 parts used)
  int budget;  // max hits enumerated per strand
  int reverse_order = 0;  // testing only: verify clusters least-supported first
};

struct Result {
  int64_t ref_pos;              // 0-based start in concatenated reference
  uint32_t cigar[MAX_CIGAR];    // BAM encoding: (len << 4) | op
  uint8_t n_cigar;
  uint8_t certified;            // 1 if the certificate holds
  uint8_t tier;                 // Tier: S0 exact, S1 1..R edits, SR repeat
  uint8_t strand;               // 0 forward, 1 reverse complement
  uint8_t reason;               // Reason when not certified
  int8_t radius;                // certified radius R (-1 for SR / uncertified)
  int8_t d1;                    // best edit distance
  int8_t d2;                    // second-best distance within R, -1 if none
  uint16_t n_best;              // loci tied at d1 (a lower bound for SR)
  uint16_t n_clusters;          // loci verified (diagnostic)
  uint8_t parts;                // parts the read was split into
  uint8_t used;                 // parts enumerated per strand (|S|)
};

struct Cluster {
  int64_t lo;       // lowest diagonal (reference start) in band
  uint8_t width;    // number of diagonals in band
  uint8_t strand;
  int8_t dist;      // exact distance, or radius + 1 if it cannot beat the best two
  uint16_t beg, end;  // member diagonals: ws.cand[strand][beg, end)
};

struct Workspace {
  uint8_t seq[2][LMAX];
  uint64_t rlo[2][PMAX][SMAX], rhi[2][PMAX][SMAX];  // index ranges per part/shift
  uint32_t count[2][PMAX];                          // hits per part (BIG if 'N')
  uint8_t order[2][PMAX];                           // parts, rarest first
  int64_t cand[2][MAX_CAND];
  int ncand[2];
  Cluster clusters[MAX_CLUST];
  uint16_t corder[MAX_CLUST];  // evaluation order of clusters
  int row_a[BMAX], row_b[BMAX];
  uint8_t tb[(LMAX + 1) * BMAX];
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
  uint64_t b = key >> (2 * ix.q - ix.dir_bits);
  uint64_t l = ix.dir[b], r = ix.dir[b + 1];
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

CERTA_HD inline void sort_i64(int64_t* a, int n) {
  for (int i = 1; i < n; ++i) {
    int64_t v = a[i];
    int j = i - 1;
    while (j >= 0 && a[j] > v) { a[j + 1] = a[j]; --j; }
    a[j + 1] = v;
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
        ws.cand[st][n++] = (int64_t)ix.pos[e] - (off + t);
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
  r.d2 = 0;
  r.n_best = (uint16_t)found;
  return true;
}

// Map one read. `codes` holds L bases coded 0..3 (4 = N). `name_hash`
// breaks ties between equally good loci deterministically.
CERTA_HD inline void process_read(const IndexView& ix, const Params& p,
                                  const uint8_t* codes, int L,
                                  uint64_t name_hash, Workspace& ws,
                                  Result& r) {
  r.ref_pos = -1; r.n_cigar = 0; r.certified = 0; r.tier = kTierS0; r.strand = 0;
  r.reason = kOk; r.radius = -1; r.d1 = -1; r.d2 = -1; r.n_best = 0;
  r.n_clusters = 0; r.parts = 0; r.used = 0;

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
  const int m0 = parts_that_fit(ws, 0, P, p.k + 1, p.budget);
  const int m1 = parts_that_fit(ws, 1, P, p.k + 1, p.budget);
  const int m = m0 < m1 ? m0 : m1;
  r.parts = (uint8_t)P;
  r.used = (uint8_t)m;
  const int R = m - 1;
  r.radius = (int8_t)R;
  if (R < 0) {
    if (!certify_repeat(ix, L, P, name_hash, ws, r)) r.reason = kRadiusNegative;
    return;
  }
  enumerate_parts(ix, L, P, m, ws, 0);
  enumerate_parts(ix, L, P, m, ws, 1);

  // Cluster sorted diagonals whose bands overlap.
  const int pad = 2 * R;  // indels shift the start diagonal by <= R
  int nclust = 0;
  for (int st = 0; st < 2; ++st) {
    int64_t* c = ws.cand[st];
    int n = ws.ncand[st];
    sort_i64(c, n);
    int i = 0;
    while (i < n) {
      int64_t dmin = c[i], dmax = c[i];
      int j = i + 1;
      while (j < n && c[j] - dmax <= 2 * pad) { dmax = c[j]; ++j; }
      int64_t width = dmax - dmin + 2 * pad + 1;
      if (width > BMAX) { r.reason = kClusterTooWide; return; }
      Cluster& cl = ws.clusters[nclust++];
      cl.lo = dmin - pad;
      cl.width = (uint8_t)width;
      cl.strand = (uint8_t)st;
      cl.beg = (uint16_t)i;
      cl.end = (uint16_t)j;
      i = j;
    }
  }
  r.n_clusters = (uint16_t)nclust;

  // Verify clusters, best-supported first, tracking the two smallest
  // distances b1 <= b2. Only a distance below b2 can change (b1, b2), so each
  // cluster is evaluated with limit b2 - 1: d1 = b1 is exact, and when b1 < b2
  // the second-best distance b2 (<= R) is exact too. Once two loci tie at b1,
  // the rest are skipped and n_best becomes a lower bound (MAPQ 0 either way).
  // A Hamming check at the member diagonals gives an upper bound first, so
  // exact matches never reach the DP.
  for (int x = 0; x < nclust; ++x) ws.corder[x] = (uint16_t)x;
  for (int x = 1; x < nclust; ++x) {  // stable sort by member count, descending
    uint16_t v = ws.corder[x];
    int sv = ws.clusters[v].end - ws.clusters[v].beg, y = x - 1;
    while (y >= 0 && ws.clusters[ws.corder[y]].end - ws.clusters[ws.corder[y]].beg < sv) {
      ws.corder[y + 1] = ws.corder[y];
      --y;
    }
    ws.corder[y + 1] = v;
  }
  if (p.reverse_order)  // results must not depend on the order (tested)
    for (int x = 0, y = nclust - 1; x < y; ++x, --y) {
      uint16_t t = ws.corder[x];
      ws.corder[x] = ws.corder[y];
      ws.corder[y] = t;
    }
  int b1 = R + 1, b2 = R + 1;
  for (int x = 0; x < nclust; ++x) {
    Cluster& cl = ws.clusters[ws.corder[x]];
    cl.dist = (int8_t)(R + 1);  // "cannot improve (b1, b2)"
    const int limit = b2 - 1;
    if (limit < 0) continue;
    const uint8_t* rd = ws.seq[cl.strand];
    int u = limit + 1;  // Hamming upper bound, only tracked below limit + 1
    for (int e = cl.beg; e < cl.end && u > 0; ++e) {
      int h = hamming(ix, rd, L, ws.cand[cl.strand][e], u - 1);
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

  const int d1 = b1;
  if (d1 > R) { r.reason = kNotFound; return; }
  int n_best = 0;
  for (int x = 0; x < nclust; ++x) n_best += ws.clusters[x].dist == d1;
  const int d2 = n_best > 1 ? d1 : (b2 <= R ? b2 : -1);
  // Deterministic, name-seeded choice among tied loci.
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
  r.certified = 1;
  r.tier = d1 == 0 ? kTierS0 : kTierS1;
  r.strand = cl.strand;
  r.d1 = (int8_t)d1;
  r.d2 = (int8_t)d2;
  r.n_best = (uint16_t)n_best;
}

}  // namespace certa
