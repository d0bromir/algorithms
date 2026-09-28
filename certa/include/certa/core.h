// CERTA certified short-read fast path: per-read core.
//
// This header is compiled both by the host C++ compiler (CPU mode, any
// architecture) and by nvcc (GPU mode). Every function is CERTA_HD, so the
// CPU and GPU back-ends run the *same* code and must produce bit-identical
// results.
//
// Certificate (see RESEARCH_ROADMAP.md, section 3.3):
//   A read of length L is split into P = k + 1 non-overlapping parts. Under
//   unit-cost edit distance, k edits touch at most k parts, so every locus
//   within k edits contains at least one part *exactly* (pigeonhole).
//   The index stores the q-mer at every s-th reference position. A part of
//   length >= q + s - 1 is queried at all s shifts, so an exact occurrence of
//   the part always contains a sampled position => its locus is found.
//   A part whose total hit count exceeds `cap` is not enumerated ("capped").
//   With c capped parts, a locus with d edits has >= k + 1 - d exact parts,
//   so it is still found whenever d <= k - c. The certified radius is
//   R = k - max(c_forward, c_reverse). All loci within R edits on both strands
//   are enumerated, and each is verified by exact banded edit distance.
//   An 'N' never matches (it costs 1, even against 'N'), so a part containing
//   'N' can never be an exact witness and is skipped without loss.
#pragma once
#include <cstdint>

#ifdef __CUDACC__
#define CERTA_HD __host__ __device__
#else
#define CERTA_HD
#endif

namespace certa {

constexpr int KMAX = 5;                         // max edits (k) supported
constexpr int CAP_MAX = 32;                     // max hits enumerated per part
constexpr int SMAX = 16;                        // max index sampling step
constexpr int LMAX = 320;                       // max read length handled
constexpr int BMAX = 48;                        // max band width (diagonals)
constexpr int MAX_CAND = (KMAX + 1) * CAP_MAX;  // candidates per strand
constexpr int MAX_CLUST = 2 * MAX_CAND;         // clusters over both strands
constexpr int MAX_CIGAR = 2 * KMAX + 2;
constexpr int BIG = 1 << 20;

enum Reason : uint8_t {
  kOk = 0,
  kBadLength = 1,        // read shorter than parts need, or longer than LMAX
  kRadiusNegative = 2,   // too many capped (repetitive) parts
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
  int k;    // requested edit budget (certificate radius before capping)
  int cap;  // max hits enumerated per part
};

struct Result {
  int64_t ref_pos;              // 0-based start in concatenated reference
  uint32_t cigar[MAX_CIGAR];    // BAM encoding: (len << 4) | op
  uint8_t n_cigar;
  uint8_t certified;            // 1 if the certificate holds
  uint8_t tier;                 // 0 = S0 exact, 1 = S1 certified <= R edits
  uint8_t strand;               // 0 forward, 1 reverse complement
  uint8_t reason;               // Reason when not certified
  int8_t radius;                // certified radius R (may be < 0)
  int8_t d1;                    // best edit distance
  int8_t d2;                    // second-best distance within R, -1 if none
  uint16_t n_best;              // number of loci tied at d1
  uint16_t n_clusters;          // loci verified (diagnostic)
  uint8_t capped_fwd, capped_rev;
};

struct Cluster {
  int64_t lo;       // lowest diagonal (reference start) in band
  uint8_t width;    // number of diagonals in band
  uint8_t strand;
  int8_t dist;      // banded edit distance, or > radius if none within radius
};

struct Workspace {
  uint8_t seq[2][LMAX];
  int64_t cand[2][MAX_CAND];
  int ncand[2];
  Cluster clusters[MAX_CLUST];
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

// Minimum read length for which every part is long enough (q + s - 1).
CERTA_HD inline int min_read_length(int k, int q, int s) {
  return (k + 1) * (q + s - 1);
}

CERTA_HD inline void sort_i64(int64_t* a, int n) {
  for (int i = 1; i < n; ++i) {
    int64_t v = a[i];
    int j = i - 1;
    while (j >= 0 && a[j] > v) { a[j + 1] = a[j]; --j; }
    a[j + 1] = v;
  }
}

// Collect candidate diagonals (reference start of the read) for one strand.
// Returns the number of capped parts.
CERTA_HD inline int collect_candidates(const IndexView& ix, const Params& p,
                                       const uint8_t* seq, int L, int64_t* cand,
                                       int* ncand) {
  const int P = p.k + 1;
  int capped = 0, n = 0;
  uint64_t rlo[SMAX], rhi[SMAX];
  for (int j = 0; j < P; ++j) {
    int off, len;
    part_geometry(L, P, j, &off, &len);
    bool has_n = false;
    for (int i = off; i < off + len; ++i) has_n |= (seq[i] > 3);
    if (has_n) continue;  // cannot be an exact witness; no loss
    uint64_t total = 0;
    for (int t = 0; t < ix.s; ++t) {
      uint64_t key = 0;
      kmer_key(seq + off + t, ix.q, &key);
      lookup(ix, key, &rlo[t], &rhi[t]);
      total += rhi[t] - rlo[t];
    }
    if (total > (uint64_t)p.cap) { ++capped; continue; }
    for (int t = 0; t < ix.s; ++t)
      for (uint64_t e = rlo[t]; e < rhi[t]; ++e)
        cand[n++] = (int64_t)ix.pos[e] - (off + t);
  }
  *ncand = n;
  return capped;
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

// Map one read. `codes` holds L bases coded 0..3 (4 = N). `name_hash`
// breaks ties between equally good loci deterministically.
CERTA_HD inline void process_read(const IndexView& ix, const Params& p,
                                  const uint8_t* codes, int L,
                                  uint64_t name_hash, Workspace& ws,
                                  Result& r) {
  r.ref_pos = -1; r.n_cigar = 0; r.certified = 0; r.tier = 0; r.strand = 0;
  r.reason = kOk; r.radius = -1; r.d1 = -1; r.d2 = -1; r.n_best = 0;
  r.n_clusters = 0; r.capped_fwd = 0; r.capped_rev = 0;

  if (L > LMAX || L < min_read_length(p.k, ix.q, ix.s)) {
    r.reason = kBadLength;
    return;
  }
  for (int i = 0; i < L; ++i) {
    ws.seq[0][i] = codes[i];
    uint8_t c = codes[L - 1 - i];
    ws.seq[1][i] = c < 4 ? (uint8_t)(3 - c) : (uint8_t)4;
  }
  int capped[2];
  for (int st = 0; st < 2; ++st)
    capped[st] = collect_candidates(ix, p, ws.seq[st], L, ws.cand[st],
                                    &ws.ncand[st]);
  r.capped_fwd = (uint8_t)capped[0];
  r.capped_rev = (uint8_t)capped[1];
  const int R = p.k - (capped[0] > capped[1] ? capped[0] : capped[1]);
  r.radius = (int8_t)R;
  if (R < 0) { r.reason = kRadiusNegative; return; }

  // Cluster sorted diagonals whose bands overlap, then verify each cluster.
  const int pad = 2 * p.k;  // indels shift the start diagonal by <= k
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
      int d = band_distance(ix, ws.seq[st], L, cl.lo, cl.width, R, ws.row_a,
                            ws.row_b);
      cl.dist = (int8_t)(d > R ? R + 1 : d);
      i = j;
    }
  }
  r.n_clusters = (uint16_t)nclust;

  int d1 = R + 1, n_best = 0;
  for (int x = 0; x < nclust; ++x) {
    if (ws.clusters[x].dist < d1) { d1 = ws.clusters[x].dist; n_best = 1; }
    else if (ws.clusters[x].dist == d1) ++n_best;
  }
  if (d1 > R) { r.reason = kNotFound; return; }
  int d2 = -1;
  if (n_best > 1) {
    d2 = d1;
  } else {
    for (int x = 0; x < nclust; ++x) {
      int d = ws.clusters[x].dist;
      if (d > d1 && d <= R && (d2 < 0 || d < d2)) d2 = d;
    }
  }
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
  r.tier = d1 == 0 ? 0 : 1;
  r.strand = cl.strand;
  r.d1 = (int8_t)d1;
  r.d2 = (int8_t)d2;
  r.n_best = (uint16_t)n_best;
}

}  // namespace certa
