#include "sam_output.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>

namespace certa {

namespace {

int mapq_of(const Result& r, int L, int s2_limit);

// BWA-MEM's single-end MAPQ, 30 (1 - S2/S1) ln(L), with alignment scores
// approximated from edit distances (each edit costs match + mismatch = 5).
int bwa_mapq(int L, int d1, int d2) {
  const double s1 = L * kMatch - d1 * (kMatch + kMismatch);
  const double s2 = L * kMatch - d2 * (kMatch + kMismatch);
  if (s1 <= 0 || s2 >= s1) return 0;
  const int q = static_cast<int>(30.0 * (1.0 - s2 / s1) * std::log(static_cast<double>(L)) + 0.499);
  return std::max(0, std::min(60, q));
}

// BWA-MEM's pairing: q_pe from the pair score margin; a mate's MAPQ is raised
// towards it (by at most 40 over its single-end MAPQ).
int mapq_pair(const Result& r, const Result& mate, const Result& se, int L, int s2_limit) {
  if (r.tier != kTierPR) return mapq_of(r, L, s2_limit);
  int q_se = se.certified == 1 ? mapq_of(se, L, s2_limit) : 0;
  const int q_pe = r.n_best == 1 ? mapq_margin(r.score + mate.score - r.sub_score) : 0;
  if (q_pe > q_se) q_se = std::min(q_pe, q_se + 40);
  return q_se;
}

// Interval mode, both mates certified: the single-end bounds of mate r (its
// single-end result se for a PR mate) and XQ. For a PR mate, XQ uses the
// certified pair margin R1 - S2 the way mapq_pair uses it.
Placement certified_placement(const Result& r, const Result& mate, const Result& se, int L, int mapq,
                              int s2_limit) {
  if (r.tier == kTierPR)
    return pr_placement(r, mate, se, L, mapq, se.certified == 1 ? mapq_of(se, L, s2_limit) : 0);
  return place_single(r, mapq, nullptr, L);
}

bool within_contig(const Result& r, const Reference& ref) {
  const int c = ref.contig_of(r.ref_pos);
  if (c < 0) return false;
  const int64_t end = r.ref_pos + ref_span_of(r) - 1;  // its last base: in contig c too
  return end >= static_cast<int64_t>(ref.offsets[c]) && end < static_cast<int64_t>(ref.offsets[c] + ref.lengths[c]);
}

// An alignment that runs off its contig (e.g. across the end of the
// circular chrM into the N padding) is not reported; neither are the claims
// that rest on it: XL (unless another alignment attains it), XH, XQ.
void keep_within_contigs(Placement& p, const Reference& ref) {
  if (p.xl_sep && !within_contig(p.xl_r, ref)) {
    p.xl_sep = false;
    p.xl = p.mapped ? p.r.score : kNoSub;
    if (p.pr_xh) {  // a PR mate's single-end term of XQ rested on it
      p.pr_xh = false;
      p.xq = std::min(p.mapq, std::min(p.q_pe, 40));
    }
  }
  if (!p.mapped || within_contig(p.r, ref)) return;
  p.mapped = false;
  p.mapq = p.xq = 0;
  if (!p.xl_sep) p.xl = kNoSub;
}

// TLEN of mate x (0 or 1) aligned at a, with its mate at b on the same contig.
int64_t tlen(const Result& a, const Result& b, int x) {
  const int64_t lo = std::min(a.ref_pos, b.ref_pos);
  const int64_t hi = std::max(a.ref_pos + ref_span_of(a), b.ref_pos + ref_span_of(b));
  const bool left = a.ref_pos < b.ref_pos || (a.ref_pos == b.ref_pos && x == 0);
  return left ? hi - lo : -(hi - lo);
}

// The complement of an ASCII base (A, C, G, T in either case; anything else
// gives N).
struct Complement {
  char c[256];
  Complement() {
    for (int i = 0; i < 256; ++i) c[i] = 'N';
    c['A'] = c['a'] = 'T';
    c['C'] = c['c'] = 'G';
    c['G'] = c['g'] = 'C';
    c['T'] = c['t'] = 'A';
  }
  char operator()(char b) const { return c[static_cast<uint8_t>(b)]; }
};
const Complement comp;

void fastq_record(const FastqRecord& rec, std::string& fq) {
  fq += '@'; fq += rec.name; fq += '\n';
  fq += rec.seq; fq += "\n+\n";
  fq += rec.qual; fq += '\n';
}

// Mate fields of a SAM line (single-end: none).
struct MateInfo {
  int flag = 0;
  std::string rnext = "*";
  int64_t pnext = 0, tlen = 0;
};

// Decimal digits of v at p (as printf's %d / %lld / %u); returns the end.
char* put_int(char* p, int64_t v) {
  uint64_t u = v < 0 ? 0 - static_cast<uint64_t>(v) : static_cast<uint64_t>(v);
  if (v < 0) *p++ = '-';
  char t[20];
  int n = 0;
  do {
    t[n++] = static_cast<char>('0' + u % 10);
    u /= 10;
  } while (u);
  while (n) *p++ = t[--n];
  return p;
}

char* put_str(char* p, const char* s, size_t n) {
  std::memcpy(p, s, n);
  return p + n;
}

// The CIGAR of r at p (at most 12 characters per operation); returns the end.
char* put_cigar(char* p, const Result& r) {
  for (int i = 0; i < r.n_cigar; ++i) {
    p = put_int(p, r.cigar[i] >> 4);
    *p++ = "MIDNS"[r.cigar[i] & 0xF];
  }
  return p;
}

// One record, formatted into a buffer and appended once (SAM formatting is
// a tenth of the CPU pass otherwise). `extra`: further tags, or null.
void sam_line(const Result& r, const FastqRecord& rec, const Reference& ref, int mapq, const MateInfo& mi,
              std::string& sam, const char* extra = nullptr) {
  const int c = ref.contig_of(r.ref_pos);
  const size_t L = rec.seq.size();
  const size_t nx = extra ? std::strlen(extra) : 0;
  const size_t need = rec.name.size() + ref.names[c].size() + mi.rnext.size() + L + rec.qual.size() + nx +
                      static_cast<size_t>(r.n_cigar) * 12 + 256;
  char local[2048];
  std::vector<char> big;
  char* const buf = need <= sizeof local ? local : (big.resize(need), big.data());
  char* p = buf;
  p = put_str(p, rec.name.data(), rec.name.size());
  *p++ = '\t';
  p = put_int(p, mi.flag | (r.strand ? 16 : 0));
  *p++ = '\t';
  p = put_str(p, ref.names[c].data(), ref.names[c].size());
  *p++ = '\t';
  p = put_int(p, r.ref_pos - static_cast<int64_t>(ref.offsets[c]) + 1);
  *p++ = '\t';
  p = put_int(p, mapq);
  *p++ = '\t';
  p = put_cigar(p, r);
  *p++ = '\t';
  p = put_str(p, mi.rnext.data(), mi.rnext.size());
  *p++ = '\t';
  p = put_int(p, mi.pnext);
  *p++ = '\t';
  p = put_int(p, mi.tlen);
  *p++ = '\t';
  if (r.strand) {
    for (size_t i = L; i-- > 0;) *p++ = comp(rec.seq[i]);
    *p++ = '\t';
    p = std::reverse_copy(rec.qual.begin(), rec.qual.end(), p);
  } else {
    p = put_str(p, rec.seq.data(), L);
    *p++ = '\t';
    p = put_str(p, rec.qual.data(), rec.qual.size());
  }
  static const char* const kTierNames[] = {"S0", "S1", "SR", "S2", "SL", "PR", "R", "IV", "MR"};
  // NM: edits in the reported (affine, possibly clipped) alignment.
  // XE: the certified minimum end-to-end edit distance over the reference
  // (-1 for SL). SL: AS (clip penalties included) is the certified maximum
  // local score over the reference, and XF the threshold it exceeds.
  auto tag = [&](const char* name, int64_t v) {
    p = put_str(p, name, 6);  // "\tXX:i:"
    p = put_int(p, v);
  };
  tag("\tNM:i:", r.nm);
  tag("\tAS:i:", r.score);
  tag("\tXE:i:", r.d1);
  p = put_str(p, "\tXT:Z:", 6);
  p = put_str(p, kTierNames[r.tier], std::strlen(kTierNames[r.tier]));
  tag("\tXR:i:", r.radius);
  tag("\tXD:i:", r.d2);
  tag("\tXB:i:", r.n_best);
  if (r.tier == kTierSL) tag("\tXF:i:", r.floor);
  // PR: XF is the threshold the certified pair score exceeds, XS a bound on
  // every other proper pair.
  if (r.tier == kTierPR) {
    tag("\tXF:i:", r.floor);
    tag("\tXS:i:", r.sub_score);
  }
  if (nx) p = put_str(p, extra, nx);
  *p++ = '\n';
  sam.append(buf, static_cast<size_t>(p - buf));
}

// An unmapped record (SAM: placed at the mate's position when it has one).
void sam_unmapped(const FastqRecord& rec, const Reference& ref, const MateInfo& mi, int64_t pos, int contig,
                  const char* extra, std::string& sam) {
  char buf[96];
  sam += rec.name;
  std::snprintf(buf, sizeof buf, "\t%d\t%s\t%lld\t0\t*\t%s\t%lld\t0\t", mi.flag | 4,
                contig >= 0 ? ref.names[contig].c_str() : "*",
                static_cast<long long>(contig >= 0 ? pos - ref.offsets[contig] + 1 : 0),
                contig >= 0 ? "=" : "*", static_cast<long long>(contig >= 0 ? pos - ref.offsets[contig] + 1 : 0));
  sam += buf;
  sam += rec.seq;
  sam += '\t';
  sam += rec.qual;
  sam += "\tXT:Z:UN";
  if (extra) sam += extra;
  sam += '\n';
}

// Interval mode: the tags every record carries. XL/XU: OPT in [XL, XU]
// (XL absent when no alignment was found); XQ: certified MAPQ lower bound.
// --witness adds what certa check needs to verify the claims (Theorem V):
//   XW:Z:t;o-e,o-e,...;o-e,...  the support t and the probes per strand
//                               (strand 1: offsets in the reverse complement);
//   XH:i:  the proven bound on SUB that XQ was computed from (single-end XQ;
//          for a PR mate, the single-end term of XQ, about loci not
//          overlapping its single-end alignment: XV if present, else its own);
//   XV:Z:contig,pos,strand,CIGAR  the alignment attaining XL, when it is not
//                               the record's own (PR mate: its single-end
//                               alignment, when it is not the record's own);
//   XY:Z:D;t;o-e,...;o-e,...    PR mates: the proper-pair distance D and the
//                               pair certificate's support and parts, which
//                               prove XS (the bound on every other proper pair).
std::string interval_tags(const Placement& p, const Reference& ref, const OutputOptions& opt) {
  const int xv = opt.witness && p.xl_sep ? ref.contig_of(p.xl_r.ref_pos) : -1;  // XV's contig
  // Tags of 6 + 20 characters; two witnesses of up to 2 PMAX probes each,
  // "o-e," <= 14 characters; XV's contig name and CIGAR.
  const size_t need = 16 * 26 + 2 * (2 * PMAX * 14 + 16) + (xv >= 0 ? ref.names[xv].size() : 0) + MAX_CIGAR * 12;
  char local[2048];
  std::vector<char> big;
  char* const buf = need <= sizeof local ? local : (big.resize(need), big.data());
  char* q = buf;
  auto tag = [&](const char* name, int64_t v) {
    q = put_str(q, name, 6);  // "\tXX:i:"
    q = put_int(q, v);
  };
  if (p.xl != kNoSub) tag("\tXL:i:", p.xl);
  tag("\tXU:i:", p.xu);
  tag("\tXQ:i:", p.xq);
  if (opt.witness) {
    if (p.xq > 0 && (p.r.tier != kTierPR || p.pr_xh)) tag("\tXH:i:", p.sub_hi);
    auto probes = [&](const Witness& w) {
      for (int st = 0; st < 2; ++st) {
        *q++ = ';';
        for (int j = 0; j < w.n[st]; ++j) {
          if (j) *q++ = ',';
          q = put_int(q, w.off[st][j]);
          *q++ = '-';
          q = put_int(q, w.end[st][j]);
        }
      }
    };
    if (p.w.t > 0) {
      tag("\tXW:Z:", p.w.t);
      probes(p.w);
    }
    if (p.mapped && p.pw.t > 0) {
      tag("\tXY:Z:", opt.max_dist);
      *q++ = ';';
      q = put_int(q, p.pw.t);
      probes(p.pw);
    }
    if (xv >= 0) {
      const Result& a = p.xl_r;
      q = put_str(q, "\tXV:Z:", 6);
      q = put_str(q, ref.names[xv].data(), ref.names[xv].size());
      *q++ = ',';
      q = put_int(q, a.ref_pos - static_cast<int64_t>(ref.offsets[xv]) + 1);
      *q++ = ',';
      *q++ = a.strand ? '-' : '+';
      *q++ = ',';
      q = put_cigar(q, a);
    }
  }
  return std::string(buf, q);
}

int mapq_of(const Result& r, int L, int s2_limit) {
  if (r.n_best > 1) return 0;
  if (r.tier == kTierS2) {
    // Second best: exact among the evaluated candidates, else at least D + 1.
    int q = bwa_mapq(L, r.d1, r.d2 >= 0 ? r.d2 : s2_limit + 1);
    // Clusters beyond the S2_TOP best-supported ones were not evaluated, so
    // an unseen tie is possible: keep MAPQ below the usual caller cut-off.
    if (r.n_clusters > S2_TOP) q = std::min(q, 10);
    return q;
  }
  // Certified: the second locus is exact within R, or beyond R the best other
  // candidate found (d2x), as bwa-mem uses its suboptimal hit. MAPQ follows
  // bwa-mem's default single-end formula (mem_approx_mapq_se with
  // mapQ_coef_len = 50) on affine scores, with sub = 19 (min_seed_len x a)
  // when no second locus was found.
  if (r.score > 0) {
    int sub = r.sub_score != kNoSub ? r.sub_score : 19;
    if (sub < 19) sub = 19;
    if (sub >= r.score) return 0;
    int q = mapq_se_raw(r.score, sub, L);
    // bwa subtracts 4.343 ln(sub_n + 1) for near-optimal alternative hits;
    // lower bound: one when the second locus is within one edit of the best.
    const int second = r.d2 >= 0 ? r.d2 : r.d2x;
    if (second >= 0 && second <= r.d1 + 1) q -= static_cast<int>(4.343 * std::log(2.0) + 0.499);
    return std::max(0, std::min(60, q));
  }
  const int second = r.d2 >= 0 ? r.d2 : r.d2x;
  return second >= 0 ? bwa_mapq(L, r.d1, second) : 60;
}

// A mate's result as paired output reports it: tier S2 is not reported for
// mates (kNotFound), and the contig check applies.
Result mate_as_reported(const Result& r0, const Reference& ref) {
  Result r = r0;
  if (r.certified == 2) {
    r.certified = 0;
    r.reason = kNotFound;
  }
  check_contig(r, ref);
  return r;
}

// Proper pair (as certified): one contig, opposite strands, and the reverse
// mate's first aligned base 0..max_dist after the forward mate's.
bool is_proper(const Result& a, const Result& b, const Reference& ref, int max_dist) {
  if (a.strand == b.strand || ref.contig_of(a.ref_pos) != ref.contig_of(b.ref_pos)) return false;
  const Result& f = a.strand ? b : a;
  const Result& r = a.strand ? a : b;
  return r.ref_pos >= f.ref_pos && r.ref_pos - f.ref_pos <= max_dist;
}

// The reason an interval-mode read is counted under: kOk if certified, else
// why not (a read placed by tier S2 is counted as not found).
int interval_reason(const Result& r) {
  return r.certified == 1 ? int(kOk) : (r.reason == kOk ? int(kNotFound) : int(r.reason));
}

}  // namespace

void check_contig(Result& r, const Reference& ref) {
  if (!r.certified) return;
  if (!within_contig(r, ref)) {
    r.certified = 0;
    r.reason = kCrossContig;
  }
}

int RecordWriter::single(const Result& r0, const FastqRecord& rec, Records& out) const {
  Result r = r0;
  check_contig(r, ref_);
  if (!r.certified) {
    if (opt_.unc_tags) {
      FastqRecord t = rec;
      t.name += " parts=" + std::to_string(r.parts) + " used=" + std::to_string(r.used) +
                " reason=" + std::to_string(r.reason);
      fastq_record(t, out.fq);
    } else {
      fastq_record(rec, out.fq);
    }
    return r.reason;
  }
  sam_line(r, rec, ref_, mapq_of(r, static_cast<int>(rec.seq.size()), opt_.s2_limit), MateInfo(), out.sam);
  return kOk;
}

// Both mates certified (r: as reported): their SAM lines, with mate fields;
// with `placed` (interval mode) also the interval tags. Returns whether the
// pair is proper.
bool RecordWriter::certified_pair(const Result* r, const Result* se, const FastqRecord* rec, Records& out,
                                  Placement* placed) const {
  const bool proper = is_proper(r[0], r[1], ref_, opt_.max_dist);
  const int c0 = ref_.contig_of(r[0].ref_pos), c1 = ref_.contig_of(r[1].ref_pos);
  for (int x = 0; x < 2; ++x) {
    const Result& a = r[x];
    const Result& b = r[1 - x];
    MateInfo mi;
    mi.flag = 1 | (proper ? 2 : 0) | (b.strand ? 32 : 0) | (x == 0 ? 64 : 128);
    const int cb = x == 0 ? c1 : c0;
    mi.rnext = c0 == c1 ? "=" : ref_.names[cb];
    mi.pnext = b.ref_pos - static_cast<int64_t>(ref_.offsets[cb]) + 1;
    if (c0 == c1) mi.tlen = tlen(a, b, x);
    const int L = static_cast<int>(rec[x].seq.size());
    const int mapq = mapq_pair(a, b, se[x], L, opt_.s2_limit);
    if (placed) {
      Placement pl = certified_placement(a, b, se[x], L, mapq, opt_.s2_limit);
      keep_within_contigs(pl, ref_);
      if (opt_.mapq_cons) pl.mapq = pl.xq;
      sam_line(a, rec[x], ref_, pl.mapq, mi, out.sam, interval_tags(pl, ref_, opt_).c_str());
      placed[x] = pl;
    } else {
      sam_line(a, rec[x], ref_, mapq, mi, out.sam);
    }
  }
  return proper;
}

bool RecordWriter::pair(const Result* r0, const Result* se, const FastqRecord* rec, Records& out, int* reason) const {
  const Result r[2] = {mate_as_reported(r0[0], ref_), mate_as_reported(r0[1], ref_)};
  if (r[0].certified && r[1].certified) {
    reason[0] = reason[1] = kOk;
    return certified_pair(r, se, rec, out, nullptr);
  }
  for (int x = 0; x < 2; ++x) reason[x] = r[x].reason;
  if (!r[0].certified && !r[1].certified) {
    fastq_record(rec[0], out.fq);
    fastq_record(rec[1], out.fq2);
    return false;
  }
  // One mate (k) keeps its single-end certificate.
  const int k = r[0].certified ? 0 : 1;
  reason[k] = kOk;
  const int mapq = mapq_of(r[k], static_cast<int>(rec[k].seq.size()), opt_.s2_limit);
  MateInfo mi;
  if (opt_.single_mates) {
    // The other mate is mapped single-end by the fallback, so both records
    // are written unpaired.
    sam_line(r[k], rec[k], ref_, mapq, mi, out.sam);
    fastq_record(rec[1 - k], out.single);
  } else {
    // Both mates still go to the fallback (so it can pair them), named
    // "name~k" where k is the certified mate, whose fallback record is
    // dropped downstream (bench/drop_certified_mates.awk). Mate fields are
    // unknown here.
    mi.flag = 1 | (k == 0 ? 64 : 128);
    sam_line(r[k], rec[k], ref_, mapq, mi, out.sam);
    FastqRecord m0 = rec[0], m1 = rec[1];
    m0.name += k == 0 ? "~1" : "~2";
    m1.name = m0.name;
    fastq_record(m0, out.fq);
    fastq_record(m1, out.fq2);
  }
  return false;
}

int RecordWriter::single_interval(const Result& r0, const RepeatOut* iv, const FastqRecord& rec, Records& out,
                                  Placement* placed) const {
  Result r = r0;
  check_contig(r, ref_);
  const int L = static_cast<int>(rec.seq.size());
  Placement p = place_single(r, r.certified == 1 ? mapq_of(r, L, opt_.s2_limit) : 0, iv, L);
  keep_within_contigs(p, ref_);
  if (opt_.mapq_cons) p.mapq = p.xq;
  const std::string t = interval_tags(p, ref_, opt_);
  if (p.mapped) sam_line(p.r, rec, ref_, p.mapq, MateInfo(), out.sam, t.c_str());
  else sam_unmapped(rec, ref_, MateInfo(), 0, -1, t.c_str(), out.sam);
  *placed = p;
  return interval_reason(r);
}

bool RecordWriter::pair_interval(const Result* r0, const Result* se_res, const RepeatOut* const* iv,
                                 const uint8_t* const* codes, const FastqRecord* rec, PairScratch& sc, Records& out,
                                 int* reason, Placement* placed) const {
  const Result r[2] = {mate_as_reported(r0[0], ref_), mate_as_reported(r0[1], ref_)};
  if (r[0].certified == 1 && r[1].certified == 1) {
    reason[0] = reason[1] = kOk;
    return certified_pair(r, se_res, rec, out, placed);
  }
  // A mate without a certificate: both placed from their single-end
  // placements and tier R outputs, with mate rescue in the pair window.
  Placement se[2];
  int L[2];
  for (int x = 0; x < 2; ++x) {
    reason[x] = interval_reason(r[x]);
    L[x] = static_cast<int>(rec[x].seq.size());
    se[x] = place_single(r[x], r[x].certified == 1 ? mapq_of(r[x], L[x], opt_.s2_limit) : 0, iv[x], L[x]);
    keep_within_contigs(se[x], ref_);
  }
  auto contig_of = [&](int64_t pos) { return ref_.contig_of(pos); };
  const bool proper = place_pair(view_, opt_.max_dist, 12, se, iv, codes, L, sc, contig_of, placed);
  if (opt_.mapq_cons)
    for (int x = 0; x < 2; ++x) placed[x].mapq = placed[x].xq;
  const int c0 = placed[0].mapped ? ref_.contig_of(placed[0].r.ref_pos) : -1;
  const int c1 = placed[1].mapped ? ref_.contig_of(placed[1].r.ref_pos) : -1;
  for (int x = 0; x < 2; ++x) {
    const Placement& a = placed[x];
    const Placement& b = placed[1 - x];
    const int ca = x == 0 ? c0 : c1, cb = x == 0 ? c1 : c0;
    MateInfo mi;
    mi.flag = 1 | (proper ? 2 : 0) | (x == 0 ? 64 : 128) | (b.mapped ? (b.r.strand ? 32 : 0) : 8);
    if (b.mapped) {
      mi.rnext = a.mapped && ca == cb ? "=" : ref_.names[cb];
      mi.pnext = b.r.ref_pos - static_cast<int64_t>(ref_.offsets[cb]) + 1;
      if (!a.mapped) mi.rnext = "=";
    } else if (a.mapped) {
      mi.rnext = "=";
      mi.pnext = a.r.ref_pos - static_cast<int64_t>(ref_.offsets[ca]) + 1;
    }
    if (a.mapped && b.mapped && ca == cb) mi.tlen = tlen(a.r, b.r, x);
    const std::string t = interval_tags(a, ref_, opt_);
    if (a.mapped) sam_line(a.r, rec[x], ref_, a.mapq, mi, out.sam, t.c_str());
    else sam_unmapped(rec[x], ref_, mi, b.mapped ? b.r.ref_pos : 0, b.mapped ? cb : -1, t.c_str(), out.sam);
  }
  return proper;
}

}  // namespace certa
