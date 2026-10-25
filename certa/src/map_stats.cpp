#include "map_stats.h"

#include <cstdio>
#include <stdexcept>

namespace certa {

namespace {

const char* const kReasonNames[kNumReasons] = {
    "certified", "bad_length", "radius_negative", "not_found",
    "cluster_too_wide", "cross_contig"};

unsigned long long ull(uint64_t v) { return static_cast<unsigned long long>(v); }

}  // namespace

void Stats::count_read(const Result& r, size_t n_bases, int final_reason) {
  ++reads;
  cands += r.ncand;
  bases += n_bases;
  ++reason[final_reason];
  if (final_reason != kOk) {
    ++unc_by_used[r.used];
    return;
  }
  if (r.certified == 2) {  // placed without certificate
    ++s2;
    return;
  }
  ++certified;
  switch (r.tier) {
    case kTierS0: ++s0; break;
    case kTierS1: ++s1; break;
    case kTierSR: ++sr; break;
    case kTierSL: ++sl; break;
    default: ++pr; break;
  }
  ties += r.n_best > 1;
  if (r.radius >= 0 && r.radius <= VMAX + 1) ++by_radius[r.radius];
  if (r.d1 >= 0 && r.d1 <= VMAX) ++by_d1[r.d1];
}

void Stats::count_iv(const Placement& p) {
  if (!p.mapped) {
    ++iv_un;
    iv_unalignable += p.xu < kMinScore;
    return;
  }
  iv_r += p.r.tier == kTierR;
  iv_iv += p.r.tier == kTierIV;
  iv_mr += p.r.tier == kTierMR;
  const int w = p.xl == kNoSub ? BIG : p.xu - p.xl;
  ++iv_width[w == 0 ? 0 : w <= 5 ? 1 : w <= 10 ? 2 : w <= 25 ? 3 : w <= 50 ? 4 : 5];
  iv_point += w == 0;
  iv_mapq20 += p.mapq >= 20;
  iv_xq20 += p.xq >= 20;
}

void Stats::add(const Stats& o) {
  reads += o.reads; bases += o.bases; certified += o.certified;
  s0 += o.s0; s1 += o.s1; sr += o.sr; sl += o.sl; s2 += o.s2; ties += o.ties;
  for (int i = 0; i < kNumReasons; ++i) reason[i] += o.reason[i];
  for (int i = 0; i < VMAX + 2; ++i) by_radius[i] += o.by_radius[i];
  for (int i = 0; i < VMAX + 1; ++i) by_d1[i] += o.by_d1[i];
  for (int i = 0; i <= PMAX; ++i) unc_by_used[i] += o.unc_by_used[i];
  cands += o.cands;
  pr += o.pr; pairs += o.pairs; pairs_cert += o.pairs_cert; pairs_proper += o.pairs_proper;
  iv_r += o.iv_r; iv_iv += o.iv_iv; iv_mr += o.iv_mr; iv_un += o.iv_un; iv_unalignable += o.iv_unalignable;
  iv_point += o.iv_point; iv_proper += o.iv_proper; iv_mapq20 += o.iv_mapq20; iv_xq20 += o.iv_xq20;
  for (int i = 0; i < 6; ++i) iv_width[i] += o.iv_width[i];
}

void PassCounts::add(const PassCounts& o) {
  pass2 += o.pass2; pass3 += o.pass3; pass_v += o.pass_v; pass_v_certified += o.pass_v_certified;
  pass_l += o.pass_l; index2 += o.index2; index2_certified += o.index2_certified;
  pairs += o.pairs; tier_r += o.tier_r;
}

void print_summary(const MapConfig& c, const Stats& s, const PassCounts& n, const Timings& t) {
  auto pct = [&](uint64_t x) { return s.reads ? 100.0 * x / s.reads : 0.0; };
  std::fprintf(stderr,
               "[map] %llu reads: certified %.2f%% (S0 exact %.2f%%, S1 <=R edits %.2f%%, "
               "SR repeat %.2f%%, SL local %.2f%%), ties %.2f%%\n",
               ull(s.reads), pct(s.certified), pct(s.s0), pct(s.s1), pct(s.sr), pct(s.sl), pct(s.ties));
  for (int i = 1; i < kNumReasons; ++i)
    std::fprintf(stderr, "[map]   uncertified/%s: %.2f%%\n", kReasonNames[i], pct(s.reason[i]));
  std::fprintf(stderr, "[map] placed without certificate (S2, <= %d edits): %.2f%%; to fallback: %.2f%%\n",
               c.p2.s2_limit, pct(s.s2), pct(s.reads - s.certified - s.s2));
  std::fprintf(stderr, "[map] second pass (budget %d, S2 <= %d): %.2f%% of reads; third pass (CPU, budget %d): %.2f%%\n",
               c.p2.budget, c.p2.s2_limit, pct(n.pass2), c.p3.budget, pct(n.pass3));
  std::fprintf(stderr, "[map] local pass (CPU, support %d): %.2f%% of reads\n", c.pl.local, pct(n.pass_l));
  if (c.pass_v())
    std::fprintf(stderr, "[map] pass V (CPU, adaptive parts): %.2f%% of reads tried, %.2f%% certified there\n",
                 pct(n.pass_v), pct(n.pass_v_certified));
  if (c.index2())
    std::fprintf(stderr, "[map] second index (CPU, every pass): %.2f%% of reads tried, %.2f%% certified there\n",
                 pct(n.index2), pct(n.index2_certified));
  if (c.paired) {
    const double np = s.pairs ? 100.0 / s.pairs : 0.0;
    std::fprintf(stderr,
                 "[map] pairs %llu: both mates certified %.2f%% (proper %.2f%%); PR pair certificate "
                 "%.2f%% of reads; pair pass tried %.2f%% of pairs\n",
                 ull(s.pairs), np * s.pairs_cert, np * s.pairs_proper, pct(s.pr), np * n.pairs);
  }
  if (c.intervals) {
    const uint64_t* wd = s.iv_width;
    std::fprintf(stderr,
                 "[map] intervals: tier R on %.2f%% of reads (tau %d); records: tier R exact %.2f%%, interval "
                 "%.2f%%, pair window %.2f%%, unmapped %.2f%% (proven unalignable %.2f%%); exact (XL = XU) "
                 "%.2f%%; XU - XL 1-5 %.2f%%, 6-10 %.2f%%, 11-25 %.2f%%, 26-50 %.2f%%, >50 %.2f%%; "
                 "MAPQ >= 20 %.2f%%, XQ >= 20 %.2f%%%s\n",
                 pct(n.tier_r), c.rp.tau, pct(s.iv_r), pct(s.iv_iv), pct(s.iv_mr), pct(s.iv_un),
                 pct(s.iv_unalignable), pct(wd[0]), pct(wd[1]), pct(wd[2]), pct(wd[3]), pct(wd[4]),
                 pct(wd[5]), pct(s.iv_mapq20), pct(s.iv_xq20), c.out.mapq_cons ? " (MAPQ = XQ)" : "");
    if (c.paired)
      std::fprintf(stderr, "[map] intervals: proper pairs %.2f%%\n", s.pairs ? 100.0 * s.iv_proper / s.pairs : 0.0);
  }
  std::fprintf(stderr, "[map] candidates enumerated per read (last pass): %.1f\n",
               s.reads ? double(s.cands) / s.reads : 0.0);
  std::fprintf(stderr, "[map] uncertified by parts searched |S|:");
  for (int i = 0; i <= PMAX; ++i) std::fprintf(stderr, " %d:%.2f%%", i, pct(s.unc_by_used[i]));
  std::fprintf(stderr, "\n");
  const double map = t.map();
  std::fprintf(stderr,
               "[map] time: load %.2f s, upload %.2f s, read+encode %.2f s (overlapped; "
               "waited %.2f s), map %.2f s (%.0f reads/s), output %.2f s, total %.2f s "
               "(%.0f reads/s end to end)\n",
               t.load, t.upload, t.read_encode, t.read_wait, map, map > 0 ? s.reads / map : 0.0, t.output, t.total,
               t.total > 0 ? s.reads / t.total : 0.0);
}

void write_stats_json(const std::string& path, const MapConfig& c, const Index& ix, const std::string& device,
                      const Stats& s, const PassCounts& n, const Timings& t) {
  std::FILE* js = std::fopen(path.c_str(), "wb");
  if (!js) throw std::runtime_error("cannot open stats output");
  std::fprintf(js, "{\n  \"mode\": \"%s\",\n  \"device\": \"%s\",\n  \"threads\": %d,\n", c.gpu ? "gpu" : "cpu",
               device.c_str(), c.threads);
  std::fprintf(js, "  \"k\": %d, \"k1\": %d, \"budget\": %d, \"budget2\": %d, \"second_pass\": %llu, \"q\": %d, \"s\": %d,\n",
               c.p2.k, c.p1.k, c.p1.budget, c.p2.budget, ull(n.pass2), ix.q, ix.s);
  std::fprintf(js, "  \"reads\": %llu, \"bases\": %llu, \"certified\": %llu, \"s0\": %llu, \"s1\": %llu, \"sr\": %llu, \"sl\": %llu, \"s2\": %llu, \"s2_limit\": %d, \"ties\": %llu,\n",
               ull(s.reads), ull(s.bases), ull(s.certified), ull(s.s0), ull(s.s1), ull(s.sr), ull(s.sl), ull(s.s2),
               c.p2.s2_limit, ull(s.ties));
  std::fprintf(js, "  \"uncertified\": {");
  for (int i = 1; i < kNumReasons; ++i)
    std::fprintf(js, "%s\"%s\": %llu", i > 1 ? ", " : "", kReasonNames[i], ull(s.reason[i]));
  // Up to KMAX, or further when pass V certified larger distances or radii.
  int top = KMAX;
  for (int i = KMAX + 1; i <= VMAX; ++i)
    if (s.by_d1[i] || s.by_radius[i]) top = i;
  std::fprintf(js, "},\n  \"certified_by_d1\": [");
  for (int i = 0; i <= top; ++i) std::fprintf(js, "%s%llu", i ? ", " : "", ull(s.by_d1[i]));
  std::fprintf(js, "],\n  \"certified_by_radius\": [");
  for (int i = 0; i <= top; ++i) std::fprintf(js, "%s%llu", i ? ", " : "", ull(s.by_radius[i]));
  if (c.intervals) {
    std::fprintf(js, "],\n  \"intervals\": {\"tier_r_reads\": %llu, \"tier_r_exact\": %llu, \"interval\": %llu, "
                     "\"pair_window\": %llu, \"unmapped\": %llu, \"unalignable\": %llu, \"proper_pairs\": %llu, "
                     "\"mapq20\": %llu, \"xq20\": %llu},\n  \"iv_width\": [",
                 ull(n.tier_r), ull(s.iv_r), ull(s.iv_iv), ull(s.iv_mr), ull(s.iv_un), ull(s.iv_unalignable),
                 ull(s.iv_proper), ull(s.iv_mapq20), ull(s.iv_xq20));
    for (int i = 0; i < 6; ++i) std::fprintf(js, "%s%llu", i ? ", " : "", ull(s.iv_width[i]));
  }
  std::fprintf(js, "],\n  \"seconds\": {\"load\": %.3f, \"upload\": %.3f, \"read_encode\": %.3f, "
                   "\"read_wait\": %.3f, \"map\": %.3f, \"output\": %.3f, \"total\": %.3f}\n}\n",
               t.load, t.upload, t.read_encode, t.read_wait, t.map(), t.output, t.total);
  std::fclose(js);
}

}  // namespace certa
