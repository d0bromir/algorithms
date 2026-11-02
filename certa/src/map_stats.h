// What certa map counts and how long it took: the summary on stderr and the
// --stats JSON.
#pragma once
#include <cstdint>
#include <string>

#include "certa/interval.h"
#include "index.h"
#include "map_config.h"

namespace certa {

// Per record, summed over the records written.
struct Stats {
  uint64_t reads = 0, bases = 0, certified = 0, s0 = 0, s1 = 0, sr = 0, sl = 0, s2 = 0, ties = 0;
  uint64_t reason[kNumReasons] = {};
  uint64_t by_radius[VMAX + 2] = {};  // certified reads per radius R (pass V reaches VMAX)
  uint64_t by_d1[VMAX + 1] = {};
  uint64_t unc_by_used[PMAX + 1] = {};  // uncertified reads per parts searched |S|
  uint64_t pr = 0;                      // reads certified by the pair certificate (PR)
  uint64_t pairs = 0, pairs_cert = 0, pairs_proper = 0;
  uint64_t cands = 0;                   // candidate diagonals enumerated (last pass per read)
  // Interval mode, per record: placed by tier R exact / as an interval / by
  // the pair window (MR) / unmapped (of which proven unalignable, XU < T);
  // XU - XL of mapped records: 0, 1-5, 6-10, 11-25, 26-50, > 50.
  uint64_t iv_r = 0, iv_iv = 0, iv_mr = 0, iv_un = 0, iv_unalignable = 0, iv_proper = 0;
  uint64_t iv_width[6] = {};
  uint64_t iv_mapq20 = 0, iv_xq20 = 0;
  // --bio, per record: certified ML placements (BM), of which elsewhere than
  // the edit certificate's best locus, placed jointly with the mate, and
  // moved there by the mate; MAPQ >= 20; reads whose parts were enumerated;
  // placed only by the mate (pair rescue).
  uint64_t bm = 0, bm_moved = 0, bm_pair = 0, bm_pair_moved = 0, bm_mapq20 = 0, bio_tried = 0, bm_rescued = 0;

  // A read of n_bases bases with its final result r, written with
  // final_reason (kOk: certified, or placed by tier S2 when r.certified == 2).
  void count_read(const Result& r, size_t n_bases, int final_reason);
  // An interval-mode record.
  void count_iv(const Placement& p);
  // A --bio record (after count_read).
  void count_bio(const Result& r, const BioOut& o, int final_reason);
  void add(const Stats& o);
};

// Reads each pass was given (and, where a pass keeps only certificates, how
// many it certified), summed over batches.
struct PassCounts {
  uint64_t pass2 = 0;       // re-run with --budget2
  uint64_t pass3 = 0;       // re-run on the CPU with --budget3
  uint64_t pass_v = 0, pass_v_certified = 0;
  uint64_t pass_l = 0;      // re-run on the CPU for tier SL
  uint64_t index2 = 0, index2_certified = 0;
  uint64_t pairs = 0;       // pairs given to the pair certificate
  uint64_t tier_r = 0;      // reads given to tier R (interval mode)
  void add(const PassCounts& o);
};

// Wall time in seconds. Mapping is summed over the first stage (passes 1
// and 2) and the rest (with a GPU they overlap; fused, the rest includes
// formatting).
struct Timings {
  double load = 0, upload = 0, read_encode = 0, read_wait = 0, map_first = 0, map_rest = 0, output = 0, total = 0;
  double map() const { return map_first + map_rest; }
};

// The summary certa map prints on stderr.
void print_summary(const MapConfig& c, const Stats& s, const PassCounts& n, const Timings& t);

// The --stats JSON, to c.stats_path; `device` is the GPU's name, or empty.
void write_stats_json(const MapConfig& c, const Index& ix, const std::string& device, const Stats& s,
                      const PassCounts& n, const Timings& t);

}  // namespace certa
