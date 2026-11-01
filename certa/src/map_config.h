// The options of certa map, parsed and checked before anything is loaded.
#pragma once
#include <cstddef>
#include <string>

#include "certa/bio.h"
#include "certa/pair.h"
#include "certa/repeat.h"
#include "cli.h"
#include "sam_output.h"

namespace certa {

// A read goes through the passes in order, each re-running only the reads the
// passes before it left uncertified:
//   1  every read, small budget and radius (--budget, --k1; GPU or CPU);
//   2  larger budget and radius (--budget2, -k), and tier S2 (--s2);
//   3  CPU, a larger budget for reads whose parts did not all fit (--budget3);
//   V  CPU, parts chosen per read on an s = 1 index (--adaptive-parts);
//   L  CPU, certified local alignment (--local);
// then, on the CPU, every pass again on a second index (--index2), the pair
// certificate (paired input) and tier R (--intervals).
// With --bio a read instead goes through one CPU pass (pb), which certifies
// its maximum-likelihood placement (certa/bio.h), and pairs are placed
// jointly from both mates' loci.
struct MapConfig {
  // Inputs and outputs. Empty: not given (SAM then goes to stdout).
  std::string index_path, reads_path, mates_path, index2_path;  // mates_path empty: single-end
  std::string sam_path, fq_path, fq2_path, single_path, stats_path;
  bool intervals = false;  // --intervals: every read to SAM, no fallback
  bool full_cert = false;  // --full-cert
  bool bio = false;        // --bio: certified maximum-likelihood placement
  bool bio_table = true;   // --bio-run-table: the reference's tandem runs tabulated once (half a byte a base)
  bool gpu = false;
  int device = 0;
  bool serial_stages = false;
  int threads = 1, io_threads = 1;
  size_t batch_size = 0;

  Params p1, p2, p3, pv, pl;  // passes 1, 2, 3, V and L
  PairParams pp;              // the pair certificate
  RepeatParams rp;            // tier R (interval mode)
  Params pb;                  // --bio: the one pass
  BioParams bp;               // --bio: the model
  OutputOptions out;

  bool paired() const { return !mates_path.empty(); }
  bool pass2() const { return p2.budget > p1.budget || p2.s2_limit > 0 || p2.k > p1.k; }
  bool pass3() const { return p3.budget > p2.budget; }
  bool pass_v() const { return pv.vparts > 0; }
  bool pass_l() const { return pl.local > 0; }
  bool index2() const { return !index2_path.empty(); }

  // Throws std::runtime_error, with a usage message, on a bad option.
  static MapConfig from_args(const Args& a);
  // The loaded index's sampling step s: pass V needs s = 1, and is turned off
  // (with a note on stderr) otherwise.
  void fit_index(int s);
  // The @PG command line: the passes and their parameters.
  std::string command_line() const;
};

}  // namespace certa
