// SAM records and fallback FASTQ records of certa map.
#pragma once
#include <string>

#include "certa/bio.h"
#include "certa/interval.h"
#include "index.h"
#include "seqio.h"

namespace certa {

// The certa map options that change what a record says.
struct OutputOptions {
  int s2_limit = 0;           // --s2: tier S2's edit limit (its second best lies beyond it)
  int max_dist = 850;         // --max-dist: proper pairs start 0..max_dist apart
  bool unc_tags = false;      // --unc-tags: uncertified reads carry parts, used and reason in a FASTQ comment
  bool single_mates = false;  // --fallback-mates single: a lone uncertified mate goes to the fallback alone
  bool witness = false;       // --witness: interval records carry what certa check needs (Theorem V)
  bool mapq_cons = false;     // --iv-mapq cons: in interval mode MAPQ is the certified lower bound XQ
};

// The records of some reads, one string per output.
struct Records {
  std::string sam;
  std::string fq, fq2;  // for the fallback: reads (mate 1) and their mates
  std::string single;   // --fallback-mates single: lone mates
};

// Host-side check: a certified alignment must lie within one contig; one that
// does not becomes uncertified (kCrossContig).
void check_contig(Result& r, const Reference& ref);

// Appends the records of reads in the order they are written. Each call
// returns how the read (or pair) ended up: its final reason (kOk:
// certified, or placed by tier S2), after the contig check.
class RecordWriter {
 public:
  RecordWriter(const Reference& ref, const IndexView& view, const OutputOptions& opt)
      : ref_(ref), view_(view), opt_(opt) {}

  // Single-end: a SAM line if certified, else the FASTQ record for the
  // fallback. Returns the final reason.
  int single(const Result& r, const FastqRecord& rec, Records& out) const;

  // A pair (r, se and rec point at both mates; se: the single-end results,
  // for MAPQ): both mates as SAM lines if both are certified, else to the
  // fallback. Writes each mate's final reason to reason[0..1] and returns
  // whether both are certified as a proper pair.
  bool pair(const Result* r, const Result* se, const FastqRecord* rec, Records& out, int* reason) const;

  // Interval mode, single-end: every read gets a SAM record; its placement
  // goes to *placed. Returns the final reason (kNotFound for tier S2).
  int single_interval(const Result& r, const RepeatOut* iv, const FastqRecord& rec, Records& out,
                      Placement* placed) const;

  // Interval mode, a pair: both records (mates without a certificate placed
  // by tier R and the pair window, certa/interval.h), their placements to
  // placed[0..1] and final reasons to reason[0..1]. Returns whether the pair
  // is proper.
  bool pair_interval(const Result* r, const Result* se, const RepeatOut* const* iv, const uint8_t* const* codes,
                     const FastqRecord* rec, PairScratch& sc, Records& out, int* reason, Placement* placed) const;

  // --bio, single-end: a SAM line for a certified maximum-likelihood
  // placement (tier BM, with the Y tags) or a proven repeat (SR), else the
  // FASTQ record for the fallback. Returns the final reason.
  int bio_single(const Result& r, const BioOut& o, const FastqRecord& rec, Records& out) const;

  // --bio, a pair (r, o and rec point at both mates): both mates as SAM
  // lines if both are reportable (BM or SR), placed jointly when o[x].pair_ok;
  // one reportable mate as in pair(). Returns whether the pair is proper.
  bool bio_pair(const Result* r, const BioOut* o, const FastqRecord* rec, Records& out, int* reason) const;

 private:
  bool certified_pair(const Result* r, const Result* se, const FastqRecord* rec, Records& out,
                      Placement* placed) const;

  const Reference& ref_;
  const IndexView& view_;
  OutputOptions opt_;
};

}  // namespace certa
