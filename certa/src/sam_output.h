// SAM records and fallback FASTQ records of certa map.
#pragma once
#include <string>

#include "certa/interval.h"
#include "index.h"
#include "seqio.h"

namespace certa {

// The certa map options that change what a record says.
struct OutputOptions {
  int s2_limit = 0;        // --s2: tier S2's edit limit (its second best lies beyond it)
  int max_dist = 850;      // --max-dist: proper pairs start 0..max_dist apart
  bool unc_tags = false;   // --unc-tags: uncertified reads carry parts, used and reason in a FASTQ comment
  bool witness = false;    // --witness: interval records carry what certa check needs (Theorem V)
  bool mapq_cons = false;  // --iv-mapq cons: in interval mode MAPQ is the certified lower bound XQ
};

// MAPQ of a single-end result for a read of length L (BWA-MEM's formulas).
int mapq_of(const Result& r, int L, int s2_limit);

// Host-side check: a certified alignment must lie within one contig; one that
// does not becomes uncertified (kCrossContig).
void check_contig(Result& r, const Reference& ref);

// A mate's result as paired output reports it: tier S2 is not reported for
// mates (kNotFound), and the contig check applies.
Result mate_as_reported(const Result& r, const Reference& ref);

// Proper pair (as certified): one contig, opposite strands, and the reverse
// mate's first aligned base 0..max_dist after the forward mate's.
bool is_proper(const Result& a, const Result& b, const Reference& ref, int max_dist);

// Appends the records of reads to output strings (SAM, and the FASTQ files
// the fallback aligner reads), in the order they are written.
class RecordWriter {
 public:
  RecordWriter(const Reference& ref, const IndexView& view, const OutputOptions& opt)
      : ref_(ref), view_(view), opt_(opt) {}

  // Single-end: a SAM line if certified, else the FASTQ record for the
  // fallback. Returns the read's final reason (after the contig check).
  int single(const Result& r, const FastqRecord& rec, std::string& sam, std::string& fq) const;

  // A pair (r, se and rec point at both mates; se: the single-end results,
  // for MAPQ): both mates as SAM lines if both are certified, else to the
  // fallback. Writes each mate's final reason to reason[0..1]. With
  // `single_fq` (--fallback-mates single), the uncertified mate of a pair with
  // one certified mate goes there alone. With `placed` (interval mode, both
  // mates certified), the records carry the interval tags and their
  // placements go to placed[0..1].
  void pair(const Result* r, const Result* se, const FastqRecord* rec, std::string& sam, std::string& fq1,
            std::string& fq2, std::string* single_fq, int* reason, Placement* placed = nullptr) const;

  // Interval mode, single-end: every read gets a SAM record; its placement
  // goes to *out.
  void single_interval(const Result& r, const RepeatOut* iv, const FastqRecord& rec, std::string& sam,
                       Placement* out) const;

  // Interval mode, a pair with at least one uncertified mate
  // (certa/interval.h): both records, their placements to out[0..1]. Returns
  // whether the pair is proper.
  bool pair_interval(const Result* r, const RepeatOut* const* iv, const uint8_t* const* codes, const FastqRecord* rec,
                     PairScratch& sc, std::string& sam, Placement* out) const;

 private:
  const Reference& ref_;
  const IndexView& view_;
  OutputOptions opt_;
};

}  // namespace certa
