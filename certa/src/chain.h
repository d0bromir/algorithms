// A UCSC chain file (e.g. hs1ToHg38.over.chain.gz) as a position lift from
// its source assembly (the chain's target, e.g. T2T-CHM13) to a reference
// (the chain's query, e.g. GRCh38), for certa map --t2t. The rule is
// bench/lift_sam.py's, which scored the T2T second opinion
// (bench/VALIDATION_RESULTS.md, "A MAPQ from T2T-CHM13").
#pragma once
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include "certa/pair_model.h"
#include "index.h"

namespace certa {

// --t2t: a mate's placement on T2T-CHM13, lifted to the reference
// (MapPipeline::t2t_chunk).
struct T2tLift {
  int state = 0;      // 0: unmapped there; 1: lifted; 3: mapped but not lifted
  int contig = -1;    // lifted: the reference contig
  int64_t left = 0;   // lifted: the unclipped leftmost base, as a reference (concatenated) position
  int strand = 0;     // lifted: the strand on the reference
  int tier = -1, mapq = 0, xe = -1, xd = -1, xo = -1;  // the T2T record's tier, MAPQ, XE, XD and XO
};

// What the pair model reads of a lift, relative to the record p placed on
// the reference (T2tView::cls): within 150 bp of its unclipped leftmost
// base, on its contig and strand, is "here".
T2tView t2t_view(const T2tLift& l, const Placement& p, const Reference& ref);

// A reference contig name as galaxy's Ensembl GRCh38 spells it: UCSC's chrN
// becomes N, chrM MT, chrUn_KI270302v1 and chr1_KI270706v1_random
// KI270302.1 and KI270706.1 (bench/truth_lift.py ensembl).
std::string ensembl_name(const std::string& ucsc);

class Chain {
 public:
  // Loads `path` (gzip or plain), resolving the chain's query contigs among
  // ref's names (as given, else by ensembl_name); blocks on other contigs
  // lift to nothing.
  static Chain load(const std::string& path, const Reference& ref);

  // A 0-based position on source contig `src`: the reference contig (index
  // in ref, or -1: the block's contig is not in ref), the 0-based forward
  // position there and whether the block is reversed. False when no block
  // holds it (the block with the largest start <= pos, ties in file order).
  bool lift(const std::string& src, int64_t pos, int* contig, int64_t* out, bool* reversed) const;

  // A record on source contig `src` with unclipped leftmost base `left`
  // (0-based) and unclipped span `span` (M, D, N, S, H, = and X): the
  // reference contig and the unclipped leftmost base there, and whether the
  // strand flips. Its first base on the source is lifted; if that falls in a
  // gap, the other end, offset by the span. False when neither end lifts,
  // or the end that lifts lies on a contig not in ref.
  bool lift_record(const std::string& src, int64_t left, int64_t span, int* contig, int64_t* out,
                   bool* flip) const;

  size_t blocks() const { return n_blocks_; }

 private:
  struct Block {
    int64_t start, end;  // on the source
    int contig;          // in ref, or -1
    bool reversed;
    int64_t qsize, qstart;
  };
  struct Contig {
    std::vector<int64_t> starts;
    std::vector<Block> b;
  };
  std::unordered_map<std::string, Contig> by_src_;
  size_t n_blocks_ = 0;
};

}  // namespace certa
