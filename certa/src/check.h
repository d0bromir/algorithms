// certa check: an independent checker of CERTA's interval output (Theorem V,
// EXACT_ALIGNMENT_PROPOSAL.md §4.7, Aim 4).
//
// The checker shares no alignment, bound or lookup code with the mapper. It
// shares the index file, which validate_index proves to be exactly the
// N-free q-mers of the reference at every s-th position (and the reference
// to be the FASTA), and the scoring scheme (BWA-MEM's defaults, -5 per
// clipped end). For each SAM record of `certa map --intervals --witness`:
//  - AS: the CIGAR, rescored against the reference, scores AS and stays
//    within one contig;
//  - XL: XL <= XU, and an alignment scoring XL is shown (the record's own,
//    or XV);
//  - XU: OPT <= XU, from the witness (XW): disjoint read intervals
//    ("probes") per strand, each at least q + s - 1 long, and a support t.
//      (i)   An alignment with fewer than t exactly occurring probes scores
//            at most U_t (Lemma L, recomputed here). Required: U_t <= XU.
//      (ii)  An alignment scoring above XU has gap total G <= pad = L - XU
//            - 7 (a gap costs >= 6 + its length). Its >= t exact probes
//            then have diagonals within G of each other, so they fall into
//            one chain of probe occurrences (all occurrences, enumerated
//            here, sorted by diagonal and split where consecutive ones are
//            more than pad apart), and the alignment lies in that chain's
//            reference window widened by pad.
//      (iii) It scores at most the chain's Lemma L bound with the chain's
//            probes free. Every chain whose bound exceeds XU is searched by
//            a local DP over the diagonals [first - pad, last + pad] of its
//            occurrences, which hold every cell of such an alignment; none
//            may score above XU.
//  - XH: SUB <= XH (the bound XQ was computed from), the same way, counting
//    only alignments that do not overlap the record's: on its strand the DP
//    runs left of its reference interval and right of it, separately;
//  - XQ: XQ <= MAPQ and XQ <= BWA-MEM's MAPQ formula at (AS, XH).
// The XQ of PR records rests on the pair certificate, which this checker
// does not re-derive; it is counted as unchecked.
#pragma once
#include <climits>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include "index.h"

namespace certa {
namespace check {

constexpr int kNone = INT_MIN;

struct Aln {
  int64_t pos = 0;  // 0-based, in the concatenated reference
  int strand = 0;
  std::vector<uint32_t> cigar;  // (length << 4) | op; op M = 0, I = 1, D = 2, S = 4
};

struct Probe {
  int off, end;
};

// The claims of one record.
struct Claim {
  std::vector<uint8_t> read;  // as sequenced, codes 0..3, 4 = N
  bool mapped = false;
  Aln aln;
  int as = kNone, mapq = 0;
  int xl = kNone, xu = kNone, xh = kNone, xq = 0;
  bool pair_xq = false;  // XQ rests on the pair certificate (PR)
  bool has_alt = false;
  Aln alt;               // XV: the alignment attaining XL
  int t = 0;             // XW: the support; 0 = no witness
  std::vector<Probe> probes[2];
};

enum Verdict : int {
  kAccepted = 0,
  kMalformed,   // unparsable record, CIGAR not covering the read, crossing contigs
  kWrongScore,  // AS differs from the rescored CIGAR
  kWrongXL,     // XL > XU, or no alignment shown that scores XL
  kXUUnproven,  // the witness does not prove OPT <= XU
  kXHUnproven,  // the witness does not prove SUB <= XH (or XH missing)
  kWrongXQ,     // XQ above MAPQ or above the formula at (AS, XH)
  kNumVerdicts
};
extern const char* const kVerdictNames[kNumVerdicts];

// Work done per record (for the cost of checking).
struct Work {
  uint64_t occurrences = 0, chains = 0, windows = 0, cells = 0;
  bool xu_checked = false, xh_checked = false, pair_xq = false;
};

class Checker {
 public:
  Checker(const Reference& ref, const Index& ix);
  Verdict check(const Claim& c, Work* w) const;
  // BWA-MEM local score of an alignment, or kNone if it is malformed, does
  // not cover the read or leaves its contig.
  int rescore(const std::vector<uint8_t>& read, const Aln& a) const;
  // A SAM record of certa map --intervals --witness; false if malformed.
  bool parse_sam(const std::string& line, Claim* c) const;

 private:
  int contig_of(int64_t pos) const;
  void range(uint64_t key, uint64_t* lo, uint64_t* hi) const;
  void occurrences(const uint8_t* rd, const Probe& p, std::vector<int64_t>& out) const;
  int best_band(const uint8_t* rd, int L, int64_t dlo, int64_t dhi, int64_t j0, int64_t j1, Work* w) const;
  bool proven(const Claim& c, const std::vector<uint8_t>* rd, int bound, bool sub, Work* w) const;

  const uint8_t* ref_;
  int64_t ref_len_;
  const Reference& names_;
  std::unordered_map<std::string, int> contig_;
  int q_, s_, dir_bits_;
  const uint64_t* keys_;
  const uint16_t* keys16_;
  const uint32_t* pos_;
  const uint64_t* dir_;
};

// Proves that the index lists exactly the N-free q-mers of its reference at
// every s-th position (each entry's key matches the reference at its
// position, positions are distinct, and their number equals the number of
// such q-mers), that its directory and key order make range lookups exact,
// and, given a FASTA path, that the stored reference equals the FASTA
// (contig names, lengths, bases; N between contigs). False with a reason.
bool validate_index(const Reference& ref, const Index& ix, const std::string& fasta, int threads,
                    std::string* why);

}  // namespace check
}  // namespace certa
