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
// PR mates (XT:Z:PR) follow BWA-MEM's pairing rule instead:
//   XQ <= max(q_se, min(q_pe, q_se + 40)), which rises with both terms, so
//   proven lower bounds of them suffice:
//  - q_se: the single-end formula at (score of the single-end alignment, XH),
//    XH proven as above about loci not overlapping that alignment (XV if
//    present, else the record's own); 0 without XH;
//  - q_pe: BWA-MEM's raw MAPQ of the pair margin AS + AS(mate) - XS, when the
//    two records form a proper pair (one contig, opposite strands, the
//    reverse mate's first aligned base 0..D after the forward mate's) and
//    every proper pair that differs from it in some mate (that mate's
//    alignment does not overlap the record's) scores <= XS. The pair
//    witness (XY: D, a support t, each mate's parts per strand) proves it:
//      (P1) a pair in which neither mate has t exactly occurring parts
//           scores <= U_t(mate 1) + U_t(mate 2) (Lemma L per mate);
//      (P2) a pair scoring above XS has gap total <= pad = L1 + L2 - XS - 7
//           in each mate. If mate x has t exact parts, they lie in one chain
//           of its part occurrences (split where diagonals are > pad apart),
//           x starts in that chain's start range, and the partner y, on the
//           other strand, starts within D of x's first aligned base;
//      (P3) per chain, the exact best proper pair through it must not exceed
//           XS: the best x per start (exhaustive DP over the chain's segment)
//           with the best y per start within D of it (exhaustive DP over the
//           window), the shown pair excluded (x wholly left or right of the
//           record's alignment, or y of its mate's). A chain is skipped when
//           its Lemma L bound plus a bound on y there (U_t(y), or the bound
//           of a chain of y whose start range meets the window) is <= XS.
//   Every pair (P3) counts is a real proper pair (or one across a contig
//   boundary), so an honest XS is proven, up to the work caps.
//   check() alone cannot see the mate: a PR record whose XQ exceeds q_se is
//   kXSUnproven there, and check_pair() decides it from both records.
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
  std::string qname;
  std::vector<uint8_t> read;  // as sequenced, codes 0..3, 4 = N
  bool mapped = false;
  Aln aln;
  int as = kNone, mapq = 0;
  int xl = kNone, xu = kNone, xh = kNone, xq = 0;
  bool pair_xq = false;  // a PR mate: XQ follows the pairing rule
  int xs = kNone;        // XS: the bound on every other proper pair
  bool has_alt = false;
  Aln alt;               // XV: the alignment attaining XL (PR: the single-end alignment)
  int t = 0;             // XW: the support; 0 = no witness
  std::vector<Probe> probes[2];
  int pd = -1, pt = 0;   // XY: proper-pair distance D and support (pt = 0: no pair witness)
  std::vector<Probe> pprobes[2];
};

enum Verdict : int {
  kAccepted = 0,
  kMalformed,   // unparsable record, CIGAR not covering the read, crossing contigs
  kWrongScore,  // AS differs from the rescored CIGAR
  kWrongXL,     // XL > XU, or no alignment shown that scores XL
  kXUUnproven,  // the witness does not prove OPT <= XU
  kXHUnproven,  // the witness does not prove SUB <= XH (or XH missing)
  kWrongXQ,     // XQ above MAPQ or above the formula at (AS, XH)
  kXSUnproven,  // a PR mate's XQ needs XS, which its pair does not prove
  kNumVerdicts
};
extern const char* const kVerdictNames[kNumVerdicts];

// Work done per record (for the cost of checking).
struct Work {
  uint64_t occurrences = 0, chains = 0, windows = 0, cells = 0;
  bool xu_checked = false, xh_checked = false;
  bool pair_xq = false;     // a PR mate with XQ > 0
  bool needs_pair = false;  // its XQ exceeds q_se: it rests on XS
  bool xs_checked = false;  // XS proven from the pair witness
  int q_se = 0;
};

class Checker {
 public:
  Checker(const Reference& ref, const Index& ix);
  // One record. A PR mate whose XQ exceeds its single-end term is
  // kXSUnproven here: that claim needs its mate (check_pair).
  Verdict check(const Claim& c, Work* w) const;
  // Both records of a pair (the same QNAME): the verdicts of check(), with
  // the pair certificate decided for PR mates.
  void check_pair(const Claim& a, const Claim& b, Work* wa, Work* wb, Verdict* va, Verdict* vb) const;
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
  // excl: null for OPT <= bound; else SUB <= bound about loci not overlapping excl.
  bool proven(const Claim& c, const std::vector<uint8_t>* rd, int bound, const Aln* excl, Work* w) const;
  // Per reference position j in [j0, j1): the best score of the read over
  // alignments starting at j with all cells in [j0, j1).
  void starts(const uint8_t* rd, int L, int64_t j0, int64_t j1, std::vector<int>& best, Work* w) const;
  bool pair_proven(const Claim& a, const Claim& b, int bound, Work* w) const;

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
