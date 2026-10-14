// Lemma L bounds for long reads (tools/lr_ceiling, EXACT_ALIGNMENT_PROPOSAL.md
// Aim 6). Parts are disjoint read intervals [off, end), each >= 3 long, not
// necessarily tiling the read. BWA-MEM's scores: +1, -4, gaps -(6 + k), -5
// per clipped read end; an inexact part inside an alignment costs >= 5.
#pragma once
#include <algorithm>
#include <climits>
#include <cstdint>
#include <vector>

#include "certa/core.h"

namespace certa {

struct ReadPart {
  int off, end;
  uint32_t occ;  // exact occurrences (0: absent)
};

// Lemma L, t = 1, with the parts flagged in `free` uncharged (null: none):
// every local alignment with no exact occurrence of a charged part scores
// <= max over spans [x, y) of (y - x) - 5 (charged parts inside) - 5
// (clipped ends). Span candidates x in {0, off + 1}, y in {L, end - 1}
// (Lemma S). O(P).
inline long lemma_l_span(int L, const std::vector<ReadPart>& p, const std::vector<char>* free = nullptr) {
  const int P = (int)p.size();
  const long c5 = kMismatch + kMatch;
  long G = 0, best = LONG_MIN, C = 0;  // C: charged parts among 0..j-1
  for (int j = 0; j <= P; ++j) {
    const long y = j == P ? L : p[j].end - 1;
    best = std::max(best, y - c5 * C - (y < L ? kClip : 0) + G);
    if (j < P) {
      C += !(free && (*free)[j]);
      G = std::max(G, -(long)(p[j].off + 1) + c5 * C - kClip);
    }
  }
  return best;
}

// Jump alignments (split alignments of a read across structural variants):
// pieces A_1 .. A_k, each a gapped local alignment of a read interval
// [a_j, b_j) on its own strand and locus, in read order (b_j <= a_{j+1});
// score sum sigma(A_j) - P (k - 1) - 5 [a_1 > 0] - 5 [b_k < L], read bases
// outside the pieces scoring 0. k = 1 is the local objective.
//
// Lemma L-J: if no piece holds an exact occurrence of a part of its strand,
// each piece scores <= (b_j - a_j) - 5 (parts of its strand inside it)
// (Lemma L), so the alignment scores <= U_J, the best such decomposition:
// returned here, exactly, by one pass over the read.
//
// part[s]: the parts of strand s in forward read coordinates (a strand-1
// part [o, e) of the reverse complement is [L - e, L - o) here), sorted.
// States: inside a piece on strand s, with c = 0 when the part holding the
// current base has been covered from its first base (it is charged 5 at its
// last base), or between pieces (a new piece pays P, the first one 5 if it
// starts after base 0).
inline long lemma_l_jump(int L, const std::vector<ReadPart>* part, long P) {
  const long NEG = LONG_MIN / 4, c5 = kMismatch + kMatch;
  std::vector<int> first[2], last[2], owner_off[2];
  for (int s = 0; s < 2; ++s) {
    first[s].assign(L, 0);
    last[s].assign(L, 0);
    owner_off[s].assign(L, -1);
    for (const ReadPart& q : part[s]) {
      first[s][q.off] = 1;
      last[s][q.end - 1] = 1;
      for (int i = q.off; i < q.end; ++i) owner_off[s][i] = q.off;
    }
  }
  long in[2][2] = {{NEG, NEG}, {NEG, NEG}}, mid = NEG;
  for (int i = 0; i < L; ++i) {
    long nx[2][2] = {{NEG, NEG}, {NEG, NEG}};
    const long start = std::max(i > 0 ? -(long)kClip : 0L, mid - P);
    for (int s = 0; s < 2; ++s) {
      // Continue a piece of strand s through base i.
      const long cont = std::max(in[s][0], in[s][1]);
      if (first[s][i]) {
        nx[s][0] = std::max(nx[s][0], cont);  // enters a new part at its first base
      } else if (owner_off[s][i] >= 0) {
        nx[s][0] = std::max(nx[s][0], in[s][0]);
        nx[s][1] = std::max(nx[s][1], in[s][1]);
      } else {
        nx[s][1] = std::max(nx[s][1], cont);  // no part here
      }
      // Start a piece at base i (inside a part past its first base: that
      // part is not charged).
      const int c = (owner_off[s][i] >= 0 && owner_off[s][i] == i) ? 0 : 1;
      nx[s][c] = std::max(nx[s][c], start);
      for (int k = 0; k < 2; ++k) {
        if (nx[s][k] <= NEG) continue;
        nx[s][k] += kMatch;
        if (k == 0 && last[s][i]) nx[s][k] -= c5;
      }
    }
    for (int s = 0; s < 2; ++s)
      for (int k = 0; k < 2; ++k) {
        in[s][k] = nx[s][k];
        mid = std::max(mid, in[s][k]);  // a piece may end after base i
      }
  }
  long best = mid - kClip;  // the last piece ends before L
  for (int s = 0; s < 2; ++s)
    for (int k = 0; k < 2; ++k) best = std::max(best, in[s][k]);  // or at L
  return best;
}

}  // namespace certa
