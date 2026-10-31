// bio_score: the --bio model's penalty (certa/bio.h, bio_rescore) of every
// end-to-end alignment in a SAM file. The soundness check of
// bench/VALIDATION_PLAN.md: no aligner may report an alignment of a read that
// is cheaper than the optimum `certa map --bio` certified for it.
//   bio_score INDEX.cidx < aligner.sam > scores.tsv
// One line per primary, mapped record: QNAME, mate (0 single, 1, 2) and the
// penalty, or -1 when the record is not end to end (clips, a CIGAR op other
// than M/=/X/I/D, a leading or trailing gap), leaves its contig, or has no
// SEQ/QUAL, and -2 when an insertion is next to a deletion: the model's
// alignments exclude that (THEORY.md section 12), so such a record is outside
// what --bio certifies. The model is --bio's default (only the pair options
// differ, and they do not change a single read's penalty). SEQ and QUAL are
// taken as aligned, so every record is scored on strand 0.
#include <cstdio>
#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>

#include "certa/bio.h"
#include "index.h"

using namespace certa;

namespace {

// Splits a SAM line in place into its first `n` tab-separated fields.
int split(char* s, char** f, int n) {
  int k = 0;
  f[k++] = s;
  for (char* p = s; *p && k < n; ++p)
    if (*p == '\t') {
      *p = 0;
      f[k++] = p + 1;
    }
  return k;
}

// The record's CIGAR in BAM encoding, false if it is not end to end.
bool parse_cigar(const char* s, Result& r, int64_t& span, int& qlen) {
  r.n_cigar = 0;
  span = 0;
  qlen = 0;
  while (*s) {
    char* end;
    const long len = std::strtol(s, &end, 10);
    if (end == s || len <= 0 || r.n_cigar == MAX_CIGAR) return false;
    uint32_t op;
    switch (*end) {
      case 'M': case '=': case 'X': op = kOpM; span += len; qlen += len; break;
      case 'I': op = kOpI; qlen += len; break;
      case 'D': op = kOpD; span += len; break;
      default: return false;  // S, H, N, P
    }
    if (r.n_cigar > 0 && (r.cigar[r.n_cigar - 1] & 0xF) == op)
      r.cigar[r.n_cigar - 1] += static_cast<uint32_t>(len) << 4;  // M runs split by =/X
    else
      r.cigar[r.n_cigar++] = (static_cast<uint32_t>(len) << 4) | op;
    s = end + 1;
  }
  return r.n_cigar > 0;
}

// An insertion next to a deletion (outside the model's alignments).
bool adjacent_gaps(const Result& r) {
  for (int c = 1; c < r.n_cigar; ++c) {
    const uint32_t a = r.cigar[c - 1] & 0xF, b = r.cigar[c] & 0xF;
    if ((a == kOpI && b == kOpD) || (a == kOpD && b == kOpI)) return true;
  }
  return false;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 2) {
    std::fprintf(stderr, "usage: bio_score INDEX.cidx < aligner.sam > scores.tsv\n");
    return 2;
  }
  Reference ref;
  Index index;
  load_index(argv[1], ref, index);
  const IndexView ix = index.view(ref);
  std::unordered_map<std::string, int> contig;
  for (size_t c = 0; c < ref.names.size(); ++c) contig[ref.names[c]] = static_cast<int>(c);
  const BioModel model;
  BioRead br;
  std::vector<uint8_t> codes(LMAX);
  std::vector<char> line(1 << 20);
  long n = 0, scored = 0, adjacent = 0;
  char* f[11];
  while (std::fgets(line.data(), static_cast<int>(line.size()), stdin)) {
    if (line[0] == '@') continue;
    if (split(line.data(), f, 11) < 11) continue;
    const int flag = std::atoi(f[1]);
    if (flag & (4 | 256 | 2048)) continue;
    char* nl = std::strpbrk(f[10], "\t\n");
    if (nl) *nl = 0;
    const int mate = (flag & 64) ? 1 : (flag & 128) ? 2 : 0;
    ++n;
    int pen = -1;
    Result r{};
    int64_t span = 0;
    int qlen = 0;
    const auto c = contig.find(f[2]);
    const int L = static_cast<int>(std::strlen(f[9]));
    // A missing SEQ or QUAL is "*"; a quality string may start with '*' (Q9).
    if (c != contig.end() && std::strcmp(f[9], "*") != 0 && std::strcmp(f[10], "*") != 0 && L <= LMAX &&
        static_cast<int>(std::strlen(f[10])) == L && parse_cigar(f[5], r, span, qlen) && qlen == L) {
      const int64_t pos = std::atoll(f[3]) - 1;
      if (adjacent_gaps(r)) {
        pen = -2;
        ++adjacent;
      } else if (pos >= 0 && pos + span <= static_cast<int64_t>(ref.lengths[c->second])) {
        for (int i = 0; i < L; ++i) codes[i] = encode_base(f[9][i]);
        bio_read(model, codes.data(), f[10], L, br);
        r.ref_pos = static_cast<int64_t>(ref.offsets[c->second]) + pos;
        r.strand = 0;
        pen = bio_rescore(ix, model, br, r);
      }
    }
    scored += pen >= 0;
    std::printf("%s\t%d\t%d\n", f[0], mate, pen);
  }
  std::fprintf(stderr, "[bio_score] %ld primary mapped records, %ld end to end and scored, %ld with an insertion next to a deletion\n",
               n, scored, adjacent);
  return 0;
}
