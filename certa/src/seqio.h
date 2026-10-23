// Buffered line reading for plain or gzip-compressed FASTA/FASTQ.
#pragma once
#include <cstdio>
#include <string>
#include <vector>

namespace certa {

class LineReader {
 public:
  explicit LineReader(const std::string& path);
  ~LineReader();
  LineReader(const LineReader&) = delete;
  LineReader& operator=(const LineReader&) = delete;
  // Reads one line without the trailing newline; false at end of file.
  bool getline(std::string& line);
  // Appends one line plus '\n' to `out`; returns the line length, or -1 at
  // end of file.
  long append_line(std::string& out);

 private:
  bool fill();
  void* gz_ = nullptr;     // gzFile when built with zlib
  std::FILE* fp_ = nullptr;
  std::vector<char> buf_;
  size_t beg_ = 0, end_ = 0;
  bool eof_ = false;
};

struct FastqRecord {
  std::string name, seq, qual;  // name without '@' and without comment
};

// The records of one batch as text: record i is base[starts[i], starts[i + 1]).
struct FastqSpans {
  const char* base = nullptr;
  std::vector<size_t> starts;
  std::string block;  // the records' text, when read from a stream
  size_t size() const { return starts.empty() ? 0 : starts.size() - 1; }
};

// Builds the records of `sp` into out[0], out[stride], out[2 * stride], ...
// on several threads. Each record's strings are assigned, so storage left by
// an earlier batch is reused.
void parse_spans(const FastqSpans& sp, FastqRecord* out, size_t stride, int threads);

// Splits up to `max_records` records off a stream (serial); returns how many.
size_t read_fastq_spans(LineReader& in, FastqSpans& sp, size_t max_records);

// Reads up to `max_records` records into `out` (resized to that many); returns
// the number read. Line splitting is serial; building the records is spread
// over `threads` threads.
size_t read_fastq_batch(LineReader& in, std::vector<FastqRecord>& out,
                        size_t max_records, int threads = 1);

// Uncompressed FASTQ read in place from a memory mapping. Record boundaries
// are found on several threads where that is exact (see next_spans), and
// records are built on several threads. The spans of one next_spans call are
// to be parsed before the next call, which drops the pages before them from
// memory (so memory does not grow with the input).
class MappedFastq {
 public:
  explicit MappedFastq(const std::string& path);
  ~MappedFastq();
  MappedFastq(const MappedFastq&) = delete;
  MappedFastq& operator=(const MappedFastq&) = delete;
  // True for a regular file that is not gzip-named (pipes and .gz use LineReader).
  static bool usable(const std::string& path);
  size_t next_spans(FastqSpans& sp, size_t max_records, int threads);
  size_t next_batch(std::vector<FastqRecord>& out, size_t max_records, int threads);

 private:
  const char* base_ = nullptr;
  size_t size_ = 0, pos_ = 0;
  size_t released_ = 0;          // the pages before this offset were dropped
  double bytes_per_record_ = 0;  // in the previous batch
};

}  // namespace certa
