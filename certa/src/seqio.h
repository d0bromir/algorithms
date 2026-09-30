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

// Reads up to `max_records` records; returns the number read. Line splitting
// is serial; building the records is spread over `threads` threads.
size_t read_fastq_batch(LineReader& in, std::vector<FastqRecord>& out,
                        size_t max_records, int threads = 1);

// Uncompressed FASTQ read in place from a memory mapping: one thread only
// finds record boundaries (memchr); records are built on several threads.
class MappedFastq {
 public:
  explicit MappedFastq(const std::string& path);
  ~MappedFastq();
  MappedFastq(const MappedFastq&) = delete;
  MappedFastq& operator=(const MappedFastq&) = delete;
  // True for a regular file that is not gzip-named (pipes and .gz use LineReader).
  static bool usable(const std::string& path);
  size_t next_batch(std::vector<FastqRecord>& out, size_t max_records, int threads);

 private:
  const char* base_ = nullptr;
  size_t size_ = 0, pos_ = 0;
};

}  // namespace certa
