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

}  // namespace certa
