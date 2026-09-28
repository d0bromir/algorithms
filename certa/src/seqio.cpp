#include "seqio.h"

#include <algorithm>
#include <cstring>
#include <exception>
#include <stdexcept>
#include <thread>

#ifdef CERTA_WITH_ZLIB
#include <zlib.h>
#endif

namespace certa {

LineReader::LineReader(const std::string& path) : buf_(1 << 22) {
#ifdef CERTA_WITH_ZLIB
  gz_ = gzopen(path.c_str(), "rb");  // transparently reads plain files too
  if (!gz_) throw std::runtime_error("cannot open " + path);
  gzbuffer(static_cast<gzFile>(gz_), 1 << 20);
#else
  if (path.size() > 3 && path.compare(path.size() - 3, 3, ".gz") == 0)
    throw std::runtime_error("built without zlib; cannot read " + path);
  fp_ = std::fopen(path.c_str(), "rb");
  if (!fp_) throw std::runtime_error("cannot open " + path);
#endif
}

LineReader::~LineReader() {
#ifdef CERTA_WITH_ZLIB
  if (gz_) gzclose(static_cast<gzFile>(gz_));
#endif
  if (fp_) std::fclose(fp_);
}

bool LineReader::fill() {
  if (eof_) return false;
  if (beg_ > 0) {
    std::memmove(buf_.data(), buf_.data() + beg_, end_ - beg_);
    end_ -= beg_;
    beg_ = 0;
  }
  if (end_ == buf_.size()) buf_.resize(buf_.size() * 2);
  long n;
#ifdef CERTA_WITH_ZLIB
  n = gzread(static_cast<gzFile>(gz_), buf_.data() + end_,
             static_cast<unsigned>(buf_.size() - end_));
  if (n < 0) throw std::runtime_error("gzip read error");
#else
  n = static_cast<long>(std::fread(buf_.data() + end_, 1, buf_.size() - end_, fp_));
#endif
  if (n == 0) eof_ = true;
  end_ += static_cast<size_t>(n);
  return n > 0;
}

bool LineReader::getline(std::string& line) {
  for (;;) {
    char* start = buf_.data() + beg_;
    char* nl = static_cast<char*>(std::memchr(start, '\n', end_ - beg_));
    if (nl) {
      size_t len = static_cast<size_t>(nl - start);
      if (len > 0 && start[len - 1] == '\r') --len;
      line.assign(start, len);
      beg_ += static_cast<size_t>(nl - start) + 1;
      return true;
    }
    if (!fill()) {
      if (beg_ == end_) return false;
      size_t len = end_ - beg_;
      if (buf_[end_ - 1] == '\r') --len;
      line.assign(buf_.data() + beg_, len);
      beg_ = end_;
      return true;
    }
  }
}

namespace {

// Parses one 4-line record from [p, end) (lines end with '\n').
void parse_record(const char* p, const char* end, FastqRecord& r) {
  const char* lines[5];
  lines[0] = p;
  for (int i = 1; i <= 4; ++i) {
    const char* nl = static_cast<const char*>(std::memchr(lines[i - 1], '\n', end - lines[i - 1]));
    lines[i] = nl + 1;
  }
  auto line = [&](int i) { return std::string(lines[i], lines[i + 1] - 1); };
  std::string h = line(0);
  if (h.empty() || h[0] != '@') throw std::runtime_error("malformed FASTQ header: " + h);
  size_t sp = h.find_first_of(" \t");
  r.name = h.substr(1, sp == std::string::npos ? std::string::npos : sp - 1);
  r.seq = line(1);
  r.qual = line(3);
  if (r.qual.size() != r.seq.size())
    throw std::runtime_error("SEQ/QUAL length mismatch: " + r.name);
}

}  // namespace

size_t read_fastq_batch(LineReader& in, std::vector<FastqRecord>& out,
                        size_t max_records, int threads) {
  // Serial part: slice whole records into one block (no per-record allocation).
  std::string block, line;
  std::vector<size_t> starts;
  block.reserve(max_records * 64);
  while (starts.size() < max_records) {
    size_t start = block.size();
    int n = 0;
    while (n < 4 && in.getline(line)) {
      if (n == 0 && line.empty()) continue;  // tolerate blank lines between records
      block.append(line);
      block.push_back('\n');
      ++n;
    }
    if (n == 0) break;
    if (n < 4) throw std::runtime_error("truncated FASTQ record at end of input");
    starts.push_back(start);
  }
  starts.push_back(block.size());

  // Parallel part: build records.
  const size_t nrec = starts.size() - 1;
  out.clear();
  out.resize(nrec);
  threads = std::max(1, std::min<int>(threads, static_cast<int>(nrec / 4096) + 1));
  std::vector<std::exception_ptr> errors(threads);
  std::vector<std::thread> pool;
  for (int t = 0; t < threads; ++t) {
    pool.emplace_back([&, t] {
      try {
        const char* base = block.data();
        for (size_t i = nrec * t / threads; i < nrec * (t + 1) / threads; ++i)
          parse_record(base + starts[i], base + starts[i + 1], out[i]);
      } catch (...) {
        errors[t] = std::current_exception();
      }
    });
  }
  for (auto& th : pool) th.join();
  for (auto& e : errors)
    if (e) std::rethrow_exception(e);
  return nrec;
}

}  // namespace certa
