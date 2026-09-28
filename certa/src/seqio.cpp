#include "seqio.h"

#include <cstring>
#include <stdexcept>

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

size_t read_fastq_batch(LineReader& in, std::vector<FastqRecord>& out,
                        size_t max_records) {
  out.clear();
  std::string h, plus;
  while (out.size() < max_records && in.getline(h)) {
    if (h.empty()) continue;
    if (h[0] != '@') throw std::runtime_error("malformed FASTQ header: " + h);
    FastqRecord r;
    size_t sp = h.find_first_of(" \t");
    r.name = h.substr(1, sp == std::string::npos ? std::string::npos : sp - 1);
    if (!in.getline(r.seq) || !in.getline(plus) || !in.getline(r.qual))
      throw std::runtime_error("truncated FASTQ record: " + r.name);
    if (r.qual.size() != r.seq.size())
      throw std::runtime_error("SEQ/QUAL length mismatch: " + r.name);
    out.push_back(std::move(r));
  }
  return out.size();
}

}  // namespace certa
