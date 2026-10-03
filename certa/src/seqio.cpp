#include "seqio.h"

#include <algorithm>
#include <cstring>
#include <exception>
#include <stdexcept>
#include <thread>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

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

long LineReader::append_line(std::string& out) {
  for (;;) {
    char* start = buf_.data() + beg_;
    char* nl = static_cast<char*>(std::memchr(start, '\n', end_ - beg_));
    size_t len;
    if (nl) {
      len = static_cast<size_t>(nl - start);
      beg_ += len + 1;
    } else if (!fill()) {  // fill() may move or reallocate the buffer
      if (beg_ == end_) return -1;
      start = buf_.data() + beg_;
      len = end_ - beg_;
      beg_ = end_;
    } else {
      continue;
    }
    if (len > 0 && start[len - 1] == '\r') --len;
    out.append(start, len);
    out.push_back('\n');
    return static_cast<long>(len);
  }
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

// Parses one 4-line record from [p, end). Lines end with '\n' (optionally
// "\r\n"); the last line may end at `end` without a newline.
void parse_record(const char* p, const char* end, FastqRecord& r) {
  const char* beg[4];
  const char* stop[4];
  const char* at = p;
  for (int i = 0; i < 4; ++i) {
    const char* nl = static_cast<const char*>(std::memchr(at, '\n', end - at));
    const char* e = nl ? nl : end;
    beg[i] = at;
    stop[i] = (e > at && e[-1] == '\r') ? e - 1 : e;
    at = nl ? nl + 1 : end;
  }
  auto line = [&](int i) { return std::string(beg[i], stop[i]); };
  std::string h = line(0);
  if (h.empty() || h[0] != '@') throw std::runtime_error("malformed FASTQ header: " + h);
  size_t sp = h.find_first_of(" \t");
  r.name = h.substr(1, sp == std::string::npos ? std::string::npos : sp - 1);
  r.seq = line(1);
  r.qual = line(3);
  if (r.qual.size() != r.seq.size())
    throw std::runtime_error("SEQ/QUAL length mismatch: " + r.name);
}

// Builds records [starts[i], starts[i+1]) from `base` on several threads.
void parse_records(const char* base, const std::vector<size_t>& starts,
                   std::vector<FastqRecord>& out, int threads) {
  const size_t nrec = starts.size() - 1;
  out.clear();
  out.resize(nrec);
  threads = std::max(1, std::min<int>(threads, static_cast<int>(nrec / 4096) + 1));
  std::vector<std::exception_ptr> errors(threads);
  std::vector<std::thread> pool;
  for (int t = 0; t < threads; ++t) {
    pool.emplace_back([&, t] {
      try {
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
}

}  // namespace

MappedFastq::MappedFastq(const std::string& path) {
  int fd = ::open(path.c_str(), O_RDONLY);
  if (fd < 0) throw std::runtime_error("cannot open " + path);
  struct stat st;
  if (fstat(fd, &st) != 0) {
    ::close(fd);
    throw std::runtime_error("cannot stat " + path);
  }
  size_ = static_cast<size_t>(st.st_size);
  if (size_ > 0) {
    void* m = mmap(nullptr, size_, PROT_READ, MAP_PRIVATE, fd, 0);
    if (m == MAP_FAILED) {
      ::close(fd);
      throw std::runtime_error("cannot mmap " + path);
    }
    madvise(m, size_, MADV_SEQUENTIAL);
    base_ = static_cast<const char*>(m);
  }
  ::close(fd);
}

MappedFastq::~MappedFastq() {
  if (base_) munmap(const_cast<char*>(base_), size_);
}

bool MappedFastq::usable(const std::string& path) {
  if (path.size() > 3 && path.compare(path.size() - 3, 3, ".gz") == 0) return false;
  struct stat st;
  return ::stat(path.c_str(), &st) == 0 && S_ISREG(st.st_mode);
}

size_t MappedFastq::next_batch(std::vector<FastqRecord>& out, size_t max_records, int threads) {
  // Serial part: only find record boundaries (memchr over the mapping).
  std::vector<size_t> starts;
  starts.reserve(max_records + 1);
  while (starts.size() < max_records && pos_ < size_) {
    // Skip blank lines between records.
    while (pos_ < size_ && (base_[pos_] == '\n' || base_[pos_] == '\r')) ++pos_;
    if (pos_ >= size_) break;
    const size_t start = pos_;
    int n = 0;
    while (n < 4 && pos_ < size_) {
      const char* nl = static_cast<const char*>(std::memchr(base_ + pos_, '\n', size_ - pos_));
      pos_ = nl ? static_cast<size_t>(nl - base_) + 1 : size_;
      ++n;
    }
    if (n < 4) throw std::runtime_error("truncated FASTQ record at end of input");
    starts.push_back(start);
  }
  starts.push_back(pos_);
  parse_records(base_, starts, out, threads);
  return out.size();
}

size_t read_fastq_batch(LineReader& in, std::vector<FastqRecord>& out,
                        size_t max_records, int threads) {
  // Serial part: slice whole records into one block (no per-record allocation).
  std::string block;
  std::vector<size_t> starts;
  block.reserve(max_records * 400);  // ~2 x 150 bp + name, typical short reads
  starts.reserve(max_records + 1);
  while (starts.size() < max_records) {
    size_t start = block.size();
    int n = 0;
    long len;
    while (n < 4 && (len = in.append_line(block)) >= 0) {
      if (n == 0 && len == 0) {  // tolerate blank lines between records
        block.resize(start);
        continue;
      }
      ++n;
    }
    if (n == 0) break;
    if (n < 4) throw std::runtime_error("truncated FASTQ record at end of input");
    starts.push_back(start);
  }
  starts.push_back(block.size());
  parse_records(block.data(), starts, out, threads);  // parallel part
  return out.size();
}

}  // namespace certa
