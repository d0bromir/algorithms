// certa-qc: fastp's QC, trimming and filtering and certa map in one process.
//
//   certa-qc ref.cidx [certa map options] --fastp [fastp options]
//
// fastp reads the input (its -i/-I), and the reads it keeps go straight to
// certa map's pipeline: nothing is compressed, written or read back in
// between. The SAM, the uncertified FASTQ and both reports are what
//   fastp -o trimmed_R1.fq -O trimmed_R2.fq [fastp options]
//   certa map ref.cidx trimmed_R1.fq trimmed_R2.fq [certa map options]
// would write (tests/qc_equivalence.sh checks this). Reads fastp does not
// keep go to --failed_out and --unpaired1/2, as in fastp.
//
// fastp runs on its own worker threads (its -w), ahead of the mapper (-t),
// so trimming overlaps loading the index and mapping.

#include <condition_variable>
#include <cstdio>
#include <cstring>
#include <deque>
#include <exception>
#include <map>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#ifdef CERTA_WITH_ZLIB
#include <zlib.h>
#endif

#include "cli.h"
#include "fastp_core.h"
#include "map_config.h"
#include "map_pipeline.h"

namespace {

const char* kUsage =
    "usage: certa-qc ref.cidx [certa map options] --fastp [fastp options]\n"
    "  The input is fastp's -i (and -I, or --interleaved_in). certa map's options are as in\n"
    "  `certa map`, without the FASTQ inputs. fastp's options are as in fastp, without -o/-O,\n"
    "  --stdout, --merge and --split: the reads fastp keeps are mapped instead of written.\n";

bool ends_with(const std::string& s, const char* suffix) {
  const size_t n = std::strlen(suffix);
  return s.size() >= n && s.compare(s.size() - n, n, suffix) == 0;
}

// A FASTQ file in fastp's format: plain, or gzip when the name ends in .gz.
class FastqFile {
 public:
  explicit FastqFile(const std::string& path) {
    if (ends_with(path, ".gz")) {
#ifdef CERTA_WITH_ZLIB
      gz_ = gzopen(path.c_str(), "wb4");
      if (!gz_) throw std::runtime_error("cannot open " + path);
      return;
#else
      throw std::runtime_error(path + ": gzip output needs certa built with zlib");
#endif
    }
    fp_ = certa::open_output(path, path);
  }
  ~FastqFile() { close(); }
  void write(const fastpcore::Record& r, bool with_reason) {
    buf_.clear();
    buf_ += r.name;
    if (with_reason) {
      buf_ += ' ';
      buf_ += r.reason;
    }
    buf_ += '\n';
    buf_ += r.seq;
    buf_ += '\n';
    buf_ += r.strand;
    buf_ += '\n';
    buf_ += r.qual;
    buf_ += '\n';
    put(buf_);
  }
  void close() {
#ifdef CERTA_WITH_ZLIB
    if (gz_ && gzclose(static_cast<gzFile>(gz_)) != Z_OK) failed_ = true;
    gz_ = nullptr;
#endif
    if (fp_ && std::fclose(fp_) != 0) failed_ = true;
    fp_ = nullptr;
  }
  bool failed() const { return failed_; }

 private:
  void put(const std::string& s) {
#ifdef CERTA_WITH_ZLIB
    if (gz_) {
      if (gzwrite(static_cast<gzFile>(gz_), s.data(), static_cast<unsigned>(s.size())) != static_cast<int>(s.size()))
        failed_ = true;
      return;
    }
#endif
    if (std::fwrite(s.data(), 1, s.size(), fp_) != s.size()) failed_ = true;
  }
  std::FILE* fp_ = nullptr;
  void* gz_ = nullptr;
  bool failed_ = false;
  std::string buf_;
};

// The files fastp would write for the reads it does not keep.
struct ExtraOutputs {
  std::unique_ptr<FastqFile> failed, unpaired1, unpaired2, overlapped;
  void write(const fastpcore::BatchResult& r) {
    if (failed)
      for (const auto& x : r.failed) failed->write(x, true);
    if (unpaired1)
      for (const auto& x : r.unpaired1) unpaired1->write(x, false);
    if (unpaired2)
      for (const auto& x : r.unpaired2) unpaired2->write(x, false);
    if (overlapped)
      for (const auto& x : r.overlapped) overlapped->write(x, false);
  }
  bool close() {
    bool ok = true;
    for (FastqFile* f : {failed.get(), unpaired1.get(), unpaired2.get(), overlapped.get()})
      if (f) {
        f->close();
        ok = ok && !f->failed();
      }
    return ok;
  }
};

// certa's record name: the FASTQ name without '@' and without its comment
// (seqio.cpp parse_spans).
void set_name(std::string& out, const std::string& fastq_name) {
  size_t b = !fastq_name.empty() && fastq_name[0] == '@' ? 1 : 0;
  size_t e = b;
  while (e < fastq_name.size() && fastq_name[e] != ' ' && fastq_name[e] != '\t') ++e;
  out.assign(fastq_name, b, e - b);
}

// fastp's reads, in input order, for certa map. One thread reads batches,
// fastp's worker slots process them, and next() hands the kept reads over
// batch after batch. At most kAhead batches are read ahead of the mapper.
class FastpSource : public certa::RecordSource {
 public:
  FastpSource(fastpcore::Session& s, ExtraOutputs& extra) : s_(s), extra_(extra) {
    ahead_ = 2 * s_.workers() + 2;
    reader_ = std::thread([this] { guard([this] { read_loop(); }); });
    for (int w = 0; w < s_.workers(); ++w) workers_.emplace_back([this, w] { guard([this, w] { work_loop(w); }); });
  }
  ~FastpSource() { stop(); }

  size_t next(std::vector<certa::FastqRecord>& out, size_t max_records) override {
    const bool paired = s_.paired();
    if (paired) max_records -= max_records % 2;
    size_t n = 0;
    while (n < max_records) {
      if (pos_ >= cur_.passed.size() && !take_next_batch()) break;
      const size_t take = std::min(max_records - n, cur_.passed.size() - pos_);
      if (out.size() < n + take) out.resize(n + take);
      for (size_t i = 0; i < take; ++i) {
        const fastpcore::Record& r = cur_.passed[pos_ + i];
        certa::FastqRecord& f = out[n + i];
        set_name(f.name, r.name);
        f.seq = r.seq;
        f.qual = r.qual;
      }
      pos_ += take;
      n += take;
    }
    out.resize(n);
    return n;
  }

  // After the mapper has taken every read: joins the threads and writes
  // fastp's summary and reports.
  void finish() {
    stop();
    rethrow();
    s_.finish();
  }

 private:
  // A thread's exception is kept for the mapper's thread, and stops the others.
  template <class F>
  void guard(F&& f) {
    try {
      f();
    } catch (...) {
      std::lock_guard<std::mutex> lock(m_);
      if (!error_) error_ = std::current_exception();
      stop_ = true;
      cv_.notify_all();
    }
  }
  void rethrow() {
    std::lock_guard<std::mutex> lock(m_);
    if (error_) std::rethrow_exception(error_);
  }

  void read_loop() {
    while (true) {
      {
        std::unique_lock<std::mutex> lock(m_);
        cv_.wait(lock, [&] { return stop_ || read_ - delivered_ < ahead_; });
        if (stop_) return;
      }
      fastpcore::Batch* b = s_.readBatch();
      std::lock_guard<std::mutex> lock(m_);
      if (!b) {
        reading_done_ = true;
        cv_.notify_all();
        return;
      }
      todo_.push_back(b);
      ++read_;
      cv_.notify_all();
    }
  }

  void work_loop(int w) {
    while (true) {
      fastpcore::Batch* b = nullptr;
      {
        std::unique_lock<std::mutex> lock(m_);
        cv_.wait(lock, [&] { return stop_ || !todo_.empty() || reading_done_; });
        if (stop_ || todo_.empty()) return;
        b = todo_.front();
        todo_.pop_front();
      }
      const long number = fastpcore::Session::batchNumber(b);
      fastpcore::BatchResult r;
      s_.process(b, w, r);
      std::lock_guard<std::mutex> lock(m_);
      done_[number] = std::move(r);
      cv_.notify_all();
    }
  }

  // Moves the next batch's result into cur_; false at the end of the input.
  bool take_next_batch() {
    std::unique_lock<std::mutex> lock(m_);
    cv_.wait(lock, [&] { return stop_ || done_.count(delivered_) || (reading_done_ && delivered_ == read_); });
    if (error_) std::rethrow_exception(error_);
    if (stop_ || !done_.count(delivered_)) return false;
    cur_ = std::move(done_[delivered_]);
    done_.erase(delivered_);
    ++delivered_;
    pos_ = 0;
    cv_.notify_all();
    lock.unlock();
    extra_.write(cur_);  // in input order, on the mapper's reading thread
    return true;
  }

  void stop() {
    {
      std::lock_guard<std::mutex> lock(m_);
      if (joined_) return;
      // Normal end: every batch was delivered, the threads are already done.
      // Early end (an error): wake and stop them.
      if (!(reading_done_ && delivered_ == read_)) stop_ = true;
      cv_.notify_all();
    }
    if (reader_.joinable()) reader_.join();
    {
      std::lock_guard<std::mutex> lock(m_);
      reading_done_ = true;
      cv_.notify_all();
    }
    for (auto& t : workers_)
      if (t.joinable()) t.join();
    std::lock_guard<std::mutex> lock(m_);
    joined_ = true;
  }

  fastpcore::Session& s_;
  ExtraOutputs& extra_;
  long ahead_ = 4;
  std::mutex m_;
  std::condition_variable cv_;
  std::deque<fastpcore::Batch*> todo_;
  std::map<long, fastpcore::BatchResult> done_;
  long read_ = 0, delivered_ = 0;
  bool reading_done_ = false, stop_ = false, joined_ = false;
  std::exception_ptr error_;
  fastpcore::BatchResult cur_;
  size_t pos_ = 0;
  std::thread reader_;
  std::vector<std::thread> workers_;
};

// fastp's options that certa-qc reads for itself.
struct FastpArgs {
  std::vector<std::string> args;
  std::string in1, in2, failed, unpaired1, unpaired2, overlapped;
  bool interleaved = false;

  explicit FastpArgs(const std::vector<std::string>& a) : args(a) {
    // fastp also takes "--option=value".
    std::string inline_value;
    bool has_inline = false;
    auto value = [&](size_t& i) -> std::string {
      if (has_inline) return inline_value;
      if (i + 1 >= args.size()) throw std::runtime_error("fastp option " + args[i] + " needs a value");
      return args[++i];
    };
    for (size_t i = 0; i < args.size(); ++i) {
      std::string s = args[i];
      has_inline = false;
      const size_t eq = s.find('=');
      if (s.compare(0, 2, "--") == 0 && eq != std::string::npos) {
        inline_value = s.substr(eq + 1);
        s.resize(eq);
        has_inline = true;
      }
      if (s == "-o" || s == "--out1" || s == "-O" || s == "--out2" || s == "--stdout")
        throw std::runtime_error("fastp option " + s + ": the reads fastp keeps are mapped, not written");
      if (s == "-m" || s == "--merge" || s == "--merged_out" || s == "--split" || s == "--split_by_lines" ||
          s == "-s" || s == "-S")
        throw std::runtime_error("fastp option " + s + " is not supported by certa-qc");
      if (s == "-i" || s == "--in1") in1 = value(i);
      else if (s == "-I" || s == "--in2") in2 = value(i);
      else if (s == "--failed_out") failed = value(i);
      else if (s == "--unpaired1") unpaired1 = value(i);
      else if (s == "--unpaired2") unpaired2 = value(i);
      else if (s == "--overlapped_out") overlapped = value(i);
      else if (s == "--interleaved_in") interleaved = true;
    }
    if (in1.empty()) throw std::runtime_error("fastp needs an input: -i reads.fq[.gz]");
  }
};

int run(int argc, char** argv) {
  int sep = 1;
  while (sep < argc && std::strcmp(argv[sep], "--fastp") != 0) ++sep;
  if (argc < 2 || sep >= argc) {
    std::fputs(kUsage, stderr);
    return 1;
  }
  FastpArgs fa(std::vector<std::string>(argv + sep + 1, argv + argc));

  // fastp parses its own options and evaluates the input (adapters, read
  // length, two-colour chemistry) before anything is mapped.
  std::vector<std::string> fastp_argv_s = {"fastp"};
  fastp_argv_s.insert(fastp_argv_s.end(), fa.args.begin(), fa.args.end());
  std::vector<char*> fastp_argv;
  for (auto& s : fastp_argv_s) fastp_argv.push_back(&s[0]);
  fastpcore::Session session(static_cast<int>(fastp_argv.size()), fastp_argv.data());

  // certa map's options, with fastp's inputs standing for the FASTQ files
  // (they only tell certa map whether the input is paired).
  certa::Args a = certa::Args::parse(sep, argv, 1);
  if (a.pos.size() != 1) throw std::runtime_error(std::string("certa map takes only the index here\n") + kUsage);
  a.pos.push_back(fa.in1);
  if (session.paired()) a.pos.push_back(fa.in2.empty() ? fa.in1 : fa.in2);
  const certa::MapConfig c = certa::MapConfig::from_args(a);

  ExtraOutputs extra;
  if (!fa.failed.empty()) extra.failed.reset(new FastqFile(fa.failed));
  if (!fa.unpaired1.empty()) extra.unpaired1.reset(new FastqFile(fa.unpaired1));
  if (!fa.unpaired2.empty() && fa.unpaired2 != fa.unpaired1) extra.unpaired2.reset(new FastqFile(fa.unpaired2));
  if (!fa.overlapped.empty()) extra.overlapped.reset(new FastqFile(fa.overlapped));

  FastpSource source(session, extra);
  const int rc = certa::run_map(c, &source);
  source.finish();
  if (!extra.close()) throw std::runtime_error("writing a fastp output failed");
  return rc;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    return run(argc, argv);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "certa-qc: error: %s\n", e.what());
    return 1;
  }
}
