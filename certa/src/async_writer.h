// An output file written from its own thread.
#pragma once
#include <condition_variable>
#include <cstdio>
#include <deque>
#include <mutex>
#include <string>
#include <thread>

namespace certa {

// Writes a FILE from its own thread, in the order chunks are put, so mapping
// never waits on a slow reader: a fallback aligner reading a pipe only a
// batch ahead would otherwise stall the mapper on every full pipe. Up to
// `cap` bytes queue in memory; past that, put() waits.
class AsyncWriter {
 public:
  AsyncWriter(std::FILE* f, size_t cap) : f_(f), cap_(cap) { th_ = std::thread([this] { run(); }); }
  AsyncWriter(const AsyncWriter&) = delete;
  AsyncWriter& operator=(const AsyncWriter&) = delete;
  ~AsyncWriter() { finish(); }
  void put(std::string&& s) {
    if (s.empty()) return;
    std::unique_lock<std::mutex> lk(m_);
    room_.wait(lk, [&] { return queued_ < cap_ || q_.empty(); });
    queued_ += s.size();
    q_.push_back(std::move(s));
    ready_.notify_one();
  }
  // Writes everything queued and stops; false if a write failed.
  bool finish() {
    {
      std::lock_guard<std::mutex> lk(m_);
      done_ = true;
    }
    ready_.notify_one();
    if (th_.joinable()) th_.join();
    return !failed_;
  }

 private:
  void run() {
    for (;;) {
      std::string s;
      bool drained;  // nothing else queued: flush, so a reader gets everything so far
      {
        std::unique_lock<std::mutex> lk(m_);
        ready_.wait(lk, [&] { return !q_.empty() || done_; });
        if (q_.empty()) return;
        s = std::move(q_.front());
        q_.pop_front();
        drained = q_.empty();
      }
      if (std::fwrite(s.data(), 1, s.size(), f_) != s.size()) failed_ = true;
      if (drained) std::fflush(f_);
      {
        std::lock_guard<std::mutex> lk(m_);
        queued_ -= s.size();
      }
      room_.notify_all();
    }
  }
  std::FILE* f_;
  size_t cap_;
  size_t queued_ = 0;
  bool done_ = false;
  bool failed_ = false;
  std::deque<std::string> q_;
  std::mutex m_;
  std::condition_variable ready_, room_;
  std::thread th_;
};

}  // namespace certa
