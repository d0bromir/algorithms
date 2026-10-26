// Spreading a loop over threads.
#pragma once
#include <algorithm>
#include <atomic>
#include <cstddef>
#include <exception>
#include <thread>
#include <vector>

namespace certa {

// Runs fn(i) for every i in [0, n) on up to `threads` threads, the caller's
// among them. The threads take the items in turn, so a thread that is not
// scheduled (certa map's readers share the cores with its mapping threads)
// holds up one item, not a fixed share of the loop. Rethrows the error of the
// lowest item that failed, as a serial loop would.
template <class F>
void parallel_for(size_t n, int threads, F&& fn) {
  const size_t T = std::min(n, static_cast<size_t>(std::max(1, threads)));
  if (T <= 1) {
    for (size_t i = 0; i < n; ++i) fn(i);
    return;
  }
  std::atomic<size_t> next{0};
  std::vector<std::exception_ptr> errors(n);
  auto work = [&] {
    for (size_t i; (i = next.fetch_add(1, std::memory_order_relaxed)) < n;) {
      try {
        fn(i);
      } catch (...) {
        errors[i] = std::current_exception();
      }
    }
  };
  std::vector<std::thread> pool;
  for (size_t t = 1; t < T; ++t) pool.emplace_back(work);
  work();
  for (auto& th : pool) th.join();
  for (auto& e : errors)
    if (e) std::rethrow_exception(e);
}

}  // namespace certa
