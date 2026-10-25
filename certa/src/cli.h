// What the certa subcommands share: arguments, timing, output files, threads.
#pragma once
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <iterator>
#include <map>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace certa {

using Clock = std::chrono::steady_clock;

// Seconds since `a`.
inline double secs(Clock::time_point a) {
  return std::chrono::duration<double>(Clock::now() - a).count();
}

inline int default_threads() {
  unsigned n = std::thread::hardware_concurrency();
  return n ? static_cast<int>(n) : 1;
}

// Calls work(w) for w in [0, n): w = 0 on this thread, the others on threads
// of their own.
template <class F>
void run_workers(int n, F&& work) {
  std::vector<std::thread> pool;
  for (int w = 1; w < n; ++w) pool.emplace_back([&work, w] { work(w); });
  work(0);
  for (auto& th : pool) th.join();
}

// fopen(path, "wb"), or "cannot open <what>".
inline std::FILE* open_output(const std::string& path, const std::string& what) {
  std::FILE* f = std::fopen(path.c_str(), "wb");
  if (!f) throw std::runtime_error("cannot open " + what);
  return f;
}

// Positional arguments, and options given as "-x value" (flags as "-x").
struct Args {
  std::vector<std::string> pos;
  std::map<std::string, std::string> opt;
  bool has(const std::string& k) const { return opt.count(k) > 0; }
  std::string get(const std::string& k, const std::string& d) const {
    auto it = opt.find(k);
    return it == opt.end() ? d : it->second;
  }
  int geti(const std::string& k, int d) const {
    return has(k) ? std::atoi(opt.at(k).c_str()) : d;
  }

  // argv[start..argc). Options take the next argument as their value, except
  // the flags, which take none.
  static Args parse(int argc, char** argv, int start) {
    static const char* flags[] = {"--gpu", "--unc-tags", "--intervals", "--witness", "--serial-stages", "--full-cert"};
    Args a;
    for (int i = start; i < argc; ++i) {
      std::string s = argv[i];
      if (s.size() > 1 && s[0] == '-') {
        bool flag = std::find(std::begin(flags), std::end(flags), s) != std::end(flags);
        if (flag) { a.opt[s] = "1"; continue; }
        if (i + 1 >= argc) throw std::runtime_error("missing value for " + s);
        a.opt[s] = argv[++i];
      } else {
        a.pos.push_back(s);
      }
    }
    return a;
  }
};

// The subcommands (main.cpp; certa map in cmd_map.cpp).
int cmd_index(const Args& a);
int cmd_map(const Args& a);
int cmd_repeat(const Args& a);
int cmd_check(const Args& a);
int cmd_check_index(const Args& a);

}  // namespace certa
