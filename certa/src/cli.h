// Command-line arguments and timing, shared by the certa subcommands.
#pragma once
#include <algorithm>
#include <chrono>
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
