#ifndef MGCVST_MEMORY_PROBE_H
#define MGCVST_MEMORY_PROBE_H

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

namespace mgcvst_memory {

inline double file_number(const std::string& path) {
  std::ifstream in(path.c_str());
  std::string token;
  if (!(in >> token)) return NAN;
  if (token == "max") return std::numeric_limits<double>::infinity();
  char* end = NULL;
  const double value = std::strtod(token.c_str(), &end);
  return (end && *end == '\0' && std::isfinite(value) && value >= 0) ? value : NAN;
}

inline double meminfo_available() {
  std::ifstream in("/proc/meminfo");
  std::string line;
  while (std::getline(in, line)) {
    if (line.compare(0, 13, "MemAvailable:") != 0) continue;
    std::istringstream fields(line.substr(13));
    double kb = NAN;
    std::string unit;
    if ((fields >> kb >> unit) && unit == "kB" && std::isfinite(kb) && kb >= 0) {
      return kb * 1024.0;
    }
  }
  return NAN;
}

// Smallest limit - usage over the process cgroup and its ancestors (v1 or v2).
inline double cgroup_headroom() {
  double best = std::numeric_limits<double>::infinity();
  std::ifstream in("/proc/self/cgroup");
  std::string line;
  while (std::getline(in, line)) {
    const size_t a = line.find(':');
    const size_t b = a == std::string::npos ? a : line.find(':', a + 1);
    if (b == std::string::npos) continue;
    const std::string controllers = line.substr(a + 1, b - a - 1);
    std::string path = line.substr(b + 1);
    bool v2 = controllers.empty();
    bool v1 = false;
    if (!v2) {
      std::istringstream list(controllers);
      std::string item;
      while (std::getline(list, item, ',')) if (item == "memory") v1 = true;
    }
    if (!v1 && !v2) continue;
    const std::string base = v2 ? "/sys/fs/cgroup" : "/sys/fs/cgroup/memory";
    const std::string limit_file = v2 ? "memory.max" : "memory.limit_in_bytes";
    const std::string usage_file = v2 ? "memory.current" : "memory.usage_in_bytes";
    while (path.size() > 1 && path[path.size() - 1] == '/') path.erase(path.size() - 1);
    std::string dir = path == "/" ? base : base + path;
    while (true) {
      const double limit = file_number(dir + "/" + limit_file);
      const double usage = file_number(dir + "/" + usage_file);
      if (std::isfinite(limit) && limit > 0 && limit < 1.15e18 &&
          std::isfinite(usage)) {
        best = std::min(best, std::max(0.0, limit - usage));
      }
      if (dir.size() <= base.size()) break;
      dir = dir.substr(0, dir.find_last_of('/'));
      if (dir.size() < base.size()) break;
    }
  }
  return best;
}

// Available physical memory in bytes; NaN when no signal is readable.
inline double available_physical_memory() {
#ifdef _WIN32
  MEMORYSTATUSEX status;
  status.dwLength = sizeof(status);
  if (GlobalMemoryStatusEx(&status)) return static_cast<double>(status.ullAvailPhys);
  return NAN;
#else
  double value = meminfo_available();
  const double cgroup = cgroup_headroom();
  if (std::isfinite(cgroup)) value = std::isfinite(value) ? std::min(value, cgroup) : cgroup;
  return value;
#endif
}

}  // namespace mgcvst_memory

#endif
