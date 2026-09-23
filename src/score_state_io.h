#ifndef MGCVST_SCORE_STATE_IO_H
#define MGCVST_SCORE_STATE_IO_H

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>
#include <atomic>

namespace mgcvst_io {

struct State {
  uint64_t q = 0;
  std::vector<double> a, upper, width;
  std::vector<std::string> width_names;
  std::string error;
};

struct Trace {
  uint64_t n_ref = 0;
  std::vector<double> a, cross;
  std::string error;
};

inline void require_little_endian() {
  const uint32_t marker = 0x01020304;
  if (*reinterpret_cast<const unsigned char*>(&marker) != 4 || sizeof(double) != 8 ||
      !std::numeric_limits<double>::is_iec559) {
    throw std::runtime_error("Native score shards require little-endian IEEE double storage.");
  }
}

template <typename T> inline void append(std::vector<char>& out, const T& x) {
  const char* p = reinterpret_cast<const char*>(&x);
  out.insert(out.end(), p, p + sizeof(T));
}
inline void append_bytes(std::vector<char>& out, const char* p, size_t n) {
  out.insert(out.end(), p, p + n);
}
inline void append_doubles(std::vector<char>& out, const std::vector<double>& x) {
  if (!x.empty()) append_bytes(out, reinterpret_cast<const char*>(x.data()), x.size() * 8);
}
inline uint64_t hash_bytes(const char* p, size_t n, uint64_t h = 14695981039346656037ULL) {
  for (size_t i = 0; i < n; ++i) {
    h ^= static_cast<unsigned char>(p[i]);
    h *= 1099511628211ULL;
  }
  return h;
}

template <typename T> inline T take(const std::vector<char>& x, size_t& at) {
  if (at > x.size() || x.size() - at < sizeof(T)) throw std::runtime_error("Truncated native score shard.");
  T z;
  std::memcpy(&z, x.data() + at, sizeof(T));
  at += sizeof(T);
  return z;
}
inline std::string take_string(const std::vector<char>& x, size_t& at, uint64_t n) {
  if (n > x.size() - at) throw std::runtime_error("Truncated native score shard string.");
  std::string z(x.data() + at, x.data() + at + static_cast<size_t>(n));
  at += static_cast<size_t>(n);
  return z;
}
inline std::vector<double> take_doubles(const std::vector<char>& x, size_t& at, uint64_t n) {
  if (n > (x.size() - at) / 8) throw std::runtime_error("Truncated native score shard values.");
  std::vector<double> z(static_cast<size_t>(n));
  if (n) std::memcpy(z.data(), x.data() + at, static_cast<size_t>(n) * 8);
  at += static_cast<size_t>(n) * 8;
  if (!std::all_of(z.begin(), z.end(), [](double v) { return std::isfinite(v); })) {
    throw std::runtime_error("Native score shard has non-finite numeric values.");
  }
  return z;
}

struct Record {
  uint32_t type = 0, flags = 0;
  uint64_t q = 0, a_len = 0, width_len = 0, n_ref = 0;
  std::vector<char> body;
};

inline std::string temporary_name(const std::string& path) {
  static std::atomic<uint64_t> serial{0};
  const auto clock = std::chrono::steady_clock::now().time_since_epoch().count();
  const auto thread = std::hash<std::thread::id>{}(std::this_thread::get_id());
  return path + ".tmp-" + std::to_string(clock) + "-" +
    std::to_string(thread) + "-" + std::to_string(serial.fetch_add(1));
}

inline void write_record(const std::string& path, const Record& r) {
  require_little_endian();
  const auto target = std::filesystem::u8path(path);
  if (std::filesystem::exists(target)) throw std::runtime_error("Native shard already exists: " + path);
  std::vector<char> header;
  append_bytes(header, "MGSTIO01", 8);
  append<uint32_t>(header, 1);
  append<uint32_t>(header, 0x01020304);
  append(header, r.type); append(header, r.flags);
  append(header, r.q); append(header, r.a_len); append(header, r.width_len);
  append(header, r.n_ref);
  // The body stores lengths before strings, then numeric payloads.
  append<uint64_t>(header, static_cast<uint64_t>(r.body.size()));
  const uint64_t checksum = hash_bytes(r.body.data(), r.body.size(),
    hash_bytes(header.data(), header.size()));
  const std::string temp = temporary_name(path);
  std::ofstream out(std::filesystem::u8path(temp), std::ios::binary | std::ios::trunc);
  if (!out) throw std::runtime_error("Could not create native shard: " + temp);
  out.write(header.data(), static_cast<std::streamsize>(header.size()));
  out.write(reinterpret_cast<const char*>(&checksum), 8);
  out.write(r.body.data(), static_cast<std::streamsize>(r.body.size()));
  out.close();
  if (!out) throw std::runtime_error("Could not finish native shard: " + temp);
  if (std::filesystem::exists(target)) throw std::runtime_error("Native shard already exists: " + path);
  std::filesystem::rename(std::filesystem::u8path(temp), target);
}

inline Record read_record(const std::string& path, uint32_t expected_type) {
  require_little_endian();
  std::ifstream in(std::filesystem::u8path(path), std::ios::binary | std::ios::ate);
  if (!in) throw std::runtime_error("Could not open native shard: " + path);
  const std::streamoff size = in.tellg();
  if (size < 72) throw std::runtime_error("Native shard is too short: " + path);
  std::vector<char> bytes(static_cast<size_t>(size));
  in.seekg(0);
  in.read(bytes.data(), size);
  if (!in) throw std::runtime_error("Could not read native shard: " + path);
  if (std::memcmp(bytes.data(), "MGSTIO01", 8) != 0) throw std::runtime_error("Native shard magic mismatch: " + path);
  size_t at = 8;
  if (take<uint32_t>(bytes, at) != 1 || take<uint32_t>(bytes, at) != 0x01020304) {
    throw std::runtime_error("Native shard version or byte order mismatch: " + path);
  }
  Record r;
  r.type = take<uint32_t>(bytes, at); r.flags = take<uint32_t>(bytes, at);
  r.q = take<uint64_t>(bytes, at); r.a_len = take<uint64_t>(bytes, at);
  r.width_len = take<uint64_t>(bytes, at); r.n_ref = take<uint64_t>(bytes, at);
  const uint64_t body_len = take<uint64_t>(bytes, at);
  const size_t hash_at = at;
  const uint64_t saved_hash = take<uint64_t>(bytes, at);
  if (r.type != expected_type || r.flags > 1 || body_len != bytes.size() - at) {
    throw std::runtime_error("Native shard type, flags, or length mismatch: " + path);
  }
  const uint64_t actual_hash = hash_bytes(bytes.data() + at, bytes.size() - at,
    hash_bytes(bytes.data(), hash_at));
  if (actual_hash != saved_hash) throw std::runtime_error("Native shard checksum mismatch: " + path);
  bytes.erase(bytes.begin(), bytes.begin() + at);
  r.body = std::move(bytes);
  return r;
}

inline void write_state(const std::string& path, const std::string& signature,
                        const std::string& feature_id, const State& s) {
  Record r; r.type = 1; r.flags = !s.error.empty();
  r.q = s.q; r.a_len = s.a.size(); r.width_len = s.width.size();
  if (signature.empty() || feature_id.empty() ||
      signature.size() > UINT32_MAX || feature_id.size() > UINT32_MAX ||
      s.error.size() > UINT32_MAX || s.width_names.size() != s.width.size()) {
    throw std::runtime_error("Invalid native score-state metadata.");
  }
  if (r.flags == 0 && (s.q == 0 || s.a.size() != s.q ||
      s.q == UINT64_MAX || s.q > UINT64_MAX / (s.q + 1) ||
      s.upper.size() != s.q * (s.q + 1) / 2)) {
    throw std::runtime_error("Invalid native score-state dimensions.");
  }
  append<uint32_t>(r.body, static_cast<uint32_t>(signature.size()));
  append<uint32_t>(r.body, static_cast<uint32_t>(feature_id.size()));
  append<uint32_t>(r.body, static_cast<uint32_t>(s.error.size()));
  append_bytes(r.body, signature.data(), signature.size());
  append_bytes(r.body, feature_id.data(), feature_id.size());
  append_bytes(r.body, s.error.data(), s.error.size());
  append_doubles(r.body, s.a); append_doubles(r.body, s.width);
  for (const auto& name : s.width_names) {
    if (name.size() > UINT32_MAX) throw std::runtime_error("Native width name is too long.");
    append<uint32_t>(r.body, static_cast<uint32_t>(name.size()));
    append_bytes(r.body, name.data(), name.size());
  }
  append_doubles(r.body, s.upper);
  write_record(path, r);
}

inline State read_state(const std::string& path, const std::string& signature,
                       const std::string& feature_id) {
  const Record r = read_record(path, 1);
  size_t at = 0;
  const uint32_t ns = take<uint32_t>(r.body, at), ni = take<uint32_t>(r.body, at);
  const uint32_t ne = take<uint32_t>(r.body, at);
  if (take_string(r.body, at, ns) != signature ||
      take_string(r.body, at, ni) != feature_id) {
    throw std::runtime_error("Native score shard signature or feature ID mismatch: " + path);
  }
  State s; s.q = r.q; s.error = take_string(r.body, at, ne);
  if (static_cast<bool>(r.flags) != !s.error.empty()) throw std::runtime_error("Native score shard error flag mismatch.");
  if (r.flags && (r.q || r.a_len || r.width_len || r.n_ref)) throw std::runtime_error("Invalid native error state.");
  if (!r.flags && (r.q == 0 || r.a_len != r.q || r.n_ref ||
      r.q == UINT64_MAX || r.q > UINT64_MAX / (r.q + 1))) {
    throw std::runtime_error("Invalid native score shard dimensions.");
  }
  s.a = take_doubles(r.body, at, r.a_len);
  s.width = take_doubles(r.body, at, r.width_len);
  for (uint64_t k = 0; k < r.width_len; ++k) {
    const uint32_t length = take<uint32_t>(r.body, at);
    s.width_names.push_back(take_string(r.body, at, length));
  }
  if (!r.flags) s.upper = take_doubles(r.body, at, r.q * (r.q + 1) / 2);
  if (at != r.body.size()) throw std::runtime_error("Native score shard has trailing bytes.");
  return s;
}

inline void write_trace(const std::string& path, const std::string& signature,
                        const std::string& feature_id, const Trace& t) {
  Record r; r.type = 2; r.flags = !t.error.empty();
  r.a_len = t.a.size(); r.n_ref = t.n_ref;
  if (signature.empty() || feature_id.empty() || signature.size() > UINT32_MAX ||
      feature_id.size() > UINT32_MAX || t.error.size() > UINT32_MAX ||
      (!r.flags && (t.n_ref == 0 || t.n_ref > UINT64_MAX / 4 ||
                    t.cross.size() != 4 * t.n_ref))) {
    throw std::runtime_error("Invalid native landmark summary dimensions.");
  }
  append<uint32_t>(r.body, static_cast<uint32_t>(signature.size()));
  append<uint32_t>(r.body, static_cast<uint32_t>(feature_id.size()));
  append<uint32_t>(r.body, static_cast<uint32_t>(t.error.size()));
  append_bytes(r.body, signature.data(), signature.size());
  append_bytes(r.body, feature_id.data(), feature_id.size());
  append_bytes(r.body, t.error.data(), t.error.size());
  append_doubles(r.body, t.a); append_doubles(r.body, t.cross);
  write_record(path, r);
}

inline Trace read_trace(const std::string& path, const std::string& signature,
                        const std::string& feature_id, uint64_t expected_ref) {
  const Record r = read_record(path, 2);
  size_t at = 0;
  const uint32_t ns = take<uint32_t>(r.body, at), ni = take<uint32_t>(r.body, at);
  const uint32_t ne = take<uint32_t>(r.body, at);
  if (take_string(r.body, at, ns) != signature ||
      take_string(r.body, at, ni) != feature_id) {
    throw std::runtime_error("Native trace shard signature or feature ID mismatch: " + path);
  }
  Trace t; t.n_ref = r.n_ref; t.error = take_string(r.body, at, ne);
  if (static_cast<bool>(r.flags) != !t.error.empty() || r.q || r.width_len ||
      r.n_ref != expected_ref || r.n_ref > UINT64_MAX / 4 ||
      (r.flags && r.a_len)) throw std::runtime_error("Native trace shard dimensions mismatch.");
  t.a = take_doubles(r.body, at, r.a_len);
  if (!r.flags) t.cross = take_doubles(r.body, at, 4 * r.n_ref);
  if (at != r.body.size()) throw std::runtime_error("Native trace shard has trailing bytes.");
  return t;
}

}  // namespace mgcvst_io
#endif
