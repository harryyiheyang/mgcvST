#define EIGEN_DONT_PARALLELIZE
#include "inla_sparse.h"
#include "liu_tail.h"
#include "memory_probe.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <memory>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

using namespace mgcvst_sparse;

namespace {

using FMat = Eigen::MatrixXf;
using Half = Eigen::half;

// Resident per-gene exact Liu state: projected score a (double) and the upper
// triangle of the reduced curvature M, column-major with the diagonal, stored
// as two-byte halves.
struct HalfCache {
  int n = 0;
  int r = 0;
  Mat a;
  std::vector<std::vector<Half> > m;
  std::vector<int> state;  // 0 empty, 1 ready, 2 failed
  std::vector<std::string> error;
  double bytes = 0;
};

const char kMagic[8] = {'M', 'G', 'S', 'T', 'F', '1', '6', '\0'};
const int32_t kVersion = 1;

inline size_t packed_length(int r) { return (size_t)r * (r + 1) / 2; }
inline double gene_bytes(int r) { return 2.0 * packed_length(r) + 8.0 * r; }

HalfCache* cache_pointer(SEXP pointer) {
  if (TYPEOF(pointer) != EXTPTRSXP || R_ExternalPtrAddr(pointer) == NULL) {
    Rcpp::stop("cache must be a live fp16 score-state cache.");
  }
  Rcpp::XPtr<HalfCache> cache(pointer);
  return cache.get();
}

void check_threads(int threads) {
  if (threads < 1) Rcpp::stop("threads must be positive.");
#ifndef _OPENMP
  if (threads > 1) Rcpp::stop("mgcvST was compiled without OpenMP; use threads = 1.");
#endif
}

// Coefficients, offsets and family parameters of k features.
struct Compact {
  SpMat A;
  SpMat X;
  Rcpp::NumericMatrix B;
  Rcpp::NumericMatrix C;
  Rcpp::NumericMatrix O;
  Rcpp::IntegerVector family;
  Rcpp::NumericVector size;
  Rcpp::NumericVector dispersion;
  int k = 0;

  FeatureModel model(int f) const {
    FeatureModel z;
    z.b = &B(0, f);
    z.c = C.nrow() ? &C(0, f) : NULL;
    z.offset = &O(0, O.ncol() == 1 ? 0 : f);
    z.family = family[f];
    z.size = size[f];
    z.dispersion = dispersion[f];
    return z;
  }
};

Compact parse_compact(const Eigen::MappedSparseMatrix<double>& A,
                      const Eigen::MappedSparseMatrix<double>& X,
                      const Rcpp::NumericMatrix& B, const Rcpp::NumericMatrix& C,
                      const Rcpp::NumericMatrix& O, const Rcpp::IntegerVector& family,
                      const Rcpp::NumericVector& size,
                      const Rcpp::NumericVector& dispersion) {
  Compact z;
  z.A = A;
  z.X = X;
  z.B = B;
  z.C = C;
  z.O = O;
  z.family = family;
  z.size = size;
  z.dispersion = dispersion;
  z.k = B.ncol();
  const int n = A.rows();
  if (X.rows() != n || B.nrow() != A.cols() || C.nrow() != X.cols() ||
      C.ncol() != z.k || O.nrow() != n || (O.ncol() != 1 && O.ncol() != z.k) ||
      family.size() != z.k || size.size() != z.k || dispersion.size() != z.k) {
    Rcpp::stop("The compact feature coefficients, offsets and family parameters are misaligned.");
  }
  for (int f = 0; f < z.k; ++f) {
    if (family[f] < 0 || family[f] > 2) Rcpp::stop("Unknown feature family code.");
  }
  return z;
}

Mat penalty_matrix(const Rcpp::Nullable<Rcpp::NumericMatrix>& precision, int px, int k) {
  Mat out = Mat::Zero(px, k);
  if (precision.isNotNull()) out = Rcpp::as<Mat>(precision.get());
  if (out.rows() != px || out.cols() != k || !out.allFinite() || (out.array() < 0).any()) {
    Rcpp::stop("nuisance_precision must be a finite non-negative ncol(X)-by-feature matrix.");
  }
  return out;
}

void unpack(const std::vector<Half>& x, FMat& m) {
  size_t at = 0;
  for (int j = 0; j < m.cols(); ++j) {
    for (int i = 0; i < j; ++i) m(i, j) = m(j, i) = static_cast<float>(x[at++]);
    m(j, j) = static_cast<float>(x[at++]);
  }
}

template <class T> void put(std::ofstream& out, T value) {
  out.write(reinterpret_cast<const char*>(&value), sizeof(T));
}
template <class T> T get(std::ifstream& in) {
  T value;
  in.read(reinterpret_cast<char*>(&value), sizeof(T));
  if (!in) throw std::runtime_error("truncated");
  return value;
}
void put_string(std::ofstream& out, const std::string& x) {
  put<int32_t>(out, (int32_t)x.size());
  out.write(x.data(), x.size());
}
std::string get_string(std::ifstream& in) {
  const int32_t length = get<int32_t>(in);
  if (length < 0 || length > (1 << 24)) throw std::runtime_error("invalid string");
  std::string x(length, '\0');
  in.read(&x[0], length);
  if (!in) throw std::runtime_error("truncated");
  return x;
}

double available_bytes(double budget) {
  double free = mgcvst_memory::available_physical_memory();
  if (std::isfinite(budget)) free = std::isfinite(free) ? std::min(free, budget) : budget;
  return free;
}

}  // namespace

// Shared working-state step: eta = offset + X c + A b, mu and working
// variance for each feature from its saved coefficients and family parameters.
// [[Rcpp::export]]
Rcpp::List mgcvst_inla_working_state_cpp(
    const Eigen::MappedSparseMatrix<double>& A,
    const Eigen::MappedSparseMatrix<double>& X,
    const Rcpp::NumericMatrix& B, const Rcpp::NumericMatrix& C,
    const Rcpp::NumericMatrix& O, const Rcpp::IntegerVector& family,
    const Rcpp::NumericVector& size, const Rcpp::NumericVector& dispersion,
    int threads = 1) {
  check_threads(threads);
  const Compact z = parse_compact(A, X, B, C, O, family, size, dispersion);
  const int n = A.rows();
  Rcpp::NumericMatrix eta(n, z.k), mu(n, z.k), variance(n, z.k);
  std::vector<std::string> error(z.k);
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
#endif
  for (int f = 0; f < z.k; ++f) {
    Vec e, m, v;
    try {
      recover_working(z.A, z.X, z.model(f), e, m, v);
      std::memcpy(&eta(0, f), e.data(), sizeof(double) * n);
      std::memcpy(&mu(0, f), m.data(), sizeof(double) * n);
      std::memcpy(&variance(0, f), v.data(), sizeof(double) * n);
    } catch (const std::exception& x) {
      error[f] = x.what();
      for (int i = 0; i < n; ++i) eta(i, f) = mu(i, f) = variance(i, f) = NA_REAL;
    }
  }
  Rcpp::CharacterVector error_out(z.k);
  for (int f = 0; f < z.k; ++f) {
    error_out[f] = error[f].empty() ? NA_STRING : Rcpp::String(error[f]);
  }
  return Rcpp::List::create(Rcpp::Named("eta") = eta, Rcpp::Named("mu") = mu,
                            Rcpp::Named("variance") = variance,
                            Rcpp::Named("error") = error_out);
}

// Serializable double reconstruction units from saved coefficients; a is the
// saved full-space score of each feature.
// [[Rcpp::export]]
Rcpp::List mgcvst_inla_compact_units_cpp(
    const Eigen::MappedSparseMatrix<double>& A,
    const Eigen::MappedSparseMatrix<double>& Q_map,
    const Eigen::Map<Eigen::VectorXd> constraint,
    const Eigen::MappedSparseMatrix<double>& X,
    const Rcpp::NumericMatrix& B, const Rcpp::NumericMatrix& C,
    const Rcpp::NumericMatrix& O, const Rcpp::IntegerVector& family,
    const Rcpp::NumericVector& size, const Rcpp::NumericVector& dispersion,
    const Eigen::Map<Eigen::VectorXd> tau, const Eigen::Map<Eigen::MatrixXd> a,
    int threads = 1, SEXP prepared = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericMatrix> nuisance_precision = R_NilValue) {
  check_threads(threads);
  const Compact z = parse_compact(A, X, B, C, O, family, size, dispersion);
  const SpMat Q = Q_map;
  const Vec g = constraint;
  const int m = A.cols();
  const SparsePrepared* cache = prepared_cache(prepared, m, g);
  const Mat penalty = penalty_matrix(nuisance_precision, X.cols(), z.k);
  if (Q.rows() != m || tau.size() != z.k || a.rows() != m || a.cols() != z.k) {
    Rcpp::stop("Q, tau and a must align with the compact features.");
  }
  const double constraint_norm = std::abs(cache->coordinate_constraint.norm() - 1.0);
  std::vector<ReconstructionUnit> units(z.k);
  std::vector<std::string> error(z.k);
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
#endif
  for (int f = 0; f < z.k; ++f) {
    try {
      if (!std::isfinite(tau[f]) || tau[f] <= 0) {
        throw std::runtime_error("The feature has invalid dispersion or smoothing parameters.");
      }
      if (!a.col(f).allFinite()) throw std::runtime_error("The saved score vector is non-finite.");
      Vec eta, mu, V;
      recover_working(z.A, z.X, z.model(f), eta, mu, V);
      curvature_unit(z.A, z.X, Q, g, V, tau[f], penalty.col(f).data(), units[f]);
      units[f].a = a.col(f);
      units[f].constraint_norm = constraint_norm;
    } catch (const std::exception& x) {
      error[f] = x.what();
    }
  }
  Rcpp::List out(z.k);
  for (int f = 0; f < z.k; ++f) {
    if (!error[f].empty()) {
      out[f] = Rcpp::List::create(Rcpp::Named("error") = error[f]);
      continue;
    }
    const ReconstructionUnit& u = units[f];
    out[f] = Rcpp::List::create(
      Rcpp::Named("a") = u.a, Rcpp::Named("statistic") = u.a.squaredNorm(),
      Rcpp::Named("expected_vp") = u.Vp, Rcpp::Named("tau") = u.tau,
      Rcpp::Named("K") = u.K, Rcpp::Named("U") = u.U,
      Rcpp::Named("H_L") = u.H_L, Rcpp::Named("H_D") = u.H_D,
      Rcpp::Named("H_perm") = u.H_perm, Rcpp::Named("hinv_g") = u.hinv_g,
      Rcpp::Named("hden") = u.hden, Rcpp::Named("width") = m,
      Rcpp::Named("normalization") = m - 1,
      Rcpp::Named("constraint_diagnostic") = u.constraint_norm
    );
  }
  return out;
}

// [[Rcpp::export]]
SEXP mgcvst_fp16_cache_cpp(int n, int r) {
  if (n < 1 || r < 1) Rcpp::stop("The fp16 cache needs at least one gene and one coordinate.");
  std::unique_ptr<HalfCache> cache(new HalfCache());
  cache->n = n;
  cache->r = r;
  cache->a = Mat::Constant(r, n, NA_REAL);
  cache->m.resize(n);
  cache->state.assign(n, 0);
  cache->error.resize(n);
  Rcpp::XPtr<HalfCache> pointer(cache.release(), true);
  pointer.attr("class") = "mgcvst_fp16_cache";
  return pointer;
}

// [[Rcpp::export]]
void mgcvst_fp16_cache_release_cpp(SEXP pointer) {
  if (TYPEOF(pointer) != EXTPTRSXP || R_ExternalPtrAddr(pointer) == NULL) return;
  Rcpp::XPtr<HalfCache> cache(pointer);
  cache.release();
}

// [[Rcpp::export]]
Rcpp::List mgcvst_fp16_cache_info_cpp(SEXP pointer) {
  const HalfCache* cache = cache_pointer(pointer);
  Rcpp::CharacterVector error(cache->n);
  for (int f = 0; f < cache->n; ++f) {
    error[f] = cache->error[f].empty() ? NA_STRING : Rcpp::String(cache->error[f]);
  }
  return Rcpp::List::create(
    Rcpp::Named("state") = Rcpp::wrap(cache->state), Rcpp::Named("error") = error,
    Rcpp::Named("bytes") = cache->bytes, Rcpp::Named("n") = cache->n,
    Rcpp::Named("r") = cache->r, Rcpp::Named("gene_bytes") = gene_bytes(cache->r)
  );
}

// Build exact Liu states of features into cache slots (1-based): W from the
// shared working-state step, M in double in the common projected space, then
// fp16 upper-triangle packing. budget caps the detected available memory
// (NaN: detected only); admission stops at the safe line
// 0.7 * available - (per-thread build and pair workspace + 2 GiB headroom).
// [[Rcpp::export]]
Rcpp::List mgcvst_fp16_build_cpp(
    SEXP pointer, const Rcpp::IntegerVector& slots,
    const Eigen::MappedSparseMatrix<double>& A,
    const Eigen::MappedSparseMatrix<double>& Q_map,
    const Eigen::Map<Eigen::VectorXd> constraint,
    const Eigen::MappedSparseMatrix<double>& X,
    const Rcpp::NumericMatrix& B, const Rcpp::NumericMatrix& C,
    const Rcpp::NumericMatrix& O, const Rcpp::IntegerVector& family,
    const Rcpp::NumericVector& size, const Rcpp::NumericVector& dispersion,
    const Eigen::Map<Eigen::VectorXd> tau, const Eigen::Map<Eigen::MatrixXd> a,
    const Eigen::Map<Eigen::MatrixXd> coordinate,
    const Eigen::Map<Eigen::MatrixXd> basis,
    double budget, int threads = 1,
    Rcpp::Nullable<Rcpp::NumericMatrix> nuisance_precision = R_NilValue) {
  check_threads(threads);
  HalfCache* cache = cache_pointer(pointer);
  const Compact z = parse_compact(A, X, B, C, O, family, size, dispersion);
  const SpMat Q = Q_map;
  const Vec g = constraint;
  const int n = A.rows(), m = A.cols(), r = cache->r, px = X.cols();
  const Mat penalty = penalty_matrix(nuisance_precision, px, z.k);
  if (slots.size() != z.k || Q.rows() != m || tau.size() != z.k ||
      a.rows() != m || a.cols() != z.k || coordinate.rows() != m ||
      coordinate.cols() != r || basis.rows() != m || basis.cols() != r) {
    Rcpp::stop("The fp16 build inputs are misaligned with the cache and basis.");
  }
  for (int f = 0; f < z.k; ++f) {
    if (slots[f] < 1 || slots[f] > cache->n) Rcpp::stop("A cache slot is out of range.");
    if (cache->state[slots[f] - 1] != 0) Rcpp::stop("A cache slot is already filled.");
  }
  const double per = gene_bytes(r);
  const double scratch = double(threads) * 8 *
    (8.0 * m * m + 6.0 * m * px + 4.0 * n + 4.0 * r * r) +
    double(threads) * 16.0 * r * r + 2.0 * 1024 * 1024 * 1024;
  const double free = available_bytes(budget);
  const double limit = std::isfinite(free) ?
    cache->bytes + std::max(0.0, 0.7 * free - scratch) :
    std::numeric_limits<double>::infinity();
  if (cache->bytes + per * z.k > limit) {
    Rcpp::stop("Insufficient memory: the fp16 gene cache would exceed the safe memory line "
               "(%.2f GiB needed beyond %.2f GiB resident, safe line %.2f GiB). "
               "Increase the job memory or reduce threads.",
               per * z.k / 1073741824.0, cache->bytes / 1073741824.0, limit / 1073741824.0);
  }
  const size_t L = packed_length(r);
  bool full = false;
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
#endif
  for (int f = 0; f < z.k; ++f) {
    const int slot = slots[f] - 1;
    bool admitted = false;
#ifdef _OPENMP
#pragma omp critical(mgcvst_fp16_admission)
#endif
    {
      const double now = mgcvst_memory::available_physical_memory();
      if (!full && cache->bytes + per <= limit &&
          (!std::isfinite(now) || now > 0.3 * scratch)) {
        cache->bytes += per;
        admitted = true;
      } else full = true;
    }
    if (!admitted) continue;
    try {
      if (!std::isfinite(tau[f]) || tau[f] <= 0) {
        throw std::runtime_error("The feature has invalid dispersion or smoothing parameters.");
      }
      if (!a.col(f).allFinite()) throw std::runtime_error("The saved score vector is non-finite.");
      Vec eta, mu, V;
      recover_working(z.A, z.X, z.model(f), eta, mu, V);
      ReconstructionUnit unit;
      curvature_unit(z.A, z.X, Q, g, V, tau[f], penalty.col(f).data(), unit);
      const Mat M = reduced_curvature(unit, basis, g);
      if (!M.allFinite()) throw std::runtime_error("The reduced curvature is non-finite.");
      std::vector<Half> packed(L);
      size_t at = 0;
      for (int j = 0; j < r; ++j) {
        for (int i = 0; i <= j; ++i) {
          const Half h(static_cast<float>(M(i, j)));
          if (!std::isfinite(static_cast<float>(h))) {
            throw std::runtime_error("The reduced curvature exceeds the fp16 range.");
          }
          packed[at++] = h;
        }
      }
      cache->a.col(slot) = coordinate.transpose() * a.col(f);
      cache->m[slot].swap(packed);
      cache->state[slot] = 1;
    } catch (const std::exception& x) {
      cache->error[slot] = x.what();
      cache->state[slot] = 2;
#ifdef _OPENMP
#pragma omp critical(mgcvst_fp16_admission)
#endif
      cache->bytes -= per;
    }
  }
  if (full) {
    Rcpp::stop("Insufficient memory: fp16 state admission stopped at the safe memory line "
               "(%.2f GiB resident). Increase the job memory or reduce threads.",
               cache->bytes / 1073741824.0);
  }
  Rcpp::IntegerVector state(z.k);
  for (int f = 0; f < z.k; ++f) state[f] = cache->state[slots[f] - 1];
  return Rcpp::List::create(Rcpp::Named("state") = state,
                            Rcpp::Named("bytes") = cache->bytes);
}

// Write built cache slots as one binary shard: header (r, signature, gene
// indices and ids, states, error messages), double projected scores, then the
// raw fp16 upper-triangle payloads.
// [[Rcpp::export]]
void mgcvst_fp16_write_cpp(SEXP pointer, const Rcpp::IntegerVector& slots,
                           const Rcpp::IntegerVector& feature_index,
                           const Rcpp::CharacterVector& feature_id,
                           std::string signature, std::string path) {
  const HalfCache* cache = cache_pointer(pointer);
  const int k = slots.size();
  if (feature_index.size() != k || feature_id.size() != k) {
    Rcpp::stop("Shard feature metadata is misaligned.");
  }
  for (int f = 0; f < k; ++f) {
    if (slots[f] < 1 || slots[f] > cache->n || cache->state[slots[f] - 1] == 0) {
      Rcpp::stop("Only built cache slots can be written.");
    }
  }
  std::ofstream out(path.c_str(), std::ios::binary | std::ios::trunc);
  if (!out) Rcpp::stop("Could not open the fp16 state shard for writing.");
  out.write(kMagic, 8);
  put<int32_t>(out, kVersion);
  put<int32_t>(out, cache->r);
  put<int32_t>(out, k);
  put_string(out, signature);
  for (int f = 0; f < k; ++f) {
    const int slot = slots[f] - 1;
    put<int32_t>(out, feature_index[f]);
    put_string(out, Rcpp::as<std::string>(feature_id[f]));
    put<int32_t>(out, cache->state[slot]);
    put_string(out, cache->error[slot]);
  }
  for (int f = 0; f < k; ++f) {
    const int slot = slots[f] - 1;
    if (cache->state[slot] == 1) {
      out.write(reinterpret_cast<const char*>(cache->a.col(slot).data()),
                sizeof(double) * cache->r);
    }
  }
  const size_t L = packed_length(cache->r);
  for (int f = 0; f < k; ++f) {
    const int slot = slots[f] - 1;
    if (cache->state[slot] == 1) {
      out.write(reinterpret_cast<const char*>(cache->m[slot].data()), sizeof(Half) * L);
    }
  }
  out.close();
  if (!out) Rcpp::stop("Could not write the fp16 state shard.");
}

// Load one shard directly into empty cache slots, validating its header against
// the expected signature, dimension and genes.
// [[Rcpp::export]]
double mgcvst_fp16_read_cpp(SEXP pointer, const Rcpp::IntegerVector& slots,
                            const Rcpp::IntegerVector& feature_index,
                            const Rcpp::CharacterVector& feature_id,
                            std::string signature, std::string path,
                            double budget) {
  HalfCache* cache = cache_pointer(pointer);
  const int k = slots.size();
  const int r = cache->r;
  const size_t L = packed_length(r);
  for (int f = 0; f < k; ++f) {
    if (slots[f] < 1 || slots[f] > cache->n || cache->state[slots[f] - 1] != 0) {
      Rcpp::stop("Resumed cache slots must be empty and in range.");
    }
  }
  std::ifstream in(path.c_str(), std::ios::binary);
  if (!in) Rcpp::stop("Could not open the fp16 state shard %s.", path);
  std::vector<int> state(k);
  std::vector<std::string> error(k);
  std::vector<Eigen::VectorXd> a(k);
  std::vector<std::vector<Half> > m(k);
  int ready = 0;
  try {
    char magic[8];
    in.read(magic, 8);
    if (!in || std::memcmp(magic, kMagic, 8) != 0) throw std::runtime_error("bad magic");
    if (get<int32_t>(in) != kVersion || get<int32_t>(in) != r ||
        get<int32_t>(in) != k || get_string(in) != signature) {
      throw std::runtime_error("incompatible header");
    }
    for (int f = 0; f < k; ++f) {
      if (get<int32_t>(in) != feature_index[f] ||
          get_string(in) != Rcpp::as<std::string>(feature_id[f])) {
        throw std::runtime_error("unexpected genes");
      }
      state[f] = get<int32_t>(in);
      error[f] = get_string(in);
      if (state[f] != 1 && state[f] != 2) throw std::runtime_error("invalid state");
      ready += state[f] == 1;
    }
    const double free = available_bytes(budget);
    if (std::isfinite(free) && gene_bytes(r) * ready > std::max(0.0, 0.7 * free)) {
      Rcpp::stop("Insufficient memory: resumed fp16 states exceed the safe memory line.");
    }
    for (int f = 0; f < k; ++f) {
      if (state[f] != 1) continue;
      a[f].resize(r);
      in.read(reinterpret_cast<char*>(a[f].data()), sizeof(double) * r);
    }
    for (int f = 0; f < k; ++f) {
      if (state[f] != 1) continue;
      m[f].resize(L);
      in.read(reinterpret_cast<char*>(m[f].data()), sizeof(Half) * L);
    }
    if (!in) throw std::runtime_error("truncated");
    in.peek();
    if (!in.eof()) throw std::runtime_error("trailing bytes");
  } catch (const std::runtime_error& x) {
    Rcpp::stop("The fp16 state shard %s is damaged or incompatible (%s).", path, x.what());
  }
  for (int f = 0; f < k; ++f) {
    const int slot = slots[f] - 1;
    cache->state[slot] = state[f];
    cache->error[slot] = error[f];
    if (state[f] == 1) {
      cache->a.col(slot) = a[f];
      cache->m[slot].swap(m[f]);
      cache->bytes += gene_bytes(r);
    }
  }
  return cache->bytes;
}

// Exact Liu pairs from the fp16 cache. Pairs are cache slots (1-based,
// left < right): every right slot after each left slot in
// [left_first, left_last] or, when supplied, the explicit left/right vectors
// sorted by left. Work is scheduled dynamically over left-gene segments; each
// worker expands its left matrix once to fp32, forms P = L R and P^2 in fp32 and
// accumulates tr(P^s), s = 1..4, the score and the Liu tail in double.
// Output i/j are the global feature indices used[slot].
// [[Rcpp::export]]
Rcpp::List mgcvst_fp16_pairs_cpp(SEXP pointer, const Rcpp::IntegerVector& used,
                                 int left_first, int left_last,
                                 Rcpp::Nullable<Rcpp::IntegerVector> left = R_NilValue,
                                 Rcpp::Nullable<Rcpp::IntegerVector> right = R_NilValue,
                                 int threads = 1) {
  check_threads(threads);
  const HalfCache* cache = cache_pointer(pointer);
  const int n = cache->n, r = cache->r;
  if (used.size() != n) Rcpp::stop("used must give one feature index per cache slot.");
  const bool explicit_pairs = left.isNotNull();
  if (explicit_pairs != right.isNotNull()) Rcpp::stop("Supply both left and right slots.");
  Rcpp::IntegerVector li, lj;
  struct Task { int left; R_xlen_t first, last, base; };
  std::vector<Task> tasks;
  R_xlen_t N = 0;
  const R_xlen_t grain = 32;
  if (explicit_pairs) {
    li = Rcpp::IntegerVector(left.get());
    lj = Rcpp::IntegerVector(right.get());
    if (li.size() != lj.size()) Rcpp::stop("left and right must have equal length.");
    N = li.size();
    for (R_xlen_t k = 0; k < N; ++k) {
      if (li[k] == NA_INTEGER || lj[k] == NA_INTEGER || li[k] < 1 || lj[k] > n ||
          li[k] >= lj[k]) {
        Rcpp::stop("Every pair must satisfy 1 <= left < right <= number of cached genes.");
      }
      if (k && li[k] < li[k - 1]) Rcpp::stop("Explicit pairs must be sorted by left slot.");
    }
    R_xlen_t k = 0;
    while (k < N) {
      R_xlen_t end = k;
      while (end < N && li[end] == li[k]) ++end;
      for (R_xlen_t s = k; s < end; s += grain) {
        tasks.push_back(Task{li[k] - 1, s, std::min(end, s + grain), -1});
      }
      k = end;
    }
  } else {
    if (left_first < 1 || left_last > n || left_first > left_last) {
      Rcpp::stop("left_first and left_last must define a range of cached genes.");
    }
    for (int l = left_first - 1; l < left_last; ++l) {
      const R_xlen_t count = n - 1 - l;
      for (R_xlen_t s = 0; s < count; s += grain) {
        tasks.push_back(Task{l, N + s, N + std::min(count, s + grain), N});
      }
      N += count;
    }
  }
  Rcpp::IntegerVector I(N), J(N);
  Rcpp::NumericVector S(N), P(N);
  int* ip = I.begin();
  int* jp = J.begin();
  double* sp = S.begin();
  double* pp = P.begin();
  const int* ljp = explicit_pairs ? lj.begin() : NULL;
  const int* up = used.begin();
  const long ntask = tasks.size();
#ifdef _OPENMP
#pragma omp parallel num_threads(threads)
#endif
  {
    FMat Lm(r, r), Rm(r, r), P1(r, r), P2(r, r);
    int current = -1;
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
    for (long t = 0; t < ntask; ++t) {
      const Task& task = tasks[t];
      const int l = task.left;
      const bool left_ready = cache->state[l] == 1;
      if (left_ready && current != l) {
        unpack(cache->m[l], Lm);
        current = l;
      }
      for (R_xlen_t k = task.first; k < task.last; ++k) {
        const int j = explicit_pairs ? ljp[k] - 1 : int(l + 1 + (k - task.base));
        ip[k] = up[l];
        jp[k] = up[j];
        sp[k] = pp[k] = NA_REAL;
        if (!left_ready || cache->state[j] != 1) continue;
        unpack(cache->m[j], Rm);
        P1.noalias() = Lm * Rm;
        P2.noalias() = P1 * P1;
        double v[4] = {0, 0, 0, 0};
        for (int i = 0; i < r; ++i) v[0] += double(P1(i, i));
        for (int c = 0; c < r; ++c) {
          for (int i = 0; i < r; ++i) {
            const double p = P1(i, c), pt = P1(c, i), h = P2(i, c), ht = P2(c, i);
            v[1] += p * pt;
            v[2] += h * pt;
            v[3] += h * ht;
          }
        }
        const double score = cache->a.col(l).dot(cache->a.col(j));
        sp[k] = score;
        const bool good = std::isfinite(score) && std::isfinite(v[0]) &&
          std::isfinite(v[1]) && std::isfinite(v[2]) && std::isfinite(v[3]) &&
          v[0] > 1e-10 && v[1] > 0 && v[2] > 0 && v[3] > 0;
        if (!good) continue;
        double lp[3];
        mgcvst_liu::liu_log_p(score, v[0], v[1], v[2], v[3], lp);
        if (std::isfinite(lp[0])) pp[k] = -lp[0] / M_LN10;
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("i") = I, Rcpp::Named("j") = J,
                            Rcpp::Named("score") = S, Rcpp::Named("mlog10p") = P);
}
