#define EIGEN_DONT_PARALLELIZE
#include <RcppEigen.h>
#include "pca_reduce.h"
#include "spa_pair.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

// PCAlearning saddlepoint path.
//
// Symmetric q x q matrices are stored in weighted-vech form: the upper
// triangle in column-major order (entry (i, j), i <= j, at j (j + 1) / 2 + i)
// with off-diagonal weight sqrt(2), so that inner products of weighted-vech
// vectors are Frobenius inner products. Training matrices H_j are kept in this
// form as float32 raw vectors; the basis B is returned in this form (L x r,
// L = q (q + 1) / 2) and its columns B_k are the symmetric basis matrices.
//
// A pair needs the first two trace moments t_s = tr((H_i H_j)^s), s = 1, 2. With
// X = sum_k x_k B_k and Y = sum_k y_k B_k:
//   s = 1: tr(XY) = x'y, since B is orthonormal.
//   s = 2: tr(XYXY) = sum tr(Y_ab Y_cd), Y_ab = B_a B_b, and tr(Y_ab Y_cd) =
//          <Y_ba, Y_cd>: one Gram matrix of the r^2 products. It is written as
//          kappa_2(x)' Tsym2 kappa_2(y) with the degree-2 monomials kappa_2
//          (nondecreasing index pairs in lexicographic order), which are formed
//          once for all genes and reused by every block.
// The k leading singular values of the pair spectrum come from the shared
// basis V of the training genes: R_g = chol(V' H_g V / scale_g) and
// svd(R_i R_j'). The remainder is a Satterthwaite node matching t_1 and t_2.

namespace {

using MatD = Eigen::MatrixXd;
using MatF = Eigen::MatrixXf;
using VecD = Eigen::VectorXd;
using Index = Eigen::Index;
using mgcvst_pca::unpack_symmetric;

double now() {
  return std::chrono::duration<double>(
    std::chrono::steady_clock::now().time_since_epoch()).count();
}

// Multisets of size s over {0, ..., r - 1}: nondecreasing tuples in
// lexicographic order (first coordinate slowest).
struct Multisets {
  int r = 0, s = 0, d = 0;
  std::vector<int> tuple;  // d x s, row-major
  std::vector<int> id;     // lexicographic index of a sorted tuple -> multiset id
  int at(int t, int k) const { return tuple[(size_t)t * s + k]; }
  int find(int* x) const {
    std::sort(x, x + s);
    long z = 0;
    for (int k = 0; k < s; ++k) z = z * r + x[k];
    return id[z];
  }
};

Multisets make_multisets(int r, int s) {
  Multisets m;
  m.r = r;
  m.s = s;
  long N = 1;
  for (int k = 0; k < s; ++k) N *= r;
  m.id.assign(N, -1);
  std::vector<int> x(s);
  for (long t = 0; t < N; ++t) {
    long z = t;
    for (int k = s - 1; k >= 0; --k) {
      x[k] = z % r;
      z /= r;
    }
    bool ok = true;
    for (int k = 1; k < s; ++k) if (x[k] < x[k - 1]) ok = false;
    if (!ok) continue;
    m.id[t] = m.d++;
    for (int k = 0; k < s; ++k) m.tuple.push_back(x[k]);
  }
  return m;
}

Rcpp::IntegerMatrix multiset_matrix(const Multisets& m) {
  Rcpp::IntegerMatrix out(m.d, m.s);
  for (int t = 0; t < m.d; ++t) for (int k = 0; k < m.s; ++k) out(t, k) = m.at(t, k) + 1;
  return out;
}

// Degree-2 monomials kappa_2(c_g) for all genes g (d2 x n).
MatD monomials2(const Multisets& m, const MatD& C) {
  const Index n = C.rows();
  MatD K(m.d, n);
  for (Index g = 0; g < n; ++g) {
    for (int t = 0; t < m.d; ++t) K(t, g) = C(g, m.at(t, 0)) * C(g, m.at(t, 1));
  }
  return K;
}

// G(upper tiles) += P' P for P (len x n) float; tiles run in parallel with
// float products and double accumulation.
void gram_accumulate(const float* data, Index len, Index n, MatD& G,
                     int threads, int tile) {
  Eigen::Map<const MatF> P(data, len, n);
  std::vector<std::pair<Index, Index> > tasks;
  for (Index u = 0; u < n; u += tile) for (Index v = u; v < n; v += tile) tasks.push_back({u, v});
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
#endif
  for (long k = 0; k < (long)tasks.size(); ++k) {
    const Index u = tasks[k].first, v = tasks[k].second;
    const Index nu = std::min<Index>(tile, n - u), nv = std::min<Index>(tile, n - v);
    MatF tmp = P.middleCols(u, nu).transpose() * P.middleCols(v, nv);
    G.block(u, v, nu, nv) += tmp.cast<double>();
  }
}

void symmetrize_upper(MatD& G) {
  for (Index j = 0; j < G.cols(); ++j) for (Index i = j + 1; i < G.rows(); ++i) G(i, j) = G(j, i);
}

const float* raw_float(SEXP x) { return reinterpret_cast<const float*>(RAW(x)); }

// Work space of one thread of the pair kernel.
struct PairWork {
  MatD Mk;
  Eigen::BDCSVD<MatD> svd;
  mgcvst_spa::Scratch scratch;
  long above = 0;  // pairs whose remainder node lies above the largest leading value
  explicit PairWork(int k) : Mk(k, k), svd(k, k, 0) { scratch.reserve(k + 2); }
};

// One pair: score U, the projected trace moments t1, t2 (unnormalized), the
// scales and the factors R_i, R_j (k x k, column-major). Writes the eight
// output columns U, t1, t2, three log p-values, remainder kind and status.
inline void spa_pair_row(double* out, Index rows, Index row, double U, double t1,
                         double t2, double scale_i, double scale_j,
                         const double* Ri, const double* Rj, int k, PairWork& w) {
  out[row] = U;
  out[row + rows] = t1;
  out[row + 2 * rows] = t2;
  for (int z = 3; z < 6; ++z) out[row + z * rows] = NA_REAL;
  out[row + 6 * rows] = 0;
  out[row + 7 * rows] = 1;
  if (!std::isfinite(scale_i) || scale_i <= 0 || !std::isfinite(scale_j) ||
      scale_j <= 0 || !std::isfinite(U) || !std::isfinite(t1) || !std::isfinite(t2)) {
    return;
  }
  Eigen::Map<const MatD> Mi(Ri, k, k), Mj(Rj, k, k);
  w.Mk.noalias() = Mi * Mj.transpose();
  if (!w.Mk.allFinite()) return;
  w.svd.compute(w.Mk);
  const VecD& sv = w.svd.singularValues();
  long double lead[4];
  mgcvst_spa::leading_sums(sv.data(), k, lead);
  const double units2 = scale_i * scale_j;
  const double tn[4] = {t1 / units2, t2 / (units2 * units2), 0.0, 0.0};
  const mgcvst_spa::Remainder rem = mgcvst_spa::make_remainder(tn, lead, 2);
  double lp[3] = {NA_REAL, NA_REAL, NA_REAL};
  const double x = U / (std::sqrt(scale_i) * std::sqrt(scale_j));
  const int status = mgcvst_spa::spa_pair(x, sv.data(), k, rem, w.scratch, lp);
  out[row + 7 * rows] = status;
  if (status == 0) {
    out[row + 3 * rows] = lp[0];
    out[row + 4 * rows] = lp[1];
    out[row + 5 * rows] = lp[2];
    out[row + 6 * rows] = rem.kind;
    if (mgcvst_spa::node_above_leading(rem, sv[0])) ++w.above;
  }
}

Rcpp::CharacterVector pair_columns() {
  return Rcpp::CharacterVector::create("U", "t1", "t2", "logp_two_sided",
    "logp_positive", "logp_negative", "remainder_kind", "status");
}

}  // namespace

// Gram matrix G_jk = tau_j tau_k <H_j, H_k>_F of packed float32 training
// matrices. The rows are cut into chunks of fixed size; each chunk gives its
// own partial Gram matrix, and the partial matrices are added in chunk order,
// so the result is bitwise the same for any number of threads. The partial
// matrices need (number of chunks) x n^2 doubles.
// [[Rcpp::export]]
Rcpp::NumericMatrix mgcvst_pca_gram_cpp(const Rcpp::List& packed,
                                        const Rcpp::NumericVector& tau,
                                        int threads = 1, int chunk = 8192) {
  const int n = packed.size();
  const Index L = XLENGTH(packed[0]) / 4;
  std::vector<const float*> ptr(n);
  for (int j = 0; j < n; ++j) ptr[j] = raw_float(packed[j]);
  const Index nb = (L + chunk - 1) / chunk;
  std::vector<MatD> parts(nb);
#ifdef _OPENMP
#pragma omp parallel num_threads(threads)
#endif
  {
    MatD V(chunk, n);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
    for (Index b = 0; b < nb; ++b) {
      const Index r0 = b * chunk, m = std::min<Index>(chunk, L - r0);
      for (int j = 0; j < n; ++j) {
        const float* x = ptr[j] + r0;
        for (Index i = 0; i < m; ++i) V(i, j) = x[i];
      }
      parts[b] = MatD::Zero(n, n);
      parts[b].selfadjointView<Eigen::Lower>().rankUpdate(V.topRows(m).transpose());
    }
  }
  MatD G = MatD::Zero(n, n);
  for (Index b = 0; b < nb; ++b) G += parts[b];
  for (int j = 0; j < n; ++j) for (int i = j; i < n; ++i) G(j, i) = G(i, j) = G(i, j) * tau[i] * tau[j];
  return Rcpp::wrap(G);
}

// Basis in weighted-vech form: B = V R, V = [tau_j vech_w(H_j)], R = n x r
// (eigenvectors scaled by 1 / sqrt(eigenvalues)).
// [[Rcpp::export]]
Rcpp::NumericMatrix mgcvst_pca_basis_cpp(const Rcpp::List& packed,
                                         const Rcpp::NumericVector& tau,
                                         const Eigen::Map<Eigen::MatrixXd> R,
                                         int threads = 1, int chunk = 8192) {
  const int n = packed.size(), r = R.cols();
  const Index L = XLENGTH(packed[0]) / 4;
  std::vector<const float*> ptr(n);
  for (int j = 0; j < n; ++j) ptr[j] = raw_float(packed[j]);
  MatD TR = Eigen::Map<const Eigen::VectorXd>(tau.begin(), n).asDiagonal() * R;
  Rcpp::NumericMatrix out(L, r);
  Eigen::Map<MatD> B(out.begin(), L, r);
  const Index nb = (L + chunk - 1) / chunk;
#ifdef _OPENMP
#pragma omp parallel num_threads(threads)
#endif
  {
    MatD V(chunk, n);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
    for (Index b = 0; b < nb; ++b) {
      const Index r0 = b * chunk, m = std::min<Index>(chunk, L - r0);
      for (int j = 0; j < n; ++j) {
        const float* x = ptr[j] + r0;
        for (Index i = 0; i < m; ++i) V(i, j) = x[i];
      }
      B.middleRows(r0, m).noalias() = V.topRows(m) * TR;
    }
  }
  return out;
}

// Trace table Tsym2 of a symmetric orthonormal basis B given as q^2 x r (full
// columns, symmetrized) or q (q + 1) / 2 x r (weighted vech). The products
// Y_ab = B_a B_b are stored for a <= b only (Y_ba is its transpose).
// [[Rcpp::export]]
Rcpp::List mgcvst_pca_tables_cpp(const Rcpp::NumericMatrix& B, int q,
                                 int threads = 1, int block = 32, int tile = 192) {
  const double t_start = now();
  const int r = B.ncol();
  const Index qq = (Index)q * q;
  const bool full = B.nrow() == qq;
  const Multisets m2 = make_multisets(r, 2);
  const int d2 = m2.d;
  std::vector<double> timing;
  std::vector<std::string> tnames;
  std::vector<double> memory;
  std::vector<std::string> mnames;
  auto mark = [&](const char* name, double t0) { timing.push_back(now() - t0); tnames.push_back(name); };

  // B_k (double) and Y_ab = B_a B_b (double GEMM, stored float).
  double t0 = now();
  std::vector<MatD> Bd(r);
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
#endif
  for (int k = 0; k < r; ++k) Bd[k] = unpack_symmetric(&B(0, k), q, full);
  std::vector<MatF> Yf(d2);
  std::vector<int> pid((size_t)r * r);
  for (int t = 0; t < d2; ++t) {
    const int a = m2.at(t, 0), b = m2.at(t, 1);
    pid[(size_t)a * r + b] = pid[(size_t)b * r + a] = t;
  }
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
#endif
  for (int t = 0; t < d2; ++t) {
    const int a = m2.at(t, 0), b = m2.at(t, 1);
    MatD y = Bd[a] * Bd[b];
    Yf[t] = y.cast<float>();
  }
  Bd.clear();
  mark("products", t0);
  memory.push_back(8.0 * qq * r + 4.0 * qq * d2); mnames.push_back("B_Y");

  // Level 2: Gram of the r^2 products over column blocks of full entries.
  t0 = now();
  MatD Gy = MatD::Zero(r * r, r * r);
  {
    MatF P((Index)q * block, r * r);
    memory.push_back(4.0 * P.size()); mnames.push_back("panel");
    for (int c0 = 0; c0 < q; c0 += block) {
      const int w = std::min(block, q - c0);
      const Index len = (Index)q * w;
      for (int a = 0; a < r; ++a) {
        for (int b = 0; b < r; ++b) {
          float* dst = P.data() + ((Index)a * r + b) * len;
          const float* src = Yf[pid[(size_t)a * r + b]].data();
          if (a <= b) {
            std::memcpy(dst, src + (Index)c0 * q, sizeof(float) * len);
          } else {
            // Y_ba = Y_ab': columns c0, ..., c0 + w - 1 of the transpose.
            for (Index i = 0; i < q; ++i) {
              const float* row = src + i * q + c0;
              for (int c = 0; c < w; ++c) dst[(Index)c * q + i] = row[c];
            }
          }
        }
      }
      gram_accumulate(P.data(), len, r * r, Gy, threads, tile);
    }
  }
  symmetrize_upper(Gy);
  MatD T2 = MatD::Zero(d2, d2);
  for (int a1 = 0; a1 < r; ++a1) for (int b1 = 0; b1 < r; ++b1)
    for (int a2 = 0; a2 < r; ++a2) for (int b2 = 0; b2 < r; ++b2) {
      int x[2] = {a1, a2}, y[2] = {b1, b2};
      T2(m2.find(x), m2.find(y)) += Gy(b1 * r + a1, a2 * r + b2);
    }
  mark("level2", t0);
  mark("total", t_start);
  memory.push_back(8.0 * (double)r * r * r * r); mnames.push_back("Gram2");

  Rcpp::NumericVector tv(timing.begin(), timing.end()), mv(memory.begin(), memory.end());
  tv.names() = tnames;
  mv.names() = mnames;
  return Rcpp::List::create(
    Rcpp::Named("Tsym2") = T2, Rcpp::Named("ms2") = multiset_matrix(m2),
    Rcpp::Named("r") = r, Rcpp::Named("q") = q,
    Rcpp::Named("timing") = tv, Rcpp::Named("memory") = mv);
}

// Degree-2 monomials of the basis coefficients of all genes (d2 x n), formed
// once and reused by every pair block.
// [[Rcpp::export]]
Rcpp::NumericMatrix mgcvst_pca_monomials_cpp(const Eigen::Map<Eigen::MatrixXd> C) {
  const Multisets m2 = make_multisets(C.cols(), 2);
  return Rcpp::wrap(monomials2(m2, C));
}

// Sum of the normalized training matrices, sum_g sym(H_g) / max|H_g|, from
// packed float32 weighted-vech vectors, added in gene order (so the result does
// not depend on the thread count). Its leading eigenvectors are the shared
// basis V of the PCAlearning route.
// [[Rcpp::export]]
Rcpp::NumericMatrix mgcvst_pca_packed_sum_cpp(const Rcpp::List& packed, int q) {
  MatD S = MatD::Zero(q, q);
  int used = 0;
  for (int g = 0; g < packed.size(); ++g) {
    if (packed[g] == R_NilValue) continue;
    const MatD M = mgcvst_pca::unpack_float_vech(raw_float(packed[g]), q);
    const double scale = mgcvst_pca::matrix_scale(M);
    if (!(std::isfinite(scale) && scale > 0)) continue;
    S += M / scale;
    ++used;
  }
  Rcpp::NumericMatrix out = Rcpp::wrap(S);
  out.attr("used") = used;
  return out;
}

// Scale max|H| and the factor R = chol(V' H / scale V) of packed training
// matrices. Returns R (k^2 x n), the scales and an error message per gene.
// [[Rcpp::export]]
Rcpp::List mgcvst_pca_packed_project_cpp(const Rcpp::List& packed,
                                         const Eigen::Map<Eigen::MatrixXd> V,
                                         int threads = 1) {
  const int n = packed.size(), q = V.rows(), k = V.cols();
  Rcpp::NumericMatrix Rout((Index)k * k, n);
  Rcpp::NumericVector scale(n, NA_REAL);
  std::vector<const float*> ptr(n);
  for (int g = 0; g < n; ++g) ptr[g] = packed[g] == R_NilValue ? nullptr : raw_float(packed[g]);
  std::vector<std::string> error(n);
  Eigen::Map<MatD> Rall(Rout.begin(), (Index)k * k, n);
  double* sp = scale.begin();
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
#endif
  for (int g = 0; g < n; ++g) {
    Rall.col(g).setConstant(NA_REAL);
    if (!ptr[g]) { error[g] = "No packed matrix."; continue; }
    const MatD M = mgcvst_pca::unpack_float_vech(ptr[g], q);
    const double s = mgcvst_pca::matrix_scale(M);
    sp[g] = s;
    MatD R;
    if (!(std::isfinite(s) && s > 0)) { error[g] = "The curvature matrix is zero or not finite."; continue; }
    if (!mgcvst_pca::project_factor(M, s, V, R)) {
      error[g] = "The shared-basis compression is not positive definite.";
      continue;
    }
    Rall.col(g) = Eigen::Map<const VecD>(R.data(), (Index)k * k);
  }
  Rcpp::CharacterVector error_out(n);
  for (int g = 0; g < n; ++g) error_out[g] = error[g].empty() ? NA_STRING : Rcpp::String(error[g]);
  return Rcpp::List::create(Rcpp::Named("R") = Rout, Rcpp::Named("scale") = scale,
                            Rcpp::Named("error") = error_out);
}

// PCAlearning materialization of dense score states: for every state (a, M)
// the score coordinates a, the basis coefficients c = <B_k, H>_F, ||H||_F^2,
// the scale max|H| and, with a shared basis V, the factor R; features with
// pack = TRUE also keep the float32 weighted-vech copy of H. The same outputs
// as the sparse INLA materialization, for mgcv states.
// [[Rcpp::export]]
Rcpp::List mgcvst_pca_dense_cpp(const Rcpp::List& H, const Eigen::Map<Eigen::MatrixXd> a,
                                Rcpp::Nullable<Rcpp::NumericMatrix> pca_basis,
                                const Rcpp::LogicalVector& pack,
                                Rcpp::Nullable<Rcpp::NumericMatrix> V,
                                int threads = 1) {
  const int n = H.size(), q = a.rows();
  const Index L = (Index)q * (q + 1) / 2;
  if (a.cols() != n || pack.size() != n) Rcpp::stop("H, a and pack must align.");
  Rcpp::NumericMatrix PBr = pca_basis.isNotNull() ? Rcpp::NumericMatrix(pca_basis.get())
                                                  : Rcpp::NumericMatrix(L, 0);
  const Eigen::Map<const MatD> PB(PBr.begin(), L, PBr.ncol());
  const int r = PB.cols();
  Rcpp::NumericMatrix Vr = V.isNotNull() ? Rcpp::NumericMatrix(V.get()) : Rcpp::NumericMatrix(q, 0);
  if (Vr.nrow() != q) Rcpp::stop("V must have q rows.");
  const Eigen::Map<const MatD> Vm(Vr.begin(), q, Vr.ncol());
  const int k = Vm.cols();
  std::vector<const double*> Hptr(n);
  for (int g = 0; g < n; ++g) {
    Rcpp::NumericMatrix current(H[g]);
    if (current.nrow() != q || current.ncol() != q) Rcpp::stop("Every state matrix must be q x q.");
    Hptr[g] = current.begin();
  }
  MatD C(n, r), Aout = a;
  VecD fro2(n);
  Rcpp::NumericMatrix Rout((Index)k * k, n);
  Eigen::Map<MatD> Rall(Rout.begin(), (Index)k * k, n);
  Rcpp::NumericVector scale(n, NA_REAL);
  double* sp = scale.begin();
  std::vector<std::vector<float> > packed(n);
  std::vector<std::string> error(n);
#ifdef _OPENMP
#pragma omp parallel num_threads(threads)
#endif
  {
    VecD h(L);
    VecD ccol(r), rcol((Index)k * k);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
    for (int g = 0; g < n; ++g) {
      Eigen::Map<const MatD> Mraw(Hptr[g], q, q);
      const MatD M = 0.5 * (Mraw + Mraw.transpose());
      double s = NA_REAL, f2 = NA_REAL;
      error[g] = mgcvst_pca::reduce_feature(M, q, PB, Vm, pack[g] != 0, h, s, f2,
                                            ccol, packed[g], rcol);
      sp[g] = s;
      fro2[g] = f2;
      C.row(g) = ccol.transpose();
      if (k) Rall.col(g) = rcol;
      if (!error[g].empty()) Aout.col(g).setConstant(NA_REAL);
    }
  }
  Rcpp::List packed_out = mgcvst_pca::packed_to_raw(packed, L);
  Rcpp::CharacterVector error_out(n);
  int failed = 0;
  for (int g = 0; g < n; ++g) {
    error_out[g] = error[g].empty() ? NA_STRING : Rcpp::String(error[g]);
    if (!error[g].empty()) ++failed;
  }
  return Rcpp::List::create(
    Rcpp::Named("a") = Rcpp::wrap(Aout), Rcpp::Named("C") = Rcpp::wrap(C),
    Rcpp::Named("fro2") = Rcpp::wrap(fro2), Rcpp::Named("packed") = packed_out,
    Rcpp::Named("R") = Rout, Rcpp::Named("scale") = scale,
    Rcpp::Named("error") = error_out, Rcpp::Named("failed") = failed,
    Rcpp::Named("width") = q);
}

// Saddlepoint pair values for listed pairs (1-based i, j): U = a_i' a_j,
// t1 = c_i' c_j, t2 = kappa_2(c_i)' Tsym2 kappa_2(c_j) with the monomials K2
// of all genes, and the k leading singular values of R_i R_j'. Columns U, t1,
// t2, logp_two_sided, logp_positive, logp_negative, remainder_kind, status.
// [[Rcpp::export]]
Rcpp::NumericMatrix mgcvst_pca_spa_pairs_cpp(const Eigen::Map<Eigen::MatrixXd> A,
                                             const Eigen::Map<Eigen::MatrixXd> C,
                                             const Eigen::Map<Eigen::MatrixXd> K2,
                                             const Eigen::Map<Eigen::MatrixXd> T2,
                                             const Eigen::Map<Eigen::MatrixXd> R,
                                             const Rcpp::NumericVector& scale,
                                             const Rcpp::IntegerVector& i,
                                             const Rcpp::IntegerVector& j,
                                             int threads = 1) {
  const int n = C.rows();
  const R_xlen_t np = i.size();
  const int k = (int)std::lround(std::sqrt((double)R.rows()));
  if ((Index)k * k != R.rows() || R.cols() != n || scale.size() != n ||
      K2.cols() != n || T2.rows() != K2.rows() || A.cols() != n || j.size() != np) {
    Rcpp::stop("The PCAlearning inputs are not aligned.");
  }
  std::vector<char> used(n, 0);
  for (R_xlen_t t = 0; t < np; ++t) {
    if (i[t] < 1 || i[t] > n || j[t] < 1 || j[t] > n) Rcpp::stop("Pair indices are out of range.");
    used[i[t] - 1] = 1;
  }
  std::vector<int> slot(n, -1);
  std::vector<int> ucols;
  for (int g = 0; g < n; ++g) if (used[g]) { slot[g] = ucols.size(); ucols.push_back(g); }
  MatD S2(K2.rows(), ucols.size());
  {
    MatD Ku(K2.rows(), ucols.size());
    for (size_t g = 0; g < ucols.size(); ++g) Ku.col(g) = K2.col(ucols[g]);
    const Index nu = ucols.size(), cb = 64;
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
#endif
    for (Index b = 0; b < (nu + cb - 1) / cb; ++b) {
      const Index c0 = b * cb, w = std::min(cb, nu - c0);
      S2.middleCols(c0, w).noalias() = T2 * Ku.middleCols(c0, w);
    }
  }
  Rcpp::NumericMatrix out(np, 8);
  double* o = out.begin();
  const double* sc = scale.begin();
  long nodes_above = 0;
#ifdef _OPENMP
#pragma omp parallel num_threads(threads)
#endif
  {
    PairWork work(k);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 64)
#endif
    for (R_xlen_t t = 0; t < np; ++t) {
      const int a = i[t] - 1, b = j[t] - 1;
      const double U = A.col(a).dot(A.col(b));
      const double t1 = C.row(a).dot(C.row(b));
      const double t2 = S2.col(slot[a]).dot(K2.col(b));
      spa_pair_row(o, np, t, U, t1, t2, sc[a], sc[b], &R(0, a), &R(0, b), k, work);
    }
#ifdef _OPENMP
#pragma omp atomic
#endif
    nodes_above += work.above;
  }
  Rcpp::colnames(out) = pair_columns();
  out.attr("nodes_above_leading") = (double)nodes_above;
  return out;
}

// All pairs (i, j), first <= i <= last, i < j <= n (1-based), ordered by i
// then j. Per block of genes i, t1 and t2 are GEMMs C_i C_j' and
// (Tsym2 K2_i)' K2_j; the saddlepoint of each pair follows.
// [[Rcpp::export]]
Rcpp::NumericMatrix mgcvst_pca_spa_block_cpp(const Eigen::Map<Eigen::MatrixXd> A,
                                             const Eigen::Map<Eigen::MatrixXd> C,
                                             const Eigen::Map<Eigen::MatrixXd> K2,
                                             const Eigen::Map<Eigen::MatrixXd> T2,
                                             const Eigen::Map<Eigen::MatrixXd> R,
                                             const Rcpp::NumericVector& scale,
                                             int first, int last,
                                             int threads = 1,
                                             int gene_block = 32, int pair_block = 256) {
  const int n = C.rows();
  const int k = (int)std::lround(std::sqrt((double)R.rows()));
  if ((Index)k * k != R.rows() || R.cols() != n || scale.size() != n ||
      K2.cols() != n || T2.rows() != K2.rows() || A.cols() != n) {
    Rcpp::stop("The PCAlearning inputs are not aligned.");
  }
  const int i0 = first - 1, ni = last - first + 1;
  std::vector<Index> offset(ni + 1, 0);
  for (int t = 0; t < ni; ++t) offset[t + 1] = offset[t] + (n - 1 - (i0 + t));
  const Index np = offset[ni];
  MatD S2(K2.rows(), ni);
  {
    const MatD Ki = K2.middleCols(i0, ni);
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
#endif
    for (int b = 0; b < (ni + 63) / 64; ++b) {
      const int c0 = b * 64, w = std::min(64, ni - c0);
      S2.middleCols(c0, w).noalias() = T2 * Ki.middleCols(c0, w);
    }
  }
  struct Task { int ib, nb, j0, nj; };
  std::vector<Task> tasks;
  for (int ib = 0; ib < ni; ib += gene_block) {
    const int nb = std::min(gene_block, ni - ib);
    for (int j0 = i0 + ib + 1; j0 < n; j0 += pair_block) tasks.push_back({ib, nb, j0, std::min(pair_block, n - j0)});
  }
  Rcpp::NumericMatrix out(np, 8);
  double* o = out.begin();
  const double* sc = scale.begin();
  long nodes_above = 0;
#ifdef _OPENMP
#pragma omp parallel num_threads(threads)
#endif
  {
    PairWork work(k);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
    for (long t = 0; t < (long)tasks.size(); ++t) {
      const Task& tk = tasks[t];
      const MatD Uab = A.middleCols(i0 + tk.ib, tk.nb).transpose() * A.middleCols(tk.j0, tk.nj);
      const MatD T1 = C.middleRows(i0 + tk.ib, tk.nb) * C.middleRows(tk.j0, tk.nj).transpose();
      const MatD T2ab = S2.middleCols(tk.ib, tk.nb).transpose() * K2.middleCols(tk.j0, tk.nj);
      for (int jj = 0; jj < tk.nj; ++jj) {
        const int jg = tk.j0 + jj;
        for (int ii = 0; ii < tk.nb; ++ii) {
          const int il = tk.ib + ii, ig = i0 + il;
          if (jg <= ig) continue;
          spa_pair_row(o, np, offset[il] + (jg - ig - 1), Uab(ii, jj), T1(ii, jj),
                       T2ab(ii, jj), sc[ig], sc[jg], &R(0, ig), &R(0, jg), k, work);
        }
      }
    }
#ifdef _OPENMP
#pragma omp atomic
#endif
    nodes_above += work.above;
  }
  Rcpp::colnames(out) = pair_columns();
  out.attr("nodes_above_leading") = (double)nodes_above;
  return out;
}
