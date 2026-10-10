#define EIGEN_DONT_PARALLELIZE
#include "inla_sparse.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include <cstring>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

// Shared compact working-state reconstruction for sparse INLA fits: saved
// coefficients -> eta, mu and working variance, and the double
// reconstruction units that PCAlearning materializes from. These two entry
// points are what the PCAlearning pair path (R/pair-pcalearning.R) and the
// INLA scale computation consume.

using namespace mgcvst_sparse;

namespace {

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
