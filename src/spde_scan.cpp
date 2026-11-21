#define USE_FC_LEN_T
#define EIGEN_DONT_PARALLELIZE
#include "sparse_score.h"
#include "spa_quadform.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include <R_ext/Lapack.h>
#include <string>
#include <vector>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

// spde.scan score kernel (R/spde-scan.R): the variance-component score test of
// an SPDE field against a GLM null, for a batch of features that share the
// null design X and the field A, Q.
//
// For feature f with null mean mu, overdispersion alpha (0 = Poisson) and
// working weights W = mu / (1 + alpha mu), the working residual is
// r = (y - mu) / mu and the null projection is
//   P0 = W - W X (X' W X)^{-1} X' W.
// With B B' = Q^{-1} (B = P^{-1} L^{-T} from the sparse Cholesky Q = P'LL'P),
//   a = B' A' P0 r,      T = |a|^2,      M = B' A' P0 A B,
// and T is calibrated against the positive spectrum of M by the saddlepoint
// (spa_quadform.h), feature by feature inside the parallel loop. P0 already removes
// the span of X, constant included, so the field is used as it is: no
// intercept projection, no constraint, Q stays sparse.

using namespace mgcvst_sparse;

namespace {

// Eigenvalues of the symmetric M above spectrum_tol times the largest one, in
// decreasing order, through LAPACK dsyevd.
Vec scan_spectrum(const Mat& value, double spectrum_tol) {
  Mat work_matrix = value;
  int n = static_cast<int>(work_matrix.rows());
  int info = 0, lwork = -1, liwork = -1, iwork_query = 0;
  double work_query = 0;
  Vec values(n);
  F77_CALL(dsyevd)("N", "L", &n, work_matrix.data(), &n, values.data(),
                   &work_query, &lwork, &iwork_query, &liwork, &info FCONE FCONE);
  if (info != 0) throw std::runtime_error("The score curvature eigendecomposition failed.");
  lwork = static_cast<int>(work_query);
  liwork = iwork_query;
  std::vector<double> work(std::max(1, lwork));
  std::vector<int> iwork(std::max(1, liwork));
  F77_CALL(dsyevd)("N", "L", &n, work_matrix.data(), &n, values.data(),
                   work.data(), &lwork, iwork.data(), &liwork, &info FCONE FCONE);
  if (info != 0 || !values.allFinite()) {
    throw std::runtime_error("The score curvature eigendecomposition failed.");
  }
  const double cutoff = spectrum_tol * values.maxCoeff();
  int kept = 0;
  for (int j = 0; j < values.size(); ++j) if (values[j] > 0 && values[j] > cutoff) ++kept;
  Vec out(kept);
  for (int j = values.size() - 1, k = 0; j >= 0; --j) {
    if (values[j] > 0 && values[j] > cutoff) out[k++] = values[j];
  }
  return out;
}

struct ScanResult {
  double statistic = NA_REAL;
  double p_value = NA_REAL;
  std::string error;
};

} // namespace

// Yt: n x features counts (column f = feature f); offset: n x 1 or n x
// features; beta: ncol(X) x features null coefficients; alpha: overdispersion
// per feature (0 = Poisson). A feature with a non-finite beta or alpha is
// reported as an error and skipped.
// [[Rcpp::export]]
Rcpp::List mgcvst_spde_score_cpp(
    const Eigen::MappedSparseMatrix<double>& A_map,
    const Eigen::MappedSparseMatrix<double>& Q_map,
    const Eigen::Map<Eigen::MatrixXd> X,
    const Eigen::MappedSparseMatrix<double>& Yt,
    const Eigen::Map<Eigen::MatrixXd> offset,
    const Eigen::Map<Eigen::MatrixXd> beta,
    const Eigen::Map<Eigen::VectorXd> alpha,
    int threads = 1,
    int block_size = 32,
    double spectrum_tol = 1e-12,
    double spa_centre_tol = 1e-6,
    double spa_root_tol = 1e-12,
    int spa_root_maxit = 100) {
  const SpMat A = A_map;
  const SpMat Q = Q_map;
  const int n = A.rows();
  const int m = A.cols();
  const int p = X.cols();
  const int features = Yt.cols();
  if (Q.rows() != m || Q.cols() != m) Rcpp::stop("Q must be square and aligned with A.");
  if (X.rows() != n || Yt.rows() != n) Rcpp::stop("X and Yt must have one row per row of A.");
  if (offset.rows() != n || (offset.cols() != 1 && offset.cols() != features)) {
    Rcpp::stop("offset must be n x 1 or n x features.");
  }
  if (beta.rows() != p || beta.cols() != features || alpha.size() != features) {
    Rcpp::stop("beta must be ncol(X) x features and alpha one value per feature.");
  }
  if (threads < 1) Rcpp::stop("threads must be positive.");
  if (block_size < 1) Rcpp::stop("block_size must be positive.");
  if (!(spectrum_tol >= 0)) Rcpp::stop("spectrum_tol must be non-negative.");
  mgcvst_spa::QuadformSpaControl spa;
  spa.centre_tol = spa_centre_tol;
  spa.root_tol = spa_root_tol;
  spa.root_maxit = spa_root_maxit;
#ifndef _OPENMP
  if (threads > 1) Rcpp::stop("spdeST was compiled without OpenMP; use threads = 1.");
#endif

  QFactor qfactor;
  qfactor.compute(Q);
  if (qfactor.info() != Eigen::Success) Rcpp::stop("The sparse Cholesky factorization of Q failed.");
  const Mat Xd = X;

  std::vector<ScanResult> result(features);
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
#endif
  for (int f = 0; f < features; ++f) {
    try {
      if (!beta.col(f).allFinite() || !std::isfinite(alpha[f]) || alpha[f] < 0) {
        throw std::runtime_error("The null fit has no finite estimate.");
      }
      const Vec off = offset.cols() == 1 ? Vec(offset.col(0)) : Vec(offset.col(f));
      const Vec mu = (Xd * beta.col(f) + off).array().exp();
      Vec y = Vec::Zero(n);
      for (Eigen::MappedSparseMatrix<double>::InnerIterator it(Yt, f); it; ++it) {
        y[it.row()] = it.value();
      }
      Vec W(n), Wr(n);
      for (int i = 0; i < n; ++i) {
        W[i] = mu[i] / (1.0 + alpha[f] * mu[i]);
        Wr[i] = W[i] * (y[i] - mu[i]) / mu[i];
      }
      SpMat WA = A;
      for (int col = 0; col < WA.outerSize(); ++col) {
        for (SpMat::InnerIterator it(WA, col); it; ++it) it.valueRef() *= W[it.row()];
      }
      SpMat K = A.transpose() * WA;
      K.makeCompressed();
      const Mat L = A.transpose() * (W.asDiagonal() * Xd);
      Mat J = Xd.transpose() * W.asDiagonal() * Xd;
      Eigen::LDLT<Mat> jfactor(0.5 * (J + J.transpose()));
      if (jfactor.info() != Eigen::Success || !jfactor.isPositive()) {
        throw std::runtime_error("The null information X'WX is not positive definite.");
      }
      const Mat Vp = jfactor.solve(Mat::Identity(p, p));
      const Vec h = A.transpose() * Wr - L * (Vp * (Xd.transpose() * Wr));
      const Vec a = apply_Bt(qfactor, h);
      result[f].statistic = a.squaredNorm();

      Mat M(m, m);
      for (int first = 0; first < m; first += block_size) {
        const int width = std::min(block_size, m - first);
        Mat Z = Mat::Zero(m, width);
        Z.block(first, 0, width, width).setIdentity();
        const Mat V = apply_B(qfactor, Z);
        Mat Yb = K * V;
        Yb.noalias() -= L * (Vp * (L.transpose() * V));
        M.middleCols(first, width) = apply_Bt(qfactor, Yb);
      }
      const Vec lambda = scan_spectrum(0.5 * (M + M.transpose()), spectrum_tol);
      if (!lambda.size() || !std::isfinite(result[f].statistic)) {
        throw std::runtime_error("The score has no finite positive mixture spectrum.");
      }
      const double pv = mgcvst_spa::quadform_upper(
        result[f].statistic, std::vector<double>(lambda.data(), lambda.data() + lambda.size()),
        spa);
      if (!std::isfinite(pv) || pv < 0 || pv > 1) {
        throw std::runtime_error("Invalid score p-value.");
      }
      result[f].p_value = pv;
    } catch (const std::exception& error) {
      result[f].error = error.what();
    }
  }

  Rcpp::NumericVector statistic(features), p_value(features);
  Rcpp::CharacterVector error(features);
  for (int f = 0; f < features; ++f) {
    statistic[f] = result[f].statistic;
    p_value[f] = result[f].error.empty() ? result[f].p_value : NA_REAL;
    error[f] = result[f].error.empty() ? NA_STRING : Rcpp::String(result[f].error);
  }
  return Rcpp::List::create(Rcpp::Named("statistic") = statistic,
                            Rcpp::Named("p_value") = p_value,
                            Rcpp::Named("error") = error);
}
