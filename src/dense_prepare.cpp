#define ARMA_WARN_LEVEL 0
#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <exception>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include "la_dispatch.h"

// [[Rcpp::depends(RcppArmadillo)]]

namespace {

struct DenseFeatureResult {
  arma::vec a;
  arma::mat H;
  std::string error;
};

bool finite_matrix(const arma::mat& x) {
  return x.is_finite();
}

// One feature through the SIMD kernel (la_kernel.cpp, score_state): no BLAS
// call, so the OpenMP loop below never meets a threaded BLAS.
DenseFeatureResult one_feature(const arma::mat& T0,
                               const double* variance,
                               const double* error,
                               double scale,
                               const arma::mat& X,
                               const arma::mat* Vp,
                               bool score_only) {
  DenseFeatureResult out;
  const arma::uword n = T0.n_rows;
  const arma::uword q = T0.n_cols;
  if (!std::isfinite(scale) || scale <= 0.0) {
    out.error = "scale must be one positive finite value per feature.";
    return out;
  }
  bool positive = true;
  for (arma::uword i = 0; i < n; ++i) {
    if (!std::isfinite(variance[i]) || !std::isfinite(error[i])) {
      out.error = "variance and error must be finite vectors with n observations.";
      return out;
    }
    if (variance[i] <= 0.0) positive = false;
  }
  if (!positive) {
    out.error = "variance must contain strictly positive values.";
    return out;
  }
  out.a.set_size(q);
  if (!score_only) out.H.set_size(q, q);
  const int status = mgcvst::la().score_state(
    T0.memptr(), static_cast<int>(n), static_cast<int>(q), variance, error, scale,
    X.n_cols ? X.memptr() : NULL, static_cast<int>(X.n_cols),
    Vp ? Vp->memptr() : NULL, out.a.memptr(), score_only ? NULL : out.H.memptr());
  if (status != 0) {
    out.error = status == 1 ? "Woodbury solve failed." :
      status == 2 ? "Nuisance projection produced non-finite values." :
      "Dense score result contains non-finite values.";
    out.a.reset();
    out.H.reset();
  }
  return out;
}

} // namespace

// [[Rcpp::export]]
Rcpp::List mgcvst_dense_score_batch_cpp(const arma::mat& T0,
                                        const arma::mat& variance,
                                        const arma::mat& error,
                                        const arma::vec& scale,
                                        const arma::mat& X,
                                        const Rcpp::List& nuisance,
                                        int threads = 1,
                                        bool score_only = false) {
  if (T0.n_rows == 0L || T0.n_cols == 0L) Rcpp::stop("T0 must be non-empty.");
  if (!finite_matrix(T0)) Rcpp::stop("T0 must be finite.");
  if (variance.n_rows != T0.n_rows || error.n_rows != T0.n_rows ||
      variance.n_cols != error.n_cols || scale.n_elem != variance.n_cols) {
    Rcpp::stop("T0, variance, error, and scale have incompatible dimensions.");
  }
  if (X.n_rows != T0.n_rows || !finite_matrix(X)) {
    Rcpp::stop("X must be finite and have the same number of rows as T0.");
  }
  if (threads < 1) Rcpp::stop("threads must be a positive integer.");
#ifndef _OPENMP
  if (threads > 1) Rcpp::stop("threads > 1 requires an OpenMP-enabled build.");
#endif

  const bool fixed_projection = nuisance.size() == 0L;
  std::vector<arma::mat> vp;
  std::vector<std::string> pre_errors(variance.n_cols);
  if (!fixed_projection) {
    if (static_cast<arma::uword>(nuisance.size()) != variance.n_cols) {
      Rcpp::stop("nuisance must be empty or contain one p-by-p matrix per feature.");
    }
    vp.reserve(nuisance.size());
    for (R_xlen_t j = 0; j < nuisance.size(); ++j) {
      vp.emplace_back();
      if (nuisance[j] == R_NilValue) {
        pre_errors[j] = "nuisance covariance entry is NULL.";
        continue;
      }
      try {
        vp.back() = Rcpp::as<arma::mat>(nuisance[j]);
      } catch (const std::exception& ex) {
        pre_errors[j] = std::string("invalid nuisance covariance: ") + ex.what();
        continue;
      }
      if (vp.back().n_rows != X.n_cols || vp.back().n_cols != X.n_cols ||
          !finite_matrix(vp.back())) {
        pre_errors[j] = "nuisance covariance must be a finite p-by-p matrix.";
      }
    }
  }

  std::vector<DenseFeatureResult> results(variance.n_cols);
  const int nthreads = std::max(1, std::min(threads,
                                             static_cast<int>(variance.n_cols)));
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) schedule(dynamic, 1)
#endif
  for (int jj = 0; jj < static_cast<int>(variance.n_cols); ++jj) {
    try {
      if (!pre_errors[jj].empty()) {
        results[jj].error = pre_errors[jj];
        continue;
      }
      results[jj] = one_feature(T0, variance.colptr(jj), error.colptr(jj), scale(jj), X,
                                fixed_projection ? NULL : &vp[jj], score_only);
    } catch (const std::exception& ex) {
      results[jj].error = ex.what();
    } catch (...) {
      results[jj].error = "Unknown dense score failure.";
    }
  }

  Rcpp::List ans(variance.n_cols);
  for (R_xlen_t j = 0; j < variance.n_cols; ++j) {
    if (!results[j].error.empty()) {
      ans[j] = Rcpp::List::create(Rcpp::Named("error") = results[j].error);
    } else {
      Rcpp::NumericVector avec(results[j].a.n_elem);
      std::copy(results[j].a.begin(), results[j].a.end(), avec.begin());
      ans[j] = Rcpp::List::create(
        Rcpp::Named("a") = avec,
        Rcpp::Named("H") = score_only ? SEXP(R_NilValue) : Rcpp::wrap(results[j].H),
        Rcpp::Named("error") = R_NilValue);
    }
  }
  return ans;
}
