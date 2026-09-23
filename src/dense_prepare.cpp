#define ARMA_WARN_LEVEL 0
#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <exception>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

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

bool finite_vector(const arma::vec& x) {
  return x.is_finite();
}

bool solve_checked(const arma::mat& A, const arma::mat& B, arma::mat& out) {
  return arma::solve(out, A, B);
}

bool solve_checked(const arma::mat& A, const arma::vec& b, arma::vec& out) {
  return arma::solve(out, A, b);
}

arma::mat apply_w(const arma::vec& Dinv,
                  const arma::mat& F,
                  const arma::mat& DinvF,
                  const arma::mat& M,
                  const arma::mat& rhs) {
  arma::mat Dinv_rhs = rhs.each_col() % Dinv;
  arma::mat Mrhs;
  if (!solve_checked(M, F.t() * Dinv_rhs, Mrhs)) {
    throw std::runtime_error("Woodbury solve failed.");
  }
  return Dinv_rhs - DinvF * Mrhs;
}

arma::mat apply_projection(const arma::mat& WZ,
                           const arma::mat& WX,
                           const arma::mat& XWZ,
                           const arma::mat& X,
                           const arma::mat* Vp) {
  if (X.n_cols == 0L) return WZ;
  arma::mat adjustment;
  if (Vp == NULL) {
    arma::mat XtWX = X.t() * WX;
    arma::mat pinv_XtWX = arma::pinv(XtWX, 5e-16);
    adjustment = WX * pinv_XtWX * XWZ;
  } else {
    adjustment = WX * (*Vp) * XWZ;
  }
  if (!finite_matrix(adjustment)) {
    throw std::runtime_error("Nuisance projection produced non-finite values.");
  }
  return WZ - adjustment;
}

DenseFeatureResult one_feature(const arma::mat& T0,
                               const arma::vec& variance,
                               const arma::vec& error,
                               double scale,
                               const arma::mat& X,
                               const arma::mat* Vp) {
  DenseFeatureResult out;
  const arma::uword n = T0.n_rows;
  const arma::uword q = T0.n_cols;
  if (!std::isfinite(scale) || scale <= 0.0) {
    out.error = "scale must be one positive finite value per feature.";
    return out;
  }
  if (variance.n_elem != n || error.n_elem != n ||
      !finite_vector(variance) || !finite_vector(error)) {
    out.error = "variance and error must be finite vectors with n observations.";
    return out;
  }
  if (arma::any(variance <= 0.0)) {
    out.error = "variance must contain strictly positive values.";
    return out;
  }

  arma::mat F = std::sqrt(scale) * T0;
  arma::vec Dinv_vec = 1.0 / variance;
  arma::mat DinvF = F.each_col() % Dinv_vec;
  arma::mat M = arma::eye<arma::mat>(q, q) + F.t() * DinvF;
  if (!finite_matrix(M)) {
    out.error = "Woodbury system contains non-finite values.";
    return out;
  }

  arma::mat WX;
  if (X.n_cols > 0L) {
    WX = apply_w(Dinv_vec, F, DinvF, M, X);
  }

  arma::mat rhs(n, q + 1L, arma::fill::zeros);
  rhs.cols(0L, q - 1L) = F;
  rhs.col(q) = error;
  arma::mat Wrhs = apply_w(Dinv_vec, F, DinvF, M, rhs);
  arma::mat WF = Wrhs.cols(0L, q - 1L);
  arma::mat We = Wrhs.col(q);
  arma::mat PF = apply_projection(WF, WX, X.n_cols ? X.t() * WF : arma::mat(), X, Vp);
  arma::mat Pe = apply_projection(We, WX, X.n_cols ? X.t() * We : arma::mat(), X, Vp);
  out.a = F.t() * Pe;
  out.H = F.t() * PF;
  out.H = 0.5 * (out.H + out.H.t());
  if (!finite_vector(out.a) || !finite_matrix(out.H)) {
    out.error = "Dense score result contains non-finite values.";
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
                                        int threads = 1) {
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
      arma::vec v = variance.col(jj);
      arma::vec e = error.col(jj);
      results[jj] = one_feature(T0, v, e, scale(jj), X,
                                fixed_projection ? NULL : &vp[jj]);
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
      ans[j] = Rcpp::List::create(Rcpp::Named("a") = avec,
                                   Rcpp::Named("H") = results[j].H,
                                   Rcpp::Named("error") = R_NilValue);
    }
  }
  return ans;
}
