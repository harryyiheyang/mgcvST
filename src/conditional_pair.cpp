// Conditional pairwise variances and stable Cauchy calibration.
#include <Rcpp.h>
#include <R_ext/BLAS.h>
#include <R_ext/RS.h>
#include <cmath>
#include <limits>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

struct Cotangent {
  double log_abs;
  int sign;
};

Cotangent cotangent_from_z(double z, double log_p) {
  const double log_half = -std::log(2.0);
  if (log_p == log_half) return {-std::numeric_limits<double>::infinity(), 0};
  const double log_q = z == 0.0 ? -std::numeric_limits<double>::infinity() :
    std::log(std::erf(z / std::sqrt(2.0)));
  const double low = std::min(log_p, log_q);
  const int sign = log_p < log_half ? 1 : -1;
  const double log_abs = low < -30.0 ? -std::log(M_PI) - low :
    -std::log(std::tan(M_PI * std::exp(low)));
  return {log_abs, sign};
}

double cauchy_log_p(double z1, double z2) {
  const double log_p1 = std::min(0.0, std::log(2.0) +
    R::pnorm5(z1, 0.0, 1.0, false, true));
  const double log_p2 = std::min(0.0, std::log(2.0) +
    R::pnorm5(z2, 0.0, 1.0, false, true));
  const Cotangent a = cotangent_from_z(z1, log_p1);
  const Cotangent b = cotangent_from_z(z2, log_p2);
  if (a.sign == 0 && b.sign == 0) return -std::log(2.0);
  const double high = std::max(a.log_abs, b.log_abs);
  const double low = std::min(a.log_abs, b.log_abs);
  int sign = 0;
  double log_sum = -std::numeric_limits<double>::infinity();
  if (a.sign == 0) {
    sign = b.sign;
    log_sum = b.log_abs;
  } else if (b.sign == 0) {
    sign = a.sign;
    log_sum = a.log_abs;
  } else if (a.sign == b.sign) {
    sign = a.sign;
    log_sum = std::isinf(high) ? high : high + std::log1p(std::exp(low - high));
  } else if (a.log_abs == b.log_abs) {
    return -std::log(2.0);
  } else {
    sign = a.log_abs > b.log_abs ? a.sign : b.sign;
    log_sum = high + std::log(-std::expm1(low - high));
  }
  const double log_t = log_sum - std::log(2.0);
  if (sign > 0 && log_t > std::log(1e6)) {
    return -std::log(M_PI) - log_t;
  }
  if (sign < 0 && log_t > std::log(1e6)) {
    return std::log1p(-std::exp(-std::log(M_PI) - log_t));
  }
  const double t = sign * std::exp(log_t);
  return std::log(std::atan2(1.0, t) / M_PI);
}

}  // namespace

// [[Rcpp::export]]
Rcpp::NumericMatrix mgcvst_conditional_variance_rows_cpp(
    const Rcpp::NumericMatrix& A, const Rcpp::List& matrices,
    int threads = 1, int block_size = 512) {
  const int q = A.nrow();
  const int G = A.ncol();
  const int n = matrices.size();
  if (q < 1 || G < 2 || n < 1 || threads < 1 || block_size < 1) {
    Rcpp::stop("Conditional variance inputs need positive dimensions, threads and block size.");
  }
#ifndef _OPENMP
  if (threads > 1) Rcpp::stop("mgcvST was compiled without OpenMP; use threads = 1.");
#endif
  std::vector<const double*> matrix_ptr(n);
  for (int i = 0; i < n; ++i) {
    Rcpp::NumericMatrix M = matrices[i];
    if (M.nrow() != q || M.ncol() != q) {
      Rcpp::stop("Every conditional covariance must match the score width.");
    }
    matrix_ptr[i] = M.begin();
  }
  Rcpp::NumericMatrix V(n, G);
  const double alpha = 1.0;
  const double beta = 0.0;
  const double* a = A.begin();
  double* v = V.begin();
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
#endif
  for (int i = 0; i < n; ++i) {
    std::vector<double> buffer(static_cast<size_t>(q) *
                               static_cast<size_t>(std::min(block_size, G)));
    for (int first = 0; first < G; first += block_size) {
      const int width = std::min(block_size, G - first);
      F77_CALL(dsymm)("L", "U", &q, &width, &alpha, matrix_ptr[i], &q,
                      a + static_cast<size_t>(first) * q, &q, &beta,
                      buffer.data(), &q FCONE FCONE);
      for (int col = 0; col < width; ++col) {
        const double* score = a + static_cast<size_t>(first + col) * q;
        const double* product = buffer.data() + static_cast<size_t>(col) * q;
        double value = 0.0;
        for (int row = 0; row < q; ++row) value += score[row] * product[row];
        v[i + static_cast<size_t>(first + col) * n] = value;
      }
    }
  }
  return V;
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix mgcvst_conditional_all_pairs_cpp(int n) {
  if (n < 2) Rcpp::stop("At least two features are required.");
  const double count = static_cast<double>(n) * (n - 1) / 2;
  if (count > std::numeric_limits<int>::max()) {
    Rcpp::stop("The pair count exceeds R's matrix row limit.");
  }
  Rcpp::IntegerMatrix out(static_cast<int>(count), 2);
  int at = 0;
  for (int i = 0; i < n - 1; ++i) {
    for (int j = i + 1; j < n; ++j) {
      out(at, 0) = i + 1;
      out(at, 1) = j + 1;
      ++at;
    }
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::List mgcvst_conditional_pairs_cpp(
    const Rcpp::NumericMatrix& S, const Rcpp::NumericMatrix& V,
    const Rcpp::IntegerMatrix& pairs) {
  const int n = S.nrow();
  if (S.ncol() != n || V.nrow() != n || V.ncol() != n || pairs.ncol() != 2) {
    Rcpp::stop("Conditional scores, variances and pair indices are incompatible.");
  }
  const int m = pairs.nrow();
  Rcpp::NumericVector score(m);
  Rcpp::NumericVector log_p(m);
  for (int k = 0; k < m; ++k) {
    const int i = pairs(k, 0) - 1;
    const int j = pairs(k, 1) - 1;
    if (i < 0 || j < 0 || i >= n || j >= n || i == j) {
      Rcpp::stop("Conditional pair indices must name distinct valid features.");
    }
    const double s = S(i, j);
    const double v1 = V(i, j);
    const double v2 = V(j, i);
    if (!std::isfinite(s) || !std::isfinite(v1) || !std::isfinite(v2) ||
        v1 <= 0.0 || v2 <= 0.0) {
      Rcpp::stop("A conditional pair has a non-finite score or non-positive variance.");
    }
    score[k] = s;
    log_p[k] = cauchy_log_p(std::abs(s) / std::sqrt(v1),
                           std::abs(s) / std::sqrt(v2));
  }
  return Rcpp::List::create(Rcpp::Named("score") = score,
                            Rcpp::Named("log_p") = log_p);
}
