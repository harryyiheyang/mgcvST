#include <Rcpp.h>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

// A marginal TAPS spectral-moment kernel. This is intentionally separate from
// score_pair.cpp: sum(lambda^k), NOT trace((H1 H2)^k) or squared-bilinear moments.
// R constructs the powers before entry. No R API, allocation or tail-probability
// function is called inside the parallel region. Each sum retains column order.
// [[Rcpp::export]]
Rcpp::NumericMatrix mgcvst_marginal_liu_moments_cpp(Rcpp::List powers, int threads = 1) {
  if (threads < 1) Rcpp::stop("threads must be positive.");
#ifndef _OPENMP
  if (threads > 1) Rcpp::stop("mgcvST was compiled without OpenMP; use threads = 1.");
#endif
  int n = powers.size();
  std::vector<const double*> ptr(n);
  std::vector<int> rows(n);
  for (int i = 0; i < n; ++i) {
    if (!Rf_isReal(powers[i]) || !Rf_isMatrix(powers[i]))
      Rcpp::stop("Each spectral power block must be a numeric matrix.");
    Rcpp::NumericMatrix x(powers[i]);
    if (x.ncol() != 4 || x.nrow() < 1)
      Rcpp::stop("Spectral power blocks must have four columns and at least one row.");
    ptr[i] = x.begin();
    rows[i] = x.nrow();
  }
  Rcpp::NumericMatrix out(n, 4);
  double* dest = out.begin();
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < 4; ++k) {
      long double sum = 0.0;
      for (int j = 0; j < rows[i]; ++j) sum += ptr[i][j + k * rows[i]];
      dest[i + k * n] = static_cast<double>(sum);
    }
  }
  return out;
}
