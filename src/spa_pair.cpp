#include "spa_pair.h"

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::plugins(openmp)]]

// Saddlepoint p-values of signed scores from leading singular values and the
// total power sums: U (length N), S (k x N, or k x 1 for a common spectrum),
// Tm (4 x N, or 4 x 1) with Tm[r, ] = sum_i s_i^(2r) of the full spectrum.
// Columns: log p two-sided, positive, negative, remainder kind (0 none, 1 one
// node, 2 two nodes, 3 Gaussian) and status (0 evaluated, 1 invalid input,
// 2 invalid p-value). This is the unit-test entry of the shared kernel.
// [[Rcpp::export]]
Rcpp::NumericMatrix mgcvst_spa_cpp(const Rcpp::NumericVector& U,
                                   const Rcpp::NumericMatrix& S,
                                   const Rcpp::NumericMatrix& Tm,
                                   int order = 4, int threads = 1) {
  const R_xlen_t N = U.size();
  const int k = S.nrow();
  if (Tm.nrow() != 4) Rcpp::stop("Tm must have four rows.");
  if (!(S.ncol() == 1 || S.ncol() == N) || !(Tm.ncol() == 1 || Tm.ncol() == N)) {
    Rcpp::stop("S and Tm must have one column or one column per score.");
  }
  if (order != 2 && order != 4) Rcpp::stop("order must be 2 or 4.");
  if (threads < 1) Rcpp::stop("threads must be a positive integer.");
#ifndef _OPENMP
  if (threads > 1) Rcpp::stop("mgcvST was compiled without OpenMP support; use threads = 1.");
#endif
  Rcpp::NumericMatrix out(N, 5);
  double* o = out.begin();
  const double* sp = S.begin();
  const double* tp = Tm.begin();
  const bool s_common = S.ncol() == 1, t_common = Tm.ncol() == 1;
#ifdef _OPENMP
#pragma omp parallel num_threads(threads)
#endif
  {
    mgcvst_spa::Scratch scratch;
    scratch.reserve(k + 2);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (R_xlen_t j = 0; j < N; ++j) {
      const double* s = sp + (s_common ? 0 : j * k);
      const double* t = tp + (t_common ? 0 : 4 * j);
      double lp[3] = {NA_REAL, NA_REAL, NA_REAL};
      int status = 1, kind = 0;
      bool ok = std::isfinite(U[j]);
      for (int r = 0; r < 4; ++r) ok = ok && std::isfinite(t[r]) && t[r] > 0;
      if (ok) {
        long double lead[4];
        mgcvst_spa::leading_sums(s, k, lead);
        const mgcvst_spa::Remainder rem = mgcvst_spa::make_remainder(t, lead, order);
        status = mgcvst_spa::spa_pair(U[j], s, k, rem, scratch, lp);
        if (status == 0) kind = rem.kind;
      }
      o[j] = lp[0];
      o[j + N] = lp[1];
      o[j + 2 * N] = lp[2];
      o[j + 3 * N] = kind;
      o[j + 4 * N] = status;
    }
  }
  Rcpp::colnames(out) = Rcpp::CharacterVector::create(
    "log_p_two_sided", "log_p_positive", "log_p_negative", "remainder_kind", "status");
  return out;
}
