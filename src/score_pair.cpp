#include <RcppArmadillo.h>

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat mgcvst_pair_trace_powers_cpp(const Rcpp::List& matrixList,
                                       const Rcpp::IntegerMatrix& pairs,
                                       int maxPower = 4,
                                       int threads = 1) {
  int n = matrixList.size();
  if (n < 1) {
    Rcpp::stop("matrixList must contain at least one matrix.");
  }
  if (pairs.ncol() != 2 || pairs.nrow() < 1) {
    Rcpp::stop("pairs must be a non-empty two-column integer matrix.");
  }
  if (maxPower < 1 || maxPower > 4) {
    Rcpp::stop("maxPower must be an integer from 1 to 4.");
  }
  if (threads < 1) {
    Rcpp::stop("threads must be a positive integer.");
  }
#ifndef _OPENMP
  if (threads > 1) {
    Rcpp::stop("mgcvST was compiled without OpenMP support; use threads = 1.");
  }
#endif

  Rcpp::NumericMatrix first(matrixList[0]);
  int q = first.nrow();
  if (q < 1 || first.ncol() != q) {
    Rcpp::stop("Every matrix must be non-empty and square.");
  }
  std::vector<double*> matrixPointers(n);
  matrixPointers[0] = first.begin();
  for (int i = 1; i < n; ++i) {
    Rcpp::NumericMatrix current(matrixList[i]);
    if (current.nrow() != q || current.ncol() != q) {
      Rcpp::stop("Every matrix must have the same square dimension.");
    }
    matrixPointers[i] = current.begin();
  }
  for (int k = 0; k < pairs.nrow(); ++k) {
    int i = pairs(k, 0) - 1;
    int j = pairs(k, 1) - 1;
    if (i < 0 || i >= n || j < 0 || j >= n) {
      Rcpp::stop("pairs contains an index outside matrixList.");
    }
  }

  arma::mat out(pairs.nrow(), maxPower, arma::fill::none);
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
  for (int k = 0; k < pairs.nrow(); ++k) {
    int i = pairs(k, 0) - 1;
    int j = pairs(k, 1) - 1;
    arma::mat left(matrixPointers[i], q, q, false, true);
    arma::mat right(matrixPointers[j], q, q, false, true);
    arma::mat product = left * right;
    out(k, 0) = arma::trace(product);

    if (maxPower >= 2) {
      out(k, 1) = arma::accu(product % product.t());
      if (maxPower >= 3) {
        arma::mat product2 = product * product;
        out(k, 2) = arma::accu(product2 % product.t());
        if (maxPower >= 4) {
          out(k, 3) = arma::accu(product2 % product2.t());
        }
      }
    }
  }

  return out;
}
