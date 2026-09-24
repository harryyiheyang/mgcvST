#define EIGEN_DONT_PARALLELIZE
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include "liu_tail.h"

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
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
    // Keep BLAS calls out of the outer OpenMP region. The conda pthread
    // OpenBLAS build cannot safely create its own workers here.
    Eigen::Map<const Eigen::MatrixXd> left(matrixPointers[i], q, q);
    Eigen::Map<const Eigen::MatrixXd> right(matrixPointers[j], q, q);
    Eigen::MatrixXd product = left * right;
    out(k, 0) = product.trace();

    if (maxPower >= 2) {
      out(k, 1) = (product.array() * product.transpose().array()).sum();
      if (maxPower >= 3) {
        Eigen::MatrixXd product2 = product * product;
        out(k, 2) =
          (product2.array() * product.transpose().array()).sum();
        if (maxPower >= 4) {
          out(k, 3) =
            (product2.array() * product2.transpose().array()).sum();
        }
      }
    }
  }

  return out;
}

// Fused exact Liu pair kernel: score, trace powers of H_l*H_r, and the
// log-space Liu tail, all in double precision with a full-rank H. H holds
// one q x q symmetric state matrix per local feature (state `M`); a holds
// one q-vector per local feature as columns (state `a`). left/right are
// 1-based local feature indices; pairs must be pre-sorted by `left` (the R
// caller sorts before the call and restores the original order after).
// Work is split into runs of equal `left` chunked to grain 32 and scheduled
// dynamically, mirroring mgcvst_fp16_pairs_cpp in inla_fp16.cpp.
// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
// [[Rcpp::export]]
Rcpp::List mgcvst_pair_liu_cpp(const Rcpp::List& H, const Rcpp::NumericMatrix& a,
                               const Rcpp::IntegerVector& left,
                               const Rcpp::IntegerVector& right,
                               int threads = 1) {
  const int n = H.size();
  if (n < 1) Rcpp::stop("H must contain at least one matrix.");
  const int q = a.nrow();
  if (q < 1) Rcpp::stop("a must have at least one row.");
  if (a.ncol() != n) Rcpp::stop("a must have one column per element of H.");
  if (threads < 1) Rcpp::stop("threads must be a positive integer.");
#ifndef _OPENMP
  if (threads > 1) {
    Rcpp::stop("mgcvST was compiled without OpenMP support; use threads = 1.");
  }
#endif
  std::vector<const double*> Hptr(n);
  for (int i = 0; i < n; ++i) {
    Rcpp::NumericMatrix current(H[i]);
    if (current.nrow() != q || current.ncol() != q) {
      Rcpp::stop("Every matrix in H must be square with dimension nrow(a).");
    }
    Hptr[i] = current.begin();
  }
  const R_xlen_t N = left.size();
  if (right.size() != N) Rcpp::stop("left and right must have equal length.");
  for (R_xlen_t k = 0; k < N; ++k) {
    if (left[k] == NA_INTEGER || right[k] == NA_INTEGER ||
        left[k] < 1 || left[k] > n || right[k] < 1 || right[k] > n) {
      Rcpp::stop("left and right must index elements of H.");
    }
    if (k && left[k] < left[k - 1]) {
      Rcpp::stop("left must be sorted in non-decreasing order.");
    }
  }

  struct Task { int left; R_xlen_t first, last; };
  std::vector<Task> tasks;
  const R_xlen_t grain = 32;
  {
    R_xlen_t k = 0;
    while (k < N) {
      R_xlen_t end = k;
      while (end < N && left[end] == left[k]) ++end;
      for (R_xlen_t s = k; s < end; s += grain) {
        tasks.push_back(Task{left[k] - 1, s, std::min(end, s + grain)});
      }
      k = end;
    }
  }

  Rcpp::NumericVector score(N), information(N), effective_rank(N),
    log_p_two_sided(N), log_p_positive(N), log_p_negative(N);
  Rcpp::IntegerVector status(N);
  const int* rightp = right.begin();
  Eigen::Map<const Eigen::MatrixXd> A(a.begin(), q, n);
  double* scorep = score.begin();
  double* infop = information.begin();
  double* rankp = effective_rank.begin();
  double* lp2p = log_p_two_sided.begin();
  double* lp1p = log_p_positive.begin();
  double* lp0p = log_p_negative.begin();
  int* statusp = status.begin();
  const long ntask = tasks.size();
#ifdef _OPENMP
#pragma omp parallel num_threads(threads)
#endif
  {
    Eigen::MatrixXd P1(q, q), P2(q, q);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
    for (long t = 0; t < ntask; ++t) {
      const Task& task = tasks[t];
      const int l = task.left;
      Eigen::Map<const Eigen::MatrixXd> Lm(Hptr[l], q, q);
      for (R_xlen_t k = task.first; k < task.last; ++k) {
        const int r = rightp[k] - 1;
        Eigen::Map<const Eigen::MatrixXd> Rm(Hptr[r], q, q);
        P1.noalias() = Lm * Rm;
        P2.noalias() = P1 * P1;
        double v[4] = {0, 0, 0, 0};
        for (int i = 0; i < q; ++i) v[0] += P1(i, i);
        for (int c = 0; c < q; ++c) {
          for (int i = 0; i < q; ++i) {
            const double p = P1(i, c), pt = P1(c, i), h = P2(i, c), ht = P2(c, i);
            v[1] += p * pt;
            v[2] += h * pt;
            v[3] += h * ht;
          }
        }
        const double sc = A.col(l).dot(A.col(r));
        scorep[k] = sc;
        const bool good = std::isfinite(sc) && std::isfinite(v[0]) &&
          std::isfinite(v[1]) && std::isfinite(v[2]) && std::isfinite(v[3]) &&
          v[0] > 1e-10 && v[1] > 0 && v[2] > 0 && v[3] > 0;
        if (!good) {
          infop[k] = rankp[k] = NA_REAL;
          lp2p[k] = lp1p[k] = lp0p[k] = NA_REAL;
          statusp[k] = 1;
          continue;
        }
        infop[k] = v[0];
        rankp[k] = v[0] * v[0] / v[1];
        double lp[3];
        mgcvst_liu::liu_log_p(sc, v[0], v[1], v[2], v[3], lp);
        lp2p[k] = lp[0];
        lp1p[k] = lp[1];
        lp0p[k] = lp[2];
        statusp[k] = std::isfinite(lp[0]) ? 0 : 2;
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("score") = score, Rcpp::Named("information") = information,
    Rcpp::Named("effective_rank") = effective_rank,
    Rcpp::Named("log_p_two_sided") = log_p_two_sided,
    Rcpp::Named("log_p_positive") = log_p_positive,
    Rcpp::Named("log_p_negative") = log_p_negative,
    Rcpp::Named("status") = status
  );
}

