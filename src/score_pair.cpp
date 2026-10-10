#define EIGEN_DONT_PARALLELIZE
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include "spa_pair.h"

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

// Sum of the normalized state matrices, sum_g sym(M_g) / max|M_g|, over the
// features of a batch: the matrix whose leading eigenvectors are the shared
// basis V of the exact route. The sum runs in feature order, so the result does
// not depend on the thread count; a feature with a non-finite or zero matrix
// is skipped. `init` is the running sum to continue from (NULL: zero), so that
// a sum over many batches is the same as one sum over all features.
// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
// [[Rcpp::export]]
Rcpp::NumericMatrix mgcvst_pair_basis_sum_cpp(
    const Rcpp::List& H,
    Rcpp::Nullable<Rcpp::NumericMatrix> init = R_NilValue) {
  const int n = H.size();
  if (n < 1) Rcpp::stop("H must contain at least one matrix.");
  Rcpp::NumericMatrix first(H[0]);
  const int q = first.nrow();
  if (q < 1 || first.ncol() != q) Rcpp::stop("Every matrix in H must be square.");
  Eigen::MatrixXd S = Eigen::MatrixXd::Zero(q, q);
  if (init.isNotNull()) {
    Rcpp::NumericMatrix start(init.get());
    if (start.nrow() != q || start.ncol() != q) {
      Rcpp::stop("init must be a q x q matrix.");
    }
    S = Eigen::Map<const Eigen::MatrixXd>(start.begin(), q, q);
  }
  int used = 0;
  for (int g = 0; g < n; ++g) {
    Rcpp::NumericMatrix current(H[g]);
    if (current.nrow() != q || current.ncol() != q) {
      Rcpp::stop("Every matrix in H must have the same dimension.");
    }
    Eigen::Map<const Eigen::MatrixXd> M(current.begin(), q, q);
    if (!M.allFinite()) continue;
    const double scale = M.cwiseAbs().maxCoeff();
    if (!(scale > 0)) continue;
    S += (0.5 / scale) * (M + M.transpose());
    ++used;
  }
  Rcpp::NumericMatrix out = Rcpp::wrap(S);
  out.attr("used") = used;
  return out;
}

// Pair basis G_g = H_g^{1/2} V (q x k) for every feature: H_g is symmetrized
// and normalized by max|M_g|, and H^{1/2} = E sqrt(D+) E' is the symmetric
// square root (not the one-sided factor E sqrt(D)). With V = I and k = q the
// compression is exact. A feature with a non-finite or zero matrix gets a
// matrix of NA.
// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
// [[Rcpp::export]]
Rcpp::List mgcvst_pair_basis_cpp(const Rcpp::List& H, const Rcpp::NumericMatrix& V,
                                 int threads = 1) {
  const int n = H.size();
  const int q = V.nrow(), k = V.ncol();
  if (n < 1) Rcpp::stop("H must contain at least one matrix.");
  if (q < 1 || k < 1 || k > q) Rcpp::stop("V must be a q x k matrix with 1 <= k <= q.");
  if (threads < 1) Rcpp::stop("threads must be a positive integer.");
#ifndef _OPENMP
  if (threads > 1) {
    Rcpp::stop("mgcvST was compiled without OpenMP support; use threads = 1.");
  }
#endif
  std::vector<const double*> Hptr(n);
  Rcpp::List out(n);
  std::vector<double*> Gptr(n);
  for (int g = 0; g < n; ++g) {
    Rcpp::NumericMatrix current(H[g]);
    if (current.nrow() != q || current.ncol() != q) {
      Rcpp::stop("Every matrix in H must be square with dimension nrow(V).");
    }
    Hptr[g] = current.begin();
    Rcpp::NumericMatrix G(q, k);
    out[g] = G;
    Gptr[g] = G.begin();
  }
  Eigen::Map<const Eigen::MatrixXd> Vm(V.begin(), q, k);
#ifdef _OPENMP
#pragma omp parallel num_threads(threads)
#endif
  {
    Eigen::MatrixXd Hn(q, q), W(q, k);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
    for (int g = 0; g < n; ++g) {
      Eigen::Map<const Eigen::MatrixXd> M(Hptr[g], q, q);
      Eigen::Map<Eigen::MatrixXd> G(Gptr[g], q, k);
      const double scale = M.allFinite() ? M.cwiseAbs().maxCoeff() : NA_REAL;
      if (!(std::isfinite(scale) && scale > 0)) {
        G.setConstant(NA_REAL);
        continue;
      }
      Hn = (0.5 / scale) * (M + M.transpose());
      solver.compute(Hn);
      if (solver.info() != Eigen::Success) {
        G.setConstant(NA_REAL);
        continue;
      }
      W.noalias() = solver.eigenvectors().transpose() * Vm;
      for (int i = 0; i < q; ++i) {
        const double d = solver.eigenvalues()[i];
        W.row(i) *= d > 0 ? std::sqrt(d) : 0.0;
      }
      G.noalias() = solver.eigenvectors() * W;
    }
  }
  return out;
}

// Exact-moment saddlepoint pair kernel. H holds one q x q symmetric state
// matrix per local feature (state `M`), G the matching q x k pair bases from
// mgcvst_pair_basis_cpp(), a one q-vector per local feature as columns (state
// `a`). left/right are 1-based local feature indices; pairs must be sorted by
// `left`. For each pair the exact trace moments t_s = tr((H_l H_r)^s) of the
// normalized matrices give the remainder, and the k leading singular values
// are those of G_l' G_r. order = 4 matches four remainder moments with two
// nodes, order = 2 uses the Satterthwaite node. `x` optionally replaces the
// normalized score (for validation against stored scores).
// Work is split into runs of equal `left` chunked to grain 32 and scheduled
// dynamically. The parallel region uses only Eigen products and decompositions
// of the pair matrices and R::pnorm.
// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
// [[Rcpp::export]]
Rcpp::List mgcvst_pair_spa_cpp(const Rcpp::List& H, const Rcpp::List& G,
                               const Rcpp::NumericMatrix& a,
                               const Rcpp::IntegerVector& left,
                               const Rcpp::IntegerVector& right,
                               int threads = 1, int order = 4,
                               Rcpp::Nullable<Rcpp::NumericVector> x = R_NilValue) {
  const int n = H.size();
  if (n < 1) Rcpp::stop("H must contain at least one matrix.");
  const int q = a.nrow();
  if (q < 1) Rcpp::stop("a must have at least one row.");
  if (a.ncol() != n) Rcpp::stop("a must have one column per element of H.");
  if (G.size() != n) Rcpp::stop("G must have one matrix per element of H.");
  if (threads < 1) Rcpp::stop("threads must be a positive integer.");
  if (order != 2 && order != 4) Rcpp::stop("order must be 2 or 4.");
#ifndef _OPENMP
  if (threads > 1) {
    Rcpp::stop("mgcvST was compiled without OpenMP support; use threads = 1.");
  }
#endif
  // Normalize each feature before matrix products; retain public score units.
  std::vector<const double*> Hptr(n), Gptr(n);
  Eigen::VectorXd scales(n);
  Eigen::MatrixXd normalized_a(q, n);
  int k = -1;
  for (int i = 0; i < n; ++i) {
    Rcpp::NumericMatrix current(H[i]);
    if (current.nrow() != q || current.ncol() != q) {
      Rcpp::stop("Every matrix in H must be square with dimension nrow(a).");
    }
    Hptr[i] = current.begin();
    Rcpp::NumericMatrix basis(G[i]);
    if (basis.nrow() != q || basis.ncol() < 1 || basis.ncol() > q) {
      Rcpp::stop("Every matrix in G must be q x k with 1 <= k <= q.");
    }
    if (k < 0) k = basis.ncol();
    if (basis.ncol() != k) Rcpp::stop("Every matrix in G must have the same width.");
    Gptr[i] = basis.begin();
    Eigen::Map<const Eigen::MatrixXd> M(current.begin(), q, q);
    scales[i] = M.allFinite() ? M.cwiseAbs().maxCoeff() : NA_REAL;
    if (std::isfinite(scales[i]) && scales[i] > 0) {
      for (int j = 0; j < q; ++j) normalized_a(j, i) = a(j, i) / std::sqrt(scales[i]);
    } else normalized_a.col(i).setConstant(NA_REAL);
  }
  const R_xlen_t N = left.size();
  if (right.size() != N) Rcpp::stop("left and right must have equal length.");
  for (R_xlen_t j = 0; j < N; ++j) {
    if (left[j] == NA_INTEGER || right[j] == NA_INTEGER ||
        left[j] < 1 || left[j] > n || right[j] < 1 || right[j] > n) {
      Rcpp::stop("left and right must index elements of H.");
    }
    if (j && left[j] < left[j - 1]) {
      Rcpp::stop("left must be sorted in non-decreasing order.");
    }
  }
  const bool explicit_x = x.isNotNull();
  Rcpp::NumericVector xv = explicit_x ? Rcpp::NumericVector(x.get()) : Rcpp::NumericVector(0);
  if (explicit_x && xv.size() != N) Rcpp::stop("x must have one value per pair.");

  struct Task { int left; R_xlen_t first, last; };
  std::vector<Task> tasks;
  const R_xlen_t grain = 32;
  {
    R_xlen_t j = 0;
    while (j < N) {
      R_xlen_t end = j;
      while (end < N && left[end] == left[j]) ++end;
      for (R_xlen_t s = j; s < end; s += grain) {
        tasks.push_back(Task{left[j] - 1, s, std::min(end, s + grain)});
      }
      j = end;
    }
  }

  Rcpp::NumericVector score(N), log_p_two_sided(N), log_p_positive(N),
    log_p_negative(N);
  Rcpp::IntegerVector remainder_kind(N), status(N);
  const int* rightp = right.begin();
  Eigen::Map<const Eigen::MatrixXd> A(a.begin(), q, n);
  double* scorep = score.begin();
  double* lp2p = log_p_two_sided.begin();
  double* lp1p = log_p_positive.begin();
  double* lp0p = log_p_negative.begin();
  int* kindp = remainder_kind.begin();
  int* statusp = status.begin();
  const double* xp = explicit_x ? xv.begin() : nullptr;
  const long ntask = tasks.size();
#ifdef _OPENMP
#pragma omp parallel num_threads(threads)
#endif
  {
    Eigen::MatrixXd P1(q, q), P2(q, q), Ln(q, q), Rn(q, q), Mk(k, k);
    Eigen::BDCSVD<Eigen::MatrixXd> svd(k, k, 0);
    mgcvst_spa::Scratch scratch;
    scratch.reserve(k + 2);
    int current = -1;
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
    for (long t = 0; t < ntask; ++t) {
      const Task& task = tasks[t];
      const int l = task.left;
      Eigen::Map<const Eigen::MatrixXd> Lm(Hptr[l], q, q);
      Eigen::Map<const Eigen::MatrixXd> Gl(Gptr[l], q, k);
      if (current != l && std::isfinite(scales[l]) && scales[l] > 0) {
        Ln = (0.5 / scales[l]) * (Lm + Lm.transpose());
        current = l;
      }
      for (R_xlen_t j = task.first; j < task.last; ++j) {
        const int r = rightp[j] - 1;
        scorep[j] = A.col(l).dot(A.col(r));
        lp2p[j] = lp1p[j] = lp0p[j] = NA_REAL;
        kindp[j] = 0;
        statusp[j] = 1;
        if (!std::isfinite(scales[l]) || scales[l] <= 0 ||
            !std::isfinite(scales[r]) || scales[r] <= 0) continue;
        Eigen::Map<const Eigen::MatrixXd> Rm(Hptr[r], q, q);
        Eigen::Map<const Eigen::MatrixXd> Gr(Gptr[r], q, k);
        Rn = (0.5 / scales[r]) * (Rm + Rm.transpose());
        P1.noalias() = Ln * Rn;
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
        const double units = std::sqrt(scales[l]) * std::sqrt(scales[r]);
        const double sc = explicit_x ? xp[j] : normalized_a.col(l).dot(normalized_a.col(r));
        if (explicit_x) scorep[j] = sc * units;
        const bool good = std::isfinite(sc) && std::isfinite(v[0]) &&
          std::isfinite(v[1]) && std::isfinite(v[2]) && std::isfinite(v[3]) &&
          v[0] > 0 && v[1] > 0 && v[2] > 0 && v[3] > 0;
        if (!good) continue;
        Mk.noalias() = Gl.transpose() * Gr;
        if (!Mk.allFinite()) continue;
        svd.compute(Mk);
        const Eigen::VectorXd& sv = svd.singularValues();
        long double lead[4];
        mgcvst_spa::leading_sums(sv.data(), k, lead);
        const mgcvst_spa::Remainder rem = mgcvst_spa::make_remainder(v, lead, order);
        double lp[3] = {NA_REAL, NA_REAL, NA_REAL};
        const int st = mgcvst_spa::spa_pair(sc, sv.data(), k, rem, scratch, lp);
        statusp[j] = st;
        if (st == 0) {
          lp2p[j] = lp[0];
          lp1p[j] = lp[1];
          lp0p[j] = lp[2];
          kindp[j] = rem.kind;
        }
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("score") = score,
    Rcpp::Named("log_p_two_sided") = log_p_two_sided,
    Rcpp::Named("log_p_positive") = log_p_positive,
    Rcpp::Named("log_p_negative") = log_p_negative,
    Rcpp::Named("remainder_kind") = remainder_kind,
    Rcpp::Named("status") = status
  );
}
