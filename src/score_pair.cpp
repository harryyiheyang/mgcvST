#define EIGEN_DONT_PARALLELIZE
#include <RcppEigen.h>
#include "spa_pair.h"
#include "la_dispatch.h"

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Rcpp::NumericMatrix mgcvst_pair_trace_powers_cpp(const Rcpp::List& matrixList,
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

  const R_xlen_t np = pairs.nrow();
  Rcpp::NumericMatrix out(pairs.nrow(), maxPower);
  double* outp = out.begin();
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
    outp[k] = product.trace();

    if (maxPower >= 2) {
      outp[k + np] = (product.array() * product.transpose().array()).sum();
      if (maxPower >= 3) {
        Eigen::MatrixXd product2 = product * product;
        outp[k + 2 * np] =
          (product2.array() * product.transpose().array()).sum();
        if (maxPower >= 4) {
          outp[k + 3 * np] =
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
// [[Rcpp::depends(RcppEigen)]]
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
// [[Rcpp::depends(RcppEigen)]]
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
// dynamically. The parallel region uses only the SIMD kernels of la_kernel.cpp
// (products and the eigenvalues of the pair matrices) and R::pnorm. The result also holds the number of pairs
// whose remainder has a node above the largest leading value (u > s_1^2), a
// diagnostic of the compression. With screen_z > 0 a pair whose Gaussian
// z-score is below screen_z is not evaluated (status 4, see below).
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Rcpp::List mgcvst_pair_spa_cpp(const Rcpp::List& H, const Rcpp::List& G,
                               const Rcpp::NumericMatrix& a,
                               const Rcpp::IntegerVector& left,
                               const Rcpp::IntegerVector& right,
                               int threads = 1, int order = 4,
                               Rcpp::Nullable<Rcpp::NumericVector> x = R_NilValue,
                               double screen_z = 0) {
  const int n = H.size();
  if (n < 1) Rcpp::stop("H must contain at least one matrix.");
  const int q = a.nrow();
  if (q < 1) Rcpp::stop("a must have at least one row.");
  if (a.ncol() != n) Rcpp::stop("a must have one column per element of H.");
  if (G.size() != n) Rcpp::stop("G must have one matrix per element of H.");
  if (threads < 1) Rcpp::stop("threads must be a positive integer.");
  if (order != 2 && order != 4) Rcpp::stop("order must be 2 or 4.");
  if (!(screen_z == 0 || (std::isfinite(screen_z) && screen_z >= 2))) {
    Rcpp::stop("screen_z must be 0 (no screen) or a finite value of at least 2.");
  }
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
  long nodes_above = 0, screened = 0;
  const bool full = k == q;
  const mgcvst::LaFns& la = mgcvst::la();
#ifdef _OPENMP
#pragma omp parallel num_threads(threads)
#endif
  {
    long above = 0, skipped = 0;
    mgcvst::LaWork work(la, q, k, !full);
    Eigen::VectorXd sv(k);
    mgcvst_spa::Scratch scratch;
    scratch.reserve(k + 2);
    int current = -1;
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
    for (long t = 0; t < ntask; ++t) {
      const Task& task = tasks[t];
      const int l = task.left;
      if (!full && current != l && std::isfinite(scales[l]) && scales[l] > 0) {
        la.set_left(work.p, Hptr[l], scales[l]);
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
        const double units = std::sqrt(scales[l]) * std::sqrt(scales[r]);
        const double sc = explicit_x ? xp[j] : normalized_a.col(l).dot(normalized_a.col(r));
        if (explicit_x) scorep[j] = sc * units;
        if (!std::isfinite(sc)) continue;
        // Gaussian lower-bound screen. Given one of the two score vectors, U is
        // Gaussian with a random variance of mean t_1 = tr(H_l H_r), so for
        // z = |U| / sqrt(t_1) >= 2 the two-sided p-value is at least the
        // Gaussian tail 2 Phibar(z), and it decreases in |U|. A pair with
        // z < screen_z therefore has p >= 2 Phibar(screen_z): its spectrum is
        // not computed, and the lower bound 2 Phibar(max(z, 2)) is returned
        // with status 4.
        if (screen_z > 0) {
          const double t1 = la.trace1(Hptr[l], Hptr[r], q) / (scales[l] * scales[r]);
          const double z = t1 > 0 ? std::fabs(sc) / std::sqrt(t1) : R_PosInf;
          if (std::isfinite(z) && z < screen_z) {
            const double one = R::pnorm(std::max(z, 2.0), 0.0, 1.0, 0, 1);
            const double other = std::log1p(-std::exp(one));
            lp2p[j] = std::min(0.0, one + mgcvst_spa::kLn2);
            lp1p[j] = sc >= 0 ? one : other;
            lp0p[j] = sc >= 0 ? other : one;
            statusp[j] = 4;
            ++skipped;
            continue;
          }
        }
        // Exact trace moments t_s = tr((H_l H_r)^s) for the remainder. With
        // k = q the singular values are the whole spectrum, there is no
        // remainder, and the moments are the power sums of the singular
        // values: the two q x q products are not formed.
        double v[4] = {0, 0, 0, 0};
        if (!full) {
          la.moments(work.p, Hptr[r], scales[r], v);
          const bool good = std::isfinite(v[0]) && std::isfinite(v[1]) &&
            std::isfinite(v[2]) && std::isfinite(v[3]) &&
            v[0] > 0 && v[1] > 0 && v[2] > 0 && v[3] > 0;
          if (!good) continue;
        }
        // The k leading singular values: those of G_l' G_r (la_kernel.cpp).
        if (la.spectrum_tn(work.p, Gptr[l], Gptr[r], q, sv.data()) != 0) continue;
        long double lead[4];
        mgcvst_spa::leading_sums(sv.data(), k, lead);
        // Full spectrum: the trace moments are the power sums; a pair without a
        // positive singular value has non-positive moments (status 1).
        if (full && !(lead[0] > 0 && std::isfinite(static_cast<double>(lead[3])))) continue;
        // With k = q the singular values are the full spectrum of the pair and
        // the remainder is zero; total-minus-leading moments are then rounding
        // only and are not matched (no remainder).
        const mgcvst_spa::Remainder rem = full ? mgcvst_spa::Remainder() :
          mgcvst_spa::make_remainder(v, lead, order, sv[0]);
        double lp[3] = {NA_REAL, NA_REAL, NA_REAL};
        const int st = mgcvst_spa::spa_pair(sc, sv.data(), k, rem, scratch, lp);
        statusp[j] = st;
        if (st == 0) {
          lp2p[j] = lp[0];
          lp1p[j] = lp[1];
          lp0p[j] = lp[2];
          kindp[j] = rem.kind;
          if (mgcvst_spa::node_above_leading(rem, sv[0])) ++above;
        }
      }
    }
#ifdef _OPENMP
#pragma omp atomic
#endif
    nodes_above += above;
#ifdef _OPENMP
#pragma omp atomic
#endif
    screened += skipped;
  }

  return Rcpp::List::create(
    Rcpp::Named("score") = score,
    Rcpp::Named("log_p_two_sided") = log_p_two_sided,
    Rcpp::Named("log_p_positive") = log_p_positive,
    Rcpp::Named("log_p_negative") = log_p_negative,
    Rcpp::Named("remainder_kind") = remainder_kind,
    Rcpp::Named("status") = status,
    Rcpp::Named("nodes_above_leading") = (double)nodes_above,
    Rcpp::Named("screened") = (double)screened,
    Rcpp::Named("simd") = std::string(la.isa)
  );
}

// The instruction set of the SIMD kernels chosen for this process.
// [[Rcpp::export]]
std::string mgcvst_simd_cpp() { return mgcvst::la().isa; }
