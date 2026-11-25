#define EIGEN_DONT_PARALLELIZE
#include <RcppEigen.h>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

void check_threads(int threads) {
  if (threads < 1) Rcpp::stop("threads must be a positive integer.");
#ifndef _OPENMP
  if (threads > 1) {
    Rcpp::stop("spdeST was compiled without OpenMP support; use threads = 1.");
  }
#endif
}

std::vector<const double*> square_pointers(const Rcpp::List& H, int q,
                                           const char* name) {
  std::vector<const double*> out(H.size());
  for (int g = 0; g < H.size(); ++g) {
    Rcpp::NumericMatrix current(H[g]);
    if (current.nrow() != q || current.ncol() != q) {
      Rcpp::stop("Every matrix in %s must be q x q.", name);
    }
    out[g] = current.begin();
  }
  return out;
}

}  // namespace

// Moments of K cluster scores against P features. aS (q x K) holds the summed
// score vectors of the clusters and MS their summed state matrices; a (q x P)
// and M hold the score vectors and state matrices of the features. For
// cluster k and feature c: U = aS_k' a_c (the score), V = aS_k' M_c aS_k (its
// variance given the cluster score when the feature is unrelated) and
// T = tr(MS_k M_c) (its information). Features are split over the threads.
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Rcpp::List mgcvst_cluster_moments_cpp(const Rcpp::NumericMatrix& aS,
                                      const Rcpp::List& MS,
                                      const Rcpp::NumericMatrix& a,
                                      const Rcpp::List& M, int threads = 1) {
  const int q = aS.nrow(), K = aS.ncol(), P = a.ncol();
  if (q < 1 || K < 1) Rcpp::stop("aS must be a q x K matrix with K >= 1.");
  if (a.nrow() != q) Rcpp::stop("a must have q rows.");
  if (MS.size() != K) Rcpp::stop("MS must hold one matrix per column of aS.");
  if (M.size() != P) Rcpp::stop("M must hold one matrix per column of a.");
  check_threads(threads);
  const std::vector<const double*> Sp = square_pointers(MS, q, "MS");
  const std::vector<const double*> Mp = square_pointers(M, q, "M");
  Rcpp::NumericMatrix U(K, P), V(K, P), T(K, P);
  double* Up = U.begin();
  double* Vp = V.begin();
  double* Tp = T.begin();
  Eigen::Map<const Eigen::MatrixXd> A(aS.begin(), q, K);
  Eigen::Map<const Eigen::MatrixXd> F(a.begin(), q, P);
#ifdef _OPENMP
#pragma omp parallel num_threads(threads)
#endif
  {
    Eigen::MatrixXd W(q, K);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
    for (int c = 0; c < P; ++c) {
      Eigen::Map<const Eigen::MatrixXd> Mc(Mp[c], q, q);
      W.noalias() = Mc * A;
      for (int k = 0; k < K; ++k) {
        Eigen::Map<const Eigen::MatrixXd> Sk(Sp[k], q, q);
        Up[k + static_cast<R_xlen_t>(c) * K] = A.col(k).dot(F.col(c));
        Vp[k + static_cast<R_xlen_t>(c) * K] = A.col(k).dot(W.col(k));
        Tp[k + static_cast<R_xlen_t>(c) * K] = Sk.cwiseProduct(Mc.transpose()).sum();
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("U") = U, Rcpp::Named("V") = V,
                            Rcpp::Named("T") = T);
}

// Corrections of the moments above for the members left out of a cluster
// when one feature is tested against it. Row r pairs the feature `feature[r]`
// (column of a, state M) with the cluster `cluster[r]` (column of aS) and
// leaves out the members idx[ptr[r]], ..., idx[ptr[r + 1] - 1] (columns of
// ae, states Me); indices are one-based. With e the summed score vector of
// the members left out, the moments of the remaining cluster are U - dU,
// V + dV and T - dT, where dU = e' a_c, dV = e' M_c e - aS' (M_c + M_c') e and
// dT = sum over the members left out of tr(M_g M_c).
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Rcpp::List mgcvst_cluster_exclude_cpp(const Rcpp::NumericMatrix& aS,
                                      const Rcpp::NumericMatrix& a,
                                      const Rcpp::List& M,
                                      const Rcpp::NumericMatrix& ae,
                                      const Rcpp::List& Me,
                                      const Rcpp::IntegerVector& feature,
                                      const Rcpp::IntegerVector& cluster,
                                      const Rcpp::IntegerVector& ptr,
                                      const Rcpp::IntegerVector& idx,
                                      int threads = 1) {
  const int q = aS.nrow(), K = aS.ncol(), P = a.ncol(), E = ae.ncol();
  const int R = feature.size();
  if (a.nrow() != q || ae.nrow() != q) Rcpp::stop("a and ae must have q rows.");
  if (M.size() != P || Me.size() != E) {
    Rcpp::stop("M and Me must hold one matrix per column of a and ae.");
  }
  if (cluster.size() != R || ptr.size() != R + 1 || ptr[0] != 0 ||
      ptr[R] != idx.size()) {
    Rcpp::stop("feature, cluster, ptr and idx do not describe the same rows.");
  }
  for (int r = 0; r < R; ++r) {
    if (feature[r] < 1 || feature[r] > P || cluster[r] < 1 || cluster[r] > K ||
        ptr[r + 1] < ptr[r]) {
      Rcpp::stop("Invalid row %d.", r + 1);
    }
  }
  for (int i = 0; i < idx.size(); ++i) {
    if (idx[i] < 1 || idx[i] > E) Rcpp::stop("idx holds an invalid member.");
  }
  check_threads(threads);
  const std::vector<const double*> Mp = square_pointers(M, q, "M");
  const std::vector<const double*> Ep = square_pointers(Me, q, "Me");
  Rcpp::NumericVector dU(R), dV(R), dT(R);
  double* dUp = dU.begin();
  double* dVp = dV.begin();
  double* dTp = dT.begin();
  const int* fp = feature.begin();
  const int* cp = cluster.begin();
  const int* pp = ptr.begin();
  const int* ip = idx.begin();
  Eigen::Map<const Eigen::MatrixXd> A(aS.begin(), q, K);
  Eigen::Map<const Eigen::MatrixXd> F(a.begin(), q, P);
  Eigen::Map<const Eigen::MatrixXd> X(ae.begin(), q, E);
#ifdef _OPENMP
#pragma omp parallel num_threads(threads)
#endif
  {
    Eigen::VectorXd e(q), s1(q), s2(q);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 8)
#endif
    for (int r = 0; r < R; ++r) {
      const int c = fp[r] - 1, k = cp[r] - 1;
      Eigen::Map<const Eigen::MatrixXd> Mc(Mp[c], q, q);
      e.setZero();
      double t = 0.0;
      for (int i = pp[r]; i < pp[r + 1]; ++i) {
        const int g = ip[i] - 1;
        e += X.col(g);
        Eigen::Map<const Eigen::MatrixXd> Mg(Ep[g], q, q);
        t += Mg.cwiseProduct(Mc.transpose()).sum();
      }
      s1.noalias() = Mc * e;
      s2.noalias() = Mc.transpose() * e;
      dUp[r] = e.dot(F.col(c));
      dVp[r] = e.dot(s1) - A.col(k).dot(s1 + s2);
      dTp[r] = t;
    }
  }
  return Rcpp::List::create(Rcpp::Named("dU") = dU, Rcpp::Named("dV") = dV,
                            Rcpp::Named("dT") = dT);
}
