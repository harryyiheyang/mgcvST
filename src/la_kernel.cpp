// Dense linear-algebra kernels of the pair tests, compiled once per instruction
// set and chosen at run time (la_dispatch.h).
//
// Eigen picks its SIMD width at preprocessing time from __AVX2__ / __FMA__ /
// __AVX512F__, which R's default package flags never define, so Eigen products
// in ordinary translation units run on SSE2. This file is compiled three times
// by src/Makevars (the ordinary flags, -mavx2 -mfma, and AVX-512), each copy
// with SPDEST_LA_ISA set and with the Eigen namespace renamed to Eigen_<isa>,
// so that the copies of every Eigen template are distinct symbols and the
// linker never mixes them. Only Eigen and plain arrays are used here (no std
// containers, no Rcpp): nothing in this file is an inline function shared with
// the rest of the package.
#ifndef SPDEST_LA_ISA
#error "SPDEST_LA_ISA must be defined (generic, avx2 or avx512)"
#endif
#define EIGEN_DONT_PARALLELIZE
#define SPDEST_CAT2(a, b) a##b
#define SPDEST_CAT(a, b) SPDEST_CAT2(a, b)
#define Eigen SPDEST_CAT(Eigen_, SPDEST_LA_ISA)
#include <Eigen/Dense>
#include <cmath>
#include "la_dispatch.h"

namespace {

using MatD = Eigen::MatrixXd;
using MatF = Eigen::MatrixXf;
using MapCD = Eigen::Map<const MatD>;

struct Work {
  int q, k;
  MatD P1, P2, Ln, Rn, Mk, Gram;
  Eigen::SelfAdjointEigenSolver<MatD> eig;
  Work(int q_, int k_, bool moments) : q(q_), k(k_), Mk(k_, k_), Gram(k_, k_), eig(k_) {
    if (moments) {
      P1.resize(q, q);
      P2.resize(q, q);
      Ln.resize(q, q);
      Rn.resize(q, q);
    }
  }
};

void* work_new(int q, int k, int moments) { return new Work(q, k, moments != 0); }
void work_free(void* work) { delete static_cast<Work*>(work); }

void set_left(void* work, const double* H, double scale) {
  Work& w = *static_cast<Work*>(work);
  MapCD M(H, w.q, w.q);
  w.Ln = (0.5 / scale) * (M + M.transpose());
}

void moments(void* work, const double* H, double scale, double* v) {
  Work& w = *static_cast<Work*>(work);
  const int q = w.q;
  MapCD M(H, q, q);
  w.Rn = (0.5 / scale) * (M + M.transpose());
  w.P1.noalias() = w.Ln * w.Rn;
  w.P2.noalias() = w.P1 * w.P1;
  v[0] = v[1] = v[2] = v[3] = 0;
  for (int i = 0; i < q; ++i) v[0] += w.P1(i, i);
  for (int c = 0; c < q; ++c) {
    for (int i = 0; i < q; ++i) {
      const double p = w.P1(i, c), pt = w.P1(c, i), h = w.P2(i, c), ht = w.P2(c, i);
      v[1] += p * pt;
      v[2] += h * pt;
      v[3] += h * ht;
    }
  }
}

double trace1(const double* Hl, const double* Hr, int q) {
  MapCD L(Hl, q, q), R(Hr, q, q);
  double direct = 0, crossed = 0;
  for (int c = 0; c < q; ++c) {
    direct += L.col(c).dot(R.col(c));
    crossed += L.col(c).dot(R.row(c));
  }
  return 0.5 * (direct + crossed);
}

// Singular values of Mk as the square roots of the eigenvalues of Mk' Mk
// (descending). Eigen 3.4.0's BDCSVD gave invalid values on rank-deficient
// matrices (k = q = 298 on Visium-B). Only the lower triangle of the Gram
// matrix is formed (a symmetric rank-k update) and read by the solver.
int gram_values(Work& w, double* sv) {
  if (!w.Mk.allFinite()) return 1;
  w.Gram.triangularView<Eigen::Lower>().setZero();
  w.Gram.selfadjointView<Eigen::Lower>().rankUpdate(w.Mk.transpose());
  w.eig.compute(w.Gram, Eigen::EigenvaluesOnly);
  if (w.eig.info() != Eigen::Success) return 2;
  const int k = w.k;
  for (int i = 0; i < k; ++i) {
    const double e = w.eig.eigenvalues()[k - 1 - i];
    sv[i] = e > 0 ? std::sqrt(e) : 0.0;
  }
  return 0;
}

int spectrum_tn(void* work, const double* Gl, const double* Gr, int q, double* sv) {
  Work& w = *static_cast<Work*>(work);
  MapCD L(Gl, q, w.k), R(Gr, q, w.k);
  w.Mk.noalias() = L.transpose() * R;
  return gram_values(w, sv);
}

int spectrum_nt(void* work, const double* Ri, const double* Rj, double* sv) {
  Work& w = *static_cast<Work*>(work);
  MapCD Mi(Ri, w.k, w.k), Mj(Rj, w.k, w.k);
  w.Mk.noalias() = Mi * Mj.transpose();
  return gram_values(w, sv);
}

void dgemm_nn(const double* A, int M, int K, const double* B, int N, double* Out) {
  MapCD a(A, M, K), b(B, K, N);
  Eigen::Map<MatD> o(Out, M, N);
  o.noalias() = a * b;
}

void sgemm_tn(const float* A, long K, int M, const float* B, int N, float* Out) {
  Eigen::Map<const MatF> a(A, K, M), b(B, K, N);
  Eigen::Map<MatF> o(Out, M, N);
  o.noalias() = a.transpose() * b;
}

// With A = F' D^{-1} F, M = I + A, C = F' D^{-1} X and g = F' D^{-1} e:
//   F' W F = M^{-1} A,  F' W X = M^{-1} C,  F' W e = M^{-1} g,
//   X' W X = X' D^{-1} X - C' M^{-1} C,  X' W e = X' D^{-1} e - C' M^{-1} g,
// so a = M^{-1} g - B K (X' W e) and H = M^{-1} A - B K B', B = M^{-1} C.
int score_state(const double* T0, int n, int q, const double* variance, const double* error,
                double scale, const double* X, int p, const double* Vp, double* a, double* H) {
  using VecD = Eigen::VectorXd;
  MapCD T(T0, n, q);
  Eigen::Map<const VecD> var(variance, n), err(error, n);
  const double root = std::sqrt(scale);
  const VecD dinv = var.cwiseInverse();
  const VecD de = dinv.cwiseProduct(err);

  MatD A = MatD::Zero(q, q);
  {
    const MatD S = dinv.cwiseSqrt().asDiagonal() * T;
    A.selfadjointView<Eigen::Lower>().rankUpdate(S.transpose(), scale);
  }
  A.triangularView<Eigen::StrictlyUpper>() = A.transpose();
  if (!A.allFinite()) return 1;
  MatD M = A;
  M.diagonal().array() += 1.0;
  Eigen::LLT<MatD> llt(M);
  if (llt.info() != Eigen::Success) return 1;

  VecD av = llt.solve(root * (T.transpose() * de));  // F' W e
  MatD B, KBt;
  if (p > 0) {
    MapCD Xm(X, n, p);
    const MatD DX = dinv.asDiagonal() * Xm;
    const MatD C = root * (T.transpose() * DX);
    B = llt.solve(C);  // F' W X
    const VecD xwe = Xm.transpose() * de - root * (B.transpose() * (T.transpose() * de));
    MatD K;
    if (Vp) {
      K = MapCD(Vp, p, p);
    } else {
      MatD xwx = Xm.transpose() * DX - C.transpose() * B;
      xwx = (0.5 * (xwx + xwx.transpose())).eval();
      Eigen::SelfAdjointEigenSolver<MatD> es(xwx);
      if (es.info() != Eigen::Success) return 2;
      VecD inv = es.eigenvalues();
      for (int i = 0; i < p; ++i) inv[i] = std::abs(inv[i]) > 5e-16 ? 1.0 / inv[i] : 0.0;
      K = es.eigenvectors() * inv.asDiagonal() * es.eigenvectors().transpose();
    }
    KBt = K * B.transpose();
    if (!KBt.allFinite()) return 2;
    av.noalias() -= B * (K * xwe);
  }
  if (!av.allFinite()) return 3;
  Eigen::Map<VecD>(a, q) = av;
  if (!H) return 0;

  MatD Hm = llt.solve(A);  // F' W F
  if (p > 0) Hm.noalias() -= B * KBt;
  Eigen::Map<MatD> Ho(H, q, q);
  Ho = 0.5 * (Hm + Hm.transpose());
  return Ho.allFinite() ? 0 : 3;
}

}  // namespace

#define SPDEST_LA_TABLE SPDEST_CAT(mgcvst_la_, SPDEST_LA_ISA)
#define SPDEST_STR2(x) #x
#define SPDEST_STR(x) SPDEST_STR2(x)
extern const mgcvst::LaFns SPDEST_LA_TABLE;
const mgcvst::LaFns SPDEST_LA_TABLE = {work_new, work_free, set_left, moments, trace1, spectrum_tn,
                                       spectrum_nt, dgemm_nn, sgemm_tn, score_state,
                                       SPDEST_STR(SPDEST_LA_ISA)};
