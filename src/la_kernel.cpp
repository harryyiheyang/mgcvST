// Dense linear-algebra kernels of the pair tests, compiled once per instruction
// set and chosen at run time (la_dispatch.h).
//
// Eigen picks its SIMD width at preprocessing time from __AVX2__ / __FMA__ /
// __AVX512F__, which R's default package flags never define, so Eigen products
// in ordinary translation units run on SSE2. This file is compiled three times
// by src/Makevars (the ordinary flags, -mavx2 -mfma, and AVX-512), each copy
// with MGCVST_LA_ISA set and with the Eigen namespace renamed to Eigen_<isa>,
// so that the copies of every Eigen template are distinct symbols and the
// linker never mixes them. Only Eigen and plain arrays are used here (no std
// containers, no Rcpp): nothing in this file is an inline function shared with
// the rest of the package.
#ifndef MGCVST_LA_ISA
#error "MGCVST_LA_ISA must be defined (generic, avx2 or avx512)"
#endif
#define EIGEN_DONT_PARALLELIZE
#define MGCVST_CAT2(a, b) a##b
#define MGCVST_CAT(a, b) MGCVST_CAT2(a, b)
#define Eigen MGCVST_CAT(Eigen_, MGCVST_LA_ISA)
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

}  // namespace

#define MGCVST_LA_TABLE MGCVST_CAT(mgcvst_la_, MGCVST_LA_ISA)
#define MGCVST_STR2(x) #x
#define MGCVST_STR(x) MGCVST_STR2(x)
extern const mgcvst::LaFns MGCVST_LA_TABLE;
const mgcvst::LaFns MGCVST_LA_TABLE = {work_new, work_free, set_left, moments, spectrum_tn,
                                       spectrum_nt, dgemm_nn, sgemm_tn,
                                       MGCVST_STR(MGCVST_LA_ISA)};
