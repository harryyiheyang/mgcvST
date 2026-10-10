#ifndef MGCVST_PCA_REDUCE_H
#define MGCVST_PCA_REDUCE_H

// Reduction of one symmetric q x q score-curvature matrix H to what the
// PCAlearning route keeps of a feature:
//   * its weighted-vech vector h (upper triangle, column-major, entry (i, j),
//     i <= j, at j (j + 1) / 2 + i, off-diagonal weight sqrt(2)), whose inner
//     products are Frobenius inner products;
//   * the basis coefficients c = B' h and ||h||^2;
//   * the scale max|H| and, for a shared basis V (q x k), the upper Cholesky
//     factor R of V' (H / scale) V, from which the pair kernel obtains the k
//     leading singular values of the pair spectrum as svd(R_i R_j').

#include <RcppEigen.h>
#include <cmath>
#include <limits>
#include <vector>

namespace mgcvst_pca {

using Mat = Eigen::MatrixXd;
using Vec = Eigen::VectorXd;
using Index = Eigen::Index;

// Symmetric q x q matrix from a full (q^2) or weighted-vech (q (q + 1) / 2)
// column.
inline Mat unpack_symmetric(const double* x, int q, bool full) {
  Mat M(q, q);
  if (full) {
    Eigen::Map<const Mat> F(x, q, q);
    M = 0.5 * (F + F.transpose());
    return M;
  }
  const double w = 1 / std::sqrt(2.0);
  Index p = 0;
  for (int j = 0; j < q; ++j) {
    for (int i = 0; i < j; ++i) M(i, j) = M(j, i) = w * x[p++];
    M(j, j) = x[p++];
  }
  return M;
}

inline Mat unpack_float_vech(const float* x, int q) {
  const Index L = static_cast<Index>(q) * (q + 1) / 2;
  Vec h(L);
  for (Index p = 0; p < L; ++p) h[p] = x[p];
  return unpack_symmetric(h.data(), q, false);
}

inline void pack_vech(const Mat& M, int q, Vec& h) {
  const double w = std::sqrt(2.0);
  Index p = 0;
  for (int j = 0; j < q; ++j) {
    for (int i = 0; i < j; ++i) h[p++] = w * M(i, j);
    h[p++] = M(j, j);
  }
}

// max |H|, NA for a non-finite matrix.
inline double matrix_scale(const Mat& M) {
  return M.allFinite() ? M.cwiseAbs().maxCoeff() : NA_REAL;
}

// R (k x k, upper) with R'R = V' (sym(M) / scale) V. A small diagonal jitter is
// added when the compression is numerically singular. Returns false when the
// factor cannot be formed.
inline bool project_factor(const Mat& M, double scale, const Mat& V, Mat& R) {
  Mat P = V.transpose() * ((0.5 / scale) * (M + M.transpose())) * V;
  P = 0.5 * (P + P.transpose());
  if (!P.allFinite()) return false;
  const double top = P.diagonal().maxCoeff();
  const double jitter[3] = {0.0, 1e-10, 1e-8};
  for (int attempt = 0; attempt < 3; ++attempt) {
    Mat Q = P;
    if (attempt) Q.diagonal().array() += jitter[attempt] * top;
    Eigen::LLT<Mat> llt(Q);
    if (llt.info() == Eigen::Success) {
      R = llt.matrixU();
      return true;
    }
  }
  return false;
}

}  // namespace mgcvst_pca

#endif
