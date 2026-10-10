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
#include <cstring>
#include <limits>
#include <string>
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

// Reduction of one symmetric score-curvature matrix M (q x q) to what the
// PCAlearning route keeps: the scale max|M|, ||vech_w(M)||^2, the coefficients
// c = B' vech_w(M) (empty when the basis has no columns), the float32 packed
// weighted vech when `pack` is set, and the vectorized factor R of the shared
// basis V (empty when V has no columns). `h` is a work vector of length
// q (q + 1) / 2. Returns an empty string, or the reason why the feature is
// unusable (then every output is missing and `packed` is empty).
inline std::string reduce_feature(const Mat& M, int q,
                                  const Eigen::Map<const Mat>& PB,
                                  const Eigen::Map<const Mat>& V, bool pack,
                                  Vec& h, double& scale, double& fro2,
                                  Eigen::Ref<Vec> c, std::vector<float>& packed,
                                  Eigen::Ref<Vec> Rcol) {
  const Index L = static_cast<Index>(q) * (q + 1) / 2;
  auto fail = [&](const char* why) {
    fro2 = NA_REAL;
    c.setConstant(NA_REAL);
    Rcol.setConstant(NA_REAL);
    std::vector<float>().swap(packed);
    return std::string(why);
  };
  scale = matrix_scale(M);
  if (!(std::isfinite(scale) && scale > 0)) {
    return fail("The curvature matrix is zero or not finite.");
  }
  pack_vech(M, q, h);
  fro2 = h.squaredNorm();
  if (c.size()) c = PB.transpose() * h;
  if (pack) {
    packed.resize(L);
    for (Index p = 0; p < L; ++p) packed[p] = static_cast<float>(h[p]);
  }
  if (V.cols()) {
    Mat R;
    if (!project_factor(M, scale, V, R)) {
      return fail("The shared-basis compression is not positive definite.");
    }
    Rcol = Eigen::Map<const Vec>(R.data(), V.cols() * V.cols());
  }
  return std::string();
}

// Float32 vectors of a batch as raw vectors for R (4 bytes per entry).
inline Rcpp::List packed_to_raw(std::vector<std::vector<float> >& packed, Index L) {
  Rcpp::List out(packed.size());
  for (size_t g = 0; g < packed.size(); ++g) {
    if (packed[g].empty()) continue;
    Rcpp::RawVector x(4 * L);
    std::memcpy(RAW(x), packed[g].data(), 4 * L);
    out[g] = x;
    std::vector<float>().swap(packed[g]);
  }
  return out;
}

}  // namespace mgcvst_pca

#endif
