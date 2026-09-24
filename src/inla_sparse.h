#ifndef MGCVST_INLA_SPARSE_H
#define MGCVST_INLA_SPARSE_H

#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mgcvst_sparse {

using SpMat = Eigen::SparseMatrix<double>;
using Mat = Eigen::MatrixXd;
using Vec = Eigen::VectorXd;
using QFactor = Eigen::SimplicialLLT<SpMat, Eigen::Lower,
                                    Eigen::AMDOrdering<int> >;
using HFactor = Eigen::SimplicialLDLT<SpMat, Eigen::Lower,
                                     Eigen::AMDOrdering<int> >;

template <class Factor>
inline Mat constrained_solve(const Factor& factor, const Vec& g, const Mat& rhs,
                             const Vec& inverse_g, double denominator) {
  Mat answer = factor.solve(rhs);
  if (factor.info() != Eigen::Success) {
    throw std::runtime_error("A constrained sparse solve failed.");
  }
  const Eigen::RowVectorXd multiplier =
    (g.transpose() * answer) / denominator;
  answer.noalias() -= inverse_g * multiplier;
  return answer;
}

template <class Factor>
inline Vec constrained_solve(const Factor& factor, const Vec& g, const Vec& rhs,
                             const Vec& inverse_g, double denominator) {
  Mat answer = constrained_solve(factor, g, Mat(rhs), inverse_g, denominator);
  return answer.col(0);
}

inline Mat apply_B(const QFactor& factor, const Mat& value) {
  return factor.permutationPinv() * factor.matrixU().solve(value);
}

inline Mat apply_Bt(const QFactor& factor, const Mat& value) {
  return factor.matrixL().solve(factor.permutationP() * value);
}

inline Vec apply_Bt(const QFactor& factor, const Vec& value) {
  Mat answer = factor.matrixL().solve(factor.permutationP() * Mat(value));
  return answer.col(0);
}

inline Mat project_coordinates(const Mat& value, const Vec& u) {
  return value - u * (u.transpose() * value);
}

inline Vec project_coordinates(const Vec& value, const Vec& u) {
  return value - u * u.dot(value);
}

inline Mat stored_ldlt_solve(const SpMat& L, const SpMat& U, const Vec& D,
                             const Eigen::VectorXi& indices, const Mat& rhs) {
  Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> P(indices.size());
  P.indices() = indices;
  Mat answer = P * rhs;
  L.triangularView<Eigen::UnitLower>().solveInPlace(answer);
  answer.array().colwise() /= D.array();
  U.triangularView<Eigen::UnitUpper>().solveInPlace(answer);
  return P.inverse() * answer;
}

inline Mat stored_constrained_solve(const SpMat& L, const SpMat& U, const Vec& D,
                                    const Eigen::VectorXi& indices, const Vec& g,
                                    const Mat& rhs, const Vec& inverse_g,
                                    double denominator) {
  Mat answer = stored_ldlt_solve(L, U, D, indices, rhs);
  const Eigen::RowVectorXd multiplier =
    (g.transpose() * answer) / denominator;
  answer.noalias() -= inverse_g * multiplier;
  return answer;
}

struct SparsePrepared {
  SpMat Q;
  Vec constraint;
  QFactor qfactor;
  Vec qinv_g;
  Vec coordinate_constraint;
  double qden = NA_REAL;
  bool valid = false;

  SparsePrepared(const SpMat& Q_, const Vec& constraint_) :
    Q(Q_), constraint(constraint_) {
    qfactor.compute(Q);
    if (qfactor.info() != Eigen::Success) return;
    qinv_g = qfactor.solve(constraint);
    if (qfactor.info() != Eigen::Success) return;
    qden = constraint.dot(qinv_g);
    coordinate_constraint = apply_Bt(qfactor, Vec(constraint));
    const double unorm = coordinate_constraint.norm();
    if (!std::isfinite(qden) || qden <= 0 ||
        !std::isfinite(unorm) || unorm <= 0) return;
    coordinate_constraint /= unorm;
    valid = true;
  }
};

inline const SparsePrepared* prepared_cache(SEXP prepared, int m,
                                            const Vec& constraint) {
  if (TYPEOF(prepared) != EXTPTRSXP || R_ExternalPtrAddr(prepared) == NULL) {
    Rcpp::stop("prepared must be a valid sparse INLA cache.");
  }
  Rcpp::XPtr<SparsePrepared> pointer(prepared);
  const SparsePrepared* cache = pointer.get();
  if (!cache->valid || cache->Q.rows() != m ||
      !cache->constraint.isApprox(constraint, 0.0)) {
    Rcpp::stop("prepared is not aligned with Q and constraint.");
  }
  return cache;
}

struct ReconstructionUnit {
  Vec a;
  SpMat K;
  Mat U;
  Mat Vp;
  SpMat H_L;
  Vec H_D;
  Eigen::VectorXi H_perm;
  Vec hinv_g;
  double hden;
  double tau;
  double constraint_norm;
};

inline std::vector<ReconstructionUnit> parse_units(const Rcpp::List& units, int m) {
  const int features = units.size();
  std::vector<ReconstructionUnit> parsed;
  parsed.reserve(features);
  for (int f = 0; f < features; ++f) {
    const Rcpp::List unit = units[f];
    ReconstructionUnit value;
    value.a = Rcpp::as<Vec>(unit["a"]);
    value.K = Rcpp::as<SpMat>(unit["K"]);
    value.U = Rcpp::as<Mat>(unit["U"]);
    value.Vp = Rcpp::as<Mat>(unit["expected_vp"]);
    value.H_L = Rcpp::as<SpMat>(unit["H_L"]);
    value.H_D = Rcpp::as<Vec>(unit["H_D"]);
    value.H_perm = Rcpp::as<Eigen::VectorXi>(unit["H_perm"]);
    value.hinv_g = Rcpp::as<Vec>(unit["hinv_g"]);
    value.hden = Rcpp::as<double>(unit["hden"]);
    value.tau = Rcpp::as<double>(unit["tau"]);
    value.constraint_norm = Rcpp::as<double>(unit["constraint_diagnostic"]);
    if (value.K.rows() != m || value.K.cols() != m ||
        value.H_L.rows() != m || value.H_L.cols() != m ||
        value.H_D.size() != m || value.H_perm.size() != m ||
        value.hinv_g.size() != m || value.a.size() != m ||
        !std::isfinite(value.tau) || value.tau <= 0 ||
        !std::isfinite(value.hden) || value.hden <= 0) {
      Rcpp::stop("A serialized sparse INLA unit is malformed.");
    }
    parsed.push_back(std::move(value));
  }
  return parsed;
}

// Reduced curvature C' (K - K S K - U Vp U') C / tau for the physical basis C,
// built in 32-column blocks and symmetrized.
template <class Basis>
inline Mat reduced_curvature(const ReconstructionUnit& unit, const Basis& basis,
                             const Vec& constraint) {
  const int r = basis.cols();
  const SpMat H_U = unit.H_L.transpose();
  Mat M = Mat::Zero(r, r);
  for (int first = 0; first < r; first += 32) {
    const int width = std::min(32, r - first);
    Mat V = basis.middleCols(first, width) / std::sqrt(unit.tau);
    Mat Y = unit.K * V;
    Mat SY = stored_constrained_solve(
      unit.H_L, H_U, unit.H_D, unit.H_perm, constraint, Y,
      unit.hinv_g, unit.hden
    );
    Y.noalias() -= unit.K * SY;
    if (unit.U.cols()) {
      Y.noalias() -= unit.U * (unit.Vp * (unit.U.transpose() * V));
    }
    M.middleCols(first, width) = basis.transpose() * Y /
      std::sqrt(unit.tau);
  }
  return 0.5 * (M + M.transpose());
}

// Compact per-feature working model: eta = offset + X c + A b.
struct FeatureModel {
  const double* b;
  const double* c;
  const double* offset;
  int family;  // 0 gaussian, 1 poisson, 2 negative binomial
  double size;
  double dispersion;
};

inline void recover_working(const SpMat& A, const SpMat& X, const FeatureModel& f,
                            Vec& eta, Vec& mu, Vec& V) {
  const int n = A.rows();
  eta = A * Eigen::Map<const Vec>(f.b, A.cols());
  if (X.cols()) eta.noalias() += X * Eigen::Map<const Vec>(f.c, X.cols());
  eta += Eigen::Map<const Vec>(f.offset, n);
  if (!eta.allFinite()) throw std::runtime_error("The recovered linear predictor is non-finite.");
  if (f.family == 0) {
    if (!std::isfinite(f.dispersion) || f.dispersion <= 0) {
      throw std::runtime_error("The Gaussian dispersion is invalid.");
    }
    mu = eta;
    V = Vec::Constant(n, f.dispersion);
    return;
  }
  mu = eta.array().exp().matrix();
  if (!mu.allFinite() || (mu.array() <= 0).any()) {
    throw std::runtime_error("The recovered conditional mode produced invalid means.");
  }
  V = mu.cwiseInverse();
  if (f.family == 2) {
    if (!std::isfinite(f.size) || f.size <= 0) {
      throw std::runtime_error("The negative-binomial size is invalid.");
    }
    V.array() += 1.0 / f.size;
  } else if (f.family != 1) {
    throw std::runtime_error("Unknown feature family code.");
  }
  if (!V.allFinite() || (V.array() <= 0).any()) {
    throw std::runtime_error("The recovered working variance is invalid.");
  }
}

// Expected-curvature reconstruction unit from the working variance V alone.
inline void curvature_unit(const SpMat& A, const SpMat& X, const SpMat& Q,
                           const Vec& g, const Vec& V, double tau,
                           const double* penalty, ReconstructionUnit& unit) {
  const int m = A.cols();
  const int px = X.cols();
  const Vec W = V.cwiseInverse();
  SpMat WA = A;
  for (int col = 0; col < WA.outerSize(); ++col) {
    for (SpMat::InnerIterator it(WA, col); it; ++it) it.valueRef() *= W[it.row()];
  }
  SpMat K = A.transpose() * WA;
  K.makeCompressed();
  HFactor hfactor;
  hfactor.compute(SpMat(K + tau * Q));
  if (hfactor.info() != Eigen::Success) {
    throw std::runtime_error("The feature sparse expected-curvature factorization failed.");
  }
  Vec hinv_g = hfactor.solve(g);
  const double hden = g.dot(hinv_g);
  if (!std::isfinite(hden) || hden <= 0) {
    throw std::runtime_error("The feature constraint has a non-positive H-inverse norm.");
  }
  Mat U(m, 0), Vp(0, 0);
  if (px) {
    SpMat WX = X;
    for (int col = 0; col < WX.outerSize(); ++col) {
      for (SpMat::InnerIterator it(WX, col); it; ++it) it.valueRef() *= W[it.row()];
    }
    const Mat L = Mat(SpMat(A.transpose() * WX));
    const Mat SL = constrained_solve(hfactor, g, L, hinv_g, hden);
    U = L - K * SL;
    Mat J = Mat(SpMat(X.transpose() * WX)) - L.transpose() * SL;
    for (int k = 0; k < px; ++k) J(k, k) += penalty[k];
    Eigen::LDLT<Mat> jfactor(0.5 * (J + J.transpose()));
    if (jfactor.info() != Eigen::Success || !jfactor.isPositive()) {
      throw std::runtime_error("The expected nuisance information is not positive definite.");
    }
    Vp = jfactor.solve(Mat::Identity(px, px));
  }
  unit.K = std::move(K);
  unit.U = std::move(U);
  unit.Vp = std::move(Vp);
  unit.H_L = SpMat(hfactor.matrixL());
  unit.H_D = hfactor.vectorD();
  unit.H_perm = hfactor.permutationP().indices();
  unit.hinv_g = std::move(hinv_g);
  unit.hden = hden;
  unit.tau = tau;
}

}  // namespace mgcvst_sparse

#endif
