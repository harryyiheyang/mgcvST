#define EIGEN_DONT_PARALLELIZE
#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <algorithm>
#include <cstring>
#include <cmath>
#include <stdexcept>
#include <string>
#include <memory>
#include <vector>
#include <utility>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

namespace {

using SpMat = Eigen::SparseMatrix<double>;
using Mat = Eigen::MatrixXd;
using Vec = Eigen::VectorXd;
using QFactor = Eigen::SimplicialLLT<SpMat, Eigen::Lower,
                                    Eigen::AMDOrdering<int> >;
using HFactor = Eigen::SimplicialLDLT<SpMat, Eigen::Lower,
                                     Eigen::AMDOrdering<int> >;

template <class Factor>
Mat constrained_solve(const Factor& factor, const Vec& g, const Mat& rhs,
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
Vec constrained_solve(const Factor& factor, const Vec& g, const Vec& rhs,
                      const Vec& inverse_g, double denominator) {
  Mat answer = constrained_solve(factor, g, Mat(rhs), inverse_g, denominator);
  return answer.col(0);
}

Mat apply_B(const QFactor& factor, const Mat& value) {
  return factor.permutationPinv() * factor.matrixU().solve(value);
}

Mat apply_Bt(const QFactor& factor, const Mat& value) {
  return factor.matrixL().solve(factor.permutationP() * value);
}

Vec apply_Bt(const QFactor& factor, const Vec& value) {
  Mat answer = factor.matrixL().solve(factor.permutationP() * Mat(value));
  return answer.col(0);
}

Mat project_coordinates(const Mat& value, const Vec& u) {
  return value - u * (u.transpose() * value);
}

Vec project_coordinates(const Vec& value, const Vec& u) {
  return value - u * u.dot(value);
}

Vec trace_powers(const Mat& value) {
  Vec out(4);
  const Mat square = value * value;
  out[0] = value.trace();
  out[1] = value.cwiseProduct(value).sum();
  out[2] = square.cwiseProduct(value.transpose()).sum();
  out[3] = square.cwiseProduct(square.transpose()).sum();
  return out;
}

struct FeatureResult {
  Vec a;
  Mat M;
  Mat Vp;
  Vec nuisance_score;
  SpMat K;
  Mat U;
  SpMat H_L;
  Vec H_D;
  Eigen::VectorXi H_perm;
  Vec hinv_g;
  Vec moments;
  double tau = NA_REAL;
  double hden = NA_REAL;
  double statistic = NA_REAL;
  double constraint_norm = NA_REAL;
  std::string error;
};

Mat stored_ldlt_solve(const SpMat& L, const SpMat& U, const Vec& D,
                      const Eigen::VectorXi& indices, const Mat& rhs) {
  Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> P(indices.size());
  P.indices() = indices;
  Mat answer = P * rhs;
  L.triangularView<Eigen::UnitLower>().solveInPlace(answer);
  answer.array().colwise() /= D.array();
  U.triangularView<Eigen::UnitUpper>().solveInPlace(answer);
  return P.inverse() * answer;
}

Mat stored_constrained_solve(const SpMat& L, const SpMat& U, const Vec& D,
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

std::vector<ReconstructionUnit> parse_units(const Rcpp::List& units, int m) {
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
Mat reduced_curvature(const ReconstructionUnit& unit,
                      const Eigen::Map<Eigen::MatrixXd>& basis,
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

} // namespace

// [[Rcpp::export]]
SEXP mgcvst_inla_sparse_prepare_cpp(
    const Eigen::MappedSparseMatrix<double>& Q_map,
    const Eigen::Map<Eigen::VectorXd> constraint) {
  const SpMat Q = Q_map;
  if (Q.rows() != Q.cols() || Q.rows() < 2) {
    Rcpp::stop("Q must be square with at least two coefficients.");
  }
  if (constraint.size() != Q.rows() || !constraint.allFinite() ||
      constraint.norm() == 0) {
    Rcpp::stop("constraint must be one finite nonzero value per SPDE coefficient.");
  }
  std::unique_ptr<SparsePrepared> cache(new SparsePrepared(Q, constraint));
  if (!cache->valid) {
    Rcpp::stop("The sparse SPDE precision factorization or constraint failed.");
  }
  Rcpp::XPtr<SparsePrepared> pointer(cache.release(), true);
  pointer.attr("class") = "mgcvst_inla_sparse_prepared";
  return pointer;
}

// [[Rcpp::export]]
bool mgcvst_inla_sparse_prepared_valid_cpp(SEXP prepared) {
  return TYPEOF(prepared) == EXTPTRSXP && R_ExternalPtrAddr(prepared) != NULL;
}

// Constrained observation-kernel directions.  The nonzero spectrum of
// A Q_g^{-1} A' is that of Z' B' A' A B Z, where B' Q B = I and Z spans the
// complement of the native constraint in B coordinates.
// [[Rcpp::export]]
Rcpp::List mgcvst_inla_sparse_observation_basis_cpp(
    const Eigen::MappedSparseMatrix<double>& A_map,
    const Eigen::Map<Eigen::VectorXd> constraint,
    double coverage = 0.995,
    bool full_rank = false,
    SEXP prepared = R_NilValue) {
  const SpMat A = A_map;
  const int m = A.cols();
  if (m < 2) Rcpp::stop("The constrained SPDE block must contain at least two coefficients.");
  if (constraint.size() != m) Rcpp::stop("constraint must align with A.");
  if (!std::isfinite(coverage) || coverage <= 0 || coverage > 1) {
    Rcpp::stop("coverage must be in (0, 1].");
  }
  if (!mgcvst_inla_sparse_prepared_valid_cpp(prepared)) {
    Rcpp::stop("prepared must be a valid sparse INLA cache.");
  }
  Rcpp::XPtr<SparsePrepared> pointer(prepared);
  const SparsePrepared* cache = pointer.get();
  if (cache->Q.rows() != m || cache->constraint.size() != m ||
      !cache->constraint.isApprox(constraint, 0.0)) {
    Rcpp::stop("prepared is not aligned with A and constraint.");
  }

  const Vec& u = cache->coordinate_constraint;
  Vec h = u;
  h[0] += u[0] >= 0 ? 1.0 : -1.0;
  h.normalize();
  Mat reflector = Mat::Identity(m, m) - 2.0 * h * h.transpose();
  Mat Z = reflector.rightCols(m - 1);
  Mat BZ = apply_B(cache->qfactor, Z);
  SpMat AtA = A.transpose() * A;
  AtA.makeCompressed();
  Mat gram = BZ.transpose() * AtA * BZ;
  Eigen::SelfAdjointEigenSolver<Mat> eig(gram);
  if (eig.info() != Eigen::Success) {
    Rcpp::stop("The constrained observation-kernel eigendecomposition failed.");
  }
  const Vec values = eig.eigenvalues().reverse();
  const Mat vectors = eig.eigenvectors().rowwise().reverse();
  const double total = values.sum();
  int r = m - 1;
  if (!full_rank && coverage < 1) {
    double cumulative = 0;
    for (int j = 0; j < m - 1; ++j) {
      cumulative += values[j];
      if (cumulative / total >= coverage) {
        r = j + 1;
        break;
      }
    }
  }
  const Mat R = Z * vectors.leftCols(r);
  const Mat C = BZ * vectors.leftCols(r);
  const double kept = values.head(r).sum() / total;
  return Rcpp::List::create(
    Rcpp::Named("coordinate") = R,
    Rcpp::Named("basis") = C,
    Rcpp::Named("values") = values,
    Rcpp::Named("rank") = r,
    Rcpp::Named("coverage") = coverage,
    Rcpp::Named("kept") = kept,
    Rcpp::Named("tail") = 1.0 - kept
  );
}

// Exact expected-curvature score states for a single constrained SPDE block.
// Q remains sparse. The redundant whitened coordinate system has dimension m;
// its rank and WGCNA normalization are m - 1 after removing the constraint.
// [[Rcpp::export]]
Rcpp::List mgcvst_inla_sparse_batch_cpp(
    const Eigen::MappedSparseMatrix<double>& A_map,
    const Eigen::MappedSparseMatrix<double>& Q_map,
    const Eigen::Map<Eigen::VectorXd> constraint,
    const Eigen::Map<Eigen::MatrixXd> X,
    const Eigen::Map<Eigen::MatrixXd> E,
    const Eigen::Map<Eigen::MatrixXd> D,
    const Eigen::Map<Eigen::VectorXd> tau,
    int threads = 1,
    bool score_only = false,
    bool null_target = false,
    int block_size = 32,
    SEXP prepared = R_NilValue,
    bool unit_only = false,
    Rcpp::Nullable<Rcpp::NumericMatrix> nuisance_precision = R_NilValue) {
  const SpMat A = A_map;
  const SpMat Q = Q_map;
  const int n = A.rows();
  const int m = A.cols();
  const int px = X.cols();
  const int features = E.cols();
  Mat nuisance_penalty = Mat::Zero(px, features);
  if (nuisance_precision.isNotNull()) {
    nuisance_penalty = Rcpp::as<Mat>(nuisance_precision.get());
  }

  if (Q.rows() != m || Q.cols() != m) Rcpp::stop("Q must be square and aligned with A.");
  if (m < 2) Rcpp::stop("The constrained SPDE block must contain at least two coefficients.");
  if (constraint.size() != m || !constraint.allFinite() || constraint.norm() == 0) {
    Rcpp::stop("constraint must be one finite nonzero value per SPDE coefficient.");
  }
  if (X.rows() != n || E.rows() != n || D.rows() != n || D.cols() != features) {
    Rcpp::stop("X, E, and D must have rows aligned with A and matching feature columns.");
  }
  if (tau.size() != features || !tau.allFinite() || (tau.array() <= 0).any()) {
    Rcpp::stop("tau must contain one positive finite value per feature.");
  }
  if (nuisance_penalty.rows() != px ||
      nuisance_penalty.cols() != features ||
      !nuisance_penalty.allFinite() ||
      (nuisance_penalty.array() < 0).any()) {
    Rcpp::stop(
      "nuisance_precision must be a finite non-negative ncol(X)-by-feature matrix."
    );
  }
  if (threads < 1) Rcpp::stop("threads must be positive.");
  if (block_size < 1) Rcpp::stop("block_size must be positive.");
  if (unit_only && (score_only || null_target)) {
    Rcpp::stop("unit_only is available only for non-null pair states.");
  }
#ifndef _OPENMP
  if (threads > 1) Rcpp::stop("mgcvST was compiled without OpenMP; use threads = 1.");
#endif

  std::unique_ptr<SparsePrepared> local_cache;
  const SparsePrepared* cache = NULL;
  if (prepared == R_NilValue) {
    local_cache.reset(new SparsePrepared(Q, constraint));
    cache = local_cache.get();
  } else {
    if (!mgcvst_inla_sparse_prepared_valid_cpp(prepared)) {
      Rcpp::stop("prepared must be a valid, unserialized sparse INLA cache.");
    }
    Rcpp::XPtr<SparsePrepared> pointer(prepared);
    cache = pointer.get();
    if (cache->Q.rows() != m || cache->constraint.size() != m ||
        !cache->constraint.isApprox(constraint, 0.0)) {
      Rcpp::stop("prepared is not aligned with Q and constraint.");
    }
  }
  if (!cache->valid) {
    Rcpp::stop("The sparse SPDE precision factorization or constraint failed.");
  }

  const QFactor& qfactor = cache->qfactor;
  const Vec& qinv_g = cache->qinv_g;
  const Vec& coordinate_constraint = cache->coordinate_constraint;
  const double qden = cache->qden;

  std::vector<FeatureResult> result(features);
  int failed = 0;

#ifdef _OPENMP
#pragma omp parallel num_threads(threads) reduction(+:failed)
#endif
  {
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
    for (int f = 0; f < features; ++f) {
      try {
        const Vec W = D.col(f).cwiseInverse();
        SpMat WA = A;
        for (int col = 0; col < WA.outerSize(); ++col) {
          for (SpMat::InnerIterator it(WA, col); it; ++it) it.valueRef() *= W[it.row()];
        }
        SpMat K = A.transpose() * WA;
        K.makeCompressed();
        const Mat WX = W.asDiagonal() * X;
        const Mat L = A.transpose() * WX;
        const Vec tvec = A.transpose() * (W.array() * E.col(f).array()).matrix();
        const Vec xte = X.transpose() * (W.array() * E.col(f).array()).matrix();

        Mat U;
        Mat Vp;
        Vec nuisance_score;
        Vec h;
        HFactor hfactor;
        Vec hinv_g;
        double hden = NA_REAL;

        if (null_target) {
          U = L;
          nuisance_score = xte;
          if (px) {
            Mat J = X.transpose() * WX;
            J.diagonal() += nuisance_penalty.col(f);
            Eigen::LDLT<Mat> jfactor(0.5 * (J + J.transpose()));
            if (jfactor.info() != Eigen::Success || !jfactor.isPositive()) {
              throw std::runtime_error("The expected marginal nuisance information is not positive definite.");
            }
            Vp = jfactor.solve(Mat::Identity(px, px));
            h = tvec - U * (Vp * nuisance_score);
          } else {
            Vp.resize(0, 0);
            h = tvec;
          }
        } else {
          SpMat H = K + tau[f] * Q;
          hfactor.compute(H);
          if (hfactor.info() != Eigen::Success) {
            throw std::runtime_error("The feature sparse expected-curvature factorization failed.");
          }
          hinv_g = hfactor.solve(constraint);
          hden = constraint.dot(hinv_g);
          if (!std::isfinite(hden) || hden <= 0) {
            throw std::runtime_error("The feature constraint has a non-positive H-inverse norm.");
          }
          Mat SL = px ? constrained_solve(hfactor, constraint, L, hinv_g, hden) : Mat(m, 0);
          Vec St = constrained_solve(hfactor, constraint, tvec, hinv_g, hden);
          U = px ? L - K * SL : Mat(m, 0);
          nuisance_score = px ? xte - L.transpose() * St : Vec(0);
          Mat J = px ? X.transpose() * WX - L.transpose() * SL : Mat(0, 0);
          if (px) {
            J.diagonal() += nuisance_penalty.col(f);
            Eigen::LDLT<Mat> jfactor(0.5 * (J + J.transpose()));
            if (jfactor.info() != Eigen::Success || !jfactor.isPositive()) {
              throw std::runtime_error("The expected nuisance information is not positive definite.");
            }
            Vp = jfactor.solve(Mat::Identity(px, px));
          } else Vp.resize(0, 0);
          h = tvec - K * St;
          if (px) h.noalias() -= U * (Vp * nuisance_score);
        }

        Vec a = project_coordinates(apply_Bt(qfactor, h), coordinate_constraint) /
          std::sqrt(tau[f]);
        result[f].a = a;
        result[f].statistic = a.squaredNorm();
        result[f].Vp = Vp;
        result[f].nuisance_score = nuisance_score;
        result[f].constraint_norm = std::abs(coordinate_constraint.norm() - 1.0);

        if (unit_only) {
          result[f].K = K;
          result[f].U = U;
          result[f].H_L = SpMat(hfactor.matrixL());
          result[f].H_D = hfactor.vectorD();
          result[f].H_perm = hfactor.permutationP().indices();
          result[f].hinv_g = hinv_g;
          result[f].hden = hden;
          result[f].tau = tau[f];
          continue;
        }

        if (!score_only) {
          Mat M = Mat::Zero(m, m);
          for (int first = 0; first < m; first += block_size) {
            const int width = std::min(block_size, m - first);
            Mat Z = Mat::Zero(m, width);
            Z.block(first, 0, width, width).setIdentity();
            Z = project_coordinates(Z, coordinate_constraint);
            Mat V = apply_B(qfactor, Z) / std::sqrt(tau[f]);
            Mat Y = K * V;
            if (!null_target) {
              Mat SY = constrained_solve(hfactor, constraint, Y, hinv_g, hden);
              Y.noalias() -= K * SY;
            }
            if (px) Y.noalias() -= U * (Vp * (U.transpose() * V));
            M.middleCols(first, width) =
              project_coordinates(apply_Bt(qfactor, Y), coordinate_constraint) /
              std::sqrt(tau[f]);
          }
          result[f].M = 0.5 * (M + M.transpose());
          if (null_target) {
            result[f].moments = trace_powers(result[f].M);
            result[f].M.resize(0, 0);
          }
        }
      } catch (const std::exception& error) {
        result[f].error = error.what();
        failed += 1;
      }
    }
  }

  Rcpp::List features_out(features);
  for (int f = 0; f < features; ++f) {
    if (!result[f].error.empty()) {
      features_out[f] = Rcpp::List::create(Rcpp::Named("error") = result[f].error);
    } else if (unit_only) {
      features_out[f] = Rcpp::List::create(
        Rcpp::Named("a") = result[f].a,
        Rcpp::Named("statistic") = result[f].statistic,
        Rcpp::Named("expected_vp") = result[f].Vp,
        Rcpp::Named("nuisance_score") = result[f].nuisance_score,
        Rcpp::Named("tau") = result[f].tau,
        Rcpp::Named("K") = result[f].K,
        Rcpp::Named("U") = result[f].U,
        Rcpp::Named("H_L") = result[f].H_L,
        Rcpp::Named("H_D") = result[f].H_D,
        Rcpp::Named("H_perm") = result[f].H_perm,
        Rcpp::Named("hinv_g") = result[f].hinv_g,
        Rcpp::Named("hden") = result[f].hden,
        Rcpp::Named("width") = m,
        Rcpp::Named("normalization") = m - 1,
        Rcpp::Named("constraint_diagnostic") = result[f].constraint_norm
      );
    } else {
      features_out[f] = Rcpp::List::create(
        Rcpp::Named("a") = result[f].a,
        Rcpp::Named("M") = (score_only || null_target) ? R_NilValue : Rcpp::wrap(result[f].M),
        Rcpp::Named("statistic") = result[f].statistic,
        Rcpp::Named("moments") = (score_only || !null_target) ? R_NilValue : Rcpp::wrap(result[f].moments),
        Rcpp::Named("expected_vp") = result[f].Vp,
        Rcpp::Named("width") = m,
        Rcpp::Named("normalization") = m - 1,
        Rcpp::Named("constraint_diagnostic") = result[f].constraint_norm
      );
    }
  }
  features_out.attr("failed") = failed;
  features_out.attr("coordinate_width") = m;
  features_out.attr("normalization") = m - 1;
  return features_out;
}

// Serializable reconstruction units for pairwise score curvature.
// [[Rcpp::export]]
Rcpp::List mgcvst_inla_sparse_units_cpp(
    const Eigen::MappedSparseMatrix<double>& A,
    const Eigen::MappedSparseMatrix<double>& Q,
    const Eigen::Map<Eigen::VectorXd> constraint,
    const Eigen::Map<Eigen::MatrixXd> X,
    const Eigen::Map<Eigen::MatrixXd> E,
    const Eigen::Map<Eigen::MatrixXd> D,
    const Eigen::Map<Eigen::VectorXd> tau,
    int threads = 1,
    SEXP prepared = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericMatrix> nuisance_precision = R_NilValue) {
  return mgcvst_inla_sparse_batch_cpp(
    A, Q, constraint, X, E, D, tau, threads, false, false, 32,
    prepared, true, nuisance_precision
  );
}

// Materialize curvature matrices temporarily from serialized units.
// [[Rcpp::export]]
Rcpp::List mgcvst_inla_sparse_materialize_cpp(
    const Rcpp::List& units,
    const Eigen::MappedSparseMatrix<double>& Q_map,
    const Eigen::Map<Eigen::VectorXd> constraint,
    int threads = 1,
    SEXP prepared = R_NilValue,
    int block_size = 32) {
  const SpMat Q = Q_map;
  const int m = Q.rows();
  if (Q.cols() != m || constraint.size() != m) {
    Rcpp::stop("Q and constraint must align with the serialized units.");
  }
  if (threads < 1 || block_size < 1) {
    Rcpp::stop("threads and block_size must be positive.");
  }
#ifndef _OPENMP
  if (threads > 1) Rcpp::stop("mgcvST was compiled without OpenMP; use threads = 1.");
#endif

  std::unique_ptr<SparsePrepared> local_cache;
  const SparsePrepared* cache = NULL;
  if (prepared == R_NilValue) {
    local_cache.reset(new SparsePrepared(Q, constraint));
    cache = local_cache.get();
  } else {
    if (!mgcvst_inla_sparse_prepared_valid_cpp(prepared)) {
      Rcpp::stop("prepared must be a valid, unserialized sparse INLA cache.");
    }
    Rcpp::XPtr<SparsePrepared> pointer(prepared);
    cache = pointer.get();
  }
  if (!cache->valid || cache->Q.rows() != m ||
      !cache->constraint.isApprox(constraint, 0.0)) {
    Rcpp::stop("prepared is not aligned with Q and constraint.");
  }

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
  std::vector<FeatureResult> result(features);
  int failed = 0;
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1) reduction(+:failed)
#endif
  for (int f = 0; f < features; ++f) {
    try {
      const ReconstructionUnit& unit = parsed[f];
      const SpMat H_U = unit.H_L.transpose();
      Mat M = Mat::Zero(m, m);
      for (int first = 0; first < m; first += block_size) {
        const int width = std::min(block_size, m - first);
        Mat Z = Mat::Zero(m, width);
        Z.block(first, 0, width, width).setIdentity();
        Z = project_coordinates(Z, cache->coordinate_constraint);
        Mat V = apply_B(cache->qfactor, Z) / std::sqrt(unit.tau);
        Mat Y = unit.K * V;
        Mat SY = stored_constrained_solve(
          unit.H_L, H_U, unit.H_D, unit.H_perm, constraint, Y,
          unit.hinv_g, unit.hden
        );
        Y.noalias() -= unit.K * SY;
        if (unit.U.cols()) {
          Y.noalias() -= unit.U * (unit.Vp * (unit.U.transpose() * V));
        }
        M.middleCols(first, width) = project_coordinates(
          apply_Bt(cache->qfactor, Y), cache->coordinate_constraint
        ) / std::sqrt(unit.tau);
      }
      result[f].a = unit.a;
      result[f].M = 0.5 * (M + M.transpose());
      result[f].statistic = unit.a.squaredNorm();
      result[f].Vp = unit.Vp;
      result[f].constraint_norm = unit.constraint_norm;
    } catch (const std::exception& error) {
      result[f].error = error.what();
      failed += 1;
    }
  }

  Rcpp::List out(features);
  for (int f = 0; f < features; ++f) {
    if (!result[f].error.empty()) {
      out[f] = Rcpp::List::create(Rcpp::Named("error") = result[f].error);
    } else {
      out[f] = Rcpp::List::create(
        Rcpp::Named("a") = result[f].a,
        Rcpp::Named("M") = result[f].M,
        Rcpp::Named("statistic") = result[f].statistic,
        Rcpp::Named("moments") = R_NilValue,
        Rcpp::Named("expected_vp") = result[f].Vp,
        Rcpp::Named("width") = m,
        Rcpp::Named("normalization") = m - 1,
        Rcpp::Named("constraint_diagnostic") = result[f].constraint_norm
      );
    }
  }
  out.attr("failed") = failed;
  out.attr("coordinate_width") = m;
  out.attr("normalization") = m - 1;
  return out;
}

// Materialize reduced curvature directly from sparse feature units.  basis is
// physical C = B R; coordinate is the matching whitened R.
// [[Rcpp::export]]
Rcpp::List mgcvst_inla_sparse_materialize_reduced_cpp(
    const Rcpp::List& units,
    const Eigen::MappedSparseMatrix<double>& Q_map,
    const Eigen::Map<Eigen::VectorXd> constraint,
    const Eigen::Map<Eigen::MatrixXd> coordinate,
    const Eigen::Map<Eigen::MatrixXd> basis,
    int threads = 1,
    SEXP prepared = R_NilValue) {
  const SpMat Q = Q_map;
  const int m = Q.rows();
  const int r = coordinate.cols();
  if (Q.cols() != m || constraint.size() != m || coordinate.rows() != m ||
      basis.rows() != m || basis.cols() != r) {
    Rcpp::stop("Q, constraint, coordinate, and basis must align.");
  }
  if (threads < 1) Rcpp::stop("threads must be positive.");
#ifndef _OPENMP
  if (threads > 1) Rcpp::stop("mgcvST was compiled without OpenMP; use threads = 1.");
#endif
  if (!mgcvst_inla_sparse_prepared_valid_cpp(prepared)) {
    Rcpp::stop("prepared must be a valid sparse INLA cache.");
  }
  Rcpp::XPtr<SparsePrepared> pointer(prepared);
  const SparsePrepared* cache = pointer.get();
  if (!cache->valid || cache->Q.rows() != m ||
      !cache->constraint.isApprox(constraint, 0.0)) {
    Rcpp::stop("prepared is not aligned with Q and constraint.");
  }

  const int features = units.size();
  const std::vector<ReconstructionUnit> parsed = parse_units(units, m);
  const Vec g = constraint;
  std::vector<FeatureResult> result(features);
  int failed = 0;
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1) reduction(+:failed)
#endif
  for (int f = 0; f < features; ++f) {
    try {
      const ReconstructionUnit& unit = parsed[f];
      result[f].a = coordinate.transpose() * unit.a;
      result[f].M = reduced_curvature(unit, basis, g);
      result[f].statistic = result[f].a.squaredNorm();
      result[f].Vp = unit.Vp;
      result[f].constraint_norm = unit.constraint_norm;
    } catch (const std::exception& error) {
      result[f].error = error.what();
      failed += 1;
    }
  }

  Rcpp::List out(features);
  for (int f = 0; f < features; ++f) {
    if (!result[f].error.empty()) {
      out[f] = Rcpp::List::create(Rcpp::Named("error") = result[f].error);
    } else {
      out[f] = Rcpp::List::create(
        Rcpp::Named("a") = result[f].a,
        Rcpp::Named("M") = result[f].M,
        Rcpp::Named("statistic") = result[f].statistic,
        Rcpp::Named("moments") = R_NilValue,
        Rcpp::Named("expected_vp") = result[f].Vp,
        Rcpp::Named("width") = r,
        Rcpp::Named("normalization") = r,
        Rcpp::Named("constraint_diagnostic") = result[f].constraint_norm
      );
    }
  }
  out.attr("failed") = failed;
  out.attr("coordinate_width") = r;
  out.attr("normalization") = r;
  return out;
}

// PCAlearning materialization: per feature, the reduced curvature H is formed
// in a thread-local buffer and reduced to the score coordinates a, the basis
// coefficients c = <B_k, H>_F, ||H||_F^2 and, for features with pack = TRUE,
// the float32 weighted-vech copy of H (upper triangle column-major, entry
// (i, j) at j (j + 1) / 2 + i, off-diagonal weight sqrt(2); the layout of
// src/pca_learning.cpp). pca_basis is the q (q + 1) / 2 x r weighted-vech basis
// (NULL: no projection). H itself is never returned.
// [[Rcpp::export]]
Rcpp::List mgcvst_inla_sparse_materialize_pca_cpp(
    const Rcpp::List& units,
    const Eigen::MappedSparseMatrix<double>& Q_map,
    const Eigen::Map<Eigen::VectorXd> constraint,
    const Eigen::Map<Eigen::MatrixXd> coordinate,
    const Eigen::Map<Eigen::MatrixXd> basis,
    Rcpp::Nullable<Rcpp::NumericMatrix> pca_basis,
    const Rcpp::LogicalVector& pack,
    int threads = 1,
    SEXP prepared = R_NilValue) {
  const SpMat Q = Q_map;
  const int m = Q.rows();
  const int q = coordinate.cols();
  const Eigen::Index L = (Eigen::Index)q * (q + 1) / 2;
  if (Q.cols() != m || constraint.size() != m || coordinate.rows() != m ||
      basis.rows() != m || basis.cols() != q) {
    Rcpp::stop("Q, constraint, coordinate, and basis must align.");
  }
  const int features = units.size();
  Rcpp::NumericMatrix PBr = pca_basis.isNotNull() ? Rcpp::NumericMatrix(pca_basis.get()) : Rcpp::NumericMatrix(L, 0);
  const Eigen::Map<Mat> PB(PBr.begin(), L, PBr.ncol());
  const int r = PB.cols();
#ifndef _OPENMP
  if (threads > 1) Rcpp::stop("mgcvST was compiled without OpenMP; use threads = 1.");
#endif
  if (!mgcvst_inla_sparse_prepared_valid_cpp(prepared)) {
    Rcpp::stop("prepared must be a valid sparse INLA cache.");
  }
  Rcpp::XPtr<SparsePrepared> pointer(prepared);
  const SparsePrepared* cache = pointer.get();
  if (!cache->valid || cache->Q.rows() != m ||
      !cache->constraint.isApprox(constraint, 0.0)) {
    Rcpp::stop("prepared is not aligned with Q and constraint.");
  }
  const std::vector<ReconstructionUnit> parsed = parse_units(units, m);
  const Vec g = constraint;
  Mat A(q, features), C(features, r);
  Vec fro2(features), statistic(features);
  std::vector<std::vector<float> > packed(features);
  std::vector<std::string> error(features);
  int failed = 0;
  const double w = std::sqrt(2.0);
#ifdef _OPENMP
#pragma omp parallel num_threads(threads) reduction(+:failed)
#endif
  {
    Vec h(L);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
    for (int f = 0; f < features; ++f) {
      try {
        const ReconstructionUnit& unit = parsed[f];
        const Mat M = reduced_curvature(unit, basis, g);
        Eigen::Index p = 0;
        for (int j = 0; j < q; ++j) {
          for (int i = 0; i < j; ++i) h[p++] = w * M(i, j);
          h[p++] = M(j, j);
        }
        A.col(f) = coordinate.transpose() * unit.a;
        statistic[f] = A.col(f).squaredNorm();
        fro2[f] = h.squaredNorm();
        if (r) C.row(f) = (PB.transpose() * h).transpose();
        if (pack[f]) {
          packed[f].resize(L);
          for (Eigen::Index k = 0; k < L; ++k) packed[f][k] = (float)h[k];
        }
      } catch (const std::exception& e) {
        error[f] = e.what();
        A.col(f).setConstant(NA_REAL);
        C.row(f).setConstant(NA_REAL);
        fro2[f] = statistic[f] = NA_REAL;
        failed += 1;
      }
    }
  }
  Rcpp::List packed_out(features);
  Rcpp::CharacterVector error_out(features);
  for (int f = 0; f < features; ++f) {
    error_out[f] = error[f].empty() ? NA_STRING : Rcpp::String(error[f]);
    if (!packed[f].empty()) {
      Rcpp::RawVector x(4 * L);
      std::memcpy(RAW(x), packed[f].data(), 4 * L);
      packed_out[f] = x;
      std::vector<float>().swap(packed[f]);
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("a") = A, Rcpp::Named("C") = C, Rcpp::Named("fro2") = fro2,
    Rcpp::Named("statistic") = statistic, Rcpp::Named("packed") = packed_out,
    Rcpp::Named("error") = error_out, Rcpp::Named("failed") = failed,
    Rcpp::Named("width") = q
  );
}
