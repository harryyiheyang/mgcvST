#define EIGEN_DONT_PARALLELIZE
#include <RcppEigen.h>
#include <Rmath.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

// Batched null fits of a log-link Poisson or negative-binomial GLM that shares
// its design X across features (spde.scan, R/spde-scan.R).
//
// The negative binomial is parameterised by its overdispersion alpha = 1/size,
// Var(y) = mu + alpha mu^2 (NB2), so that the Poisson model is the finite end
// point alpha = 0 of a bounded search rather than size -> infinity.  For each
// feature:
//   1. the Poisson IRLS gives beta(0) and the boundary score
//        dl/dalpha |_0 = sum_i [(y_i - mu_i)^2 - y_i] / 2
//      (Cameron and Trivedi's overdispersion score).  The profile likelihood
//      has this derivative at alpha = 0, so a score <= 0 is the KKT condition
//      of the end point and the maximum-likelihood fit is the Poisson fit;
//   2. otherwise alpha is found by the nested scheme of mgcv's nb(): IRLS for
//      beta at fixed alpha (concave in beta), and a safeguarded Newton step on
//      alpha for the profile likelihood, whose gradient is the partial
//      derivative at beta(alpha).
// The gamma-function terms depend on the counts only through the tail counts
// c_j = #{i : y_i > j}, so each outer step costs O(n + max y), exactly.
// The result is the joint maximum-likelihood estimate, which is the INLA
// conditional mode under the flat priors of the INLA branch up to the Laplace
// correction of the size.

namespace {

typedef Eigen::MatrixXd Mat;
typedef Eigen::VectorXd Vec;
typedef Eigen::MappedSparseMatrix<double> MSpMat;

struct NullFit {
  Vec beta;
  double alpha = NA_REAL;
  double loglik = NA_REAL;
  double score0 = NA_REAL;
  double phi = NA_REAL;
  int outer = 0;
  int inner = 0;
  bool converged = false;
  std::string error;
};

// q(x) / x^2 with q(x) = log1p(x) - x / (1 + x); series near 0.
inline double qt(double x) {
  if (x < 1e-2) {
    double s = 0, xk = 1;
    for (int k = 0; k <= 10; ++k) {
      s += (k % 2 ? -1.0 : 1.0) * (k + 1.0) / (k + 2.0) * xk;
      xk *= x;
    }
    return s;
  }
  return (std::log1p(x) - x / (1.0 + x)) / (x * x);
}

// r(x) = 1 / (x (1 + x)^2) - 2 qt(x) / x; series near 0.
inline double rt(double x) {
  if (x < 1e-2) {
    double s = 0, xk = 1;
    for (int k = 1; k <= 11; ++k) {
      s += (k % 2 ? -1.0 : 1.0) * (k + 1.0) * k / (k + 2.0) * xk;
      xk *= x;
    }
    return s;
  }
  return 1.0 / (x * (1.0 + x) * (1.0 + x)) - 2.0 * qt(x) / x;
}

// log1p(x) / x; 1 at 0.
inline double l1x(double x) {
  return x < 1e-8 ? 1.0 - 0.5 * x : std::log1p(x) / x;
}

// Deviance of the NB2 (alpha > 0) or Poisson (alpha = 0) fit.
double deviance(const Vec& y, const Vec& mu, double alpha) {
  double d = 0;
  for (int i = 0; i < y.size(); ++i) {
    const double yi = y[i], mi = mu[i];
    if (yi > 0) d += yi * std::log(yi / mi);
    d -= alpha > 0 ? (yi + 1.0 / alpha) * std::log1p(alpha * (yi - mi) / (1.0 + alpha * mi))
                   : yi - mi;
  }
  return 2.0 * d;
}

// IRLS for fixed alpha (alpha = 0 is Poisson), warm-started at beta. As in
// glm.fit(), it stops when the relative deviance change |dev - dev_old| /
// (|dev| + 0.1) is below tol; a coefficient that runs off under separation (a
// factor level with zero counts) then stops where its fitted means are
// negligible instead of failing. Returns the number of iterations, or -1 on
// failure.
// Weighted least squares min ||sw .* (X b - z)|| by column-pivoted QR, as
// glm.fit() solves it: a column whose pivot falls below rank_tol (relative)
// is aliased and its entry is 0, so a direction the data do not inform (a
// factor level whose fitted means have vanished) does not move.
Vec wls(const Mat& X, const Vec& sw, const Vec& z, double rank_tol) {
  Eigen::ColPivHouseholderQR<Mat> qr(sw.asDiagonal() * X);
  qr.setThreshold(rank_tol);
  return qr.solve((sw.array() * z.array()).matrix());
}

int irls(const Mat& X, const Vec& offset, const Vec& y, double alpha,
         Vec& beta, Vec& eta, Vec& mu, double tol, int maxit,
         int maxit_halving, double eta_max, double rank_tol) {
  const int n = X.rows();
  Vec sw(n), z(n);
  eta.noalias() = X * beta;
  eta += offset;
  mu = eta.array().exp();
  double dev_old = deviance(y, mu, alpha);
  for (int it = 1; it <= maxit; ++it) {
    // Fisher scoring: weights w = mu / (1 + alpha mu), working residual
    // (y - mu) / mu, so X'WX step = X'(y - mu) / (1 + alpha mu).
    for (int i = 0; i < n; ++i) {
      sw[i] = std::sqrt(mu[i] / (1.0 + alpha * mu[i]));
      z[i] = (y[i] - mu[i]) / mu[i];
    }
    const Vec step = wls(X, sw, z, rank_tol);
    if (!step.allFinite()) return -1;
    // Halve a step that would overflow the mean or raise the deviance (as
    // mgcv's gam.fit3 does); a rise below the stopping tolerance is accepted.
    const double dev_cap = dev_old + tol * (std::abs(dev_old) + 0.1);
    double scale = 1.0, dev = NA_REAL;
    bool moved = false;
    for (int h = 0; h < maxit_halving; ++h) {
      const Vec trial = beta + scale * step;
      const Vec eta_t = X * trial + offset;
      if ((eta_t.array() < eta_max).all()) {
        const Vec mu_t = eta_t.array().exp();
        const double dev_t = deviance(y, mu_t, alpha);
        if (std::isfinite(dev_t) && dev_t <= dev_cap) {
          beta = trial;
          eta = eta_t;
          mu = mu_t;
          dev = dev_t;
          moved = true;
          break;
        }
      }
      scale *= 0.5;
    }
    if (!moved) return -1;
    if (std::abs(dev - dev_old) / (std::abs(dev) + 0.1) < tol) return it;
    dev_old = dev;
  }
  return -1;
}

// Tail-count sums sum_j c_j j^m / (1 + alpha j)^m for m = 1, 2 and the log
// term sum_j c_j log1p(alpha j).
void tail_sums(const std::vector<double>& c, double alpha, double& s1,
               double& s2, double& sl) {
  s1 = 0; s2 = 0; sl = 0;
  for (size_t j = 1; j < c.size(); ++j) {
    if (c[j] == 0) continue;
    const double v = j / (1.0 + alpha * j);
    s1 += c[j] * v;
    s2 += c[j] * v * v;
    sl += c[j] * std::log1p(alpha * j);
  }
}

// NB2 log-likelihood at alpha >= 0 (alpha = 0 is Poisson).
double loglik(const Vec& mu, const Vec& y, const std::vector<double>& c,
              double alpha, double lfact) {
  double s1, s2, sl;
  tail_sums(c, alpha, s1, s2, sl);
  double ll = sl - lfact;
  for (int i = 0; i < mu.size(); ++i) {
    const double x = alpha * mu[i];
    ll += -mu[i] * l1x(x);
    if (y[i] > 0) ll += y[i] * (std::log(mu[i]) - std::log1p(x));
  }
  return ll;
}

// Partial derivatives of the log-likelihood in alpha at fixed beta.
void alpha_derivatives(const Vec& mu, const Vec& y,
                       const std::vector<double>& c, double alpha,
                       double& g, double& h) {
  double s1, s2, sl;
  tail_sums(c, alpha, s1, s2, sl);
  g = s1;
  h = -s2;
  for (int i = 0; i < mu.size(); ++i) {
    const double m = mu[i];
    const double x = alpha * m;
    const double d = 1.0 / (1.0 + x);
    g += m * m * qt(x) - y[i] * m * d;
    h += m * m * m * rt(x) + y[i] * m * m * d * d;
  }
}

} // namespace

// Y is features x observations, so it is passed as the transposed
// n x features dgCMatrix Yt: column f holds the counts of feature f.  offset
// is n x 1 (shared) or n x features.  poisson[f] = TRUE fits the Poisson
// model without the size search.  alpha_max bounds the overdispersion (size
// >= 1 / alpha_max); alpha_fixed >= 0 fixes it instead.
// [[Rcpp::export]]
Rcpp::List mgcvst_nb_null_batch_cpp(
    const Eigen::MappedSparseMatrix<double>& Yt,
    const Eigen::Map<Eigen::MatrixXd> X,
    const Eigen::Map<Eigen::MatrixXd> offset,
    const Rcpp::LogicalVector poisson,
    double alpha_max = 1e4,
    double alpha_fixed = -1,
    double tol = 1e-8,
    int maxit_outer = 25,
    int maxit_inner = 25,
    double alpha_min = 1e-12,
    double step_bound = 10,
    double step_fallback = 2,
    int maxit_halving = 30,
    double eta_max = 700,
    double rank_tol = 1e-11,
    int threads = 1) {
  const int n = Yt.rows();
  const int features = Yt.cols();
  const int p = X.cols();
  if (X.rows() != n) Rcpp::stop("X must have one row per observation.");
  if (offset.rows() != n || (offset.cols() != 1 && offset.cols() != features)) {
    Rcpp::stop("offset must be n x 1 or n x features.");
  }
  if (poisson.size() != features) Rcpp::stop("poisson must have one flag per feature.");
  if (p < 1) Rcpp::stop("X must have at least one column.");
  if (threads < 1) Rcpp::stop("threads must be positive.");
  if (!(alpha_max > 0) || !std::isfinite(alpha_max)) Rcpp::stop("alpha_max must be positive.");
  const Mat Xd = X;
  std::vector<int> is_pois(features);
  for (int f = 0; f < features; ++f) is_pois[f] = poisson[f] == TRUE;

  std::vector<NullFit> fit(features);

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
#endif
  for (int f = 0; f < features; ++f) {
    NullFit& out = fit[f];
    try {
      Vec y = Vec::Zero(n);
      double lfact = 0, ymax = 0;
      for (MSpMat::InnerIterator it(Yt, f); it; ++it) {
        const double v = it.value();
        if (v == 0) continue;
        if (v < 0 || v != std::floor(v)) {
          throw std::runtime_error("The counts must be non-negative integers.");
        }
        y[it.row()] = v;
        lfact += R::lgammafn(v + 1.0);
        ymax = std::max(ymax, v);
      }
      if (ymax > 1e8) throw std::runtime_error("A count is too large.");
      // No count: the Poisson MLE is at mu = 0, with no score to test.
      if (ymax == 0) throw std::runtime_error("All counts are zero.");
      // Tail counts c_j = #{i : y_i > j}, j = 0, ..., max y - 1.
      std::vector<double> c(static_cast<size_t>(ymax), 0.0);
      for (int i = 0; i < n; ++i) {
        const int yi = static_cast<int>(y[i]);
        if (yi > 0) c[yi - 1] += 1.0;
      }
      for (int j = static_cast<int>(c.size()) - 2; j >= 0; --j) c[j] += c[j + 1];

      const Vec off = offset.cols() == 1 ? Vec(offset.col(0)) : Vec(offset.col(f));
      // Start as glm.fit() does for the Poisson family: mu = y + 0.1, then
      // one weighted least-squares step on the working response.
      Vec beta;
      {
        const Vec m0 = (y.array() + 0.1).matrix();
        const Vec z = (m0.array().log() - off.array() + (y - m0).array() / m0.array()).matrix();
        beta = wls(Xd, m0.array().sqrt().matrix(), z, rank_tol);
        if (!beta.allFinite()) throw std::runtime_error("The start of the Poisson fit failed.");
      }
      Vec eta(n), mu(n);
      int it = irls(Xd, off, y, 0.0, beta, eta, mu, tol, maxit_inner, maxit_halving, eta_max, rank_tol);
      if (it < 0) throw std::runtime_error("The Poisson fit did not converge.");
      out.inner += it;
      // Boundary score and the Pearson dispersion of the Poisson fit.
      double score0 = 0, pearson = 0, mu2 = 0;
      for (int i = 0; i < n; ++i) {
        const double r = y[i] - mu[i];
        score0 += r * r - y[i];
        pearson += r * r / mu[i];
        mu2 += mu[i] * mu[i];
      }
      out.score0 = 0.5 * score0;
      out.phi = n > p ? pearson / (n - p) : NA_REAL;

      double alpha;
      bool converged = true;
      if (is_pois[f] || (alpha_fixed < 0 && score0 <= 0)) {
        alpha = 0;
      } else if (alpha_fixed >= 0) {
        alpha = alpha_fixed;
        if (alpha > 0) {
          it = irls(Xd, off, y, alpha, beta, eta, mu, tol, maxit_inner, maxit_halving, eta_max, rank_tol);
          if (it < 0) throw std::runtime_error("The fixed-size IRLS did not converge.");
          out.inner += it;
        }
      } else {
        // Moment start: sum (y - mu)^2 - y = alpha sum mu^2.
        alpha = std::min(alpha_max, std::max(alpha_min, score0 / mu2));
        converged = false;
        double ll_old = NA_REAL;
        for (int outer = 1; outer <= maxit_outer; ++outer) {
          out.outer = outer;
          it = irls(Xd, off, y, alpha, beta, eta, mu, tol, maxit_inner, maxit_halving, eta_max, rank_tol);
          if (it < 0) throw std::runtime_error("The fixed-size IRLS did not converge.");
          out.inner += it;
          // Stop on the relative change of the profile log-likelihood, the
          // same form as the IRLS deviance criterion.
          const double ll = loglik(mu, y, c, alpha, lfact);
          if (outer > 1 && std::abs(ll - ll_old) / (std::abs(ll) + 0.1) < tol) {
            converged = true; break;
          }
          ll_old = ll;
          double g, h;
          alpha_derivatives(mu, y, c, alpha, g, h);
          // Newton where the curvature is negative, else a step_fallback-fold
          // increase or decrease; each step stays within a factor step_bound
          // of the current alpha.
          double next = h < 0 ? alpha - g / h :
            (g > 0 ? step_fallback * alpha : alpha / step_fallback);
          next = std::min(step_bound * alpha, std::max(alpha / step_bound, next));
          next = std::min(alpha_max, next);
          if (alpha == alpha_max && g > 0) { converged = true; break; }
          alpha = next;
        }
        if (!converged) it = irls(Xd, off, y, alpha, beta, eta, mu, tol, maxit_inner, maxit_halving, eta_max, rank_tol);
        if (it < 0) throw std::runtime_error("The final IRLS did not converge.");
        out.inner += it;
      }
      out.beta = beta;
      out.alpha = alpha;
      out.loglik = loglik(mu, y, c, alpha, lfact);
      out.converged = converged;
    } catch (const std::exception& e) {
      out.error = e.what();
    }
  }

  Mat B = Mat::Constant(p, features, NA_REAL);
  Rcpp::NumericVector alpha(features), size(features), loglik_out(features),
    score0(features), phi(features);
  Rcpp::IntegerVector outer(features), inner(features);
  Rcpp::LogicalVector converged(features);
  Rcpp::CharacterVector error(features);
  for (int f = 0; f < features; ++f) {
    score0[f] = fit[f].score0;
    phi[f] = fit[f].phi;
    outer[f] = fit[f].outer;
    inner[f] = fit[f].inner;
    if (fit[f].error.empty()) {
      B.col(f) = fit[f].beta;
      alpha[f] = fit[f].alpha;
      size[f] = fit[f].alpha > 0 ? 1.0 / fit[f].alpha : R_PosInf;
      loglik_out[f] = fit[f].loglik;
      converged[f] = fit[f].converged;
      error[f] = NA_STRING;
    } else {
      alpha[f] = NA_REAL;
      size[f] = NA_REAL;
      loglik_out[f] = NA_REAL;
      converged[f] = false;
      error[f] = fit[f].error;
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("beta") = B, Rcpp::Named("alpha") = alpha,
    Rcpp::Named("size") = size, Rcpp::Named("loglik") = loglik_out,
    Rcpp::Named("score0") = score0, Rcpp::Named("phi") = phi,
    Rcpp::Named("converged") = converged, Rcpp::Named("outer") = outer,
    Rcpp::Named("inner") = inner, Rcpp::Named("error") = error
  );
}
