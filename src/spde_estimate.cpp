#define EIGEN_DONT_PARALLELIZE
#include <RcppEigen.h>
#include <Rmath.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <algorithm>
#include <cmath>
#include <functional>
#include <string>
#include <vector>
#include "nb_family.h"

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

// Spatial fits of a log-link Poisson or negative-binomial model with one SPDE
// field (spde.estimate, R/spde-estimate.R):
//
//   y_i ~ NB(mu_i, alpha),  log mu = offset + X beta + A u,
//   u ~ N(0, (tau Q)^{-1}) on the constraint c'u = 0,  beta flat.
//
// The negative binomial is NB2 with overdispersion alpha = 1/size, so the
// Poisson model is the finite end point alpha = 0. (tau, alpha) maximise the
// Laplace approximation of the marginal likelihood
//
//   L(tau, alpha) = l(y | beta, u, alpha) - tau u'Qu / 2 + (m - 1) log(tau) / 2
//                   - [log|H| + log(c~' H^{-1} c~)] / 2,
//
// where (beta, u) is the constrained mode at (tau, alpha), H the negative
// Hessian of the penalised log-likelihood there (observed information) and
// c~ = (0, c); the last bracket is the log-determinant of H on the constraint.
// Terms that depend on neither tau nor alpha are dropped.
//
// H is solved in blocks: S = A'WA + tau Q by a sparse Cholesky factor whose
// symbolic analysis is done once, and the p columns of X by their dense Schur
// complement G = X'WX - (A'WX)' S^{-1} (A'WX).
//
// (log tau, log alpha) are updated jointly: line searches over each in turn
// bring the point near the maximum, and Newton steps in both, with the
// gradient and the Hessian of L by central differences, finish it. The
// negative-binomial maximum is searched from the best alpha at the Poisson
// tau and, when that gains nothing, from the alpha of the null model; it is
// compared with the Poisson maximum at the end.
//
// Every stopping rule is a change of the objective: the inner Newton
// iteration stops on its predicted gain of the penalised log-likelihood, the
// line searches on the change of L over their three best points, and the
// joint Newton iteration on its predicted gain of L. The gradient of L at
// the returned point is reported.

namespace {

using namespace mgcvst_nb;
typedef Eigen::MatrixXd Mat;
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::MappedSparseMatrix<double> MSpMat;

struct Control {
  double tol, tol_inner, rank_tol, eta_max;
  double log_tau_min, log_tau_max, alpha_min, alpha_max;
  int maxit_inner, maxit_halving, maxit_line, maxit_outer;
  double fd_step, step_max;
};

struct Fit {
  Vec beta, u;
  double tau = NA_REAL, alpha = NA_REAL, laml = NA_REAL, laml_poisson = NA_REAL;
  double grad_log_tau = NA_REAL, grad_log_alpha = NA_REAL;
  int evaluations = 0, newton = 0, outer = 0;
  bool converged = false, poisson = false;
  std::string error;
};

// One feature's problem and the factor it reuses.
class Model {
public:
  Model(const SpMat& A, const SpMat& Q, const Mat& X, const Vec& c,
        const Control& ctl)
    : A_(A), At_(A.transpose()), Q_(Q), X_(X), c_(c), ctl_(ctl),
      n_(A.rows()), m_(A.cols()), p_(X.cols()) {
    // The pattern of S = A'WA + tau Q is fixed: analyse it once and record
    // where each product a_ij a_ik and each entry of Q lands in its values.
    SpMat pattern = At_ * A_;
    pattern += Q_;
    pattern.makeCompressed();
    S_ = pattern;
    llt_.analyzePattern(S_);
    auto slot = [&](int row, int col) {
      const int* inner = S_.innerIndexPtr();
      const int* outer = S_.outerIndexPtr();
      const int* hit = std::lower_bound(inner + outer[col], inner + outer[col + 1], row);
      if (hit == inner + outer[col + 1] || *hit != row) {
        throw std::runtime_error("The pattern of A'A + Q is incomplete.");
      }
      return static_cast<int>(hit - inner);
    };
    // At_ is column-major over observations: column i holds row i of A.
    for (int i = 0; i < n_; ++i) {
      for (SpMat::InnerIterator a(At_, i); a; ++a) {
        for (SpMat::InnerIterator b(At_, i); b; ++b) {
          aa_slot_.push_back(slot(a.row(), b.row()));
          aa_value_.push_back(a.value() * b.value());
        }
      }
      aa_end_.push_back(static_cast<int>(aa_slot_.size()));
    }
    for (int col = 0; col < m_; ++col) {
      for (SpMat::InnerIterator q(Q_, col); q; ++q) {
        q_slot_.push_back(slot(q.row(), col));
        q_value_.push_back(q.value());
      }
    }
  }

  void set_feature(const Vec& y, const Vec& offset) {
    y_ = y;
    off_ = offset;
    double ymax = y.size() ? y.maxCoeff() : 0.0;
    lfact_ = 0;
    tail_.assign(static_cast<size_t>(ymax), 0.0);
    for (int i = 0; i < n_; ++i) {
      const int yi = static_cast<int>(y[i]);
      if (yi > 0) tail_[yi - 1] += 1.0;
      lfact_ += R::lgammafn(y[i] + 1.0);
    }
    for (int j = static_cast<int>(tail_.size()) - 2; j >= 0; --j) tail_[j] += tail_[j + 1];
    separate();
    // Start as glm.fit() does for the Poisson family, with no field.
    const Vec m0 = (y.array() + 0.1).matrix();
    const Vec z = (m0.array().log() - off_.array() + (y - m0).array() / m0.array()).matrix();
    Vec sw = m0.array().sqrt().matrix();
    for (int i : gone_) sw[i] = 0;
    Eigen::CompleteOrthogonalDecomposition<Mat> qr;
    qr.setThreshold(ctl_.rank_tol);
    qr.compute(sw.asDiagonal() * X_);
    beta_ = qr.solve((sw.array() * z.array()).matrix());
    u_ = Vec::Zero(m_);
    evaluations = 0;
    newton = 0;
  }

  // L(log tau, alpha) at the constrained mode, warm-started at the last mode.
  // NaN when the mode is not found.
  double laml(double rho, double alpha) {
    ++evaluations;
    const double tau = std::exp(rho);
    Vec beta = beta_, u = u_, eta(n_), mu(n_), w(n_), s(n_);
    if (!linear_predictor(beta, u, eta, mu)) {
      beta = beta_start_or_zero();
      u.setZero();
      if (!linear_predictor(beta, u, eta, mu)) return NA_REAL;
    }
    double P = penalised(mu, u, tau, alpha);
    bool still = false;
    for (int it = 0; it <= ctl_.maxit_inner; ++it) {
      // Observed information and score in eta.
      for (int i = 0; i < n_; ++i) {
        const double f = 1.0 / (1.0 + alpha * mu[i]);
        w[i] = mu[i] * (1.0 + alpha * y_[i]) * f * f;
        s[i] = (y_[i] - mu[i]) * f;
      }
      if (!factor(w, tau)) return NA_REAL;
      const Vec rb = X_.transpose() * s;
      const Vec ru = At_ * s - tau * (Q_ * u);
      Vec db, du;
      solve(rb, ru, db, du);
      // Predicted gain of the penalised log-likelihood of this Newton step.
      const double gain = 0.5 * (rb.dot(db) + ru.dot(du));
      if (!std::isfinite(gain)) return NA_REAL;
      if (still || gain <= ctl_.tol_inner * (std::abs(P) + 0.1)) {
        beta_ = beta;
        u_ = u;
        return loglik(mu, y_, tail_, alpha, lfact_) - 0.5 * tau * u.dot(Q_ * u) +
          0.5 * (m_ - 1) * rho - 0.5 * (logdet_S_ + logdet_K_);
      }
      if (it == ctl_.maxit_inner) break;
      ++newton;
      // Halve a step that overflows the mean or lowers the objective.
      double scale = 1.0;
      bool moved = false;
      for (int h = 0; h < ctl_.maxit_halving; ++h) {
        const Vec bt = beta + scale * db, ut = u + scale * du;
        Vec et(n_), mt(n_);
        if (linear_predictor(bt, ut, et, mt)) {
          const double Pt = penalised(mt, ut, tau, alpha);
          if (std::isfinite(Pt) && Pt >= P - ctl_.tol_inner * (std::abs(P) + 0.1)) {
            // The objective changed by less than the tolerance: the mode.
            // The predicted gain can stay above the tolerance at the
            // resolution of the objective, where only a step halved many
            // times is accepted.
            still = std::abs(Pt - P) <= ctl_.tol_inner * (std::abs(P) + 0.1);
            beta = bt; u = ut; eta = et; mu = mt; P = Pt;
            moved = true;
            break;
          }
        }
        scale *= 0.5;
      }
      if (!moved) return NA_REAL;
    }
    return NA_REAL;
  }

  // The fixed effects, moved along each separating direction until the
  // fitted means of its rows are at most rank_tol: their limit is zero.
  Vec beta() const {
    Vec beta = beta_;
    for (size_t k = 0; k < sep_.size(); ++k) {
      const Vec z = X_ * sep_[k];
      const Vec eta = off_ + A_ * u_ + X_ * beta;
      double t = 0;
      for (int i : sep_rows_[k]) t = std::max(t, (eta[i] - std::log(ctl_.rank_tol)) / -z[i]);
      beta += t * sep_[k];
    }
    return beta;
  }
  const Vec& u() const { return u_; }
  int evaluations = 0, newton = 0;

private:
  Vec beta_start_or_zero() const { return Vec::Zero(p_); }

  bool linear_predictor(const Vec& beta, const Vec& u, Vec& eta, Vec& mu) const {
    eta = off_ + A_ * u;
    if (p_) eta.noalias() += X_ * beta;
    for (int i : gone_) eta[i] = 0;
    if (!eta.allFinite() || !(eta.array() < ctl_.eta_max).all()) return false;
    mu = eta.array().exp();
    for (int i : gone_) mu[i] = 0;
    return true;
  }

  // Separation: a direction v of the fixed effects with Xv = 0 on every row
  // with a count and Xv <= 0 elsewhere (a factor level without counts). The
  // likelihood rises along v for ever and the fitted means of the rows with
  // Xv < 0 go to zero. Those rows take their limit, mean zero, so that the
  // criterion is the one of the model without them and without v. The
  // directions tried are the null vectors of the rows with a count, one per
  // column a pivoted QR leaves out, each tried on the rows still in.
  void separate() {
    gone_.clear();
    sep_.clear();
    sep_rows_.clear();
    if (!p_ || !(y_.array() > 0).any()) return;
    std::vector<char> out(n_, 0);
    {
      std::vector<int> pos;
      for (int i = 0; i < n_; ++i) if (y_[i] > 0) pos.push_back(i);
      Mat Xp(pos.size(), p_);
      for (size_t k = 0; k < pos.size(); ++k) Xp.row(k) = X_.row(pos[k]);
      Eigen::ColPivHouseholderQR<Mat> qr;
      qr.setThreshold(ctl_.rank_tol);
      qr.compute(Xp);
      const int r = qr.rank();
      if (r == p_) return;
      const Mat R = qr.matrixR().topRows(r);
      const Mat R1 = R.leftCols(r).triangularView<Eigen::Upper>();
      for (int j = r; j < p_; ++j) {
        Vec v = Vec::Zero(p_);
        if (r) v.head(r) = -R1.triangularView<Eigen::Upper>().solve(R.col(j));
        v[j] = 1;
        v = qr.colsPermutation() * v;
        const Vec z = X_ * v;
        const double small = ctl_.rank_tol * X_.cwiseAbs().maxCoeff() * v.cwiseAbs().maxCoeff();
        double lo = 0, hi = 0;
        for (int i = 0; i < n_; ++i) {
          if (out[i] || y_[i] > 0) continue;
          lo = std::min(lo, z[i]);
          hi = std::max(hi, z[i]);
        }
        if ((lo < -small) == (hi > small)) continue;  // aliased, or two signs
        if (hi > small) v = -v;
        std::vector<int> rows;
        for (int i = 0; i < n_; ++i) {
          if (!out[i] && !(y_[i] > 0) && std::abs(z[i]) > small) { out[i] = 1; rows.push_back(i); }
        }
        gone_.insert(gone_.end(), rows.begin(), rows.end());
        sep_.push_back(v);
        sep_rows_.push_back(rows);
      }
    }
  }

  double penalised(const Vec& mu, const Vec& u, double tau, double alpha) const {
    return loglik(mu, y_, tail_, alpha, lfact_) - 0.5 * tau * u.dot(Q_ * u);
  }

  // Factor S = A'WA + tau Q and the Schur complement of X.
  bool factor(const Vec& w, double tau) {
    double* value = S_.valuePtr();
    std::fill(value, value + S_.nonZeros(), 0.0);
    int k = 0;
    for (int i = 0; i < n_; ++i) {
      for (const int end = aa_end_[i]; k < end; ++k) value[aa_slot_[k]] += w[i] * aa_value_[k];
    }
    for (size_t j = 0; j < q_slot_.size(); ++j) value[q_slot_[j]] += tau * q_value_[j];
    llt_.factorize(S_);
    if (llt_.info() != Eigen::Success) return false;
    const Vec d = llt_.vectorD();
    if (!(d.array() > 0).all()) return false;
    logdet_S_ = d.array().log().sum();
    // The Newton system with the constraint c'u = 0, after the field is
    // eliminated: K (db, nu) = (rb - B'S^{-1}ru, -c'S^{-1}ru), with
    //   K = [ G, -b ; -b', -d ],  G = X'WX - B'S^{-1}B,  b = B'S^{-1}c,
    //   d = c'S^{-1}c,
    // and log|G| + log(c'H^{-1}c) = log|det K|. G alone is nearly singular
    // along the intercept when the field can take a constant up, which the
    // constraint forbids; K is not.
    Sic_ = llt_.solve(c_);
    const double dcc = c_.dot(Sic_);
    if (!(dcc > 0) || !std::isfinite(dcc)) return false;
    Mat K(p_ + 1, p_ + 1);
    if (p_) {
      const Mat WX = w.asDiagonal() * X_;
      B_ = At_ * WX;
      SiB_ = llt_.solve(B_);
      const Mat G = X_.transpose() * WX - B_.transpose() * SiB_;
      K.topLeftCorner(p_, p_) = 0.5 * (G + G.transpose());
      const Vec b = B_.transpose() * Sic_;
      K.topRightCorner(p_, 1) = -b;
      K.bottomLeftCorner(1, p_) = -b.transpose();
    }
    K(p_, p_) = -dcc;
    // An aliased direction (a factor level whose fitted means have
    // vanished) is left out of the solve and of the determinant.
    Eigen::SelfAdjointEigenSolver<Mat> eig(K);
    if (eig.info() != Eigen::Success) return false;
    const Vec ev = eig.eigenvalues();
    const double top = ev.array().abs().maxCoeff();
    if (!(top > 0) || !std::isfinite(top)) return false;
    Vec inv = Vec::Zero(p_ + 1);
    logdet_K_ = 0;
    for (int j = 0; j <= p_; ++j) {
      if (std::abs(ev[j]) > ctl_.rank_tol * top) {
        inv[j] = 1.0 / ev[j];
        logdet_K_ += std::log(std::abs(ev[j]));
      }
    }
    Kinv_ = eig.eigenvectors() * inv.asDiagonal() * eig.eigenvectors().transpose();
    return true;
  }

  // The Newton step of (beta, u) under c'u = 0, from the scores rb and ru.
  void solve(const Vec& rb, const Vec& ru, Vec& db, Vec& du) const {
    const Vec t = llt_.solve(ru);
    Vec rhs(p_ + 1);
    if (p_) rhs.head(p_) = rb - B_.transpose() * t;
    rhs[p_] = -c_.dot(t);
    const Vec z = Kinv_ * rhs;
    db = z.head(p_);
    du = t - z[p_] * Sic_;
    if (p_) du.noalias() -= SiB_ * db;
  }

  const SpMat& A_;
  const SpMat At_;
  const SpMat& Q_;
  const Mat& X_;
  const Vec& c_;
  const Control& ctl_;
  const int n_, m_, p_;
  Vec y_, off_, beta_, u_;
  std::vector<int> gone_;
  std::vector<Vec> sep_;
  std::vector<std::vector<int>> sep_rows_;
  std::vector<double> tail_;
  double lfact_ = 0, logdet_S_ = 0, logdet_K_ = 0;
  Eigen::SimplicialLDLT<SpMat> llt_;
  SpMat S_;
  std::vector<int> aa_slot_, aa_end_, q_slot_;
  std::vector<double> aa_value_, q_value_;
  Mat B_, SiB_, Kinv_;
  Vec Sic_;
};

// Brent's maximiser of f on [a, b]. It stops when the objective at its three
// best points agrees within tol (relative, as glm.fit() measures a deviance
// change), or after maxit evaluations. A NaN objective is the worst value.
struct LineResult { double x, f; bool converged; };

LineResult maximise(const std::function<double(double)>& f, double a, double b,
                    double tol, int maxit, double start = NA_REAL) {
  const double golden = 0.3819660112501051;
  auto value = [&](double x) {
    const double v = f(x);
    return std::isfinite(v) ? -v : R_PosInf;
  };
  // The first point: `start` when it lies inside (a, b), else the golden point.
  double x = std::isfinite(start) && start > a && start < b ? start : a + golden * (b - a);
  double w = x, v = x;
  double fx = value(x), fw = fx, fv = fx;
  double d = 0, e = 0;
  bool converged = false;
  for (int it = 1; it < maxit; ++it) {
    const double xm = 0.5 * (a + b);
    const double eps = std::sqrt(std::numeric_limits<double>::epsilon()) * (std::abs(x) + 1.0);
    if (v != w && w != x && v != x && std::isfinite(fx) &&
        std::max(std::abs(fw - fx), std::abs(fv - fx)) <= tol * (std::abs(fx) + 0.1)) {
      converged = true;
      break;
    }
    if (std::abs(x - xm) <= 2 * eps - 0.5 * (b - a)) { converged = true; break; }
    bool golden_step = true;
    if (std::abs(e) > eps && std::isfinite(fx) && std::isfinite(fw) && std::isfinite(fv)) {
      double r = (x - w) * (fx - fv), q = (x - v) * (fx - fw);
      double pp = (x - v) * q - (x - w) * r;
      q = 2 * (q - r);
      if (q > 0) pp = -pp;
      q = std::abs(q);
      const double e_old = e;
      e = d;
      if (std::abs(pp) < std::abs(0.5 * q * e_old) && pp > q * (a - x) && pp < q * (b - x)) {
        d = pp / q;
        const double uu = x + d;
        if (uu - a < 2 * eps || b - uu < 2 * eps) d = xm >= x ? eps : -eps;
        golden_step = false;
      }
    }
    if (golden_step) {
      e = (x >= xm ? a : b) - x;
      d = golden * e;
    }
    const double uu = std::abs(d) >= eps ? x + d : x + (d > 0 ? eps : -eps);
    const double fu = value(uu);
    if (fu <= fx) {
      if (uu >= x) a = x; else b = x;
      v = w; fv = fw; w = x; fw = fx; x = uu; fx = fu;
    } else {
      if (uu < x) a = uu; else b = uu;
      if (fu <= fw || w == x) { v = w; fv = fw; w = uu; fw = fu; }
      else if (fu <= fv || v == x || v == w) { v = uu; fv = fu; }
    }
  }
  return {x, -fx, converged};
}

// A point of the search in (rho, t) = (log tau, log alpha) and the gradient of
// L there, when it was computed.
struct Point {
  double rho = NA_REAL, t = NA_REAL, f = NA_REAL;
  double g_rho = NA_REAL, g_t = NA_REAL;
  bool converged = false, failed = false, at_zero = false;
};

// Climb L from (rho, t) with value f: line searches over rho and t in turn
// until a cycle gains nothing, then Newton steps in (rho, t) jointly until the
// predicted gain is below the tolerance. at_zero marks a point that reached
// the lower end of alpha, where the model is the Poisson model.
Point climb(Model& model, const Control& ctl, double rho, double t, double f,
            int& cycles) {
  Point z;
  z.rho = rho; z.t = t; z.f = f;
  if (!std::isfinite(f)) { z.failed = true; return z; }
  const double lo = std::log(ctl.alpha_min), hi = std::log(ctl.alpha_max);
  const double h = ctl.fd_step;
  auto threshold = [&](double v) { return ctl.tol * (std::abs(v) + 0.1); };
  bool lines = true;
  for (int cycle = 1; cycle <= ctl.maxit_outer; ++cycle) {
    ++cycles;
    const double before = z.f;
    const LineResult rt = maximise(
      [&](double x) { return model.laml(x, std::exp(z.t)); },
      ctl.log_tau_min, ctl.log_tau_max, ctl.tol, ctl.maxit_line);
    lines = rt.converged;
    if (rt.f > z.f) { z.rho = rt.x; z.f = rt.f; }
    const LineResult ra = maximise(
      [&](double x) { return model.laml(z.rho, std::exp(x)); }, lo, hi, ctl.tol,
      ctl.maxit_line);
    lines = lines && ra.converged;
    if (ra.f > z.f) { z.t = ra.x; z.f = ra.f; }
    if (z.f - before <= threshold(z.f)) break;
  }
  z.converged = lines;
  for (int it = 0; it < ctl.maxit_outer; ++it) {
    if (z.t - h <= lo) { z.at_zero = true; return z; }
    // No room for the differences at an end of the interval: the point of the
    // line searches stands.
    if (z.t + h >= hi || z.rho - h <= ctl.log_tau_min || z.rho + h >= ctl.log_tau_max) return z;
    auto value = [&](double dr, double dt) { return model.laml(z.rho + dr, std::exp(z.t + dt)); };
    const double fr1 = value(h, 0), fr0 = value(-h, 0), ft1 = value(0, h), ft0 = value(0, -h);
    const double f11 = value(h, h), f00 = value(-h, -h);
    if (!std::isfinite(fr1 + fr0 + ft1 + ft0 + f11 + f00)) { z.converged = false; return z; }
    const double gr = (fr1 - fr0) / (2 * h), gt = (ft1 - ft0) / (2 * h);
    const double hrr = (fr1 - 2 * z.f + fr0) / (h * h), htt = (ft1 - 2 * z.f + ft0) / (h * h);
    const double hrt = (f11 + f00 - 2 * z.f) / (2 * h * h) - 0.5 * (hrr + htt);
    z.g_rho = gr; z.g_t = gt;
    const double det = hrr * htt - hrt * hrt;
    double dr, dt, gain;
    if (hrr < 0 && det > 0) {
      // Newton step of a concave quadratic and its gain.
      dr = -(htt * gr - hrt * gt) / det;
      dt = -(hrr * gt - hrt * gr) / det;
      gain = 0.5 * (gr * dr + gt * dt);
    } else {
      // Not concave here: a gradient step of length step_max.
      const double norm = std::sqrt(gr * gr + gt * gt);
      if (!(norm > 0)) { z.converged = true; return z; }
      dr = ctl.step_max * gr / norm;
      dt = ctl.step_max * gt / norm;
      gain = gr * dr + gt * dt;
    }
    if (!std::isfinite(gain)) { z.converged = false; return z; }
    if (gain <= threshold(z.f)) { z.converged = true; return z; }
    const double longest = std::max(std::abs(dr), std::abs(dt));
    double scale = longest > ctl.step_max ? ctl.step_max / longest : 1.0;
    bool moved = false;
    for (int k = 0; k < ctl.maxit_halving; ++k) {
      const double r1 = std::min(ctl.log_tau_max, std::max(ctl.log_tau_min, z.rho + scale * dr));
      const double t1 = std::min(hi, std::max(lo, z.t + scale * dt));
      const double f1 = model.laml(r1, std::exp(t1));
      if (std::isfinite(f1) && f1 > z.f) {
        z.rho = r1; z.t = t1; z.f = f1;
        z.g_rho = z.g_t = NA_REAL;
        moved = true;
        break;
      }
      scale *= 0.5;
    }
    // No step raises L although a gain was predicted: the differences are at
    // the resolution of L.
    if (!moved) { z.converged = false; return z; }
  }
  z.converged = false;
  return z;
}

// alpha_start is the alpha of the null model (1 / size), NaN when there is
// none.
Fit fit_feature(Model& model, const Vec& y, const Vec& offset, bool poisson,
                const Control& ctl, double alpha_start) {
  Fit out;
  if (!(y.array() > 0).any()) {
    out.error = "All counts are zero.";
    return out;
  }
  model.set_feature(y, offset);
  // Poisson: one line search over log tau.
  const LineResult r = maximise([&](double rho) { return model.laml(rho, 0.0); },
                                ctl.log_tau_min, ctl.log_tau_max, ctl.tol, ctl.maxit_line);
  if (!std::isfinite(r.f)) {
    out.error = "The Poisson spatial fit failed.";
    return out;
  }
  double rho = r.x, alpha = 0;
  bool converged = r.converged;
  out.laml_poisson = r.f;
  out.poisson = true;
  if (!poisson) {
    const double lo = std::log(ctl.alpha_min), hi = std::log(ctl.alpha_max);
    const double need = r.f + ctl.tol * (std::abs(r.f) + 0.1);
    Point best;
    bool have = false, searched = false;
    // From the best alpha at the Poisson tau.
    const LineResult ra = maximise(
      [&](double t) { return model.laml(r.x, std::exp(t)); }, lo, hi, ctl.tol,
      ctl.maxit_line, std::log(alpha_start));
    if (std::isfinite(ra.f)) {
      searched = true;
      const Point z = climb(model, ctl, r.x, ra.x, ra.f, out.outer);
      if (!z.failed) { best = z; have = true; }
    }
    // From the alpha of the null model, when the first search ends on the
    // Poisson model.
    if ((!have || best.at_zero || !(best.f > need)) && std::isfinite(alpha_start)) {
      const double t0 = std::min(hi, std::max(lo, std::log(alpha_start)));
      const Point z = climb(model, ctl, r.x, t0, model.laml(r.x, std::exp(t0)), out.outer);
      if (!z.failed) {
        searched = true;
        if (!have || z.f > best.f) { best = z; have = true; }
      }
    }
    if (!searched) {
      out.error = "The negative-binomial spatial fit failed at every overdispersion.";
      return out;
    }
    if (have && !best.at_zero && best.f > need) {
      rho = best.rho;
      alpha = std::exp(best.t);
      converged = best.converged;
      out.poisson = false;
      out.grad_log_tau = best.g_rho;
      out.grad_log_alpha = best.g_t;
    }
  }
  // The mode at the optimum (the searches end on other points).
  const double last = model.laml(rho, alpha);
  if (!std::isfinite(last)) {
    out.error = "The spatial fit failed at its optimum.";
    return out;
  }
  out.beta = model.beta();
  out.u = model.u();
  out.tau = std::exp(rho);
  out.alpha = alpha;
  out.laml = last;
  // The gradient at the returned point, where the Newton iteration did not
  // leave one (the Poisson fit, or a last accepted step).
  const double h = ctl.fd_step;
  const Vec beta = out.beta, u = out.u;
  if (rho - h > ctl.log_tau_min && rho + h < ctl.log_tau_max &&
      (out.poisson || !std::isfinite(out.grad_log_tau))) {
    out.grad_log_tau = (model.laml(rho + h, alpha) - model.laml(rho - h, alpha)) / (2 * h);
    if (!out.poisson) {
      const double t = std::log(alpha);
      out.grad_log_alpha = (model.laml(rho, std::exp(t + h)) - model.laml(rho, std::exp(t - h))) / (2 * h);
    }
  }
  out.converged = converged;
  out.evaluations = model.evaluations;
  out.newton = model.newton;
  return out;
}

} // namespace

// Yt is the transposed n x features count matrix; offset is n x 1 (shared) or
// n x features; Q is the prior precision of the field at tau = 1 and
// constraint the vector c of c'u = 0. poisson[f] = TRUE fits the Poisson
// model only.
// [[Rcpp::export]]
Rcpp::List mgcvst_spde_estimate_cpp(
    const Eigen::MappedSparseMatrix<double>& Yt,
    const Eigen::Map<Eigen::MatrixXd> X,
    const Eigen::MappedSparseMatrix<double>& A,
    const Eigen::MappedSparseMatrix<double>& Q,
    const Eigen::Map<Eigen::VectorXd> constraint,
    const Eigen::Map<Eigen::MatrixXd> offset,
    const Rcpp::LogicalVector poisson,
    const Rcpp::NumericVector alpha_start,
    double tol = 1e-8,
    double tol_inner = 1e-13,
    int maxit_inner = 50,
    int maxit_halving = 30,
    int maxit_line = 100,
    int maxit_outer = 25,
    double log_tau_min = -12,
    double log_tau_max = 12,
    double alpha_min = 1e-6,
    double alpha_max = 1e4,
    double eta_max = 700,
    double rank_tol = 1e-11,
    double fd_step = 0.01,
    double step_max = 2,
    int threads = 1) {
  const int n = Yt.rows(), features = Yt.cols(), p = X.cols(), m = A.cols();
  if (X.rows() != n || A.rows() != n) Rcpp::stop("X and A must have one row per observation.");
  if (Q.rows() != m || Q.cols() != m) Rcpp::stop("Q must be square with ncol(A) columns.");
  if (constraint.size() != m) Rcpp::stop("constraint must have ncol(A) entries.");
  if (offset.rows() != n || (offset.cols() != 1 && offset.cols() != features)) {
    Rcpp::stop("offset must be n x 1 or n x features.");
  }
  if (poisson.size() != features) Rcpp::stop("poisson must have one flag per feature.");
  if (alpha_start.size() != features) Rcpp::stop("alpha_start must have one value per feature.");
  if (threads < 1) Rcpp::stop("threads must be positive.");
  if (!(fd_step > 0) || !(step_max > 0)) Rcpp::stop("fd_step and step_max must be positive.");
  if (!(log_tau_min < log_tau_max) || !(alpha_min > 0) || !(alpha_min < alpha_max)) {
    Rcpp::stop("The search bounds are not ordered.");
  }
  const Control ctl{tol, tol_inner, rank_tol, eta_max, log_tau_min, log_tau_max,
                    alpha_min, alpha_max, maxit_inner, maxit_halving, maxit_line, maxit_outer,
                    fd_step, step_max};
  const SpMat As = A, Qs = Q;
  const Mat Xd = X;
  const Vec c = constraint;
  std::vector<int> is_pois(features);
  std::vector<double> start(features);
  for (int f = 0; f < features; ++f) {
    is_pois[f] = poisson[f] == TRUE;
    start[f] = alpha_start[f] > 0 && R_finite(alpha_start[f]) ? alpha_start[f] : NA_REAL;
  }
  std::vector<Fit> fit(features);
  // An error of the set-up (the sparse analysis) is raised here, outside the
  // parallel region, where an exception cannot be thrown.
  { Model check(As, Qs, Xd, c, ctl); }

#ifdef _OPENMP
#pragma omp parallel num_threads(threads)
#endif
  {
    Model model(As, Qs, Xd, c, ctl);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
    for (int f = 0; f < features; ++f) {
      try {
        Vec y = Vec::Zero(n);
        for (MSpMat::InnerIterator it(Yt, f); it; ++it) {
          const double v = it.value();
          if (v < 0 || v != std::floor(v)) throw std::runtime_error("The counts must be non-negative integers.");
          y[it.row()] = v;
        }
        const Vec off = offset.cols() == 1 ? Vec(offset.col(0)) : Vec(offset.col(f));
        fit[f] = fit_feature(model, y, off, is_pois[f], ctl, start[f]);
      } catch (const std::exception& e) {
        fit[f] = Fit();
        fit[f].error = e.what();
      }
    }
  }

  Mat B = Mat::Constant(p, features, NA_REAL), U = Mat::Constant(m, features, NA_REAL);
  Rcpp::NumericVector tau(features), alpha(features), laml(features), laml0(features),
    grad_tau(features), grad_alpha(features);
  Rcpp::IntegerVector evaluations(features), newton(features), outer(features);
  Rcpp::LogicalVector converged(features), pois(features);
  Rcpp::CharacterVector error(features);
  for (int f = 0; f < features; ++f) {
    const Fit& z = fit[f];
    const bool ok = z.error.empty();
    if (ok) { B.col(f) = z.beta; U.col(f) = z.u; }
    tau[f] = z.tau; alpha[f] = z.alpha; laml[f] = z.laml; laml0[f] = z.laml_poisson;
    grad_tau[f] = z.grad_log_tau; grad_alpha[f] = z.grad_log_alpha;
    evaluations[f] = z.evaluations; newton[f] = z.newton; outer[f] = z.outer;
    converged[f] = ok && z.converged;
    pois[f] = ok ? (z.poisson ? TRUE : FALSE) : NA_LOGICAL;
    if (ok) error[f] = NA_STRING; else error[f] = z.error;
  }
  return Rcpp::List::create(
    Rcpp::_["beta"] = B, Rcpp::_["u"] = U, Rcpp::_["tau"] = tau,
    Rcpp::_["alpha"] = alpha, Rcpp::_["laml"] = laml,
    Rcpp::_["laml_poisson"] = laml0, Rcpp::_["poisson"] = pois,
    Rcpp::_["converged"] = converged,
    Rcpp::_["grad_log_tau"] = grad_tau, Rcpp::_["grad_log_alpha"] = grad_alpha,
    Rcpp::_["evaluations"] = evaluations, Rcpp::_["newton"] = newton,
    Rcpp::_["outer"] = outer, Rcpp::_["error"] = error);
}

// L(tau, alpha) of one feature at the pairs (log_tau[k], alpha[k]), each at
// its own constrained mode (spde.laml, R/spde-estimate.R). alpha = 0 is the
// Poisson model. NaN where the mode is not found.
// [[Rcpp::export]]
Rcpp::NumericVector mgcvst_spde_laml_cpp(
    const Eigen::Map<Eigen::VectorXd> y,
    const Eigen::Map<Eigen::MatrixXd> X,
    const Eigen::MappedSparseMatrix<double>& A,
    const Eigen::MappedSparseMatrix<double>& Q,
    const Eigen::Map<Eigen::VectorXd> constraint,
    const Eigen::Map<Eigen::VectorXd> offset,
    const Eigen::Map<Eigen::VectorXd> log_tau,
    const Eigen::Map<Eigen::VectorXd> alpha,
    double tol_inner = 1e-13,
    int maxit_inner = 50,
    int maxit_halving = 30,
    double eta_max = 700,
    double rank_tol = 1e-11) {
  const int n = y.size(), m = A.cols();
  if (X.rows() != n || A.rows() != n || offset.size() != n) Rcpp::stop("y, X, A and offset must have one row per observation.");
  if (Q.rows() != m || Q.cols() != m || constraint.size() != m) Rcpp::stop("Q and constraint must match ncol(A).");
  if (log_tau.size() != alpha.size()) Rcpp::stop("log_tau and alpha must have one value per point.");
  Control ctl{};
  ctl.tol_inner = tol_inner; ctl.rank_tol = rank_tol; ctl.eta_max = eta_max;
  ctl.maxit_inner = maxit_inner; ctl.maxit_halving = maxit_halving;
  const SpMat As = A, Qs = Q;
  const Mat Xd = X;
  const Vec c = constraint, yv = y, off = offset;
  Model model(As, Qs, Xd, c, ctl);
  model.set_feature(yv, off);
  Rcpp::NumericVector out(log_tau.size());
  for (int k = 0; k < log_tau.size(); ++k) {
    if (!(alpha[k] >= 0) || !std::isfinite(log_tau[k])) Rcpp::stop("alpha must be non-negative and log_tau finite.");
    out[k] = model.laml(log_tau[k], alpha[k]);
  }
  return out;
}
