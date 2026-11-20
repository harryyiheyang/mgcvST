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
// Every stopping rule is a change of the objective: the inner Newton
// iteration stops on its predicted gain of the penalised log-likelihood, the
// line searches on the change of L over their three best points, and the
// alternation over (log tau, alpha) on the gain of L of a cycle.

namespace {

using namespace mgcvst_nb;
typedef Eigen::MatrixXd Mat;
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::MappedSparseMatrix<double> MSpMat;

struct Control {
  double tol, tol_inner, rank_tol, eta_max;
  double log_tau_min, log_tau_max, alpha_min, alpha_max;
  int maxit_inner, maxit_halving, maxit_line, maxit_outer;
};

struct Fit {
  Vec beta, u;
  double tau = NA_REAL, alpha = NA_REAL, laml = NA_REAL, laml_poisson = NA_REAL;
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
    // Start as glm.fit() does for the Poisson family, with no field.
    const Vec m0 = (y.array() + 0.1).matrix();
    const Vec z = (m0.array().log() - off_.array() + (y - m0).array() / m0.array()).matrix();
    const Vec sw = m0.array().sqrt().matrix();
    Eigen::ColPivHouseholderQR<Mat> qr(sw.asDiagonal() * X_);
    qr.setThreshold(ctl_.rank_tol);
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
      Vec db, du, hb, hu;
      solve(rb, ru, db, du);
      solve(Vec::Zero(p_), c_, hb, hu);
      const double chc = c_.dot(hu);
      if (!(chc > 0) || !std::isfinite(chc)) return NA_REAL;
      const double nu = c_.dot(du) / chc;
      db -= nu * hb;
      du -= nu * hu;
      // Predicted gain of the penalised log-likelihood of this Newton step.
      const double gain = 0.5 * (rb.dot(db) + ru.dot(du));
      if (!std::isfinite(gain)) return NA_REAL;
      if (gain <= ctl_.tol_inner * (std::abs(P) + 0.1)) {
        beta_ = beta;
        u_ = u;
        return loglik(mu, y_, tail_, alpha, lfact_) - 0.5 * tau * u.dot(Q_ * u) +
          0.5 * (m_ - 1) * rho - 0.5 * (logdet_S_ + logdet_G_ + std::log(chc));
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

  const Vec& beta() const { return beta_; }
  const Vec& u() const { return u_; }
  int evaluations = 0, newton = 0;

private:
  Vec beta_start_or_zero() const { return Vec::Zero(p_); }

  bool linear_predictor(const Vec& beta, const Vec& u, Vec& eta, Vec& mu) const {
    eta = off_ + A_ * u;
    if (p_) eta.noalias() += X_ * beta;
    if (!eta.allFinite() || !(eta.array() < ctl_.eta_max).all()) return false;
    mu = eta.array().exp();
    return true;
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
    logdet_G_ = 0;
    if (p_) {
      const Mat WX = w.asDiagonal() * X_;
      B_ = At_ * WX;
      SiB_ = llt_.solve(B_);
      const Mat G = X_.transpose() * WX - B_.transpose() * SiB_;
      // An aliased direction (a factor level whose fitted means have
      // vanished) is left out of the solve and of the determinant.
      Eigen::SelfAdjointEigenSolver<Mat> eig(0.5 * (G + G.transpose()));
      if (eig.info() != Eigen::Success) return false;
      const Vec ev = eig.eigenvalues();
      const double top = ev.maxCoeff();
      if (!(top > 0)) return false;
      Vec inv = Vec::Zero(p_);
      for (int j = 0; j < p_; ++j) {
        if (ev[j] > ctl_.rank_tol * top) {
          inv[j] = 1.0 / ev[j];
          logdet_G_ += std::log(ev[j]);
        }
      }
      Ginv_ = eig.eigenvectors() * inv.asDiagonal() * eig.eigenvectors().transpose();
    }
    return true;
  }

  // H^{-1} (rb, ru) by blocks.
  void solve(const Vec& rb, const Vec& ru, Vec& db, Vec& du) const {
    const Vec t = llt_.solve(ru);
    if (p_) {
      db = Ginv_ * (rb - B_.transpose() * t);
      du = t - SiB_ * db;
    } else {
      db = Vec::Zero(0);
      du = t;
    }
  }

  const SpMat& A_;
  const SpMat At_;
  const SpMat& Q_;
  const Mat& X_;
  const Vec& c_;
  const Control& ctl_;
  const int n_, m_, p_;
  Vec y_, off_, beta_, u_;
  std::vector<double> tail_;
  double lfact_ = 0, logdet_S_ = 0, logdet_G_ = 0;
  Eigen::SimplicialLDLT<SpMat> llt_;
  SpMat S_;
  std::vector<int> aa_slot_, aa_end_, q_slot_;
  std::vector<double> aa_value_, q_value_;
  Mat B_, SiB_, Ginv_;
};

// Brent's maximiser of f on [a, b]. It stops when the objective at its three
// best points agrees within tol (relative, as glm.fit() measures a deviance
// change), or after maxit evaluations. A NaN objective is the worst value.
struct LineResult { double x, f; bool converged; };

LineResult maximise(const std::function<double(double)>& f, double a, double b,
                    double tol, int maxit) {
  const double golden = 0.3819660112501051;
  auto value = [&](double x) {
    const double v = f(x);
    return std::isfinite(v) ? -v : R_PosInf;
  };
  double x = a + golden * (b - a), w = x, v = x;
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

Fit fit_feature(Model& model, const Vec& y, const Vec& offset, bool poisson,
                const Control& ctl) {
  Fit out;
  if (!(y.array() > 0).any()) {
    out.error = "All counts are zero.";
    return out;
  }
  model.set_feature(y, offset);
  // Poisson: one line search over log tau.
  LineResult r = maximise([&](double rho) { return model.laml(rho, 0.0); },
                          ctl.log_tau_min, ctl.log_tau_max, ctl.tol, ctl.maxit_line);
  if (!std::isfinite(r.f)) {
    out.error = "The Poisson spatial fit failed.";
    return out;
  }
  double rho = r.x, alpha = 0, best = r.f;
  bool converged = r.converged;
  out.laml_poisson = best;
  out.poisson = true;
  if (!poisson) {
    // The Poisson fit is the maximum unless L rises along alpha from it.
    const double lo = std::log(ctl.alpha_min), hi = std::log(ctl.alpha_max);
    converged = false;
    for (int outer = 1; outer <= ctl.maxit_outer; ++outer) {
      out.outer = outer;
      const LineResult ra = maximise(
        [&](double t) { return model.laml(rho, std::exp(t)); }, lo, hi, ctl.tol, ctl.maxit_line);
      if (!(ra.f > best + ctl.tol * (std::abs(best) + 0.1))) {
        converged = true;
        break;
      }
      alpha = std::exp(ra.x);
      out.poisson = false;
      const double before = best;
      best = ra.f;
      const LineResult rt = maximise(
        [&](double x) { return model.laml(x, alpha); },
        ctl.log_tau_min, ctl.log_tau_max, ctl.tol, ctl.maxit_line);
      if (rt.f > best) { rho = rt.x; best = rt.f; }
      if (best - before <= ctl.tol * (std::abs(best) + 0.1)) {
        converged = true;
        break;
      }
    }
  }
  // The mode at the optimum (the line searches end on other points).
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
    double tol = 1e-8,
    double tol_inner = 1e-10,
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
    int threads = 1) {
  const int n = Yt.rows(), features = Yt.cols(), p = X.cols(), m = A.cols();
  if (X.rows() != n || A.rows() != n) Rcpp::stop("X and A must have one row per observation.");
  if (Q.rows() != m || Q.cols() != m) Rcpp::stop("Q must be square with ncol(A) columns.");
  if (constraint.size() != m) Rcpp::stop("constraint must have ncol(A) entries.");
  if (offset.rows() != n || (offset.cols() != 1 && offset.cols() != features)) {
    Rcpp::stop("offset must be n x 1 or n x features.");
  }
  if (poisson.size() != features) Rcpp::stop("poisson must have one flag per feature.");
  if (threads < 1) Rcpp::stop("threads must be positive.");
  if (!(log_tau_min < log_tau_max) || !(alpha_min > 0) || !(alpha_min < alpha_max)) {
    Rcpp::stop("The search bounds are not ordered.");
  }
  const Control ctl{tol, tol_inner, rank_tol, eta_max, log_tau_min, log_tau_max,
                    alpha_min, alpha_max, maxit_inner, maxit_halving, maxit_line, maxit_outer};
  const SpMat As = A, Qs = Q;
  const Mat Xd = X;
  const Vec c = constraint;
  std::vector<int> is_pois(features);
  for (int f = 0; f < features; ++f) is_pois[f] = poisson[f] == TRUE;
  std::vector<Fit> fit(features);

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
        fit[f] = fit_feature(model, y, off, is_pois[f], ctl);
      } catch (const std::exception& e) {
        fit[f] = Fit();
        fit[f].error = e.what();
      }
    }
  }

  Mat B = Mat::Constant(p, features, NA_REAL), U = Mat::Constant(m, features, NA_REAL);
  Rcpp::NumericVector tau(features), alpha(features), laml(features), laml0(features);
  Rcpp::IntegerVector evaluations(features), newton(features), outer(features);
  Rcpp::LogicalVector converged(features), pois(features);
  Rcpp::CharacterVector error(features);
  for (int f = 0; f < features; ++f) {
    const Fit& z = fit[f];
    const bool ok = z.error.empty();
    if (ok) { B.col(f) = z.beta; U.col(f) = z.u; }
    tau[f] = z.tau; alpha[f] = z.alpha; laml[f] = z.laml; laml0[f] = z.laml_poisson;
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
    Rcpp::_["evaluations"] = evaluations, Rcpp::_["newton"] = newton,
    Rcpp::_["outer"] = outer, Rcpp::_["error"] = error);
}
