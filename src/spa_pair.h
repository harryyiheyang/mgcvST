#ifndef MGCVST_SPA_PAIR_H
#define MGCVST_SPA_PAIR_H

// Saddlepoint calibration of the signed cross-gene score, shared by the exact
// and the PCAlearning pair routes.
//
// Under the null the score is U = sum_i s_i x_i y_i with x_i, y_i iid N(0, 1),
// where s_i are the singular values of H_1^{1/2} H_2^{1/2}. The spectrum is
// represented by k leading values s_1, ..., s_k (from the k x k compression of
// the product on a shared basis) and by a remainder that matches the remaining
// power sums mu_r = sum_{i > k} s_i^(2r) = t_r - sum_{i <= k} s_i^(2r):
//   * two nodes   matches mu_1..mu_4 (exact route, four trace moments),
//   * one node    matches mu_1, mu_2 (Satterthwaite type; PCAlearning route and
//                 fallback of the two-node form),
//   * Gaussian    variance mu_1, when mu_2 is not positive,
//   * none        when mu_1 is not positive.
// A node (u, m) stands for m singular values s = sqrt(u).
//
// Cumulant generating function in t:
//   K(t) = -1/2 sum_i log(1 - s_i^2 t^2) - 1/2 sum_n m_n log(1 - u_n t^2)
//          + g t^2 / 2,
// saddlepoint K'(t) = x, w = sqrt(2 (t x - K)), v = t sqrt(K''), and the
// Lugannani-Rice tail P(U > x) = Phibar(w) + phi(w) (1 / v - 1 / w), computed in
// log space (the r* form when the bracket is not positive). The tail is
// evaluated at x = |U|; the two one-sided tails follow from the symmetry of U.
//
// Everything is scaled by the largest singular value s* (tau = t s*,
// D_i = 1 - s_i'^2 tau^2 = delta_i + s_i'^2 eta, eta = 1 - tau^2 = eps (1 + tau)),
// so that the root of K' = x is found without cancellation near tau = 1.
//
// This header is thread safe: it uses only arithmetic, R::pnorm and R::dnorm.

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace mgcvst_spa {

constexpr double kLn2 = 0.693147180559945309417232121458;
constexpr double kLogSqrt2Pi = 0.918938533204672741780329736406;

enum RemainderKind { REM_NONE = 0, REM_ONE = 1, REM_TWO = 2, REM_GAUSS = 3 };

struct Remainder {
  int kind = REM_NONE;
  int nodes = 0;
  double u[2] = {0.0, 0.0};  // node values u = s_r^2
  double m[2] = {0.0, 0.0};  // node counts
  double g = 0.0;            // Gaussian variance
};

// Power sums sum_i s_i^(2r), r = 1..4, of the leading values in long double.
inline void leading_sums(const double* s, int k, long double* lead) {
  lead[0] = lead[1] = lead[2] = lead[3] = 0.0L;
  for (int i = 0; i < k; ++i) {
    if (!(s[i] > 0)) continue;
    const long double u = static_cast<long double>(s[i]) * static_cast<long double>(s[i]);
    const long double u2 = u * u;
    lead[0] += u;
    lead[1] += u2;
    lead[2] += u2 * u;
    lead[3] += u2 * u2;
  }
}

// Remainder from the total power sums t[0..3] (t_r = tr((H_1 H_2)^r) in the
// units of s) and the leading sums. order = 4 tries two nodes first; order = 2
// starts at one node. A remainder moment smaller than 1e-10 t_r in absolute
// value is zero.
inline Remainder make_remainder(const double* t, const long double* lead, int order) {
  long double mu[4];
  for (int r = 0; r < 4; ++r) {
    mu[r] = static_cast<long double>(t[r]) - lead[r];
    if (std::fabs(mu[r]) < 1e-10L * std::fabs(static_cast<long double>(t[r]))) mu[r] = 0.0L;
  }
  Remainder R;
  if (!(mu[0] > 0)) return R;
  if (order >= 4 && mu[1] > 0 && mu[2] > 0 && mu[3] > 0) {
    // Gauss nodes of the measure with moments mu_1..mu_4: the roots of
    // det x^2 + A1 x + A0, weights from sum w = mu_1, sum w x = mu_2.
    const long double det = mu[0] * mu[2] - mu[1] * mu[1];
    const long double A1 = -(mu[0] * mu[3] - mu[1] * mu[2]);
    const long double A0 = mu[1] * mu[3] - mu[2] * mu[2];
    const long double disc = A1 * A1 - 4.0L * det * A0;
    if (det > 0 && disc >= 0) {
      const long double sq = std::sqrt(disc);
      const long double qv = -0.5L * (A1 + (A1 >= 0 ? sq : -sq));
      if (qv != 0) {
        long double x1 = qv / det, x2 = A0 / qv;
        if (x1 > x2) std::swap(x1, x2);
        if (x1 > 0 && x2 > 0 && (x2 - x1) > 1e-8L * x2) {
          const long double w2 = (mu[1] - x1 * mu[0]) / (x2 - x1);
          const long double w1 = mu[0] - w2;
          if (w1 > 0 && w2 > 0) {
            R.kind = REM_TWO;
            R.nodes = 2;
            R.u[0] = static_cast<double>(x1);
            R.u[1] = static_cast<double>(x2);
            R.m[0] = static_cast<double>(w1 / x1);
            R.m[1] = static_cast<double>(w2 / x2);
            return R;
          }
        }
      }
    }
  }
  if (mu[1] > 0) {
    R.kind = REM_ONE;
    R.nodes = 1;
    R.u[0] = static_cast<double>(mu[1] / mu[0]);
    R.m[0] = static_cast<double>(mu[0] * mu[0] / mu[1]);
    return R;
  }
  R.kind = REM_GAUSS;
  R.g = static_cast<double>(mu[0]);
  return R;
}

// TRUE when a remainder node lies above the largest leading value, u > s_1^2:
// the compression has then not captured the top of the spectrum. The kernels
// count such pairs as a diagnostic.
inline bool node_above_leading(const Remainder& R, double s1) {
  for (int n = 0; n < R.nodes; ++n) if (R.u[n] > s1 * s1) return true;
  return false;
}

// Per-thread scratch for the normalized spectrum.
struct Scratch {
  std::vector<double> a, d, mult;
  void reserve(int n) {
    if (static_cast<int>(a.size()) < n + 2) {
      a.assign(n + 2, 0.0);
      d.assign(n + 2, 1.0);
      mult.assign(n + 2, 1.0);
    }
  }
};

struct Eval {
  double f = 0, fp = 0, K = 0, K2 = 0;
};

// f(tau) = K'(tau) - x' and its tau-derivative, in the scaled variables. a, d:
// squared normalized values and 1 - a; the nodes follow in the same arrays.
// eps = 1 - tau is carried separately so that 1 - tau^2 = eps (1 + tau) keeps
// its relative accuracy.
inline void eval_point(const double* a, const double* d, const double* mult, int n,
                       double gp, double xp, double tau, double eps, bool full,
                       Eval& e) {
  const double eta = eps * (1.0 + tau);
  double sumA = 0, sumA2 = 0, sum2 = 0, sumLog = 0;
  for (int i = 0; i < n; ++i) {
    const double ai = a[i];
    if (ai <= 0) continue;
    const double D = d[i] + ai * eta;
    const double r = mult[i] * (ai / D);
    sumA += r;
    sumA2 += r * (ai / D);
    if (full) {
      sum2 += r * (2.0 - D) / D;
      sumLog += mult[i] * std::log(D);
    }
  }
  const double B = sumA + gp;
  e.f = tau * B - xp;
  e.fp = B + 2.0 * tau * tau * sumA2;
  if (full) {
    e.K = -0.5 * sumLog + 0.5 * gp * tau * tau;
    e.K2 = sum2 + gp;
  }
}

// log P(U > x), x > 0. Returns false for invalid input. `scratch` is resized on
// first use only.
inline bool spa_log_tail(double x, const double* s, int k, const Remainder& R,
                         Scratch& scratch, double* log_p) {
  if (!(x > 0) || !std::isfinite(x)) return false;
  // Largest squared singular value.
  double s2max = 0;
  for (int i = 0; i < k; ++i) {
    if (!std::isfinite(s[i])) return false;
    if (s[i] > 0) s2max = std::max(s2max, s[i] * s[i]);
  }
  for (int n = 0; n < R.nodes; ++n) s2max = std::max(s2max, R.u[n]);
  if (!(s2max > 0)) {
    // No singular part: the score is Gaussian with variance g.
    if (!(R.g > 0)) return false;
    *log_p = R::pnorm(x / std::sqrt(R.g), 0.0, 1.0, 0, 1);
    return std::isfinite(*log_p);
  }
  const double smax = std::sqrt(s2max);
  const int n = k + R.nodes;
  scratch.reserve(n);
  double* a = scratch.a.data();
  double* d = scratch.d.data();
  double* mult = scratch.mult.data();
  for (int i = 0; i < k; ++i) {
    const double si = s[i] > 0 ? s[i] : 0.0;
    a[i] = si * si / s2max;
    d[i] = (s2max - si * si) / s2max;
    mult[i] = 1.0;
  }
  for (int j = 0; j < R.nodes; ++j) {
    a[k + j] = R.u[j] / s2max;
    d[k + j] = (s2max - R.u[j]) / s2max;
    mult[k + j] = R.m[j];
  }
  const double gp = R.g / s2max;
  const double xp = x / smax;
  if (!std::isfinite(xp)) return false;

  double B0 = gp;
  for (int i = 0; i < n; ++i) B0 += mult[i] * a[i];
  Eval e;
  double tau, eps;
  // Which half of (0, 1) holds the root?
  eval_point(a, d, mult, n, gp, xp, 0.5, 0.5, false, e);
  if (e.f >= 0) {
    // tau in (0, 0.5]: Newton in tau from the right, monotone (f is convex).
    double lo = 0.0, hi = std::min(xp / B0, 0.5);
    tau = hi;
    for (int it = 0; it < 200; ++it) {
      eps = 1.0 - tau;
      eval_point(a, d, mult, n, gp, xp, tau, eps, false, e);
      if (e.f > 0) hi = tau; else if (e.f < 0) lo = tau; else break;
      double next = tau - e.f / e.fp;
      if (!(next > lo && next < hi)) next = 0.5 * (lo + hi);
      if (std::fabs(next - tau) <= 1e-15 * std::max(tau, 1e-300)) { tau = next; break; }
      tau = next;
    }
    eps = 1.0 - tau;
  } else {
    // tau in (0.5, 1): iterate on eps = 1 - tau. The terms that attain the
    // largest squared value (a = 1, total multiplicity m*) give
    // K' >= m* tau / (eps (2 - eps)), hence the lower bracket
    // eps_L = 2 m* / (2 x' + m* + sqrt(4 x'^2 + m*^2)). A remainder node with
    // a multiplicity below 1 can be the only term with a = 1.
    double mstar = 0.0;
    for (int i = 0; i < n; ++i) if (a[i] == 1.0) mstar += mult[i];
    if (!(mstar > 0)) mstar = 1.0;
    double lo = 2.0 * mstar / (2.0 * xp + mstar + std::sqrt(4.0 * xp * xp + mstar * mstar)),
           hi = 0.5;
    if (!(lo > 0)) return false;
    eps = lo;
    for (int it = 0; it < 200; ++it) {
      tau = 1.0 - eps;
      eval_point(a, d, mult, n, gp, xp, tau, eps, false, e);
      if (e.f > 0) lo = eps; else if (e.f < 0) hi = eps; else break;
      double next = eps + e.f / e.fp;
      if (!(next > lo && next < hi)) next = (hi > 4.0 * lo) ? std::sqrt(lo * hi) : 0.5 * (lo + hi);
      if (std::fabs(next - eps) <= 1e-15 * eps) { eps = next; break; }
      eps = next;
    }
    tau = 1.0 - eps;
  }
  eval_point(a, d, mult, n, gp, xp, tau, eps, true, e);
  // The root must satisfy the saddlepoint equation: a bracket or iteration
  // that failed returns status 2 rather than a wrong p-value.
  if (!(std::fabs(e.f) <= 1e-10 * xp)) return false;
  const double w2 = 2.0 * (tau * xp - e.K);
  const double w = w2 > 0 ? std::sqrt(w2) : 0.0;
  if (!(e.K2 > 0) || !std::isfinite(e.K2) || !std::isfinite(w)) return false;
  if (w < 1e-3) {
    *log_p = -kLn2;  // two-sided p = 1
    return true;
  }
  const double v = tau * std::sqrt(e.K2);
  const double lPhi = R::pnorm(w, 0.0, 1.0, 0, 1);
  const double lphi = -0.5 * w * w - kLogSqrt2Pi;
  const double term = 1.0 / v - 1.0 / w;
  const double bracket = 1.0 + std::exp(lphi - lPhi) * term;
  double lp;
  if (bracket > 0 && std::isfinite(bracket)) {
    lp = lPhi + std::log(bracket);
  } else {
    // Barndorff-Nielsen r* form.
    const double rstar = w + std::log(v / w) / w;
    lp = R::pnorm(rstar, 0.0, 1.0, 0, 1);
  }
  if (!std::isfinite(lp)) return false;
  *log_p = std::min(lp, -kLn2);
  return true;
}

// Natural-log p-values (two-sided, positive, negative) of the signed score U.
// Returns the status: 0 evaluated, 2 invalid p-value.
inline int spa_pair(double U, const double* s, int k, const Remainder& R,
                    Scratch& scratch, double* lp) {
  const double x = std::fabs(U);
  if (!std::isfinite(U)) return 2;
  if (x == 0) {
    lp[0] = 0.0;
    lp[1] = lp[2] = -kLn2;
    return 0;
  }
  double one;
  if (!spa_log_tail(x, s, k, R, scratch, &one)) return 2;
  lp[0] = std::min(0.0, one + kLn2);
  const double other = std::log1p(-std::exp(one));
  lp[1] = U >= 0 ? one : other;
  lp[2] = U <= 0 ? one : other;
  return 0;
}

}  // namespace mgcvst_spa

#endif
