#ifndef MGCVST_SPA_QUADFORM_H
#define MGCVST_SPA_QUADFORM_H

#include <Rcpp.h>
#include <cmath>
#include <vector>

// Saddlepoint upper tail of Q = sum_k lambda_k chi^2_1 (Kuonen 1999,
// Biometrika 86:929, Barndorff-Nielsen form): the Stage 1 p-value of
// spde.scan. K(t) = -1/2 sum log(1 - 2 lambda t) on t < 1 / (2 max lambda);
// the root of K'(t) = q is found by Newton steps kept inside a bisection
// bracket. Near the mean (|q - mean| <= centre_tol * mean) the formula is
// 0/0 and its limit 1/2 - rho3 / (6 sqrt(2 pi)) is used. Returns NA when
// lambda has no positive entry or q is not finite.
namespace mgcvst_spa {

struct QuadformSpaControl {
  double centre_tol = 1e-6;
  double root_tol = 1e-12;
  int root_maxit = 100;
};

inline double quadform_upper(double q, const std::vector<double>& lambda_in,
                             const QuadformSpaControl& ctl) {
  std::vector<double> lambda;
  double lmax = 0, k1 = 0, k2 = 0, k3 = 0;
  for (double l : lambda_in) {
    if (l > 0 && std::isfinite(l)) {
      lambda.push_back(l);
      lmax = std::max(lmax, l);
      k1 += l;
      k2 += l * l;
      k3 += l * l * l;
    }
  }
  if (lambda.empty() || !std::isfinite(q)) return NA_REAL;
  if (q <= 0) return 1.0;
  const double centre =
    0.5 - 8.0 * k3 / std::pow(2.0 * k2, 1.5) / (6.0 * std::sqrt(2.0 * M_PI));
  if (std::abs(q - k1) <= ctl.centre_tol * k1) return centre;
  auto K1 = [&](double t) {
    double s = 0;
    for (double l : lambda) s += l / (1.0 - 2.0 * l * t);
    return s;
  };
  auto K2 = [&](double t) {
    double s = 0;
    for (double l : lambda) {
      const double d = 1.0 - 2.0 * l * t;
      s += 2.0 * l * l / (d * d);
    }
    return s;
  };
  // Bracket [lo, hi] with K'(lo) <= q <= K'(hi); K' increases in t.
  double lo, hi;
  if (q > k1) {
    lo = 0;
    hi = 1.0 / (2.0 * lmax);
  } else {
    hi = 0;
    lo = -1.0 / (2.0 * lmax);
    while (K1(lo) > q) lo *= 2.0;
  }
  double t = 0.5 * (lo + hi);
  for (int it = 0; it < ctl.root_maxit; ++it) {
    const double g = K1(t) - q;
    if (std::abs(g) <= ctl.root_tol * q) break;
    if (g > 0) hi = t; else lo = t;
    const double next = t - g / K2(t);
    t = (next > lo && next < hi) ? next : 0.5 * (lo + hi);
  }
  double K = 0;
  for (double l : lambda) K += -0.5 * std::log1p(-2.0 * l * t);
  const double w = (t > 0 ? 1.0 : -1.0) * std::sqrt(std::max(2.0 * (t * q - K), 0.0));
  const double v = t * std::sqrt(K2(t));
  if (!std::isfinite(w) || !std::isfinite(v) || w == 0 || v / w <= 0) return centre;
  return R::pnorm(w + std::log(v / w) / w, 0.0, 1.0, 0, 0);
}

} // namespace mgcvst_spa

#endif
