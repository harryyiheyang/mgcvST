#ifndef SPDEST_NB_FAMILY_H
#define SPDEST_NB_FAMILY_H

// NB2 (Var = mu + alpha mu^2) log-likelihood pieces shared by the null fits
// (nb_null.cpp) and the spatial fits (spde_estimate.cpp). alpha = 0 is the
// Poisson model.

#include <RcppEigen.h>
#include <cmath>
#include <vector>

namespace mgcvst_nb {

typedef Eigen::VectorXd Vec;

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

// Tail-count sums sum_j c_j j^m / (1 + alpha j)^m for m = 1, 2 and the log
// term sum_j c_j log1p(alpha j).
inline void tail_sums(const std::vector<double>& c, double alpha, double& s1,
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
inline double loglik(const Vec& mu, const Vec& y, const std::vector<double>& c,
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


} // namespace mgcvst_nb

#endif
