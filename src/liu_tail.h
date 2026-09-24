#ifndef MGCVST_LIU_TAIL_H
#define MGCVST_LIU_TAIL_H

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>

namespace mgcvst_liu {

inline double log_add(double a, double b) {
  if (a < b) std::swap(a, b);
  if (b == -std::numeric_limits<double>::infinity()) return a;
  return a + std::log1p(std::exp(b - a));
}

// log P(chi2'(df, ncp) > x). Central tails use R's log pchisq. Noncentral
// tails are the Poisson mixture sum_i w_i Q(df / 2 + i, x / 2) in log space:
// Q increases in i by Q(a + 1, y) = Q(a, y) + y^a e^-y / Gamma(a + 1) (a gamma
// density at y, R's dgamma), so one log pgamma call at
// i0 = max(0, floor(lambda - 12 sqrt(lambda))) starts an
// upward recursion. Terms below i0 are bounded by Q(i0) P(Pois(lambda) < i0);
// the log-concave terms are summed past their mode until they fall below
// exp(-45) of the running sum.
inline double log_upper_nchisq(double x, double df, double ncp) {
  if (std::isnan(x) || std::isnan(df) || std::isnan(ncp)) return x + df + ncp;
  if (!(df > 0) || !(ncp >= 0)) return NAN;
  if (ncp == 0) return R::pchisq(x, df, 0, 1);
  if (x <= 0) return 0;
  const double lambda = 0.5 * ncp, y = 0.5 * x, h = 0.5 * df;
  const double start = std::floor(lambda - 12 * std::sqrt(lambda));
  const long i0 = start > 0 ? (long)start : 0L;
  double lQ = R::pgamma(y, h + i0, 1.0, 0, 1);
  double lsum = R::dpois((double)i0, lambda, 1) + lQ, prev = lsum;
  for (long i = i0 + 1; i < i0 + 10000000L; ++i) {
    lQ = log_add(lQ, R::dgamma(y, h + i, 1.0, 1));
    const double term = R::dpois((double)i, lambda, 1) + lQ;
    lsum = log_add(lsum, term);
    if (i > lambda && term < prev && term < lsum - 45) break;
    prev = term;
  }
  return lsum;
}

// Liu moment matching of the squared bilinear score (R/score-test.R) with a
// log-space tail; writes log p two-sided, positive, negative.
inline void liu_log_p(double U, double A, double B, double C, double D,
                      double* out) {
  const double c1 = A;
  const double c2 = A * A + 3 * B;
  const double c3 = R_pow(A, 3.0) + 9 * A * B + 15 * C;
  const double c4 = R_pow(A, 4.0) + 18 * (A * A) * B + 60 * A * C +
    24 * (B * B) + 105 * D;
  const double s1 = c3 / R_pow(c2, 1.5);
  const double s2 = c4 / (c2 * c2);
  const double tstar = (U * U - c1) / std::sqrt(2 * c2);
  double a = 1 / s1, delta = 0, df = 1 / (s1 * s1);
  if (s1 * s1 > s2) {
    a = 1 / (s1 - std::sqrt(s1 * s1 - s2));
    delta = s1 * R_pow(a, 3.0) - a * a;
    df = a * a - 2 * delta;
  }
  const double x = tstar * (std::sqrt(2.0) * a) + (df + delta);
  double lp = log_upper_nchisq(x, df, delta);
  if (lp > 0) lp = 0;
  const double half = lp - M_LN2;
  const double other = std::log1p(-0.5 * std::exp(lp));
  out[0] = lp;
  out[1] = U >= 0 ? half : other;
  out[2] = U <= 0 ? half : other;
  if (std::isnan(lp)) out[1] = out[2] = lp;
}

}  // namespace mgcvst_liu

#endif
