#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

// Multiple-testing adjustment of natural-log p-values, entirely in log space.
// A tail below the double range (log p = -5000, say) keeps its ordering and
// its adjusted value; nothing is exponentiated.
//
//   BH     Benjamini-Hochberg step-up: q_(i) = min_{k >= i} min(1, m p_(k) / k)
//   BY     Benjamini-Yekutieli (2001) under arbitrary dependence: BH with m
//          replaced by m c(m), c(m) = sum_{i = 1}^m 1 / i. This is
//          stats::p.adjust(p, "BY"), whose c(m) is summed in the same order.
//   Sidak  single-step Sidak: q = 1 - (1 - p)^m
//   none   q = p
//
// m counts the finite log p-values (-Inf, a p-value of zero, included).
// NA, NaN and +Inf are returned as NA and do not enter the ranking. A log p
// above 0 is a rounding artefact and is treated as 0.

namespace {

inline double log_one_minus_exp(double x) {  // log(1 - exp(x)) for x <= 0
  return x > -M_LN2 ? std::log(-std::expm1(x)) : std::log1p(-std::exp(x));
}

inline bool usable(double x) { return !std::isnan(x) && x < INFINITY; }

template <class Index>
void step_up(const double* lp, double* out, R_xlen_t n, double log_c) {
  std::vector<Index> order;
  order.reserve(static_cast<size_t>(n));
  for (R_xlen_t i = 0; i < n; ++i) {
    if (usable(lp[i])) order.push_back(static_cast<Index>(i));
  }
  const double m = static_cast<double>(order.size());
  std::sort(order.begin(), order.end(), [lp](Index a, Index b) {
    const double x = std::min(lp[a], 0.0), y = std::min(lp[b], 0.0);
    return x < y;
  });
  double run = 0.0;
  for (R_xlen_t k = static_cast<R_xlen_t>(order.size()); k >= 1; --k) {
    const double x = std::min(lp[order[k - 1]], 0.0);
    const double raw = x + std::log(m / static_cast<double>(k)) + log_c;
    if (raw < run) run = raw;
    out[order[k - 1]] = run;
  }
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List mgcvst_log_adjust_cpp(const Rcpp::NumericVector& log_p,
                                 const std::string& method) {
  const bool bh = method == "BH", by = method == "BY";
  const bool sidak = method == "Sidak", none = method == "none";
  if (!bh && !by && !sidak && !none) {
    Rcpp::stop("method must be one of \"BY\", \"BH\", \"Sidak\" or \"none\".");
  }
  const R_xlen_t n = log_p.size();
  const double* lp = log_p.begin();
  Rcpp::NumericVector out(n, NA_REAL);
  double* q = out.begin();
  R_xlen_t m = 0;
  for (R_xlen_t i = 0; i < n; ++i) m += usable(lp[i]);

  if (none) {
    for (R_xlen_t i = 0; i < n; ++i) {
      if (usable(lp[i])) q[i] = std::min(lp[i], 0.0);
    }
  } else if (sidak) {
    for (R_xlen_t i = 0; i < n; ++i) {
      if (!usable(lp[i])) continue;
      const double x = std::min(lp[i], 0.0);
      // Below exp(-700) the product m p is below 1e-292 for any m that fits
      // in memory, so 1 - (1 - p)^m = m p to the working precision, whereas
      // log(1 - p) already loses digits in the subnormal range just above.
      if (x < -700.0) {
        q[i] = x + std::log(static_cast<double>(m));
        continue;
      }
      const double a = log_one_minus_exp(x);  // log(1 - p)
      q[i] = log_one_minus_exp(static_cast<double>(m) * a);
    }
  } else if (m > 0) {
    double log_c = 0.0;
    if (by) {
      long double c = 0.0L;
      for (R_xlen_t i = 1; i <= m; ++i) c += static_cast<double>(1.0 / static_cast<double>(i));
      log_c = std::log(static_cast<double>(c));
    }
    if (n <= static_cast<R_xlen_t>(std::numeric_limits<int32_t>::max())) {
      step_up<int32_t>(lp, q, n, log_c);
    } else {
      step_up<int64_t>(lp, q, n, log_c);
    }
  }
  return Rcpp::List::create(Rcpp::Named("log_q") = out,
                            Rcpp::Named("n") = static_cast<double>(m));
}
