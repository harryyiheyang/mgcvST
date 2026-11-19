#include "spa_quadform.h"

// Unit-test entry of the Stage 1 saddlepoint (spa_quadform.h): the upper tail
// of q against sum(lambda_k chi^2_1).
// [[Rcpp::export]]
double mgcvst_spa_quadform_cpp(double q, const std::vector<double>& lambda,
                               double centre_tol = 1e-6, double root_tol = 1e-12,
                               int root_maxit = 100) {
  mgcvst_spa::QuadformSpaControl ctl;
  ctl.centre_tol = centre_tol;
  ctl.root_tol = root_tol;
  ctl.root_maxit = root_maxit;
  return mgcvst_spa::quadform_upper(q, lambda, ctl);
}
