#' SPDE smoothing and covariance testing for spatial features (mgcv and INLA implementations)
#'
#' `spdeST` combines topology-aware boundary learning and mesh construction, a proper
#' Matérn SPDE smooth for `mgcv`, compact marginal estimation, and finite-rank
#' RKHS covariance score testing.
#'
#' @importFrom mgcv Predict.matrix smooth.construct
#' @importFrom Rcpp evalCpp
#' @importFrom utils object.size
#' @useDynLib spdeST, .registration = TRUE
#' @keywords internal
"_PACKAGE"
