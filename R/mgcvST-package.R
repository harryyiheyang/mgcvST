#' SPDE smoothing and covariance testing for spatial features
#'
#' `mgcvST` combines topology-aware boundary learning and mesh construction, a proper
#' Matérn SPDE smooth for `mgcv`, compact marginal estimation, and finite-rank
#' RKHS covariance score testing.
#'
#' @importFrom mgcv Predict.matrix smooth.construct
#' @importFrom Rcpp evalCpp
#' @useDynLib mgcvST, .registration = TRUE
#' @keywords internal
"_PACKAGE"
