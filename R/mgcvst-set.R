#' Prepare a shared design for all mgcvST genes
#'
#' Supply a complete mgcv formula, including the SPDE term, covariates,
#' nuisance smoothers and any shared `offset()`. Parsing is delegated to mgcv.
#' Alternatively supply an externally constructed `mgcv::gam(..., fit = FALSE)`
#' setup as `G`. No gene is fitted here. The formal prediction matrix is
#' computed once and saved as `L`; estimate and test reuse this design.
#'
#' Coordinates, mesh, factor coding, covariates, smooths and observation order
#' are shared across genes. Gene-specific covariates are unsupported; additional
#' gene-specific offsets can be supplied to [mgcvST.estimate()]. To change the
#' design, construct a new setting. Prepared SPDE bases are supplied through
#' `s(x, y, bs = "spde", xt = basis)` or `bs = "spdePC"`.
#'
#' @param formula Complete two-sided mgcv formula. The response is a label and
#'   need not exist in `data`.
#' @param data One shared data frame with one row per observation.
#' @param family `mgcv::nb()` for raw counts or `stats::gaussian()` for SCT data.
#' @param G Optional external setup returned by `mgcv::gam(..., fit = FALSE)`.
#'   Supply either `formula` and `data`, or `G`.
#' @param ... Setup arguments passed to `mgcv::gam(..., fit = FALSE)`.
#' @return An `mgcvST_model` with shared `G`, `L`, geometry, and `timing` in
#'   seconds (`setup_seconds`, `lpmatrix_seconds`, `elapsed`). Each SPDE smooth
#'   also contains a `timing` environment recording `basis_seconds`,
#'   `prediction_seconds`, and corresponding call counts. Prediction calls
#'   update this environment; these counters are local to the R process.
#' @export
mgcvST.set <- function(formula = NULL, data = NULL, family = mgcv::nb(),
                       G = NULL, ...) {
  t0 <- proc.time()[["elapsed"]]
  if (is.null(G)) {
    if (!inherits(formula, "formula") || length(formula) != 3L ||
        !is.symbol(formula[[2L]])) stop("Supply a two-sided formula with a single response name.")
    if (as.character(formula[[2L]]) %in% all.vars(formula[[3L]])) {
      stop("The response cannot also be a covariate or offset.")
    }
    if (!is.data.frame(data) || !nrow(data)) stop("data must be one shared non-empty data frame.")
    data[[as.character(formula[[2L]])]] <- numeric(nrow(data))
    env <- new.env(parent = environment(formula))
    env$s <- mgcv::s
    environment(formula) <- env
    G <- mgcv::gam(formula, data = data, family = family, fit = FALSE,
                   na.action = stats::na.fail, ...)
  } else if (!is.null(formula) || !is.null(data) || !missing(family) || length(list(...))) {
    stop("Supply G alone, or formula, data, family and setup arguments.")
  }
  if (!is.list(G) || is.null(G$X) || is.null(G$mf) || is.null(G$terms) ||
      is.null(G$family) || is.null(G$smooth) ||
      nrow(G$mf) != nrow(G$X) || length(G$y) != nrow(G$X)) {
    stop("G must be a complete mgcv::gam(..., fit = FALSE) setup.")
  }
  response <- attr(G$terms, "response")
  if (length(response) != 1L || response < 1L || response > ncol(G$mf)) {
    stop("G must contain a single response column.")
  }
  if (names(G$mf)[response] %in% all.vars(G$formula[[3L]])) {
    stop("The response cannot also be a covariate or offset.")
  }
  if (!(G$family$family == "gaussian" || grepl("^negative binomial", tolower(G$family$family)))) {
    stop("mgcvST.set() supports negative binomial and Gaussian families.")
  }
  if (length(G$paraPen) || isTRUE(G$n.paraPen > 0) ||
      (!is.null(G$H) && any(G$H != 0))) {
    stop("mgcvST.set() does not support cross-penalties or extra coefficient penalties.")
  }
  if (length(G$term.names) != ncol(G$X) || anyNA(G$term.names) ||
      any(!nzchar(G$term.names)) || anyDuplicated(G$term.names)) {
    stop("G must contain one unique coefficient name per design column.")
  }
  for (j in seq_along(G$smooth)) {
    sm <- G$smooth[[j]]
    if (inherits(sm, "spde.smooth") || inherits(sm, "spdePC.smooth")) {
      if (is.null(sm$score.component)) sm$score.component <- sm$component
      G$smooth[[j]] <- sm
    }
  }
  pseudo <- G
  pseudo$model <- G$mf
  pseudo$coefficients <- stats::setNames(numeric(ncol(G$X)), G$term.names)
  pseudo$linear.predictors <- numeric(nrow(G$X))
  class(pseudo) <- c("gam", "glm", "lm")
  setup_seconds <- proc.time()[["elapsed"]] - t0
  t1 <- proc.time()[["elapsed"]]
  L <- .gam_training_lpmatrix(pseudo)
  if (!identical(dim(L), dim(G$X)) || !identical(colnames(L), G$term.names) ||
      !isTRUE(all.equal(as.numeric(L), as.numeric(G$X), tolerance = 1e-10))) {
    stop("The prepared design and formal lpmatrix disagree; rebuild G at the supplied coordinates.")
  }
  lpmatrix_seconds <- proc.time()[["elapsed"]] - t1
  geometry <- .mgcvst_model_geometry(pseudo, L)
  components <- geometry$score_components
  structure(list(
    G = G, L = L, geometry = geometry, shared_design = TRUE,
    setting = if (length(components) == 1L) "global" else "global_local",
    components = components, formula = G$formula, internal_formula = G$formula,
    response = names(G$mf)[response],
    kappa = stats::setNames(vapply(geometry$target, function(j) G$smooth[[j]]$kappa,
                                  numeric(1L)), components),
    offset = geometry$offset,
    timing = list(setup_seconds = setup_seconds, lpmatrix_seconds = lpmatrix_seconds,
                  elapsed = proc.time()[["elapsed"]] - t0)
  ), class = "mgcvST_model")
}
