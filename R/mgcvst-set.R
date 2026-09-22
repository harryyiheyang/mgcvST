#' Prepare a shared design for all mgcvST genes
#'
#' Supplies a complete mgcv formula or an external mgcv setup. Shared full and
#' null GAM setups are prepared with a zero response; each feature is later
#' fitted independently with BAM.
#'
#' @param formula Complete two-sided mgcv formula. The response is a label and
#'   need not exist in `data`; mgcv expands factors, interactions and contrasts.
#' @param data One shared data frame with one row per observation.
#' @param family `mgcv::nb()` for raw counts or `stats::gaussian()` for SCT data.
#' @param G Optional external setup returned by `mgcv::gam(..., fit = FALSE)`.
#'   Supply either `formula` and `data`, or `G`.
#' @param ... Setup arguments passed to `mgcv::gam(..., fit = FALSE)`.
#' @return An `mgcvST_model` containing the shared full and null GAM setups,
#'   their training designs, the target spatial geometry and setup timing.
#' @export
mgcvST.set <- function(formula = NULL, data = NULL, family = mgcv::nb(),
                       G = NULL, ...) {
  args <- c(list(formula = formula, data = data, G = G), list(...))
  if (!missing(family)) args$family <- family
  do.call(.mgcvst_set_prepare, args)
}

.mgcvst_formula_terms <- function(x) {
  if (is.call(x) && identical(x[[1L]], as.name("+"))) {
    c(.mgcvst_formula_terms(x[[2L]]), .mgcvst_formula_terms(x[[3L]]))
  } else list(x)
}

.mgcvst_rebuild_formula <- function(formula, terms) {
  rhs <- if (length(terms)) Reduce(function(x, y) call("+", x, y), terms) else 1
  stats::as.formula(call("~", formula[[2L]], rhs), env = environment(formula))
}

.mgcvst_null_formula <- function(formula, target_index) {
  smooth_index <- 0L
  keep <- Filter(function(x) {
    smooth <- is.call(x) && as.character(x[[1L]]) %in% c("s", "te", "ti", "t2")
    if (smooth) smooth_index <<- smooth_index + 1L
    bs <- if (smooth) x[["bs"]] else NULL
    !(smooth && smooth_index == target_index && !is.null(bs) &&
      as.character(bs) %in% c("spde", "spdePC"))
  }, .mgcvst_formula_terms(formula[[3L]]))
  .mgcvst_rebuild_formula(formula, keep)
}

.mgcvst_external_spec <- function(G) {
  data <- as.data.frame(G$mf)
  terms <- .mgcvst_formula_terms(G$formula[[3L]])
  has_offset <- vapply(terms, function(x) {
    is.call(x) && identical(x[[1L]], as.name("offset"))
  }, logical(1L))
  data$.mgcvST_setup_offset <- if (is.null(G$offset)) numeric(nrow(data)) else G$offset
  terms <- c(terms[!has_offset], list(call("offset", as.name(".mgcvST_setup_offset"))))
  formula <- .mgcvst_rebuild_formula(G$formula, terms)
  env <- new.env(parent = environment(formula))
  env$s <- mgcv::s
  calls <- Filter(function(x) is.call(x) && identical(x[[1L]], as.name("s")), terms)
  for (j in seq_len(min(length(calls), length(G$smooth)))) {
    xt <- calls[[j]][["xt"]]
    if (is.symbol(xt)) assign(as.character(xt), G$smooth[[j]]$xt, envir = env)
  }
  environment(formula) <- env
  list(formula = formula, data = data)
}

.mgcvst_mark_global <- function(G) {
  for (j in seq_along(G$smooth)) {
    sm <- G$smooth[[j]]
    if (inherits(sm, "spde.smooth") || inherits(sm, "spdePC.smooth")) {
      if (is.null(sm$score.component)) sm$score.component <- sm$component
      G$smooth[[j]] <- sm
    }
  }
  G
}

.mgcvst_set_prepare <- function(formula = NULL, data = NULL, family = mgcv::nb(),
                               G = NULL, ..., .allow_poisson = FALSE) {
  t0 <- proc.time()[["elapsed"]]
  external <- !is.null(G)
  if (external) {
    if (!is.null(formula) || !is.null(data) || !missing(family) || length(list(...))) {
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
    if (!(G$family$family == "gaussian" ||
          (.allow_poisson && G$family$family == "poisson") ||
          grepl("^negative binomial", tolower(G$family$family)))) {
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
    public_formula <- G$formula
    G <- .mgcvst_mark_global(G)
    spec <- .mgcvst_external_spec(G)
    formula <- spec$formula
    data <- spec$data
    family <- G$family
  } else {
    if (!inherits(formula, "formula") || length(formula) != 3L ||
        !is.symbol(formula[[2L]])) stop("Supply a two-sided formula with a single response name.")
    if (as.character(formula[[2L]]) %in% all.vars(formula[[3L]])) {
      stop("The response cannot also be a covariate or offset.")
    }
    if (!is.data.frame(data) || !nrow(data)) stop("data must be one shared non-empty data frame.")
    public_formula <- formula
    data <- as.data.frame(data)
    data[[as.character(formula[[2L]])]] <- numeric(nrow(data))
  }
  full_G <- mgcv::gam(formula, data = data, family = family, fit = FALSE,
                      na.action = stats::na.fail, ...)
  full_G <- .mgcvst_mark_global(full_G)
  target_index <- which(vapply(full_G$smooth,
    function(s) identical(s$score.component, "global"), logical(1L)))
  null_formula <- .mgcvst_null_formula(formula, target_index)
  null_G <- mgcv::gam(null_formula, data = data, family = family, fit = FALSE,
                      na.action = stats::na.fail, ...)
  L <- as.matrix(full_G$X)
  attr(full_G, "training_X") <- L
  pseudo <- full_G
  pseudo$model <- full_G$mf
  pseudo$coefficients <- stats::setNames(numeric(ncol(full_G$X)), full_G$term.names)
  pseudo$linear.predictors <- numeric(nrow(L))
  class(pseudo) <- c("gam", "glm", "lm")
  geometry <- .mgcvst_model_geometry(pseudo, L)
  response <- attr(full_G$terms, "response")
  structure(list(
    G = full_G, L = L, geometry = geometry, shared_design = TRUE,
    setting = "global", components = geometry$score_components,
    formula = public_formula, internal_formula = formula,
    full_formula = formula, full_data = data,
    response = names(full_G$mf)[response], null_formula = null_formula,
    null_data = data, null_X = as.matrix(null_G$X),
    null_response = names(full_G$mf)[response],
    kappa = stats::setNames(vapply(geometry$target, function(j) full_G$smooth[[j]]$kappa,
                                  numeric(1L)), geometry$score_components),
    offset = geometry$offset,
    timing = list(setup_seconds = proc.time()[["elapsed"]] - t0,
                  lpmatrix_seconds = 0, elapsed = proc.time()[["elapsed"]] - t0)
  ), class = "mgcvST_model")
}
