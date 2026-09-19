# Add one language object to the right-hand side of a formula.
.mgcvst_formula_add <- function(formula, term, response) {
  rhs <- formula[[3L]]
  out <- call("~", response, call("+", rhs, term))
  stats::as.formula(out, env = environment(formula))
}

# Construct one prepared-basis SPDE smooth call.
.mgcvst_spde_call <- function(coordinates, basis_name) {
  as.call(list(
    as.name("s"), as.name(coordinates[1L]), as.name(coordinates[2L]),
    bs = "spde",
    xt = as.name(basis_name)
  ))
}

#' Construct a reusable mgcvST model setting
#'
#' Builds the reusable `mgcv::gam(..., fit = FALSE)` setup consumed by
#' [mgcvST.estimate()]. The response name on the left-hand side is a label: it
#' need not be a column of `data`, because the helper replaces it with an
#' internal response bridge before constructing the mgcv setup. Linear terms,
#' factors, interactions, offsets, and additional mgcv smoothers on the
#' right-hand side are retained.
#'
#' `basis` is one object returned by [spde_basis()], and the model carries one
#' global spatial score process. Other user smoothers contribute nuisance
#' covariance to the marginal working model but are not score directions.
#'
#' @param formula A two-sided mgcv-style formula.
#' @param data Model data with one row per spatial observation.
#' @param basis A prepared [spde_basis()] object.
#' @param family An mgcv family object. The default is `mgcv::nb()`.
#' @param setting The model uses one `"global"` spatial score process.
#'   `setting = "global_local"` is not supported.
#' @param coordinates Character vector naming the two coordinate columns.
#' @param ... Additional arguments passed to `mgcv::gam(..., fit = FALSE)`.
#' @return An object of class `mgcvST_model` for [mgcvST.estimate()].
#' @export
model.set <- function(
    formula, data, basis, family = mgcv::nb(),
    setting = "global", coordinates = c("x", "y"), ...) {
  if (!identical(setting, "global")) {
    stop("setting must be \"global\". The second \"local\" geographic process ",
         "(setting = \"global_local\") was removed from mgcvST; supply one ",
         "spde_basis() object.")
  }
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop("formula must be a two-sided formula.")
  }
  response <- formula[[2L]]
  if (!is.symbol(response)) {
    stop("The formula response must be a single symbol.")
  }
  data <- as.data.frame(data)
  if (!nrow(data)) stop("data must contain at least one row.")
  coordinates <- as.character(coordinates)
  if (length(coordinates) != 2L || anyNA(coordinates) ||
      any(!nzchar(coordinates)) || anyDuplicated(coordinates)) {
    stop("coordinates must contain two different non-empty column names.")
  }
  if (!all(coordinates %in% names(data))) {
    stop("Both coordinate columns must be present in data.")
  }
  xy <- as.matrix(data[, coordinates, drop = FALSE])
  storage.mode(xy) <- "double"
  if (any(!is.finite(xy))) stop("Coordinate columns must be finite numeric values.")
  if (!inherits(basis, "mgcvST_spde_basis")) {
    stop("basis must be an object returned by spde_basis().")
  }
  basis <- list(global = basis)
  if (any(vapply(basis, function(x) is.null(x$kappa), logical(1L)))) {
    stop("model.set() requires fixed-kappa bases for covariance score testing.")
  }
  for (component in names(basis)) {
    .spde_basis_validate(basis[[component]], xy)
    basis[[component]]$component <- component
    basis[[component]]$score.component <- component
  }

  bridge <- ".mgcvST_response"
  if (bridge %in% names(data)) {
    stop("data must not contain the reserved column .mgcvST_response.")
  }
  env <- new.env(parent = environment(formula))
  env$.mgcvST_basis_global <- basis$global
  env$s <- mgcv::s
  formula0 <- formula
  environment(formula0) <- env
  formula0[[2L]] <- as.name(bridge)
  formula0 <- .mgcvst_formula_add(
    formula0,
    .mgcvst_spde_call(coordinates, ".mgcvST_basis_global"),
    as.name(bridge)
  )
  environment(formula0) <- env
  data[[bridge]] <- numeric(nrow(data))
  G <- mgcv::gam(
    formula0, data = data, family = family, fit = FALSE,
    na.action = stats::na.fail, ...
  )
  components <- "global"
  structure(
    list(
      G = G,
      setting = "global",
      components = components,
      formula = formula,
      internal_formula = formula0,
      response = as.character(response),
      coordinates = coordinates,
      kappa = stats::setNames(
        vapply(basis[components], function(x) {
          if (is.null(x$kappa)) NA_real_ else x$kappa
        }, numeric(1)), components
      ),
      offset = if (is.null(G$offset)) numeric(nrow(data)) else as.numeric(G$offset)
    ),
    class = "mgcvST_model"
  )
}

#' Print an mgcvST model setting
#'
#' @param x An object returned by `model.set()`.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.mgcvST_model <- function(x, ...) {
  cat("mgcvST model setting\n")
  cat("  setting:", x$setting, "\n")
  cat("  response bridge:", x$response, "-> .mgcvST_response\n")
  cat("  components:", paste(x$components, collapse = ", "), "\n")
  cat("  kappa:", paste(names(x$kappa), format(x$kappa), sep = "=", collapse = ", "), "\n")
  invisible(x)
}
