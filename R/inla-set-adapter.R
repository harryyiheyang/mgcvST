# Replace only SPDE xt arguments before delegating the complete formula to
# the same mgcv parser/design freezer used by mgcvST.set().
.inlast_prepare_formula <- function(formula, data) {
  if (!inherits(formula, "formula") || length(formula) != 3L ||
      !is.symbol(formula[[2L]])) {
    stop("Supply a two-sided formula with a single response name.")
  }
  env <- new.env(parent = environment(formula))
  env$s <- mgcv::s
  evaluation <- list2env(as.list(data), parent = env)
  prepared <- list()
  walk <- function(expr) {
    if (!is.call(expr)) return(expr)
    head <- expr[[1L]]
    is_s <- identical(head, as.name("s")) ||
      (is.call(head) && identical(head[[1L]], as.name("::")) &&
       identical(as.character(head[[2L]]), "mgcv") &&
       identical(as.character(head[[3L]]), "s"))
    if (is_s) {
      spec <- eval(expr, envir = env)
      if (inherits(spec, "spdePC.smooth.spec")) {
        stop("inlaST.set() requires a full SPDE basis; bs = 'spdePC' is not supported.")
      }
      if (inherits(spec, "spde.smooth.spec")) {
        if (length(spec$term) != 2L || !identical(spec$by, "NA") ||
            !is.null(spec$id)) {
          stop("INLA SPDE terms require two coordinates and do not support by or linked id terms.")
        }
        xy <- lapply(spec$term, function(term) eval(str2lang(term), evaluation))
        if (any(vapply(xy, length, integer(1L)) != nrow(data)) ||
            any(!vapply(xy, is.numeric, logical(1L)))) {
          stop("Each SPDE coordinate must be a numeric observation-length vector.")
        }
        value <- .inlast_prepare_basis(spec$xt, do.call(cbind, xy))
        component <- .spde_basis_component(spec)$component
        if (component %in% names(prepared)) {
          stop("Supply only one SPDE term for each global/local component.")
        }
        value$basis$component <- value$basis$score.component <- component
        prepared[[component]] <<- value
        key <- paste0(".inlaST_prepared_", component)
        assign(key, value$basis, env)
        expr$xt <- as.name(key)
        return(expr)
      }
      return(expr)
    }
    if (length(expr) > 1L) {
      for (i in seq.int(2L, length(expr))) expr[[i]] <- walk(expr[[i]])
    }
    expr
  }
  formula[[3L]] <- walk(formula[[3L]])
  environment(formula) <- env
  if (!length(prepared)) stop("The complete formula must contain a full SPDE term.")
  list(formula = formula, prepared = prepared)
}

# External G is already a frozen model. Preserve its coefficient coordinates
# and check that they really represent the required observation constraint.
.inlast_prepare_frozen_components <- function(base) {
  prepared <- list()
  for (component in base$components) {
    j <- base$geometry$target[[component]]
    sm <- base$G$smooth[[j]]
    if (!inherits(sm, "spde.smooth") || inherits(sm, "spdePC.smooth")) {
      stop("inlaST.set(G=...) requires a full SPDE basis; spdePC is unsupported.")
    }
    if (!identical(sm$by, "NA") || !is.null(sm$id) || length(sm$term) != 2L ||
        !all(sm$term %in% names(base$G$mf))) {
      stop("The external G must retain two SPDE coordinate columns without by/id terms.")
    }
    xy <- as.matrix(base$G$mf[, sm$term, drop = FALSE])
    value <- .inlast_prepare_basis(sm$basis, xy)
    Z <- as.matrix(sm$basis$projection)
    m <- ncol(value$raw$A)
    valid <- identical(dim(Z), c(m, m - 1L)) && all(is.finite(Z))
    if (valid) valid <- max(abs(crossprod(Z) - diag(ncol(Z)))) < 1e-8 &&
      max(abs(crossprod(value$raw$constraint, Z))) < 1e-10
    if (!valid) {
      stop("The external G does not enforce observation mean-zero; rebuild with inlaST.set(formula, data, ...).")
    }
    B <- as.matrix(value$raw$A %*% Z)
    Q <- crossprod(Z, as.matrix(value$raw$Q %*% Z))
    geometry <- base$geometry$smooth[[j]]
    if (!isTRUE(all.equal(unname(B), unname(geometry$B), tolerance = 1e-9)) ||
        length(geometry$penalties) != 1L ||
        !isTRUE(all.equal(unname(Q), unname(geometry$penalties[[1L]]),
                         tolerance = 1e-9))) {
      stop("The external G SPDE design/penalty differs from its raw constrained precision.")
    }
    value$raw$projection <- Z
    prepared[[component]] <- value
  }
  prepared
}

.inlast_family_control <- function(model, control) {
  fixed_size <- model$inla_spec$nb_size_fixed
  if (!is.null(fixed_size)) {
    if (!is.null(control[["nb_size", exact = TRUE]]) &&
        !isTRUE(all.equal(control[["nb_size", exact = TRUE]], fixed_size, tolerance = 1e-12))) {
      stop("control$nb_size conflicts with the fixed theta in the model family.")
    }
    control$nb_size <- fixed_size
  }
  control
}
