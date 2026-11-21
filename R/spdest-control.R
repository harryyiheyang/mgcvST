# The controls and the model specification of Stage 2 (spde.estimate()).

.spdest_merge_control <- function(base = list(), override = list()) {
  validate <- function(x, label) {
    if (is.null(x)) x <- list()
    if (!is.list(x)) stop(label, " must be a list.")
    nm <- names(x)
    if (length(x) && (is.null(nm) || anyNA(nm) || any(!nzchar(nm)) ||
                      anyDuplicated(nm))) {
      stop(label, " must have unique, non-empty names.")
    }
    x
  }
  out <- validate(base, "base control")
  override <- validate(override, "override control")
  for (nm in names(override)) out[nm] <- list(override[[nm]])
  out
}

# control = list(spde = <named spde.estimate() controls>).
.spdest_control <- function(control = list()) {
  if (is.null(control)) control <- list()
  if (!is.list(control)) stop("control must be a list.")
  .mgcvst_refuse_screen_control(control)
  out <- .spdest_merge_control(list(spde = list()), control)
  unknown <- setdiff(names(out), "spde")
  if (length(unknown)) {
    stop("Unknown control field", if (length(unknown) > 1L) "s" else "",
         ": ", paste(unknown, collapse = ", "),
         ". control takes spde, a named list of spde.estimate() controls.")
  }
  spde_allowed <- setdiff(names(formals(spde.estimate)), c(
    "Y", "X", "A", "Q", "offset", "poisson", "size_start", "variance_start", "variance_se",
    "constraint", "feature_id",
    "threads", "chunk_size"))
  if (!is.list(out$spde) || (length(out$spde) && (is.null(names(out$spde)) ||
      length(setdiff(names(out$spde), spde_allowed))))) {
    stop("control$spde takes named spde.estimate() controls: ",
         paste(spde_allowed, collapse = ", "), ".")
  }
  out
}

.spdest_validate_spec <- function(spec, y, offset) {
  if (!is.list(spec) || is.null(spec$fixed) || is.null(spec$random)) {
    stop("spec must contain fixed and random model blocks.")
  }
  X <- as.matrix(spec$fixed$X)
  storage.mode(X) <- "double"
  n <- length(y)
  if (length(dim(X)) != 2L || nrow(X) != n) {
    stop("spec$fixed$X must have one row per observation.")
  }
  xnames <- spec$fixed$names
  if (is.null(xnames)) xnames <- colnames(X)
  if (is.null(xnames)) {
    xnames <- if (ncol(X)) paste0("fixed", seq_len(ncol(X))) else character()
  }
  if (length(xnames) != ncol(X) || anyNA(xnames) ||
      any(!nzchar(xnames)) || anyDuplicated(xnames)) {
    stop("spec$fixed$names must uniquely name every fixed-effect column.")
  }
  if (!is.list(spec$random)) {
    stop("spec$random must be a list of latent blocks.")
  }
  rnames <- vapply(spec$random, function(z) {
    if (is.null(z$name)) "" else as.character(z$name)[1L]
  }, character(1L))
  if (any(!nzchar(rnames)) || anyNA(rnames) || anyDuplicated(rnames)) {
    stop("Every random block must have a unique non-empty name.")
  }
  random <- vector("list", length(spec$random))
  constraints <- vector("list", length(spec$random))
  for (j in seq_along(spec$random)) {
    z <- spec$random[[j]]
    A <- methods::as(z$A, "CsparseMatrix")
    Q <- methods::as(z$Q, "CsparseMatrix")
    if (nrow(A) != n || !ncol(A) || nrow(Q) != ncol(A) ||
        ncol(Q) != ncol(A)) {
      stop("Random block '", rnames[j], "' has incompatible A and Q matrices.")
    }
    kind <- if (is.null(z$kind)) "nuisance" else as.character(z$kind)[1L]
    if (!(kind %in% c("spde", "nuisance"))) {
      stop("Random block kind must be 'spde' or 'nuisance'.")
    }
    # The rank deficiency is known where Q is constructed (SPDE precisions and
    # iid identities are full rank), so it travels with the block spec.
    rankdef <- if (is.null(z$rankdef)) 0L else as.integer(z$rankdef)[1L]
    constraint <- NULL
    if (identical(kind, "spde")) {
      # This is deliberately recomputed, rather than merely trusting adapter
      # metadata.  Thus every spatial fit uses the observation mean constraint.
      constraint <- as.numeric(Matrix::crossprod(A, rep.int(1 / n, n)))
      if (!any(abs(constraint) > 0) || any(!is.finite(constraint))) {
        stop("SPDE block '", rnames[j], "' has an invalid observation constraint.")
      }
      if (!is.null(z$constraint) && !isTRUE(all.equal(
          as.numeric(z$constraint), constraint, tolerance = 1e-10,
          check.attributes = FALSE))) {
        stop("SPDE block '", rnames[j],
             "' constraint does not equal crossprod(A, 1 / n).")
      }
    } else if (!is.null(z$constraint)) {
      constraint <- as.numeric(z$constraint)
      if (length(constraint) != ncol(A) || any(!is.finite(constraint)) ||
          !any(abs(constraint) > 0)) {
        stop("Random block '", rnames[j], "' has an invalid constraint.")
      }
    }
    z$A <- A
    z$Q <- Matrix::forceSymmetric(Q)
    z$name <- rnames[j]
    z$kind <- kind
    z$rankdef <- rankdef
    z$target <- isTRUE(z$target)
    random[[j]] <- z
    constraints[j] <- list(constraint)
  }
  family <- tolower(as.character(spec$family)[1L])
  if (family %in% c("nb", "negative binomial", "negative_binomial")) {
    family <- "negative_binomial"
  }
  if (!(family %in% c("poisson", "negative_binomial"))) {
    stop("spec$family must be 'poisson' or 'negative_binomial'.")
  }
  if ((any(y < 0) || any(abs(y - round(y)) > 1e-8))) {
    stop("Poisson and negative-binomial responses must be non-negative counts.")
  }
  list(
    n = n, X = X, xnames = xnames, random = random,
    constraints = constraints, family = family, offset = offset
  )
}

.spdest_null_spec <- function(spec) {
  target <- which(vapply(spec$random, function(z) isTRUE(z$target), logical(1L)))
  random <- spec$random[-target]
  spec$random <- random
  spec$nuisance_index <- seq_len(ncol(spec$fixed$X) + sum(vapply(
    random, function(z) ncol(z$A), integer(1L)
  )))
  spec$nuisance_map <- NULL
  spec
}
