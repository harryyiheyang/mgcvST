# Post hoc covariance sensitivity for saved formal 3D null states; no refitting.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
library(mgcvST)
stopifnot(as.character(packageVersion("mgcvST")) == "0.0.1.9006")
package_path <- normalizePath(find.package("mgcvST"), winslash = "/", mustWork = TRUE)
package_version <- as.character(packageVersion("mgcvST"))
options(warn = 2)

a <- commandArgs(TRUE)
case_arg <- a[grepl("^--case=", a)]
first_arg <- a[grepl("^--first=", a)]
last_arg <- a[grepl("^--last=", a)]
partial <- "--partial" %in% a
if (length(case_arg) != 1L) stop("Supply one --case=independent-0.3, independent-3, joint-0.3, or joint-3.")
case <- sub("^--case=", "", case_arg)
if (!case %in% c("independent-0.3", "independent-3", "joint-0.3", "joint-3")) {
  stop("Unknown 3D null case: ", case)
}
first <- if (length(first_arg)) as.integer(sub("^--first=", "", first_arg)) else 1L
last <- if (length(last_arg)) as.integer(sub("^--last=", "", last_arg)) else 500L
if (length(first) != 1L || length(last) != 1L || is.na(first) || is.na(last) ||
    first < 1L || last > 500L || first > last) stop("Invalid replicate interval.")
if ((first != 1L || last != 500L) && !partial) {
  stop("A subset requires the explicit --partial flag.")
}

src <- file.path("artifacts/inla-stress-calibration/null3d-paired", case)
files <- file.path(src, sprintf("rep-%04d.rds", first:last))
missing <- (first:last)[!file.exists(files)]
if (length(missing) && !partial) {
  stop("The complete 500-replicate case is not available; first missing replicate: ", missing[1L])
}
files <- files[file.exists(files)]
if (!length(files)) stop("No saved formal replicates are available in the requested interval.")
label <- if (partial) sprintf("%s-partial-%04d-%04d", case, first, last) else case
out <- file.path("artifacts/inla-stress-calibration/null-sensitivity", label)
dir.create(out, recursive = TRUE, showWarnings = FALSE)

R <- list()
Fdiag <- list()
k <- 0L
kd <- 0L
for (file in files) {
  z <- readRDS(file)
  rep <- as.integer(sub("rep-([0-9]+)[.]rds", "\\1", basename(file)))
  if (length(z$states) != 2L || ncol(z$Y) < 2L) {
    stop("Saved formal replicate does not contain two complete states: ", file)
  }
  A <- list(native = vector("list", 2L), sandwich = vector("list", 2L),
    expected = vector("list", 2L), observed = vector("list", 2L))
  M <- list(native = vector("list", 2L), sandwich = vector("list", 2L),
    expected = vector("list", 2L), observed = vector("list", 2L))
  reason <- c(native = NA_character_, sandwich = NA_character_,
    expected = NA_character_, observed = NA_character_)
  leak_abs <- leak_rel <- leak_score_rel <- rep(NA_real_, 2L)

  for (j in 1:2) {
    s <- z$states[[j]]
    required <- c("a", "M", "Vp", "Vp_expected", "G", "w", "F", "D", "mu", "e")
    if (!all(required %in% names(s))) stop("Incomplete saved state in ", file, ", feature ", j, ".")
    F <- as.matrix(s$F)
    D <- as.numeric(s$D)
    mu <- as.numeric(s$mu)
    e <- as.numeric(s$e)
    Vp <- as.numeric(s$Vp)
    A$native[[j]] <- as.numeric(s$a)
    M$native[[j]] <- as.matrix(s$M)
    base_ok <- all(is.finite(F)) && all(is.finite(D)) && all(D > 0) &&
      all(is.finite(mu)) && all(mu > 0) && all(is.finite(e)) &&
      length(Vp) == 1L && is.finite(Vp) && Vp >= 0
    if (!base_ok) {
      reason[c("sandwich", "expected", "observed")] <- "nonfinite_or_invalid_saved_working_state"
      next
    }

    V <- diag(D) + tcrossprod(F)
    Vi <- chol2inv(chol(V))
    vx <- rowSums(Vi)
    sx <- sum(vx)
    P <- Vi - tcrossprod(vx) * Vp
    A$sandwich[[j]] <- A$native[[j]]
    M$sandwich[[j]] <- crossprod(P %*% F, V %*% (P %*% F))
    leak <- as.numeric(P %*% rep(1, nrow(P)))
    leak0 <- as.numeric(Vi %*% rep(1, nrow(P)))
    leak_abs[j] <- max(abs(leak))
    leak_rel[j] <- sqrt(sum(leak^2)) / sqrt(sum(leak0^2))
    leak_score_rel[j] <- sqrt(sum(crossprod(F, leak)^2)) /
      max(sqrt(sum(crossprod(F, leak0)^2)), .Machine$double.eps)

    expected_ok <- is.finite(sx) && sx > 0 && is.finite(1 / sx)
    if (expected_ok) {
      P0 <- Vi - tcrossprod(vx) / sx
      A$expected[[j]] <- as.numeric(crossprod(F, P0 %*% e))
      M$expected[[j]] <- crossprod(F, P0 %*% F)
    } else {
      reason[["expected"]] <- "nonfinite_expected_intercept_covariance"
    }

    h <- if (nrow(z$hyper) == 1L) 1L else j
    log_size <- z$hyper$log_size[h]
    y <- as.numeric(z$Y[, j])
    inv_size_saved <- exp(-log_size)
    if (is.finite(inv_size_saved) &&
        max(abs(D - (1 / mu + inv_size_saved))) > 1e-10) {
      stop("Saved expected NB working variance is inconsistent: ", file, ", feature ", j, ".")
    }
    observed_ok <- length(log_size) == 1L && is.finite(log_size) &&
      all(is.finite(y)) && all(y >= 0)
    if (observed_ok) {
      inv_size <- exp(-log_size)
      Wo <- mu * (1 + y * inv_size) / (1 + mu * inv_size)^2
      observed_ok <- all(is.finite(Wo)) && all(Wo > 0)
    }
    if (observed_ok) {
      Vo <- diag(1 / Wo) + tcrossprod(F)
      Voi <- chol2inv(chol(Vo))
      vxo <- rowSums(Voi)
      sxo <- sum(vxo)
      observed_ok <- is.finite(sxo) && sxo > 0 && is.finite(1 / sxo)
    }
    if (observed_ok) {
      Po <- Voi - tcrossprod(vxo) / sxo
      eta_offset <- e - (y - mu) / mu
      score_observed <- (y - mu) / (1 + mu * inv_size)
      e_observed <- eta_offset + score_observed / Wo
      observed_ok <- all(is.finite(e_observed))
    }
    if (observed_ok) {
      A$observed[[j]] <- as.numeric(crossprod(F, Po %*% e_observed))
      M$observed[[j]] <- crossprod(F, Po %*% F)
    } else {
      reason[["observed"]] <- "nonfinite_observed_hessian_or_working_response"
    }

    kd <- kd + 1L
    Fdiag[[kd]] <- data.frame(case = case, replicate = rep, feature = j,
      native_Vp = Vp, expected_Vp = if (expected_ok) 1 / sx else NA_real_,
      maximum_absolute_mean_leakage = leak_abs[j],
      relative_mean_leakage = leak_rel[j],
      relative_score_mean_leakage = leak_score_rel[j],
      expected_available = expected_ok, observed_available = observed_ok,
      stringsAsFactors = FALSE)
  }

  U <- c(native = NA_real_, sandwich = NA_real_, expected = NA_real_, observed = NA_real_)
  for (route in names(U)) {
    ready <- all(vapply(A[[route]], function(x) length(x) && all(is.finite(x)), logical(1L))) &&
      all(vapply(M[[route]], function(x) length(x) && all(is.finite(x)), logical(1L)))
    if (ready) U[[route]] <- sum(A[[route]][[1L]] * A[[route]][[2L]])
    if (!ready && is.na(reason[[route]])) reason[[route]] <- "incomplete_route_state"
  }
  original_U <- unique(z$rows$statistic)
  if (length(original_U) != 1L || !is.finite(original_U) ||
      abs(U[["native"]] - original_U) > 1e-10) stop("Native score reconstruction failed: ", file)

  diag_file <- file.path("artifacts/inla-stress-calibration/score-spectra", case,
    sprintf("rep-%04d.rds", rep))
  diagnostic_match <- NA
  if (file.exists(diag_file)) {
    dz <- readRDS(diag_file)
    diagnostic_match <- identical(z$Y[, 1:2, drop = FALSE], dz$Y[, 1:2, drop = FALSE])
    if (!diagnostic_match) stop("Formal and spectrum-diagnostic responses differ: replicate ", rep)
  }

  for (cal in c("davies", "liu")) {
    original <- z$rows[z$rows$calibration == cal, , drop = FALSE]
    if (nrow(original) != 1L) stop("Original calibration row is missing: ", file)
    for (route in names(U)) {
      if (route == "native") {
        p <- original$p_value
        info <- original$information
        route_reason <- if (is.finite(p)) NA_character_ else if (!is.na(original$error))
          original$error else "native_nonfinite_without_message"
      } else if (is.finite(U[[route]])) {
        ans <- rkhs_score_calibrate(U[[route]], M[[route]][[1L]], M[[route]][[2L]], method = cal)
        p <- ans$p_two_sided
        info <- ans$information
        route_reason <- if (is.finite(p)) NA_character_ else "calibration_returned_nonfinite"
      } else {
        p <- info <- NA_real_
        route_reason <- reason[[route]]
      }
      k <- k + 1L
      R[[k]] <- data.frame(case = case, replicate = rep, seed = original$seed,
        calibration = cal, route = route, p_value = p, statistic = U[[route]],
        information = info, available = is.finite(p), unavailable_reason = route_reason,
        native_reference_p = original$p_value,
        absolute_p_difference = abs(p - original$p_value),
        relative_p_difference = abs(p - original$p_value) /
          max(abs(original$p_value), .Machine$double.eps),
        native_reference_statistic = original_U,
        relative_statistic_difference = abs(U[[route]] - original_U) /
          max(abs(original_U), .Machine$double.eps),
        maximum_relative_mean_leakage = max(leak_rel, na.rm = TRUE),
        formal_diagnostic_match = diagnostic_match, stringsAsFactors = FALSE)
    }
  }
}

R <- do.call(rbind, R)
Fdiag <- do.call(rbind, Fdiag)
write.csv(R, file.path(out, "replicate-results.csv"), row.names = FALSE)
write.csv(Fdiag, file.path(out, "feature-diagnostics.csv"), row.names = FALSE)

S <- list()
ks <- 0L
planned <- 500L
for (cal in c("davies", "liu")) for (route in c("native", "sandwich", "expected", "observed")) {
  d <- R[R$calibration == cal & R$route == route, ]
  available <- sum(d$available)
  unavailable <- nrow(d) - available
  absent <- planned - nrow(d)
  for (alpha in c(0.05, 0.01)) {
    rejected <- sum(d$p_value < alpha, na.rm = TRUE)
    ci <- if (available) binom.test(rejected, available)$conf.int else c(NA_real_, NA_real_)
    ks <- ks + 1L
    S[[ks]] <- data.frame(case = case, analysis_status = if (partial) "partial" else "complete",
      calibration = cal, route = route, alpha = alpha, planned_attempts = planned,
      requested_replicates = last - first + 1L, read_replicates = nrow(d),
      missing_in_requested_interval = length(missing), available = available,
      unavailable = unavailable, absent_from_500 = absent, rejected = rejected,
      rejection_rate_among_available = if (available) rejected / available else NA_real_,
      binomial_CI_lower = ci[1L], binomial_CI_upper = ci[2L],
      all_attempt_lower = rejected / planned,
      all_attempt_upper = (rejected + unavailable + absent) / planned,
      stringsAsFactors = FALSE)
  }
}
S <- do.call(rbind, S)
write.csv(S, file.path(out, "calibration-summary.csv"), row.names = FALSE)
write.csv(as.data.frame(table(R$route, R$unavailable_reason, useNA = "ifany")),
  file.path(out, "unavailable-reasons.csv"), row.names = FALSE)
write.csv(data.frame(package_path = package_path, package_version = package_version),
  file.path(out, "package-provenance.csv"), row.names = FALSE)
writeLines(c(
  "Post hoc sensitivity prompted by the expected-working/native-Vp curvature mismatch.",
  "The prespecified 500-replicate native experiment remains primary; no original p-value is replaced.",
  "Sandwich retains native U and uses F' P V P F; reported mean leakage remains unresolved.",
  "Expected recomputes P, M, feature scores and U with the expected NB working state.",
  "Observed recomputes observed-Hessian weights, Newton working responses, P, M, feature scores and U.",
  "Every unavailable route is retained with a reason. Partial runs retain a 500-attempt uncertainty bound and are labelled partial.",
  "Exact binomial intervals use available replicates; all-attempt bounds include unavailable and not-yet-read replicates."
), file.path(out, "notes.txt"))
