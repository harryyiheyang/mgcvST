# Score-covariance sensitivity from saved 3D spectrum states; no refitting.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
library(mgcvST)
stopifnot(as.character(packageVersion("mgcvST")) == "0.0.1.9006")
package_path <- normalizePath(find.package("mgcvST"), winslash = "/", mustWork = TRUE)
package_version <- as.character(packageVersion("mgcvST"))
options(warn = 2)

src <- "artifacts/inla-stress-calibration/score-spectra/independent-0.3"
out <- "artifacts/inla-stress-calibration/score-spectrum-sensitivity"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
files <- file.path(src, sprintf("rep-%04d.rds", 1:10))
if (!all(file.exists(files))) stop("All ten saved spectrum diagnostics are required.")

R <- list()
D <- list()
k <- 0L
kd <- 0L
for (r in seq_along(files)) {
  z <- readRDS(files[r])
  if (length(z$spectra) != 2L || ncol(z$Y) != 2L) {
    stop("Each saved replicate must contain two score states and responses.")
  }
  a_native <- a_sandwich <- a_expected <- a_observed <- vector("list", 2L)
  M_native <- M_sandwich <- M_expected <- M_observed <- vector("list", 2L)
  leakage_max <- leakage_relative <- leakage_score_relative <- numeric(2L)

  for (j in 1:2) {
    s <- z$spectra[[j]]
    F <- s$F
    V <- s$V
    Vi <- chol2inv(chol(V))
    vx <- rowSums(Vi)
    one <- rep(1, nrow(V))
    P <- s$P
    P0 <- Vi - tcrossprod(vx) / sum(vx)

    size <- exp(z$hyper$log_size[j])
    ew <- as.numeric(s$expected_weights)
    mu <- ew * size / (size - ew)
    y <- z$Y[, j]
    eta_offset <- s$e - (y - mu) / mu
    Wo <- as.numeric(s$observed_weights)
    Vo <- diag(1 / Wo) + tcrossprod(F)
    Voi <- chol2inv(chol(Vo))
    vxo <- rowSums(Voi)
    Po <- Voi - tcrossprod(vxo) / sum(vxo)
    score_observed <- size * (y - mu) / (size + mu)
    e_observed <- eta_offset + score_observed / Wo

    a_native[[j]] <- as.numeric(crossprod(F, P %*% s$e))
    a_sandwich[[j]] <- a_native[[j]]
    a_expected[[j]] <- as.numeric(crossprod(F, P0 %*% s$e))
    a_observed[[j]] <- as.numeric(crossprod(F, Po %*% e_observed))
    M_native[[j]] <- s$M
    M_sandwich[[j]] <- crossprod(P %*% F, V %*% (P %*% F))
    M_expected[[j]] <- crossprod(F, P0 %*% F)
    M_observed[[j]] <- crossprod(F, Po %*% F)
    stopifnot(max(abs(M_sandwich[[j]] - s$M_sandwich)) < 1e-10)
    stopifnot(max(abs(M_expected[[j]] - s$M_expected)) < 1e-10)
    stopifnot(max(abs(M_observed[[j]] - s$M_observed)) < 1e-10)

    leak <- as.numeric(P %*% one)
    leak0 <- as.numeric(Vi %*% one)
    leakage_max[j] <- max(abs(leak))
    leakage_relative[j] <- sqrt(sum(leak^2)) / sqrt(sum(leak0^2))
    leakage_score_relative[j] <- sqrt(sum(crossprod(F, leak)^2)) /
      max(sqrt(sum(crossprod(F, leak0)^2)), .Machine$double.eps)

    ev <- eigen((s$M + t(s$M)) / 2, symmetric = TRUE)
    neg <- ev$values < 0
    kd <- kd + 1L
    D[[kd]] <- data.frame(
      replicate = r, feature = j,
      minimum_native_eigenvalue = min(ev$values),
      maximum_absolute_native_eigenvalue = max(abs(ev$values)),
      negative_count = sum(neg),
      removed_negative_trace = -sum(ev$values[neg]),
      removed_negative_trace_fraction = -sum(ev$values[neg]) / sum(abs(ev$values)),
      removed_negative_frobenius = sqrt(sum(ev$values[neg]^2)),
      removed_negative_frobenius_fraction = sqrt(sum(ev$values[neg]^2)) /
        sqrt(sum(ev$values^2)), stringsAsFactors = FALSE)
  }

  U <- c(
    native = sum(a_native[[1L]] * a_native[[2L]]),
    sandwich = sum(a_sandwich[[1L]] * a_sandwich[[2L]]),
    expected = sum(a_expected[[1L]] * a_expected[[2L]]),
    observed = sum(a_observed[[1L]] * a_observed[[2L]])
  )
  stopifnot(abs(U[["native"]] - z$rows$statistic[1L]) < 1e-10)
  MM <- list(native = M_native, sandwich = M_sandwich,
    expected = M_expected, observed = M_observed)

  for (cal in c("davies", "liu")) {
    ref <- z$rows$p_value[z$rows$calibration == cal]
    if (length(ref) != 1L) stop("Saved native reference calibration is incomplete.")
    for (route in names(MM)) {
      ans <- if (route == "native") {
        list(p_two_sided = ref, information = z$rows$information[
          z$rows$calibration == cal])
      } else {
        rkhs_score_calibrate(U[[route]], MM[[route]][[1L]], MM[[route]][[2L]],
          method = cal)
      }
      k <- k + 1L
      R[[k]] <- data.frame(
        replicate = r, calibration = cal, route = route,
        p_value = ans$p_two_sided, native_reference_p = ref,
        absolute_p_difference = abs(ans$p_two_sided - ref),
        relative_p_difference = abs(ans$p_two_sided - ref) /
          max(abs(ref), .Machine$double.eps),
        statistic = U[[route]], native_reference_statistic = U[["native"]],
        absolute_statistic_difference = abs(U[[route]] - U[["native"]]),
        relative_statistic_difference = abs(U[[route]] - U[["native"]]) /
          max(abs(U[["native"]]), .Machine$double.eps),
        information = ans$information,
        maximum_absolute_mean_leakage = max(leakage_max),
        maximum_relative_mean_leakage = max(leakage_relative),
        maximum_relative_score_mean_leakage = max(leakage_score_relative),
        stringsAsFactors = FALSE)
    }
  }
}

R <- do.call(rbind, R)
D <- do.call(rbind, D)
write.csv(R, file.path(out, "score-sensitivity.csv"), row.names = FALSE)
write.csv(D, file.path(out, "psd-projection-diagnostic.csv"), row.names = FALSE)

S <- aggregate(p_value ~ calibration + route, data = R,
  FUN = function(x) sum(is.finite(x)), na.action = na.pass)
names(S)[3L] <- "finite_p_count"
S$native_reference_finite_count <- vapply(seq_len(nrow(S)), function(i) {
  z <- R$native_reference_p[R$calibration == S$calibration[i] & R$route == S$route[i]]
  sum(is.finite(z))
}, integer(1L))
S$median_absolute_p_difference <- vapply(seq_len(nrow(S)), function(i) {
  z <- R$absolute_p_difference[R$calibration == S$calibration[i] & R$route == S$route[i]]
  median(z, na.rm = TRUE)
}, numeric(1L))
S$maximum_absolute_p_difference <- vapply(seq_len(nrow(S)), function(i) {
  z <- R$absolute_p_difference[R$calibration == S$calibration[i] & R$route == S$route[i]]
  max(z, na.rm = TRUE)
}, numeric(1L))
S$median_relative_statistic_difference <- vapply(seq_len(nrow(S)), function(i) {
  z <- R$relative_statistic_difference[R$calibration == S$calibration[i] & R$route == S$route[i]]
  median(z)
}, numeric(1L))
S$maximum_relative_statistic_difference <- vapply(seq_len(nrow(S)), function(i) {
  z <- R$relative_statistic_difference[R$calibration == S$calibration[i] & R$route == S$route[i]]
  max(z)
}, numeric(1L))
write.csv(S, file.path(out, "score-sensitivity-summary.csv"), row.names = FALSE)
write.csv(data.frame(package_path = package_path, package_version = package_version),
  file.path(out, "package-provenance.csv"), row.names = FALSE)

z <- readRDS(files[5L])
H <- list()
for (j in 1:2) {
  F <- z$spectra[[j]]$F
  E <- svd(F)
  nr <- sum(E$d > max(E$d) * 1e-12)
  B <- E$u[, seq_len(nr), drop = FALSE]
  one <- rep(1, nrow(F))
  H[[j]] <- data.frame(replicate = 5L, feature = j, observations = nrow(F),
    columns = ncol(F), rank = nr, relative_rank_tolerance = 1e-12,
    smallest_singular_value = min(E$d),
    relative_constant_projection_residual = sqrt(sum((one - B %*% crossprod(B, one))^2) / sum(one^2)))
}
write.csv(do.call(rbind, H), file.path(out, "factor-geometry.csv"), row.names = FALSE)

writeLines(c(
  "Ten saved independent-0.3 spectrum replicates; no model was refitted.",
  "Native results remain the reference and are not replaced.",
  "Sandwich uses the original statistic and F' P V P F covariance.",
  "Expected recomputes P, M, feature scores and U with the expected NB working covariance.",
  "Observed recomputes the NB observed-Hessian weights, Newton working response, P, M, feature scores and U.",
  "Mean leakage is reported for the native hybrid P; sandwich covariance does not remove that leakage or nuisance-estimation uncertainty.",
  "PSD-projection fields only quantify negative spectral mass; projected matrices are not used for p-values.",
  "Ten replicates measure sensitivity magnitude and do not replace the prespecified 500-replicate type-I experiment."
), file.path(out, "notes.txt"))
