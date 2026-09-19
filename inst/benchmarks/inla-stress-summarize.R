#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
Sys.setenv(LC_ALL = "C")

root <- Sys.getenv("MGCVST_STRESS_OUTPUT",
  "artifacts/inla-stress-calibration")
out <- Sys.getenv("MGCVST_STRESS_SUMMARY_OUTPUT", file.path(root, "summary"))
max.rep <- as.integer(Sys.getenv("MGCVST_STRESS_SUMMARY_MAX_REP", "500"))
target <- 500L
if (!dir.exists(root)) stop("Stress-calibration output directory is missing.")
if (!is.finite(max.rep) || max.rep < 1L || max.rep > 500L) {
  stop("MGCVST_STRESS_SUMMARY_MAX_REP must be between 1 and 500.")
}
dir.create(out, recursive = TRUE, showWarnings = FALSE)

binom_ci <- function(x, n) {
  if (!n) return(c(NA_real_, NA_real_))
  as.numeric(binom.test(x, n, conf.level = 0.95)$conf.int)
}

cases <- data.frame(dimension = c(rep("2d", 4L), rep("3d", 4L)),
  case = c("pair-0.3", "pair-3", "marginal-0.3", "marginal-3",
    "independent-0.3", "independent-3", "joint-0.3", "joint-3"),
  stringsAsFactors = FALSE)
null.rows <- list()
hyper.rows <- list()
completion.rows <- list()
unavailable.rows <- list()
spectrum.rows <- list()
failure.rows <- list()
state.rows <- list()
eigen.rows <- list()
raw.rows <- list()

for (s in seq_len(nrow(cases))) {
  dim <- cases$dimension[s]
  case <- cases$case[s]
  folder <- file.path(root, if (dim == "2d") "null2d" else "null3d-paired",
    case)
  wanted <- sprintf("rep-%04d.rds", seq_len(max.rep))
  files <- if (dir.exists(folder)) list.files(folder, pattern = "^rep-[0-9]{4}\\.rds$",
    full.names = TRUE) else character()
  files <- files[basename(files) %in% wanted]
  rep.id <- as.integer(sub("rep-([0-9]{4})\\.rds", "\\1", basename(files)))
  files <- files[order(rep.id)]
  rep.id <- sort(rep.id)
  failures <- if (dir.exists(folder)) list.files(folder,
    pattern = "^rep-[0-9]{4}-failure\\.json$", full.names = TRUE) else character()
  failure.id <- as.integer(sub("rep-([0-9]{4})-failure\\.json", "\\1",
    basename(failures)))
  keep.failure <- failure.id <= max.rep & !failure.id %in% rep.id
  failures <- failures[keep.failure]
  failure.id <- failure.id[keep.failure]
  for (j in seq_along(failures)) {
    F <- jsonlite::read_json(failures[j], simplifyVector = TRUE)
    failure.rows[[length(failure.rows) + 1L]] <- data.frame(dimension = dim,
      case = case, replicate = failure.id[j], label = F$label,
      group = F$group, exit_code = F$exit_code, status = F$status,
      wall_seconds = F$wall_seconds, peak_rss_bytes = F$peak_rss_bytes,
      peak_private_bytes = F$peak_private_bytes)
  }
  failure.id <- unique(failure.id[failure.id <= max.rep & !failure.id %in% rep.id])
  incomplete.id <- setdiff(seq_len(target), c(rep.id, failure.id))
  completion.rows[[s]] <- data.frame(dimension = dim, case = case,
    target = target, inspected_through = max.rep, completed = length(files),
    failures = length(failure.id),
    incomplete = length(incomplete.id), final = length(files) +
      length(failure.id) == target)
  if (!length(files)) next

  rows <- list()
  hyper <- list()
  for (j in seq_along(files)) {
    z <- readRDS(files[j])
    z$rows$replicate_file <- rep.id[j]
    rows[[j]] <- z$rows
    z$hyper$replicate <- rep.id[j]
    hyper[[j]] <- z$hyper
    if (dim == "3d" && length(z$states)) {
      for (h in seq_along(z$states)) {
        st <- z$states[[h]]
        evals <- eigen(st$M, symmetric = TRUE, only.values = TRUE)$values
        expected.reason <- NA_character_
        if (length(st$Vp_expected) != 1L) expected.reason <- "Vp_expected_missing"
        if (length(st$Vp_expected) == 1L && !is.finite(st$Vp_expected)) {
          expected.reason <- "Vp_expected_nonfinite"
        }
        if (!all(is.finite(st$G))) expected.reason <- "G_nonfinite"
        if (!all(is.finite(st$w))) expected.reason <- "w_nonfinite"
        expected.available <- is.na(expected.reason)
        if (expected.available) {
          expected.M <- st$G - st$Vp_expected * tcrossprod(st$w)
          expected.available <- all(is.finite(expected.M))
          if (!expected.available) expected.reason <- "G_minus_Vp_expected_wwt_nonfinite"
        }
        if (expected.available) {
          expected.evals <- eigen(expected.M, symmetric = TRUE,
            only.values = TRUE)$values
        } else {
          expected.evals <- NA_real_
        }
        scale.M <- max(abs(evals))
        negative.mass <- sum(abs(evals[evals < 0]))
        total.mass <- sum(abs(evals))
        state.rows[[length(state.rows) + 1L]] <- data.frame(dimension = dim,
          case = case, replicate = rep.id[j], state = h, eigenvalues = length(evals),
          minimum_eigenvalue = min(evals), maximum_eigenvalue = max(evals),
          maximum_absolute_eigenvalue = scale.M,
          relative_minimum_eigenvalue = if (scale.M > 0) min(evals) / scale.M else NA_real_,
          negative_count = sum(evals < 0),
          negative_within_existing_cutoff = sum(evals >= -1e-10 & evals < 0),
          negative_below_existing_cutoff = sum(evals < -1e-10),
          negative_mass = negative.mass,
          negative_mass_fraction = if (total.mass > 0) negative.mass /
            total.mass else NA_real_,
          zero_spectrum = scale.M == 0,
          symmetry_max_abs = max(abs(st$M - t(st$M))),
          reconstruction_max_abs = max(abs(st$M -
            (st$G - st$Vp * tcrossprod(st$w)))),
          Vp = st$Vp, Vp_expected = st$Vp_expected,
          Vp_ratio = if (is.finite(st$Vp) && is.finite(st$Vp_expected))
            st$Vp / st$Vp_expected else NA_real_,
          expected_spectrum_available = expected.available,
          expected_spectrum_unavailable_reason = if (expected.available) NA_character_
            else expected.reason,
          expected_minimum_eigenvalue = if (expected.available)
            min(expected.evals) else NA_real_,
          expected_negative_count = if (expected.available)
            sum(expected.evals < 0) else NA_integer_,
          F_rows = nrow(st$F), F_columns = ncol(st$F),
          finite_D = sum(is.finite(st$D)), D_min = min(st$D), D_max = max(st$D),
          finite_mu = sum(is.finite(st$mu)), mu_min = min(st$mu), mu_max = max(st$mu),
          finite_e = sum(is.finite(st$e)), e_min = min(st$e), e_max = max(st$e))
        eigen.rows[[length(eigen.rows) + 1L]] <- data.frame(dimension = dim,
          case = case, replicate = rep.id[j], state = h,
          eigen_index_descending = seq_along(evals), eigenvalue = evals,
          below_existing_cutoff = evals < -1e-10,
          negative_within_existing_cutoff = evals >= -1e-10 & evals < 0)
      }
    }
  }
  R <- do.call(rbind, rows)
  H <- do.call(rbind, hyper)
  R$dimension <- dim
  H$dimension <- dim
  H$case <- case
  raw.rows[[s]] <- data.frame(dimension = dim, case = case,
    replicate = R$replicate, seed = R$seed, calibration = R$calibration,
    p_value = R$p_value, U = R$statistic, information = R$information,
    error = R$error)

  for (cal in c("davies", "liu")) {
    Z <- R[R$calibration == cal, , drop = FALSE]
    valid <- is.finite(Z$p_value) & Z$p_value >= 0 & Z$p_value <= 1
    unavailable <- !valid
    for (alpha in c(0.05, 0.01)) {
      reject <- sum(Z$p_value[valid] < alpha)
      ci <- binom_ci(reject, sum(valid))
      unresolved <- sum(unavailable) + length(failure.id) + length(incomplete.id)
      null.rows[[length(null.rows) + 1L]] <- data.frame(dimension = dim,
        case = case, calibration = cal, alpha = alpha, target = target,
        completed = nrow(Z), valid = sum(valid), unavailable = sum(unavailable),
        failures = length(failure.id), incomplete = length(incomplete.id),
        rejections = reject, rejection_rate_valid = reject / sum(valid),
        exact_ci_low = ci[1L], exact_ci_high = ci[2L],
        all_target_lower = reject / target,
        all_target_upper = (reject + unresolved) / target,
        final = length(incomplete.id) == 0L)
    }
    if (any(unavailable)) {
      bad <- Z[unavailable, , drop = FALSE]
      reason <- as.character(bad$error)
      spectrum.error <- !is.na(reason) & grepl("Non-PSD plug-in score matrix",
        reason, fixed = TRUE)
      reason[spectrum.error] <- "minimum_M_eigenvalue_below_existing_cutoff"
      reason[is.na(reason) & !is.finite(bad$statistic)] <- "nonfinite_statistic"
      if ("minimum_M_eigenvalue" %in% names(bad)) {
        reason[is.na(reason) & bad$minimum_M_eigenvalue < -1e-10] <-
          "below_existing_spectrum_cutoff"
      }
      reason[is.na(reason)] <- "calibration_returned_unavailable_without_worker_error"
      unavailable.rows[[length(unavailable.rows) + 1L]] <- data.frame(dimension = dim,
        case = case, calibration = cal, replicate = bad$replicate,
        p_value = bad$p_value, statistic = bad$statistic,
        information = bad$information, worker_error = bad$error,
        unavailable_reason = reason,
        minimum_M_eigenvalue = if ("minimum_M_eigenvalue" %in% names(bad))
          bad$minimum_M_eigenvalue else NA_real_)
    }
  }

  conv <- if ("converged" %in% names(H)) H$converged else H$mode_status == 0L
  tau <- if ("tau" %in% names(H)) H$tau else exp(H$log_precision)
  size <- if ("size" %in% names(H)) H$size else exp(H$log_size)
  hyper.rows[[s]] <- data.frame(dimension = dim, case = case,
    fits = nrow(H), converged = sum(conv, na.rm = TRUE),
    nonconverged = sum(!conv | is.na(conv)), finite_tau = sum(is.finite(tau)),
    tau_min = min(tau, na.rm = TRUE), tau_median = median(tau, na.rm = TRUE),
    tau_max = max(tau, na.rm = TRUE), extreme_tau = sum(is.finite(tau) &
      (tau < 1e-8 | tau > 1e8)), finite_size = sum(is.finite(size)),
    size_min = min(size, na.rm = TRUE), size_median = median(size, na.rm = TRUE),
    size_max = max(size, na.rm = TRUE), extreme_size = sum(is.finite(size) &
      (size < 1e-6 | size > 1e6)))

  if ("minimum_M_eigenvalue" %in% names(R)) {
    E <- R[!duplicated(R$replicate), c("replicate", "minimum_M_eigenvalue",
      "statistic", "error")]
    neg <- is.finite(E$minimum_M_eigenvalue) & E$minimum_M_eigenvalue < 0
    spectrum.rows[[s]] <- data.frame(dimension = dim, case = case,
      completed = nrow(E), finite_spectrum = sum(is.finite(E$minimum_M_eigenvalue)),
      below_existing_cutoff = sum(E$minimum_M_eigenvalue < -1e-10, na.rm = TRUE),
      negative_within_existing_cutoff = sum(E$minimum_M_eigenvalue >= -1e-10 &
        E$minimum_M_eigenvalue < 0, na.rm = TRUE),
      minimum_eigenvalue = min(E$minimum_M_eigenvalue, na.rm = TRUE),
      q01_eigenvalue = quantile(E$minimum_M_eigenvalue, 0.01, na.rm = TRUE,
        names = FALSE), median_eigenvalue = median(E$minimum_M_eigenvalue,
        na.rm = TRUE))
    if (any(neg)) write.csv(E[neg, ], file.path(out,
      paste0("negative-spectrum-", case, ".csv")), row.names = FALSE)
  }
}

write.csv(do.call(rbind, completion.rows), file.path(out,
  "null-completion.csv"), row.names = FALSE)
if (length(null.rows)) write.csv(do.call(rbind, null.rows), file.path(out,
  "null-rejection-summary.csv"), row.names = FALSE)
if (length(hyper.rows)) write.csv(do.call(rbind, hyper.rows), file.path(out,
  "null-hyperparameter-summary.csv"), row.names = FALSE)
if (length(unavailable.rows)) write.csv(do.call(rbind, unavailable.rows), file.path(out,
  "null-unavailable-detail.csv"), row.names = FALSE)
if (length(failure.rows)) write.csv(do.call(rbind, failure.rows), file.path(out,
  "null-failure-detail.csv"), row.names = FALSE)
if (length(spectrum.rows)) write.csv(do.call(rbind, spectrum.rows), file.path(out,
  "null-spectrum-summary.csv"), row.names = FALSE)
if (length(raw.rows)) write.csv(do.call(rbind, raw.rows), file.path(out,
  "null-p-values.csv"), row.names = FALSE)
if (length(state.rows)) {
  state.summary <- do.call(rbind, state.rows)
  write.csv(state.summary, file.path(out, "null-state-matrix-summary.csv"),
    row.names = FALSE)
  state.aggregate <- list()
  state.cases <- unique(state.summary$case)
  for (case in state.cases) {
    Z <- state.summary[state.summary$case == case, , drop = FALSE]
    ratio <- Z$Vp_ratio[is.finite(Z$Vp_ratio)]
    state.aggregate[[case]] <- data.frame(case = case, states = nrow(Z),
      zero_spectrum = sum(Z$zero_spectrum),
      states_with_negative_eigenvalue = sum(Z$negative_count > 0),
      states_below_existing_cutoff = sum(Z$negative_below_existing_cutoff > 0),
      negative_within_existing_cutoff = sum(Z$negative_within_existing_cutoff),
      negative_below_existing_cutoff = sum(Z$negative_below_existing_cutoff),
      minimum_eigenvalue = min(Z$minimum_eigenvalue),
      minimum_relative_eigenvalue = min(Z$relative_minimum_eigenvalue,
        na.rm = TRUE), median_negative_mass_fraction = median(
        Z$negative_mass_fraction, na.rm = TRUE),
      maximum_negative_mass_fraction = max(Z$negative_mass_fraction,
        na.rm = TRUE), maximum_symmetry_residual = max(Z$symmetry_max_abs),
      maximum_reconstruction_residual = max(Z$reconstruction_max_abs),
      expected_spectrum_available = sum(Z$expected_spectrum_available),
      expected_spectrum_unavailable = sum(!Z$expected_spectrum_available),
      finite_Vp_ratio = length(ratio),
      Vp_ratio_min = if (length(ratio)) min(ratio) else NA_real_,
      Vp_ratio_median = if (length(ratio)) median(ratio) else NA_real_,
      Vp_ratio_max = if (length(ratio)) max(ratio) else NA_real_)
  }
  write.csv(do.call(rbind, state.aggregate), file.path(out,
    "null-state-matrix-aggregate.csv"), row.names = FALSE)
}
if (length(eigen.rows)) {
  eigenvalues <- do.call(rbind, eigen.rows)
  write.csv(eigenvalues, file.path(out, "null-state-eigenvalues.csv"),
    row.names = FALSE)
  write.csv(eigenvalues[eigenvalues$eigenvalue < 0, ], file.path(out,
    "null-negative-eigenvalues.csv"), row.names = FALSE)
}

stress.dir <- file.path(root, "stress")
groups <- if (dir.exists(stress.dir)) list.dirs(stress.dir, recursive = FALSE,
  full.names = FALSE) else character()
job.rows <- list()
fits <- list()
for (group in groups) {
  files <- list.files(file.path(stress.dir, group), pattern = "^job-[0-9]{3}\\.rds$",
    full.names = TRUE)
  rep.id <- as.integer(sub("job-([0-9]{3})\\.rds", "\\1", basename(files)))
  files <- files[rep.id <= max.rep]
  rep.id <- rep.id[rep.id <= max.rep]
  for (j in seq_along(files)) {
    z <- readRDS(files[j])
    z$metrics$group <- group
    job.rows[[length(job.rows) + 1L]] <- z$metrics
    fits[[paste(group, rep.id[j], sep = "/")]] <- z
  }
}
if (length(job.rows)) write.csv(do.call(rbind, job.rows), file.path(out,
  "stress-fit-summary.csv"), row.names = FALSE)

configs <- c("independent-w2-t1", "independent-w4-t1",
  "independent-w8-t1", "independent-w4-t2")
agreement <- list()
for (config in configs) {
  for (r in seq_len(min(10L, max.rep))) {
    ref.key <- paste("independent-w1-t1", r, sep = "/")
    key <- paste(config, r, sep = "/")
    if (!ref.key %in% names(fits) || !key %in% names(fits)) next
    ref <- fits[[ref.key]]
    z <- fits[[key]]
    agreement[[length(agreement) + 1L]] <- data.frame(config = config,
      gene = r, max_abs_u_difference = max(abs(z$u - ref$u)),
      rmse_u_difference = sqrt(mean((z$u - ref$u)^2)),
      max_abs_fixed_difference = max(abs(z$fixed$mean - ref$fixed$mean)),
      rmse_fixed_difference = sqrt(mean((z$fixed$mean - ref$fixed$mean)^2)),
      max_abs_theta_difference = max(abs(z$theta - ref$theta)),
      rmse_theta_difference = sqrt(mean((z$theta - ref$theta)^2)),
      log_size_difference = z$theta[1L] - ref$theta[1L],
      log_precision_difference = z$theta[2L] - ref$theta[2L])
  }
}
if (length(agreement)) write.csv(do.call(rbind, agreement), file.path(out,
  "stress-independent-agreement.csv"), row.names = FALSE)

joint.groups <- groups[grepl("^joint40-", groups)]
repeatability <- list()
for (group in joint.groups) {
  ref.key <- paste(group, 1L, sep = "/")
  if (!ref.key %in% names(fits)) next
  ref <- fits[[ref.key]]
  for (r in seq_len(min(10L, max.rep))) {
    key <- paste(group, r, sep = "/")
    if (!key %in% names(fits)) next
    z <- fits[[key]]
    repeatability[[length(repeatability) + 1L]] <- data.frame(group = group,
      replicate = r, max_abs_u_difference = max(abs(z$u - ref$u)),
      rmse_u_difference = sqrt(mean((z$u - ref$u)^2)),
      max_abs_fixed_difference = max(abs(z$fixed$mean - ref$fixed$mean)),
      rmse_fixed_difference = sqrt(mean((z$fixed$mean - ref$fixed$mean)^2)),
      max_abs_theta_difference = max(abs(z$theta - ref$theta)),
      rmse_theta_difference = sqrt(mean((z$theta - ref$theta)^2)),
      log_size_difference = z$theta[1L] - ref$theta[1L],
      log_precision_difference = z$theta[2L] - ref$theta[2L])
  }
}
if (length(repeatability)) write.csv(do.call(rbind, repeatability), file.path(out,
  "stress-joint-repeatability.csv"), row.names = FALSE)

jobs.file <- file.path(root, "stress-jobs.json")
if (file.exists(jobs.file)) {
  jobs <- jsonlite::read_json(jobs.file, simplifyVector = TRUE)
  expected.groups <- c("independent-w1-t1", "independent-w2-t1",
    "independent-w4-t1", "independent-w8-t1", "independent-w4-t2",
    "joint40-n5000-t4", "joint40-n97830-t4")
  expected.labels <- unlist(lapply(expected.groups, function(x)
    paste(x, seq_len(10L), sep = "-")))
  terminal.status <- c("completed", "failed", "resource_stop",
    "not_attempted_after_resource_limit")
  stress.final <- all(expected.labels %in% jobs$label) &&
    all(jobs$status[jobs$label %in% expected.labels] %in% terminal.status)
  stress.status <- if (stress.final) "final" else "snapshot_incomplete"
  jobs$summary_status <- stress.status
  jobs$peak_rss_gib <- jobs$peak_rss_bytes / 1024^3
  jobs$peak_private_gib <- jobs$peak_private_bytes / 1024^3
  write.csv(jobs, file.path(out, "stress-supervisor-jobs.csv"), row.names = FALSE)
  stress.completion <- data.frame(lane = c("independent", "joint", "all"),
    expected = c(50L, 20L, 70L),
    recorded = c(sum(jobs$group %in% expected.groups[1:5]),
      sum(jobs$group %in% expected.groups[6:7]),
      sum(jobs$group %in% expected.groups)),
    completed = c(sum(jobs$group %in% expected.groups[1:5] &
        jobs$status == "completed"),
      sum(jobs$group %in% expected.groups[6:7] & jobs$status == "completed"),
      sum(jobs$group %in% expected.groups & jobs$status == "completed")),
    summary_status = stress.status)
  write.csv(stress.completion, file.path(out, "stress-completion.csv"),
    row.names = FALSE)
}
groups.file <- file.path(root, "stress-groups.json")
if (file.exists(groups.file)) {
  G <- jsonlite::read_json(groups.file, simplifyVector = TRUE)
  if (exists("jobs", inherits = FALSE)) {
    completed <- vapply(G$group, function(x) sum(jobs$group == x &
      jobs$status == "completed"), integer(1L))
    failed <- vapply(G$group, function(x) sum(jobs$group == x &
      jobs$status == "failed"), integer(1L))
    resource.stop <- vapply(G$group, function(x) sum(jobs$group == x &
      jobs$status == "resource_stop"), integer(1L))
    not.attempted <- vapply(G$group, function(x) sum(jobs$group == x &
      jobs$status == "not_attempted_after_resource_limit"), integer(1L))
    G$completed_jobs <- completed
    G$failed_jobs <- failed
    G$resource_stops <- resource.stop
    G$not_attempted <- not.attempted
  } else {
    G$completed_jobs <- NA_integer_
    G$failed_jobs <- NA_integer_
    G$resource_stops <- NA_integer_
    G$not_attempted <- NA_integer_
  }
  G$peak_sum_rss_gib <- G$peak_sum_rss_bytes / 1024^3
  G$peak_sum_private_gib <- G$peak_sum_private_bytes / 1024^3
  G$min_system_available_gib <- G$min_system_available_bytes / 1024^3
  G$throughput_fits_per_minute <- G$completed_jobs / G$wall_seconds * 60
  G$summary_status <- if (exists("stress.status", inherits = FALSE))
    stress.status else "snapshot_unavailable"
  write.csv(G, file.path(out, "stress-supervisor-groups.csv"), row.names = FALSE)
}
