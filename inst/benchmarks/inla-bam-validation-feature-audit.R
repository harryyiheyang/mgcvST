Sys.setenv(LC_ALL = "C")

repo <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
root <- file.path(repo, "artifacts/inla-bam-validation/components")
out <- file.path(repo, "artifacts/inla-bam-validation/summary/real-feature-audit.csv")
datasets <- c("celltype", "unadjusted")
rows <- vector("list", 4L)
ids <- vector("list", 2L)
position <- 0L

for (dataset in datasets) {
  path <- file.path(root, dataset)
  input <- readRDS(file.path(path, "input.rds"))
  bam <- readRDS(file.path(path, "bam-compact-fit.rds"))
  inla <- readRDS(file.path(path, "inla-fit.rds"))
  bam_marginal <- readRDS(file.path(path, "bam-marginal-tests.rds"))
  ids[[dataset]] <- as.character(input$feature_id)

  same_input_bam <- identical(ids[[dataset]], as.character(bam$feature_id))
  same_input_inla <- identical(ids[[dataset]], as.character(inla$feature_id))
  same_bam_inla <- identical(as.character(bam$feature_id), as.character(inla$feature_id))

  bam_p <- bam_marginal$smooth.pvalue[match(ids[[dataset]], bam_marginal$feature_id)]
  inla_p <- inla$diagnostics$marginal_p_value[
    match(ids[[dataset]], inla$diagnostics$feature_id)
  ]
  both_valid <- is.finite(bam_p) & is.finite(inla_p)
  agreement <- (bam_p < 0.05) == (inla_p < 0.05)

  bam_size <- rep(NA_real_, length(ids[[dataset]]))
  size_source <- "not retained in compact fit"
  if (dataset == "celltype") {
    payload <- readRDS(file.path(path, "bam-payloads.rds"))
    bam_size <- vapply(payload, function(x) x$W$family_parameters[1L], numeric(1L))
    size_source <- "bam-payloads.rds"
    rm(payload)
  }
  inla_size <- vapply(inla$family_parameters, function(x) x[1L], numeric(1L))

  fit_list <- list(bam = bam, inla = inla)
  size_list <- list(bam = bam_size, inla = inla_size)
  size_sources <- c(bam = size_source, inla = "inla-fit.rds family_parameters")

  for (backend in c("bam", "inla")) {
    fit <- fit_list[[backend]]
    size <- size_list[[backend]]
    diagnostic <- fit$diagnostics
    lambda <- as.numeric(fit$lambda)
    finite_lambda <- lambda[is.finite(lambda)]
    dispersion <- as.numeric(fit$dispersion)
    finite_dispersion <- dispersion[is.finite(dispersion)]
    available <- is.finite(dispersion) & dispersion > 0 &
      is.finite(lambda) & lambda > 0 &
      colSums(!is.finite(fit$working_error)) == 0L &
      colSums(!is.finite(fit$working_variance)) == 0L
    finite_size <- size[is.finite(size)]
    position <- position + 1L
    rows[[position]] <- data.frame(
      dataset = dataset,
      backend = backend,
      feature_count = length(fit$feature_id),
      same_ids_input_bam = same_input_bam,
      same_ids_input_inla = same_input_inla,
      same_ids_bam_inla = same_bam_inla,
      converged_count = sum(diagnostic$converged %in% TRUE),
      nonconverged_count = sum(!(diagnostic$converged %in% TRUE)),
      error_count = sum(!is.na(diagnostic$error_message) & nzchar(diagnostic$error_message)),
      compact_feature_available_count = sum(available),
      nonconverged_compact_available_count = sum(
        available & !(diagnostic$converged %in% TRUE)
      ),
      lambda_finite_count = length(finite_lambda),
      lambda_min = min(finite_lambda),
      lambda_median = median(finite_lambda),
      lambda_max = max(finite_lambda),
      dispersion_finite_count = length(finite_dispersion),
      dispersion_min = min(finite_dispersion),
      dispersion_median = median(finite_dispersion),
      dispersion_max = max(finite_dispersion),
      nb_size_source = size_sources[[backend]],
      nb_size_finite_count = length(finite_size),
      nb_size_min = if (length(finite_size)) min(finite_size) else NA_real_,
      nb_size_median = if (length(finite_size)) median(finite_size) else NA_real_,
      nb_size_max = if (length(finite_size)) max(finite_size) else NA_real_,
      nb_size_above_1e6 = sum(size > 1e6, na.rm = TRUE),
      nb_size_above_1e12 = sum(size > 1e12, na.rm = TRUE),
      bam_marginal_valid = sum(is.finite(bam_p)),
      inla_marginal_valid = sum(is.finite(inla_p)),
      marginal_both_valid = sum(both_valid),
      marginal_both_valid_and_converged = sum(
        both_valid & bam$diagnostics$converged & inla$diagnostics$converged
      ),
      bam_marginal_raw_rejected = sum(bam_p < 0.05, na.rm = TRUE),
      inla_marginal_raw_rejected = sum(inla_p < 0.05, na.rm = TRUE),
      marginal_raw_decision_agree = sum(agreement & both_valid),
      marginal_raw_decision_agreement_rate = mean(agreement[both_valid]),
      marginal_converged_raw_decision_agree = sum(
        agreement & both_valid & bam$diagnostics$converged & inla$diagnostics$converged
      ),
      stringsAsFactors = FALSE
    )
  }
  rm(input, bam, inla, bam_marginal, fit_list)
  gc()
}

result <- do.call(rbind, rows)
result$cross_dataset_common_features <- length(intersect(ids$celltype, ids$unadjusted))
stopifnot(all(result$same_ids_input_bam))
stopifnot(all(result$same_ids_input_inla))
stopifnot(all(result$same_ids_bam_inla))
stopifnot(all(result$cross_dataset_common_features == 443L))
write.csv(result, out, row.names = FALSE)
