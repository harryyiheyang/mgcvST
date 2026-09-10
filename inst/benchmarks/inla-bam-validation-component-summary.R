#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
Sys.setenv(LC_ALL = "C")

out <- Sys.getenv("MGCVST_COMPONENT_OUTPUT",
  "artifacts/inla-bam-validation/components")
datasets <- strsplit(Sys.getenv("MGCVST_COMPONENT_DATASETS",
  "unadjusted,celltype"), ",", fixed = TRUE)[[1L]]
for (dataset in datasets) {
  d.out <- file.path(out, dataset)
  bam.fit <- readRDS(file.path(d.out, "bam-compact-fit.rds"))
  inla.fit <- readRDS(file.path(d.out, "inla-fit.rds"))
  bam.pair <- readRDS(file.path(d.out, "bam-pair-tests.rds"))$results
  inla.pair <- readRDS(file.path(d.out, "inla-pair-tests.rds"))$results
  key <- c("pair_index", "feature1", "feature2")
  if (!identical(bam.pair[, key], inla.pair[, key])) {
    stop("Pair order differs between engines for ", dataset, ".")
  }
  bam.valid <- is.finite(bam.pair$p_two_sided) & bam.pair$p_two_sided >= 0 &
    bam.pair$p_two_sided <= 1
  inla.valid <- is.finite(inla.pair$p_two_sided) & inla.pair$p_two_sided >= 0 &
    inla.pair$p_two_sided <= 1
  valid <- bam.valid & inla.valid
  bam.conv <- setNames(bam.fit$diagnostics$converged,
    bam.fit$diagnostics$feature_id)
  inla.conv <- setNames(inla.fit$diagnostics$converged,
    inla.fit$diagnostics$feature_id)
  both <- bam.conv[bam.pair$feature1] & bam.conv[bam.pair$feature2] &
    inla.conv[bam.pair$feature1] & inla.conv[bam.pair$feature2]
  rows <- list()
  subsets <- list(all_requested = rep(TRUE, nrow(bam.pair)),
    both_engines_converged = both)
  for (label in names(subsets)) {
    ix <- subsets[[label]]
    iv <- ix & valid
    bp <- bam.pair$p_two_sided[iv]
    ip <- inla.pair$p_two_sided[iv]
    bd <- bam.pair$discovered[iv]
    id <- inla.pair$discovered[iv]
    rows[[label]] <- data.frame(subset = label, attempted = sum(ix),
      bam_valid = sum(ix & bam.valid), inla_valid = sum(ix & inla.valid),
      valid = sum(iv), bam_bh_denominator = sum(bam.valid),
      inla_bh_denominator = sum(inla.valid),
      bh_scope = "original_full_call_valid_pairs",
      bam_discoveries = sum(bd, na.rm = TRUE),
      inla_discoveries = sum(id, na.rm = TRUE),
      both_discoveries = sum(bd & id, na.rm = TRUE),
      union_discoveries = sum(bd | id, na.rm = TRUE),
      discovery_jaccard = sum(bd & id, na.rm = TRUE) /
        sum(bd | id, na.rm = TRUE),
      decision_agreement = mean(bd == id, na.rm = TRUE),
      direction_agreement = mean(sign(bam.pair$signed_score[iv]) ==
        sign(inla.pair$signed_score[iv]), na.rm = TRUE),
      signed_score_pearson = cor(bam.pair$signed_score[iv],
        inla.pair$signed_score[iv], use = "complete.obs"),
      signed_score_spearman = cor(bam.pair$signed_score[iv],
        inla.pair$signed_score[iv], method = "spearman", use = "complete.obs"),
      p_pearson = cor(bp, ip, use = "complete.obs"),
      p_spearman = cor(bp, ip, method = "spearman", use = "complete.obs"),
      median_abs_p_difference = median(abs(bp - ip)),
      max_abs_p_difference = max(abs(bp - ip)),
      raw_decision_agreement = mean((bp < 0.05) == (ip < 0.05)),
      neglog10p_pearson = cor(-log10(pmax(bp, .Machine$double.xmin)),
        -log10(pmax(ip, .Machine$double.xmin)), use = "complete.obs"))
  }
  write.csv(do.call(rbind, rows), file.path(d.out,
    "pair-engine-comparison.csv"), row.names = FALSE)

  original <- read.csv(file.path(d.out, "original-primary-component-edges.csv"))
  original.key <- paste(pmin(original$feature1, original$feature2),
    pmax(original$feature1, original$feature2), sep = "--")
  preserved <- list()
  pair.results <- list(bam = bam.pair, inla = inla.pair)
  for (engine in names(pair.results)) {
    current <- pair.results[[engine]]
    current.key <- paste(pmin(current$feature1, current$feature2),
      pmax(current$feature1, current$feature2), sep = "--")
    j <- match(original.key, current.key)
    if (anyNA(j)) stop("An original component edge is absent for ", dataset, ".")
    preserved[[engine]] <- data.frame(dataset = dataset, engine = engine,
      original_positive_edges = length(j),
      current_positive_discoveries = sum(current$discovered_positive[j],
        na.rm = TRUE),
      preservation_rate = mean(current$discovered_positive[j], na.rm = TRUE),
      direction_positive_rate = mean(current$signed_score[j] > 0, na.rm = TRUE),
      signed_score_pearson = cor(original$signed_score,
        current$signed_score[j], use = "complete.obs"))
  }
  write.csv(do.call(rbind, preserved), file.path(d.out,
    "original-edge-preservation-summary.csv"), row.names = FALSE)

  bam.marginal <- readRDS(file.path(d.out, "bam-marginal-tests.rds"))
  inla.marginal <- readRDS(file.path(d.out, "inla-marginal-tests.rds"))
  j <- match(bam.marginal$feature_id, inla.marginal$feature_id)
  if (anyNA(j)) stop("Marginal feature IDs differ for ", dataset, ".")
  bp <- bam.marginal$smooth.pvalue
  ip <- inla.marginal$p_value[j]
  valid <- is.finite(bp) & is.finite(ip)
  marginal <- data.frame(dataset = dataset, features = length(bp),
    valid = sum(valid), p_pearson = cor(bp, ip, use = "complete.obs"),
    p_spearman = cor(bp, ip, method = "spearman", use = "complete.obs"),
    neglog10p_pearson = cor(-log10(pmax(bp, .Machine$double.xmin)),
      -log10(pmax(ip, .Machine$double.xmin)), use = "complete.obs"),
    bam_bh05 = sum(p.adjust(bp, "BH") < 0.05, na.rm = TRUE),
    inla_bh05 = sum(p.adjust(ip, "BH") < 0.05, na.rm = TRUE))
  write.csv(marginal, file.path(d.out, "marginal-engine-summary.csv"),
    row.names = FALSE)

  bam.W <- readRDS(file.path(d.out, "bam-wgcna.rds"))$modules
  inla.W <- readRDS(file.path(d.out, "inla-wgcna.rds"))$modules
  j <- match(bam.W$feature_id, inla.W$feature_id)
  if (anyNA(j)) stop("WGCNA feature IDs differ for ", dataset, ".")
  tab <- table(bam.W$module, inla.W$module[j])
  ch <- function(x) x * (x - 1) / 2
  n <- sum(tab)
  expected.index <- sum(ch(rowSums(tab))) * sum(ch(colSums(tab))) / ch(n)
  ari <- (sum(ch(tab)) - expected.index) /
    (0.5 * (sum(ch(rowSums(tab))) + sum(ch(colSums(tab)))) - expected.index)
  bam.grey <- bam.W$module == 0L | bam.W$color == "grey"
  inla.grey <- inla.W$module == 0L | inla.W$color == "grey"
  modules <- data.frame(dataset = dataset, features = n,
    bam_modules = length(unique(bam.W$module[!bam.grey])),
    inla_modules = length(unique(inla.W$module[!inla.grey])),
    bam_grey = sum(bam.grey), inla_grey = sum(inla.grey),
    adjusted_rand_index = ari,
    exact_numeric_label_agreement = mean(bam.W$module == inla.W$module[j]))
  write.csv(modules, file.path(d.out, "module-engine-summary.csv"),
    row.names = FALSE)

  j <- match(bam.fit$feature_id, inla.fit$feature_id)
  size <- unlist(inla.fit$family_parameters[bam.fit$feature_id],
    use.names = FALSE)
  parameters <- data.frame(dataset = dataset,
    features = length(bam.fit$feature_id),
    lambda_pearson = cor(bam.fit$lambda, inla.fit$lambda[j]),
    lambda_spearman = cor(bam.fit$lambda, inla.fit$lambda[j],
      method = "spearman"),
    lambda_median_relative_difference = median(
      (inla.fit$lambda[j] - bam.fit$lambda) / bam.fit$lambda),
    dispersion_bam_min = min(bam.fit$dispersion),
    dispersion_bam_max = max(bam.fit$dispersion),
    dispersion_inla_min = min(inla.fit$dispersion),
    dispersion_inla_max = max(inla.fit$dispersion),
    inla_nb_size_median = median(size),
    inla_nb_size_q25 = quantile(size, 0.25, names = FALSE),
    inla_nb_size_q75 = quantile(size, 0.75, names = FALSE),
    inla_nb_size_finite = sum(is.finite(size)))
  write.csv(parameters, file.path(d.out, "fit-parameter-summary.csv"),
    row.names = FALSE)

  E <- read.csv(file.path(d.out, "original-primary-component-edges.csv"))
  ek <- paste(pmin(E$feature1, E$feature2), pmax(E$feature1, E$feature2),
    sep = "--")
  preservation <- list()
  for (engine in c("bam", "inla")) {
    P <- if (engine == "bam") bam.pair else inla.pair
    pk <- paste(pmin(P$feature1, P$feature2), pmax(P$feature1, P$feature2),
      sep = "--")
    j <- match(ek, pk)
    if (anyNA(j)) stop("An original component edge is absent from ", dataset, ".")
    preservation[[engine]] <- data.frame(dataset = dataset, engine = engine,
      original_positive_edges = nrow(E),
      current_positive_discoveries = sum(P$discovered_positive[j]),
      preservation_rate = mean(P$discovered_positive[j]),
      direction_positive_rate = mean(P$signed_score[j] > 0),
      signed_score_pearson = cor(E$signed_score, P$signed_score[j]))
  }
  write.csv(do.call(rbind, preservation),
    file.path(d.out, "original-edge-preservation-summary.csv"), row.names = FALSE)
}
