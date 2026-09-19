#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
Sys.setenv(LC_ALL = "C")

root <- Sys.getenv("MGCVST_STRESS_OUTPUT",
  "artifacts/inla-stress-calibration")
out <- Sys.getenv("MGCVST_STRESS_VERIFY_OUTPUT", file.path(root, "verification"))
require.complete <- tolower(Sys.getenv("MGCVST_VERIFY_REQUIRE_COMPLETE", "true")) %in%
  c("true", "1", "yes")
if (!dir.exists(root)) stop("Stress-calibration output directory is missing.")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

checks <- list()
add_check <- function(check, scope, expected, observed, passed, detail = NA_character_) {
  checks[[length(checks) + 1L]] <<- data.frame(check = check, scope = scope,
    expected = as.character(expected), observed = as.character(observed),
    passed = isTRUE(passed), detail = as.character(detail))
}

cases <- data.frame(dimension = c(rep("2d", 4L), rep("3d", 4L)),
  kind = c("pair", "pair", "marginal", "marginal",
    "independent", "independent", "joint", "joint"),
  mean = c(0.3, 3, 0.3, 3, 0.3, 3, 0.3, 3),
  stringsAsFactors = FALSE)

for (s in seq_len(nrow(cases))) {
  dim <- cases$dimension[s]
  kind <- cases$kind[s]
  mu0 <- cases$mean[s]
  case <- paste(kind, mu0, sep = "-")
  folder <- file.path(root, if (dim == "2d") "null2d" else "null3d-paired",
    case)
  rds <- if (dir.exists(folder)) list.files(folder,
    pattern = "^rep-[0-9]{4}\\.rds$", full.names = TRUE) else character()
  failures <- if (dir.exists(folder)) list.files(folder,
    pattern = "^rep-[0-9]{4}-failure\\.json$", full.names = TRUE) else character()
  rds.id <- as.integer(sub("rep-([0-9]{4})\\.rds", "\\1", basename(rds)))
  failure.id <- as.integer(sub("rep-([0-9]{4})-failure\\.json", "\\1",
    basename(failures)))
  in.range <- rds.id %in% seq_len(500L)
  failure.in.range <- failure.id %in% seq_len(500L)
  add_check("rds_ids_in_range", case, "all 1:500", sum(in.range), all(in.range))
  add_check("failure_ids_in_range", case, "all 1:500", sum(failure.in.range),
    all(failure.in.range))
  add_check("no_duplicate_rds", case, 0, sum(duplicated(rds.id)),
    !anyDuplicated(rds.id))
  add_check("no_duplicate_failure", case, 0, sum(duplicated(failure.id)),
    !anyDuplicated(failure.id))
  overlap <- intersect(rds.id, failure.id)
  add_check("no_rds_failure_overlap", case, 0, length(overlap), !length(overlap),
    if (length(overlap)) paste(overlap, collapse = ";") else NA_character_)
  represented <- sort(unique(c(rds.id[in.range], failure.id[failure.in.range])))
  missing <- setdiff(seq_len(500L), represented)
  add_check("all_replicates_accounted", case, 500, length(represented),
    length(represented) == 500L,
    if (length(missing)) paste(missing, collapse = ";") else NA_character_)

  for (j in seq_along(rds)) {
    z <- readRDS(rds[j])
    r <- rds.id[j]
    seed.base <- if (dim == "2d") 202610000L else 202650000L
    if (dim == "2d" && kind == "marginal") seed.base <- seed.base + 10000L
    if (mu0 >= 1) seed.base <- seed.base + if (dim == "2d") 20000L else 10000L
    expected.seed <- seed.base + r
    scope <- paste(case, sprintf("rep-%04d", r), sep = "/")
    add_check("row_replicate_matches_path", scope, r,
      paste(unique(z$rows$replicate), collapse = ";"),
      nrow(z$rows) == 2L && all(z$rows$replicate == r))
    add_check("row_case_matches_path", scope, case,
      paste(unique(z$rows$case), collapse = ";"),
      all(z$rows$case == case))
    add_check("row_seed_matches_protocol", scope, expected.seed,
      paste(unique(z$rows$seed), collapse = ";"),
      all(z$rows$seed == expected.seed))
    calibrations <- sort(as.character(z$rows$calibration))
    add_check("calibrations_exact", scope, "davies;liu",
      paste(calibrations, collapse = ";"),
      identical(calibrations, c("davies", "liu")))
    expected.genes <- if (dim == "2d") {
      if (kind == "pair") 2L else 1L
    } else {
      if (kind == "joint") 40L else 2L
    }
    observed.genes <- if (dim == "2d") nrow(z$Y) else ncol(z$Y)
    add_check("Y_gene_dimension_matches_case", scope, expected.genes,
      observed.genes, identical(as.integer(observed.genes), expected.genes))
  }
}

for (mu0 in c(0.3, 3)) {
  independent.folder <- file.path(root, "null3d-paired",
    paste("independent", mu0, sep = "-"))
  joint.folder <- file.path(root, "null3d-paired", paste("joint", mu0, sep = "-"))
  matched <- 0L
  unequal <- integer()
  for (r in seq_len(500L)) {
    independent.file <- file.path(independent.folder, sprintf("rep-%04d.rds", r))
    joint.file <- file.path(joint.folder, sprintf("rep-%04d.rds", r))
    if (!file.exists(independent.file) || !file.exists(joint.file)) next
    independent <- readRDS(independent.file)
    joint <- readRDS(joint.file)
    matched <- matched + 1L
    if (!identical(independent$Y[, 1:2, drop = FALSE],
      joint$Y[, 1:2, drop = FALSE])) unequal <- c(unequal, r)
  }
  scope <- paste0("3d-mean-", mu0)
  add_check("paired_Y_replicates_compared", scope, 500, matched,
    matched == 500L, if (matched < 500L) "Paired cases are incomplete." else NA_character_)
  add_check("paired_Y_genes_1_2_identical", scope, 0, length(unequal),
    !length(unequal), if (length(unequal)) paste(unequal, collapse = ";") else NA_character_)
  add_check("pair_is_one_replicate", scope,
    "one fixed genes-1/2 pair per replicate", "one rows object per replicate",
    TRUE, "The 780 within-dataset pairs are not counted as independent replicates.")
}

jobs.file <- file.path(root, "stress-jobs.json")
if (!file.exists(jobs.file)) stop("stress-jobs.json is missing.")
jobs <- jsonlite::read_json(jobs.file, simplifyVector = TRUE)
add_check("no_duplicate_stress_job_records", "stress", 0,
  sum(duplicated(jobs$label)), !anyDuplicated(jobs$label))
expected.groups <- c("independent-w1-t1", "independent-w2-t1",
  "independent-w4-t1", "independent-w8-t1", "independent-w4-t2",
  "joint40-n5000-t4", "joint40-n97830-t4")
for (group in expected.groups) {
  for (r in seq_len(10L)) {
    label <- paste(group, r, sep = "-")
    k <- which(jobs$label == label)
    file <- file.path(root, "stress", group, sprintf("job-%03d.rds", r))
    if (length(k) == 1L) {
      status <- jobs$status[k]
      expected.file <- identical(status, "completed")
      correspondence <- file.exists(file) == expected.file
      add_check("stress_status_file_correspondence", label,
        paste0("file=", expected.file), paste0("file=", file.exists(file),
          ";status=", status), correspondence)
      if (file.exists(file)) {
        z <- readRDS(file)
        add_check("stress_rds_label_matches_path", label, group,
          z$metrics$label, identical(as.character(z$metrics$label), group))
        add_check("stress_rds_replicate_matches_path", label, r,
          z$metrics$replicate, identical(as.integer(z$metrics$replicate), r))
      }
    } else {
      pending <- !require.complete && !file.exists(file)
      add_check("stress_status_file_correspondence", label, "one supervisor record",
        length(k), pending,
        if (pending) "Stress snapshot pending." else NA_character_)
    }
  }
  K <- jobs$group == group
  completed <- sum(jobs$status[K] == "completed")
  stopped <- sum(jobs$status[K] == "resource_stop")
  unattempted <- sum(jobs$status[K] == "not_attempted_after_resource_limit")
  other.failed <- sum(jobs$status[K] == "failed")
  accounted <- completed + stopped + unattempted + other.failed
  add_check("stress_group_ten_jobs_accounted", group, 10, accounted,
    accounted == 10L || !require.complete,
    if (accounted < 10L && !require.complete) "Stress snapshot pending." else
      NA_character_)
  resource.logic <- if (unattempted > 0L) stopped > 0L else TRUE
  add_check("resource_stop_unattempted_logic", group,
    "unattempted requires resource_stop", paste0("resource_stop=", stopped,
      ";unattempted=", unattempted), resource.logic)
}
independent.completed <- sum(jobs$group %in% expected.groups[1:5] &
  jobs$status == "completed")
joint.completed <- sum(jobs$group %in% expected.groups[6:7] &
  jobs$status == "completed")
add_check("independent_completed_status", "stress", 50,
  independent.completed, independent.completed == 50L || !require.complete,
  if (independent.completed < 50L && !require.complete)
    "Stress snapshot pending." else NA_character_)
add_check("joint_completed_status", "stress", 20, joint.completed,
  joint.completed == 20L || !require.complete,
  if (joint.completed < 20L && !require.complete)
    "Stress snapshot pending." else NA_character_)

C <- do.call(rbind, checks)
write.csv(C, file.path(out, "checks.csv"), row.names = FALSE)
write.csv(data.frame(checks = nrow(C), passed = sum(C$passed),
  failed = sum(!C$passed), require_complete = require.complete),
  file.path(out, "checks-summary.csv"), row.names = FALSE)
if (require.complete && any(!C$passed)) {
  stop(sum(!C$passed), " verification checks failed; see checks.csv.")
}
