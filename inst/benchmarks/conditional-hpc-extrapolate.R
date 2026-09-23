#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: conditional-hpc-extrapolate.R RESULTS_ROOT")
root <- args[1L]
G <- c(500L, 1000L, 2000L)
R1 <- do.call(rbind, lapply(G, function(n) {
  x <- read.csv(file.path(root, sprintf("timing-%d.csv", n)))
  if (nrow(x) != 1L || x$G != n) stop("Invalid timing file for G = ", n)
  x
}))

parse_rss <- function(x) {
  x <- trimws(x)
  x <- sub("[+].*$", "", x)
  value <- as.numeric(sub("[KMGTP]?$", "", x))
  unit <- sub("^[0-9.]+", "", x)
  factor <- stats::setNames(c(1, 1024, 1024^2, 1024^3, 1024^4, 1024^5),
                            c("", "K", "M", "G", "T", "P"))
  if (anyNA(value) || any(!unit %in% names(factor))) stop("Invalid sacct MaxRSS value.")
  value * factor[unit]
}

R1$max_rss_bytes <- vapply(G, function(n) {
  x <- read.delim(file.path(root, sprintf("sacct-%d.txt", n)),
                  sep = "|", header = TRUE, check.names = FALSE,
                  stringsAsFactors = FALSE)
  rss <- x$MaxRSS[!is.na(x$MaxRSS) & nzchar(x$MaxRSS)]
  if (!length(rss)) stop("No MaxRSS recorded for G = ", n)
  max(parse_rss(rss))
}, numeric(1L))
R1$disk_peak_bytes <- vapply(G, function(n) {
  as.numeric(readLines(file.path(root, sprintf("disk-peak-%d.txt", n)),
                       warn = FALSE)[1L])
}, numeric(1L))
if (any(!is.finite(R1$disk_peak_bytes))) stop("Invalid disk peak measurement.")

metrics <- c("total_seconds", "test_seconds", "basis_seconds", "score_seconds",
             "materialize_seconds", "variance_seconds", "remaining_test_seconds",
             "max_rss_bytes", "disk_peak_bytes")
predictions <- lapply(metrics, function(metric) {
  y <- R1[[metric]]
  model <- lm(y ~ G + I(G^2), data = data.frame(y = y, G = G))
  data.frame(metric = metric, observed_500 = y[1L], observed_1000 = y[2L],
             observed_2000 = y[3L], predicted_10000 =
               as.numeric(predict(model, newdata = data.frame(G = 10000L))),
             intercept = unname(coef(model)[1L]),
             linear = unname(coef(model)[2L]),
             quadratic = unname(coef(model)[3L]))
})
R2 <- do.call(rbind, predictions)
write.csv(R1, file.path(root, "scale-observed.csv"), row.names = FALSE)
write.csv(R2, file.path(root, "scale-extrapolated-10000.csv"), row.names = FALSE)
print(R2)
cat("Three points fit three coefficients exactly; these are planning estimates without uncertainty intervals.\n")
