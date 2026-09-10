# Small, reproducible fitted-null check. This is not a tail-calibration study.
suppressPackageStartupMessages(library(mgcvST))
reps <- as.integer(Sys.getenv("MGCVST_INLA_NULL_REPS", "30"))
out <- Sys.getenv("MGCVST_INLA_OUTPUT", "inla-null-results")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
stopifnot(is.finite(reps), reps >= 10L)
set.seed(91419)
vertices <- as.matrix(expand.grid(x = seq(0, 1, length.out = 5L),
                                 y = seq(0, 1, length.out = 5L)))
mesh <- list(loc = vertices, graph = list(tv = geometry::delaunayn(vertices)))
d <- data.frame(x = runif(100L, .01, .99), y = runif(100L, .01, .99))
basis <- spde_basis(mesh, as.matrix(d), kappa = .7, project_intercept = TRUE)
model <- inlaST.set(response ~ 1, d, basis, family = gaussian())
tau <- 2
phi <- .25
# Q = R'R, so solve(R, epsilon) has covariance Q^{-1}.
innovation <- matrix(rnorm(ncol(basis$B) * 2L * reps), ncol(basis$B))
u <- backsolve(chol(tau * basis$Q), innovation)
signal <- basis$B %*% u
Y <- t(1 + signal + matrix(rnorm(length(signal), sd = sqrt(phi)), nrow(signal)))
rownames(Y) <- paste0("feature", seq_len(nrow(Y)))
fit <- inlaST.estimate(
  Y, model, control = list(fixed_precision = tau, gaussian_precision = 1 / phi),
  marginal_args = list(method = "liu"), BPPARAM = BiocParallel::SerialParam()
)
stopifnot(all(fit$diagnostics$converged))
pairs <- matrix(seq_len(2L * reps), ncol = 2L, byrow = TRUE)
test <- mgcvST.test(fit, pairs = pairs, calibration = "davies", threads = 1L)
p <- test$results$p_two_sided
stopifnot(length(p) == reps, all(is.finite(p)))
rejected <- sum(p < .05)
ci <- stats::binom.test(rejected, reps)$conf.int
summary <- data.frame(replicates = reps, rejected_at_005 = rejected,
                      rejection_rate = rejected / reps,
                      binomial_ci_lower = ci[1L], binomial_ci_upper = ci[2L],
                      median_p = median(p), estimate_seconds = fit$timing$elapsed)
utils::write.csv(test$results, file.path(out, "null-pairs.csv"), row.names = FALSE)
utils::write.csv(summary, file.path(out, "null-summary.csv"), row.names = FALSE)
writeLines(c(
  "Independent Gaussian fields simulated and fitted with identical fixed hyperparameters.",
  "All spatial fields obey the observation mean-zero constraint.",
  "This small smoke experiment does not establish calibration after hyperparameter estimation,",
  "for negative-binomial observations, or in extreme multiple-testing tails."
), file.path(out, "null-scope.txt"))
print(summary)
