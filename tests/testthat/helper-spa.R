# Saddlepoint (Lugannani-Rice) calibration shared by the exact and the
# PCAlearning pair routes. The R reference below is the research-script form
# of the same approximation: a bracketed root of K'(t) = x and the tail
# Phibar(w) + phi(w) (1 / v - 1 / w) in log space.

.spa_ref <- function(x, s, m = numeric(0), u = numeric(0), g = 0) {
  K_parts <- function(t) {
    s2 <- s^2
    c(K = -0.5 * sum(log1p(-s2 * t^2)) - 0.5 * sum(m * log1p(-u * t^2)) + 0.5 * g * t^2,
      K1 = sum(s2 * t / (1 - s2 * t^2)) + sum(m * u * t / (1 - u * t^2)) + g * t,
      K2 = sum(s2 * (1 + s2 * t^2) / (1 - s2 * t^2)^2) +
        sum(m * u * (1 + u * t^2) / (1 - u * t^2)^2) + g)
  }
  smax <- sqrt(max(c(s^2, u)))
  th <- uniroot(function(t) K_parts(t)[["K1"]] - x, c(0, (1 / smax) * (1 - 1e-13)),
                tol = 1e-15)$root
  k <- K_parts(th)
  w <- sqrt(max(2 * (th * x - k[["K"]]), 0))
  v <- th * sqrt(k[["K2"]])
  lq <- pnorm(w, lower.tail = FALSE, log.p = TRUE)
  min(0, log(2) + lq + log1p(exp(dnorm(w, log = TRUE) - lq) * (1 / v - 1 / w)))
}

# Remainder of power sums mu_1..mu_4 as nodes (u, m) in the form of the
# research scripts: two Gauss nodes, else one Satterthwaite node, else a
# Gaussian term.
.spa_rem_one <- function(mu) {
  if (mu[1] > 0 && mu[2] > 0) list(m = mu[1]^2 / mu[2], u = mu[2] / mu[1], g = 0, kind = 1L)
  else list(m = numeric(0), u = numeric(0), g = max(mu[1], 0), kind = 3L)
}
.spa_rem_two <- function(mu) {
  H <- matrix(c(mu[1], mu[2], mu[2], mu[3]), 2L)
  if (any(mu <= 0) || det(H) <= 0) return(.spa_rem_one(mu))
  ab <- solve(H, -c(mu[3], mu[4]))
  x <- Re(polyroot(c(ab[1L], ab[2L], 1)))
  if (any(!is.finite(x)) || any(x <= 0)) return(.spa_rem_one(mu))
  w <- solve(rbind(c(1, 1), x), mu[1:2])
  if (any(w <= 0)) return(.spa_rem_one(mu))
  list(m = w / x, u = x, g = 0, kind = 2L)
}
.spa_powers <- function(s) vapply(1:4, function(r) sum(s^(2 * r)), numeric(1L))

.spa_kernel <- function(x, s, Tm, order = 4L) {
  mgcvST:::mgcvst_spa_cpp(x, matrix(s, ncol = 1L), matrix(Tm, 4L, 1L), order, 1L)
}

.spa_sqrtm <- function(M) {
  E <- eigen((M + t(M)) / 2, symmetric = TRUE)
  E$vectors %*% (sqrt(pmax(E$values, 0)) * t(E$vectors))
}

