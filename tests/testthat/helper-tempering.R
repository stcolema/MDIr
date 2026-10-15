# Exact reference posteriors for tiny categorical models at inverse temperature
# beta (see verification/tempering/exact_reference.R for the derivation).

tp_lbeta <- function(a) sum(lgamma(a)) - lgamma(sum(a))

tp_log_marg_cat <- function(Xs, alpha, beta) {
  out <- 0
  for (p in seq_along(alpha)) {
    cnt <- tabulate(Xs[, p] + 1, length(alpha[[p]]))
    out <- out + tp_lbeta(alpha[[p]] + beta * cnt) - tp_lbeta(alpha[[p]])
  }
  out
}

# Label prior of one view (L = 1): Dirichlet-multinomial with concentration
# mass / K, integrated over mass ~ Gamma(2, 0.1); depends on the counts only
tp_prior_counts <- function(Nk, K, mass_shape = 2, mass_rate = 0.1) {
  N <- sum(Nk)
  f <- function(m) {
    a <- m / K
    exp(lgamma(m) - lgamma(m + N) + sum(lgamma(a + Nk) - lgamma(a)) + dgamma(m, mass_shape, mass_rate, log = TRUE))
  }
  stats::integrate(Vectorize(f), 0, Inf, rel.tol = 1e-10, subdivisions = 2000L)$value
}

tp_exact_L1 <- function(X, K, alpha, beta) {
  N <- nrow(X)
  grid <- as.matrix(expand.grid(rep(list(0:(K - 1)), N)))
  cache <- new.env()
  lp <- numeric(nrow(grid))
  for (i in seq_len(nrow(grid))) {
    c_i <- grid[i, ]
    Nk <- tabulate(c_i + 1, K)
    key <- paste(Nk, collapse = "_")
    if (is.null(cache[[key]])) cache[[key]] <- log(tp_prior_counts(Nk, K))
    ll <- 0
    for (k in 0:(K - 1)) {
      idx <- which(c_i == k)
      if (length(idx)) ll <- ll + tp_log_marg_cat(X[idx, , drop = FALSE], alpha, beta)
    }
    lp[i] <- cache[[key]] + ll
  }
  p <- exp(lp - max(lp))
  list(grid = grid, prob = p / sum(p))
}

# Rejection rate of an exchange between temperatures b1 < b2 for the Gaussian
# path pi_b proportional to N(x; 0, 1) exp(b * l(x)), l(x) = -(x - m)^2 / (2 s2),
# by Monte Carlo over independent draws from the two tempered distributions
tp_gauss_rejection <- function(b1, b2, m = 3, s2 = 0.05, n = 4e5) {
  draw <- function(b, n) {
    prec <- 1 + b / s2
    rnorm(n, (b * m / s2) / prec, sqrt(1 / prec))
  }
  l <- function(x) -(x - m)^2 / (2 * s2)
  x1 <- draw(b1, n); x2 <- draw(b2, n)
  1 - mean(pmin(1, exp((b1 - b2) * (l(x2) - l(x1)))))
}
