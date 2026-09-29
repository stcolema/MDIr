# Marginal posterior of the GP hyperparameters for a single component. After
# integrating out the mean function mu ~ N(0, K), with xbar the mean of the n items
# and SS the within-item sum of squares,
#   p(X | a, l, s2) = (2 pi s2)^(-nP/2) (2 pi s2 / n)^(P/2) exp(-SS / (2 s2)) N(xbar; 0, K + s2 / n I)
# and each hyperparameter has a log-normal(0, 1) prior.
gp_log_marginal <- function(X, a, l, s2, jitter = 1e-8) {
  n <- nrow(X); P <- ncol(X)
  xbar <- colMeans(X)
  SS <- sum(sweep(X, 2, xbar)^2)
  D <- outer(seq_len(P), seq_len(P), function(i, j) -0.5 * (i - j)^2)
  Kmat <- a * exp(D / l) + diag(jitter * a, P)
  S <- Kmat + diag(s2 / n, P)
  ch <- chol(S)
  z <- backsolve(ch, xbar, transpose = TRUE)
  ll_bar <- -0.5 * (P * log(2 * pi) + 2 * sum(log(diag(ch))) + sum(z^2))
  -0.5 * n * P * log(2 * pi * s2) + 0.5 * P * log(2 * pi * s2 / n) - SS / (2 * s2) + ll_bar
}

test_that("GP hyperparameter sampler targets the marginal posterior", {
  skip_on_cran()
  set.seed(61)
  P <- 4; n <- 15
  mu <- as.numeric(t(chol(outer(1:P, 1:P, function(i, j) 0.8 * exp(-0.5 * (i - j)^2 / 2)))) %*% rnorm(P))
  X <- matrix(rnorm(n * P, mu, 0.5), n, P, byrow = TRUE)
  rownames(X) <- seq_len(n)

  grid <- seq(-4, 3, length.out = 36)
  eg <- expand.grid(la = grid, ll = grid, ls = grid)
  lp <- mapply(function(la, ll, ls) {
    gp_log_marginal(X, exp(la), exp(ll), exp(ls)) + dnorm(la, log = TRUE) + dnorm(ll, log = TRUE) + dnorm(ls, log = TRUE)
  }, eg$la, eg$ll, eg$ls)
  w <- exp(lp - max(lp)); w <- w / sum(w)
  grid_mean <- c(sum(w * eg$la), sum(w * eg$ll), sum(w * eg$ls))
  grid_sd <- sqrt(c(sum(w * eg$la^2), sum(w * eg$ll^2), sum(w * eg$ls^2)) - grid_mean^2)

  out <- callMDI(list(X), R = 40000, thin = 5, types = "GP", K = 1)
  h <- out$hypers[[1]]
  keep <- -(1:200)
  draws <- log(cbind(h$amplitude[keep, 1], h$length[keep, 1], h$noise[keep, 1]))
  # the hyperparameters change only every fifth sweep: drop repeated rows
  draws <- draws[c(TRUE, rowSums(abs(diff(draws))) > 0), , drop = FALSE]
  expect_gt(nrow(draws), 1000)
  ess <- apply(draws, 2, function(x) rankNormalizedRhat(cbind(x[seq(1, length(x) - 1, 2)], x[seq(2, length(x), 2)][seq_len(length(x) %/% 2)][seq_len(length(seq(1, length(x) - 1, 2)))]))$ess_bulk)
  for (j in 1:3) {
    mcse <- grid_sd[j] / sqrt(ess[j])
    expect_lt(abs(mean(draws[, j]) - grid_mean[j]), 5 * mcse + 0.05 * grid_sd[j],
              label = paste("hyperparameter", j))
    expect_equal(sd(draws[, j]), grid_sd[j], tolerance = 0.15)
  }
})
