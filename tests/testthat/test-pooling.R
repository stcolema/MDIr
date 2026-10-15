# Pooled variance scale. With one occupied component the marginal posterior of the
# scale s given the data is available in closed form:
#   diagonal Gaussian, per feature: p(s | x) proportional to Gamma(s; a, a / c) s^(nu / 2) (s + Q)^(-(nu + n) / 2),
#   MVN: p(s | x) proportional to prod_p Gamma(s_p; a, a / c_p) |Psi|^(nu / 2) |Psi + Q|^(-(nu + n) / 2), Psi = diag(s),
# where Q is the within-component scatter plus kappa n / (kappa + n) (xbar - xi)(xbar - xi)'.

fit_fixed_component <- function(X, type, K, R = 40000) {
  rownames(X) <- seq_len(nrow(X))
  # every item labelled as component 1: any further component stays empty
  callMDI(list(X), R = R, thin = 4, types = type, K = K,
          initial_labels = matrix(1, nrow(X), 1), fixed = matrix(1, nrow(X), 1))
}

test_that("pooled scale of the diagonal Gaussian matches its marginal posterior (also with an empty component)", {
  skip_on_cran()
  set.seed(81)
  n <- 12
  X <- cbind(rnorm(n, 1, 2), rnorm(n, -1, 0.4))
  hp <- mdir:::densityHyperparameters(X, 1, 0, numeric(0))
  for (K in c(1, 2)) {
    hp <- mdir:::densityHyperparameters(X, K, 0, numeric(0))
    a <- hp$scale_pooling_shape; nu <- hp$nu; kappa <- hp$kappa
    out <- fit_fixed_component(X, "G", K)
    s_draws <- out$pooled_hyperparameters[[1]][-(1:200), , drop = FALSE]
    for (p in 1:2) {
      xbar <- mean(X[, p]); Q <- sum((X[, p] - xbar)^2) + kappa * n / (kappa + n) * (xbar - hp$xi[p])^2
      c_p <- hp$scale_prior_mean[p]
      grid <- exp(seq(log(c_p) - 6, log(c_p) + 4, length.out = 4000))
      ld <- dgamma(grid, a, a / c_p, log = TRUE) + 0.5 * nu * log(grid) - 0.5 * (nu + n) * log(grid + Q) + log(grid)   # + log(grid): integrate on the log scale
      mom <- grid_moments(log(grid), ld)
      lg <- log(s_draws[, p])
      expect_lt(abs(mean(lg) - mom["mean"]), 5 * mom["sd"] / sqrt(ess_of(lg)) + 0.02, label = paste("K", K, "feature", p, "mean"))
      expect_equal(sd(lg), unname(mom["sd"]), tolerance = 0.1, label = paste("K", K, "feature", p, "sd"))
    }
  }
})

test_that("pooled scale of the MVN matches its marginal posterior", {
  skip_on_cran()
  set.seed(82)
  n <- 15; P <- 2
  X <- matrix(rnorm(n * P), n) %*% chol(matrix(c(1.5, .6, .6, 0.5), 2)) + 0.3
  hp <- mdir:::densityHyperparameters(X, 1, 1, numeric(0))
  a <- hp$scale_pooling_shape; nu <- hp$nu; kappa <- hp$kappa
  xbar <- colMeans(X); S <- crossprod(sweep(X, 2, xbar))
  Q <- S + kappa * n / (kappa + n) * tcrossprod(xbar - hp$xi)
  c_p <- hp$scale_prior_mean
  g <- seq(-6, 4, length.out = 220)
  eg <- expand.grid(l1 = log(c_p[1]) + g, l2 = log(c_p[2]) + g)
  ld <- mapply(function(l1, l2) {
    s <- exp(c(l1, l2)); Psi <- diag(s)
    sum(dgamma(s, a, a / c_p, log = TRUE)) + 0.5 * nu * sum(log(s)) -
      0.5 * (nu + n) * determinant(Psi + Q)$modulus + l1 + l2
  }, eg$l1, eg$l2)
  out <- fit_fixed_component(X, "MVN", 1)
  s_draws <- out$pooled_hyperparameters[[1]][-(1:200), , drop = FALSE]
  for (p in 1:2) {
    mom <- grid_moments(eg[[p]], ld)
    lg <- log(s_draws[, p])
    expect_lt(abs(mean(lg) - mom["mean"]), 5 * mom["sd"] / sqrt(ess_of(lg)) + 0.02, label = paste("feature", p))
    expect_equal(sd(lg), unname(mom["sd"]), tolerance = 0.1)
  }
})

test_that("pooling can be switched off and the fixed scale is the empirical value", {
  X <- matrix(rnorm(60), 30); rownames(X) <- 1:30
  fit <- callMDI(list(X), R = 40, thin = 4, types = "MVN", K = 3, density_prior = densityPrior(scale_pool_shape = 0))
  hp <- mdir:::densityHyperparameters(X, 3, 1, c(0, 1, 1, 2, 1))
  expect_true(all(fit$pooled_hyperparameters[[1]] == rep(hp$scale_prior_mean, each = nrow(fit$pooled_hyperparameters[[1]]))))
  pooled <- callMDI(list(X), R = 40, thin = 4, types = "MVN", K = 3)
  expect_gt(var(pooled$pooled_hyperparameters[[1]][, 1]), 0)
})

test_that("pooling shrinks the scales of small components toward the shared scale", {
  skip_on_cran()
  # Two large tight clusters and one tiny cluster with 3 items: with pooling the
  # small component's variance is pulled toward the population value
  set.seed(83)
  X <- cbind(c(rnorm(40, 0, 0.5), rnorm(40, 6, 0.5), rnorm(3, 12, 1.5)))
  rownames(X) <- seq_len(nrow(X))
  labels <- matrix(rep(1:3, c(40, 40, 3)), ncol = 1)
  post_var <- function(dp) {
    fit <- callMDI(list(X), R = 6000, thin = 3, types = "G", K = 3, initial_labels = labels,
                   initial_labels_as_intended = TRUE, fixed = matrix(1, nrow(X), 1), density_prior = dp)
    th <- fit$parameters[[1]][-(1:200), , drop = FALSE]     # mu (3), variances (3)
    colMeans(th[, 4:6])
  }
  pooled <- post_var(densityPrior(scale_pool_shape = 2))
  fixed <- post_var(densityPrior(scale_pool_shape = 0))
  # the tiny component's variance estimate is closer to the big components' with pooling
  expect_lt(abs(log(pooled[3]) - log(pooled[1])), abs(log(fixed[3]) - log(fixed[1])) + 1e-8)
  expect_true(all(is.finite(pooled)))
})
