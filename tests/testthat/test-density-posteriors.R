# With K = 1 and a single view the sampler only updates the component
# parameters, so their posterior is available in closed form.

run_single_component <- function(X, type, R = 6000) {
  rownames(X) <- seq_len(nrow(X))
  callMDI(list(X), R = R, thin = 1, types = type, K = 1, initial_labels = matrix(1, nrow(X), 1),
          initial_labels_as_intended = FALSE)
}

test_that("diagonal Gaussian parameters follow the normal-inverse-gamma posterior", {
  skip_on_cran()
  set.seed(11)
  N <- 25; P <- 2
  X <- cbind(rnorm(N, 2, 1.5), rnorm(N, -1, 0.5))
  hp <- mdir:::densityHyperparameters(X, 1, 0, numeric(0))
  out <- run_single_component(X, "G")
  theta <- out$parameters[[1]][-(1:200), , drop = FALSE]     # mu (P), variances (P)
  n <- N; kappa_n <- hp$kappa + n; nu_n <- hp$nu + n
  xbar <- colMeans(X)
  mu_n <- (hp$kappa * hp$xi + n * xbar) / kappa_n
  ss <- colSums(sweep(X, 2, xbar)^2)
  scale_n <- hp$scale + ss + n * hp$kappa / kappa_n * (xbar - hp$xi)^2
  for (p in 1:P) {
    ev <- scale_n[p] / (nu_n - 2)                                # mean of InvGamma(nu_n/2, scale_n/2)
    expect_equal(mean(theta[, P + p]), ev, tolerance = 0.05)
    expect_lt(abs(mean(theta[, p]) - mu_n[p]), 0.05 * sqrt(ev))
    # Var(mu) = E[sigma^2] / kappa_n; the old code used variance / kappa as an sd
    expect_equal(var(theta[, p]), ev / kappa_n, tolerance = 0.1)
  }
})

test_that("MVN parameters follow the normal-inverse-Wishart posterior", {
  skip_on_cran()
  set.seed(12)
  N <- 30; P <- 3
  Sig <- matrix(c(1, .5, .2, .5, 2, -.3, .2, -.3, .7), P)
  X <- matrix(rnorm(N * P), N) %*% chol(Sig) + rep(c(1, 0, -1), each = N)
  hp <- mdir:::densityHyperparameters(X, 1, 1, numeric(0))
  out <- run_single_component(X, "MVN")
  theta <- out$parameters[[1]][-(1:200), , drop = FALSE]
  n <- N; kappa_n <- hp$kappa + n; nu_n <- hp$nu + n
  xbar <- colMeans(X)
  mu_n <- (hp$kappa * hp$xi + n * xbar) / kappa_n
  S <- crossprod(sweep(X, 2, xbar))
  scale_n <- hp$scale + S + n * hp$kappa / kappa_n * tcrossprod(xbar - hp$xi)
  Ecov <- scale_n / (nu_n - P - 1)
  cov_draws <- theta[, P + seq_len(P * P), drop = FALSE]
  expect_equal(matrix(colMeans(cov_draws), P), Ecov, tolerance = 0.05)
  expect_true(all(abs(colMeans(theta[, 1:P]) - as.numeric(mu_n)) < 0.05 * sqrt(diag(Ecov))))
})

test_that("categorical probabilities follow the Dirichlet posterior", {
  skip_on_cran()
  set.seed(13)
  N <- 40
  X <- cbind(sample(0:2, N, TRUE, c(.6, .3, .1)), sample(0:1, N, TRUE))
  hp <- mdir:::densityHyperparameters(X, 1, 2, numeric(0))
  out <- run_single_component(X, "C")
  theta <- out$parameters[[1]][-(1:200), , drop = FALSE]
  start <- 0
  for (p in 1:2) {
    ncat <- hp$n_cat[p]
    alpha_n <- as.numeric(hp$concentration[[p]]) + tabulate(X[, p] + 1, ncat)
    expect_equal(unname(colMeans(theta[, start + seq_len(ncat)])), as.numeric(alpha_n / sum(alpha_n)), tolerance = 0.03)
    start <- start + ncat
  }
})

test_that("empirical-Bayes hyperparameters use only observed entries", {
  X <- cbind(c(1, 2, 3, 4, NA, NA), c(NA, 2, 4, 6, 8, 10))
  for (type in c(0, 1)) {
    hp <- mdir:::densityHyperparameters(X, 2, type, numeric(0))
    expect_true(all(is.finite(unlist(hp))))
    expect_equal(as.numeric(hp$xi), c(2.5, 6))
  }
  # a column that is missing in the first row and complete elsewhere (the case
  # the old complete-row detection got wrong)
  X2 <- cbind(c(NA, 1:9), rnorm(10))
  expect_true(all(is.finite(mdir:::densityHyperparameters(X2, 3, 0, numeric(0))$scale)))
})
