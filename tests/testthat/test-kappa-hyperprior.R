# Hierarchical (log-uniform) prior on the precision kappa of the component-mean prior,
#   mu_k | Sigma_k ~ N(xi, Sigma_k / kappa),  log kappa ~ Uniform(log lo, log hi).
# With one occupied component and a fixed variance scale, the marginal posterior of kappa given the data is
# p(kappa | x) proportional to the collapsed marginal likelihood m(x; kappa), which is available in closed form:
#   diagonal Gaussian: prod_p [kappa / (kappa + n)]^(1/2) (s_p + SS_p + kappa n / (kappa + n) d_p^2)^(-(nu + n) / 2)
#   MVN: [kappa / (kappa + n)]^(P / 2) |Psi + S + kappa n / (kappa + n) d d'|^(-(nu + n) / 2),  Psi = diag(s),
# with d = xbar - xi.

fit_kappa_component <- function(X, type, K, dp, R = 40000) {
  rownames(X) <- seq_len(nrow(X))
  callMDI(list(X), n_iter = R, thin = 4, types = type, K_max = K, density_prior = dp,
          initial_labels = matrix(1, nrow(X), 1), fixed = matrix(1, nrow(X), 1))
}

test_that("densityPrior validates and prints the kappa options", {
  expect_equal(unname(densityPrior()[6:8]), c(0.01, 0, 0))
  dp <- densityPrior(kappa_range = c(1e-6, 1))
  expect_equal(unname(dp[7:8]), c(1e-6, 1))
  expect_output(print(dp), "kappa sampled, log-uniform")
  expect_output(print(densityPrior(kappa = 0.5)), "kappa fixed at 0.5")
  expect_error(densityPrior(kappa = -1), "positive")
  expect_error(densityPrior(kappa_range = c(1, 0.1)), "0 < lo < hi")
  expect_error(densityPrior(kappa_range = 1), "0 < lo < hi")
  expect_error(densityPrior(kappa_range = c(0, 1)), "0 < lo < hi")
  # an old five-entry vector still means a fixed kappa of 0.01
  expect_equal(mdir:::densityHyperparameters(matrix(rnorm(20), 10), 2, 1, c(2, 1, 1, 2, 1))$kappa, 0.01)
})

test_that("the default keeps kappa fixed and the pooled hyperparameters unchanged", {
  set.seed(91)
  X <- matrix(rnorm(60), 30); rownames(X) <- 1:30
  fit <- callMDI(list(X), n_iter = 40, thin = 4, types = "MVN", K_max = 3)
  expect_equal(ncol(fit$pooled_hyperparameters[[1]]), 2)           # the two variance scales only
  fit2 <- callMDI(list(X), n_iter = 40, thin = 4, types = "MVN", K_max = 3, density_prior = densityPrior(kappa_range = c(1e-4, 1)))
  expect_equal(ncol(fit2$pooled_hyperparameters[[1]]), 3)          # plus kappa
  kap <- fit2$pooled_hyperparameters[[1]][, 3]
  expect_true(all(kap >= 1e-4 & kap <= 1))
  expect_gt(length(unique(kap)), 1)
  hp <- mdir:::densityHyperparameters(X, 3, 1, unclass(densityPrior(kappa_range = c(1e-4, 1))))
  expect_equal(hp$kappa_range, c(1e-4, 1))
})

test_that("kappa of the diagonal Gaussian matches its marginal posterior (also with an empty component)", {
  skip_on_cran()
  set.seed(92)
  n <- 14
  X <- cbind(rnorm(n, 1.5, 1), rnorm(n, -2, 0.5))
  lo <- 1e-4; hi <- 1
  dp <- densityPrior(scale_pool_shape = 0, kappa_range = c(lo, hi))
  for (K in c(1, 2)) {
    hp <- mdir:::densityHyperparameters(X, K, 0, unclass(dp))
    nu <- hp$nu; s <- hp$scale
    lk <- seq(log(lo), log(hi), length.out = 4000); kap <- exp(lk)
    ld <- rep(0, length(kap))
    for (p in 1:2) {
      xbar <- mean(X[, p]); SS <- sum((X[, p] - xbar)^2); d <- xbar - hp$xi[p]
      ld <- ld + 0.5 * (log(kap) - log(kap + n)) - 0.5 * (nu + n) * log(s[p] + SS + kap * n / (kap + n) * d^2)
    }
    mom <- grid_moments(lk, ld)
    out <- fit_kappa_component(X, "G", K, dp)
    lg <- log(out$pooled_hyperparameters[[1]][-(1:200), ncol(out$pooled_hyperparameters[[1]])])
    expect_lt(abs(mean(lg) - mom["mean"]), 5 * mom["sd"] / sqrt(ess_of(lg)) + 0.02, label = paste("K", K, "mean"))
    expect_equal(sd(lg), unname(mom["sd"]), tolerance = 0.1, label = paste("K", K, "sd"))
  }
})

test_that("kappa of the MVN matches its marginal posterior", {
  skip_on_cran()
  set.seed(93)
  n <- 16; P <- 2
  X <- matrix(rnorm(n * P), n) %*% chol(matrix(c(1.5, .6, .6, 0.5), 2)) + c(0.8, -0.5)
  lo <- 1e-4; hi <- 1
  dp <- densityPrior(scale_pool_shape = 0, kappa_range = c(lo, hi))
  hp <- mdir:::densityHyperparameters(X, 1, 1, unclass(dp))
  nu <- hp$nu; Psi <- hp$scale
  xbar <- colMeans(X); S <- crossprod(sweep(X, 2, xbar)); d <- xbar - hp$xi
  lk <- seq(log(lo), log(hi), length.out = 4000); kap <- exp(lk)
  ld <- sapply(kap, function(k) 0.5 * P * (log(k) - log(k + n)) -
                 0.5 * (nu + n) * as.numeric(determinant(Psi + S + k * n / (k + n) * tcrossprod(d))$modulus))
  mom <- grid_moments(lk, ld)
  out <- fit_kappa_component(X, "MVN", 1, dp)
  lg <- log(out$pooled_hyperparameters[[1]][-(1:200), ncol(out$pooled_hyperparameters[[1]])])
  expect_lt(abs(mean(lg) - mom["mean"]), 5 * mom["sd"] / sqrt(ess_of(lg)) + 0.02)
  expect_equal(sd(lg), unname(mom["sd"]), tolerance = 0.1)
})

test_that("kappa is carried through tempering and SMC", {
  skip_on_cran()
  set.seed(94)
  X <- cbind(c(rnorm(15, 0, 1), rnorm(15, 5, 1))); rownames(X) <- 1:30
  dp <- densityPrior(kappa_range = c(1e-4, 1))
  fit <- callMDI(list(X), n_iter = 200, thin = 5, types = "G", K_max = 4, density_prior = dp, betas = c(0.3, 0.6, 1))
  kap <- fit$pooled_hyperparameters[[1]][, ncol(fit$pooled_hyperparameters[[1]])]
  expect_true(all(kap >= 1e-4 & kap <= 1)); expect_gt(length(unique(kap)), 1)
  smc <- smcMDI(list(X), types = "G", K_max = 4, n_particles = 50, density_prior = dp)
  expect_equal(length(smc$particle_weights), 50)
})
