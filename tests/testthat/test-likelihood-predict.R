# The joint (marginal over assignments) likelihood of an item, its class
# probabilities, and prediction of new items. All are checked against
# enumeration of the K_1 x ... x K_L joint assignments.

random_log_g <- function(K) {
  g <- matrix(-Inf, max(K), length(K))
  for (l in seq_along(K)) g[seq_len(K[l]), l] <- rnorm(K[l], -3, 2)
  g
}

test_that("marginal likelihood and class probabilities equal enumeration (ragged K, excluded components)", {
  set.seed(1)
  for (K in list(c(3L, 4L), c(3L, 4L, 2L), c(2L, 3L, 3L, 2L), c(2L, 2L, 3L, 2L, 3L))) {
    for (rep in 1:3) {
      st <- random_mdi_state(K)
      log_g <- random_log_g(K)
      log_g[2, 1] <- -Inf   # a component this item cannot belong to
      ref <- bf_log_marginal(log_g, st$w, K, st$phi)
      expect_equal(mdir:::mdiLogMarginalCpp(log_g, st$w, K, st$phi), ref$log_marginal, tolerance = 1e-10)
      cp <- mdir:::mdiClassProbabilitiesCpp(log_g, st$w, K, st$phi)
      expect_equal(cp, ref$class_prob, tolerance = 1e-10)
      expect_equal(colSums(cp), rep(1, length(K)), tolerance = 1e-10)
    }
  }
})

test_that("the marginal likelihood does not underflow and reduces to a mixture for one view or phi = 0", {
  set.seed(2)
  K <- c(3L, 4L, 2L)
  st <- random_mdi_state(K)
  log_g <- random_log_g(K)
  base <- mdir:::mdiLogMarginalCpp(log_g, st$w, K, st$phi)
  expect_equal(mdir:::mdiLogMarginalCpp(log_g - 3000, st$w, K, st$phi), base - 3000 * length(K), tolerance = 1e-10)

  # phi = 0: independent mixtures with weights w / sum(w) in each view
  st0 <- st; st0$phi[] <- 0
  indep <- sum(vapply(seq_along(K), function(l) {
    k <- seq_len(K[l]); log(sum(st$w[k, l] / sum(st$w[k, l]) * exp(log_g[k, l])))
  }, numeric(1)))
  expect_equal(mdir:::mdiLogMarginalCpp(log_g, st0$w, K, st0$phi), indep, tolerance = 1e-10)

  # one view
  st1 <- random_mdi_state(4L)
  lg1 <- random_log_g(4L)
  expect_equal(mdir:::mdiLogMarginalCpp(lg1, st1$w, 4L, matrix(0, 1, 1)),
               log(sum(st1$w[, 1] / sum(st1$w[, 1]) * exp(lg1[, 1]))), tolerance = 1e-10)

  # a view with nothing possible
  lg_bad <- log_g; lg_bad[, 2] <- -Inf
  expect_equal(mdir:::mdiLogMarginalCpp(lg_bad, st$w, K, st$phi), -Inf)
})

test_that("the recorded joint log-likelihood equals enumeration at the saved parameters (with missing data)", {
  set.seed(3)
  N <- 20; K <- c(3L, 4L, 3L)
  z <- sample(2, N, TRUE)
  X <- lapply(1:3, function(v) { m <- matrix(rnorm(N * 2, 2 * z), N, 2); rownames(m) <- seq_len(N); m })
  X[[2]][3, 1] <- NA      # a missing entry
  X[[3]][5, ] <- NA       # a view with nothing observed for an item
  fit <- callMDI(X, R = 40, thin = 10, types = rep("G", 3), K = K, save_pointwise = TRUE)
  expect_equal(dim(fit$pointwise_likelihood), c(5, N))
  expect_equal(rowSums(fit$pointwise_likelihood), fit$joint_likelihood, tolerance = 1e-12)

  for (s in c(1, 3, 5)) {
    st <- bf_state(fit, s)
    ref <- vapply(seq_len(N), function(n) {
      log_g <- matrix(-Inf, max(K), 3)
      for (l in 1:3) log_g[seq_len(K[l]), l] <- bf_gaussian_loglik(fit, l, s, X[[l]][n, ])
      bf_log_marginal(log_g, st$w, K, st$phi)$log_marginal
    }, numeric(1))
    expect_equal(fit$pointwise_likelihood[s, ], ref, tolerance = 1e-9, info = paste("draw", s))
  }
})

test_that("the joint likelihood differs from the per-view observed likelihood only through phi", {
  set.seed(4)
  N <- 30
  X <- matrix(rnorm(N * 2, rep(c(0, 3), each = N / 2)), N, 2); rownames(X) <- seq_len(N)
  # One view: the two coincide (after the initial state, which the older traces leave at zero)
  m <- callMixtureModel(X, R = 30, thin = 10, type = "G", K = 4)
  expect_equal(m$joint_likelihood[-1], m$observed_likelihood[-1], tolerance = 1e-10)

  # Two views with a strong association: they do not
  Y <- X + matrix(rnorm(N * 2, sd = 0.2), N, 2)
  fit <- callMDI(list(X, Y), R = 200, thin = 10, types = c("G", "G"), K = c(4, 4))
  keep <- 12:21
  expect_gt(mean(fit$phis[keep, 1]), 1)
  expect_gt(max(abs(fit$joint_likelihood[keep] - fit$observed_likelihood[keep])), 1)
})

test_that("pointwiseLogLik stacks chains, applies the burn in and says when nothing was saved", {
  set.seed(5)
  N <- 20
  X <- lapply(1:2, function(v) { m <- matrix(rnorm(N * 2, rep(c(0, 3), each = N / 2)), N, 2); rownames(m) <- seq_len(N); m })
  chains <- runMCMCChains(X, 2, R = 50, thin = 5, types = c("G", "G"), K = c(3, 3), save_pointwise = TRUE)
  ll <- pointwiseLogLik(chains, burn = 20)
  # saved draws 0, 5, ..., 50; the initial state and the first floor(20 / 5) are dropped
  expect_equal(dim(ll), c(2 * 6, N))
  expect_equal(attr(ll, "chain_id"), rep(1:2, each = 6))
  expect_equal(colnames(ll), as.character(seq_len(N)))
  expect_equal(unname(rowSums(ll[1:6, ])), chains[[1]]$joint_likelihood[6:11], tolerance = 1e-12)

  one <- pointwiseLogLik(chains[[1]])
  expect_equal(nrow(one), 10)

  no_pointwise <- runMCMCChains(X, 1, R = 20, thin = 5, types = c("G", "G"), K = c(3, 3))
  expect_error(pointwiseLogLik(no_pointwise), "save_pointwise")
  expect_error(pointwiseLogLik(chains, burn = 1000), "No saved iterations")
})

test_that("predicting the fitted items reproduces their pointwise likelihoods (all densities, outliers, missing data)", {
  skip_on_cran()
  set.seed(6)
  N <- 24
  z <- sample(2, N, TRUE)
  X <- list(
    { m <- matrix(rnorm(N * 2, 2 * z), N, 2); m[4, 2] <- NA; m },
    { m <- matrix(rnorm(N * 3, 2 * z), N, 3); m[7, ] <- NA; m },
    { m <- matrix(rnorm(N * 2, 2 * z), N, 2); m },
    matrix(rbinom(N * 3, 2, c(0.2, 0.8)[z]), N, 3)
  )
  for (v in seq_along(X)) rownames(X[[v]]) <- seq_len(N)
  types <- c("G", "MVN", "TAGM", "C")
  fit <- callMDI(X, R = 60, thin = 20, types = types, K = c(3, 3, 3, 3), save_pointwise = TRUE)
  # the last saved draw only
  pred <- predictMDI(fit, X, X, burn = 40)
  expect_equal(dim(pred$log_likelihood), c(1, N))
  expect_equal(unname(pred$log_likelihood[1, ]), unname(fit$pointwise_likelihood[4, ]), tolerance = 1e-9)
  expect_equal(unname(pred$log_predictive), unname(pred$log_likelihood[1, ]))
})

test_that("predictMDI equals enumeration for new items: class probabilities, co-clustering, outlying new data", {
  set.seed(7)
  N <- 24; K <- c(3L, 4L)
  z <- sample(2, N, TRUE)
  X <- lapply(1:2, function(v) { m <- matrix(rnorm(N * 2, 3 * z), N, 2); rownames(m) <- seq_len(N); m })
  # View 1 semi-supervised: the first 12 items have observed classes
  fixed <- matrix(0, N, 2); fixed[1:12, 1] <- 1
  labels <- matrix(1, N, 2); labels[1:12, 1] <- z[1:12]
  fit <- callMDI(X, R = 40, thin = 10, types = c("G", "G"), K = K,
                 initial_labels = labels, fixed = fixed)
  new <- list(matrix(c(0.1, 3.2, 3, 2.9, NA, 6), 3, 2, byrow = TRUE),
              matrix(c(0.3, 0.1, 5.5, 3.1, 3.2, 3.0), 3, 2, byrow = TRUE))
  pred <- predictMDI(fit, X, new, burn = 30, coclustering = TRUE)   # draw 5 only
  s <- 5
  st <- bf_state(fit, s)
  ref <- lapply(1:3, function(j) {
    log_g <- matrix(-Inf, max(K), 2)
    for (l in 1:2) log_g[seq_len(K[l]), l] <- bf_gaussian_loglik(fit, l, s, new[[l]][j, ])
    bf_log_marginal(log_g, st$w, K, st$phi)
  })
  expect_equal(unname(pred$log_likelihood[1, ]), vapply(ref, `[[`, numeric(1), "log_marginal"), tolerance = 1e-9)

  # semi-supervised view 1 returns class probabilities, view 2 does not
  expect_null(pred$class_probability[[2]])
  ref_cp <- t(vapply(ref, function(r) r$class_prob[seq_len(K[1]), 1], numeric(K[1])))
  expect_equal(unname(pred$class_probability[[1]]), unname(ref_cp), tolerance = 1e-9)

  # co-clustering with fitted items: sum_k p(c_new = k) 1[c_i = k]
  for (l in 1:2) {
    ref_cc <- t(vapply(ref, function(r) r$class_prob[fit$allocations[s, , l] + 1, l], numeric(N)))
    expect_equal(unname(pred$coclustering[[l]]), unname(ref_cc), tolerance = 1e-9)
  }
  expect_equal(dim(pred$coclustering[[1]]), c(3, N))
  expect_equal(colnames(pred$coclustering[[1]]), as.character(seq_len(N)))
})

test_that("predictMDI pools chains, thins to n_draws and checks its input", {
  set.seed(8)
  N <- 20
  X <- lapply(1:2, function(v) { m <- matrix(rnorm(N * 2, rep(c(0, 3), each = N / 2)), N, 2); rownames(m) <- seq_len(N); m })
  chains <- runMCMCChains(X, 2, R = 50, thin = 5, types = c("G", "G"), K = c(3, 3))
  new <- lapply(X, function(m) m[1:4, ])
  all_draws <- predictMDI(chains, X, new, burn = 20)
  expect_equal(dim(all_draws$log_likelihood), c(12, 4))
  some <- predictMDI(chains, X, new, burn = 20, n_draws = 5)
  expect_equal(nrow(some$log_likelihood), 5)
  expect_true(all(is.finite(some$log_predictive)))

  expect_error(predictMDI(chains, X, list(new[[1]][, 1, drop = FALSE], new[[2]])), "columns")
  expect_error(predictMDI(chains, X, new[1]), "one matrix for each")
  expect_error(predictMDI(chains, X, list(new[[1]], new[[2]][1:2, ])), "same number of items")
  expect_error(predictMDI(chains, X, new, burn = 1000), "No saved iterations")
  expect_error(predictMDI(chains, X, new, n_draws = 0), "positive")

  processed <- processMCMCChains(chains, burn = 20)
  expect_error(predictMDI(processed, X, new), "processMCMCChain")

  no_pars <- runMCMCChains(X, 1, R = 20, thin = 5, types = c("G", "G"), K = c(3, 3), save_parameters = FALSE)
  expect_error(predictMDI(no_pars, X, new), "save_parameters")
})

test_that("predictMDI rejects categories the model has not seen", {
  set.seed(9)
  N <- 20
  X <- list(matrix(rbinom(N * 2, 1, 0.5), N, 2), matrix(rnorm(N * 2), N, 2))
  for (v in 1:2) rownames(X[[v]]) <- seq_len(N)
  fit <- callMDI(X, R = 20, thin = 10, types = c("C", "G"), K = c(2, 2))
  new <- list(matrix(c(0, 2), 1, 2), matrix(c(0, 0), 1, 2))
  expect_error(predictMDI(fit, X, new, burn = 10), "category")
})
