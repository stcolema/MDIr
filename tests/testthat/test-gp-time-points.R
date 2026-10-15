# Gaussian process views whose columns are not one unit apart: `time_points`.

gp_data <- function(N = 30, P = 6, seed = 1) {
  set.seed(seed)
  X <- matrix(rnorm(N * P, rep(c(0, 2), each = N / 2)), N, P)
  rownames(X) <- seq_len(N)
  X
}

run_gp <- function(X, ..., seed = 11) {
  set.seed(seed)
  callMDI(list(X), n_iter = 400, thin = 4, types = "GP", K_max = 3, check_prior = FALSE,
          control = mdirControl(n_burn = 200), ...)
}

test_that("the default positions are 0, 1, ..., P - 1, and a shift of them changes nothing", {
  X <- gp_data()
  base <- run_gp(X)
  explicit <- run_gp(X, time_points = list(0:5))
  shifted <- run_gp(X, time_points = list(37 + 0:5))
  expect_identical(explicit$hypers, base$hypers)
  expect_identical(explicit$allocations, base$allocations)
  expect_identical(shifted$hypers, base$hypers)
  expect_identical(base$time_points, list(numeric(0)))
  expect_equal(explicit$time_points, list(as.numeric(0:5)))
})

test_that("rescaling the positions rescales the length scales and nothing else", {
  # The model is invariant under (positions, length) -> (c positions, c length), the floor and
  # the prior's upper quantile scale with c, and the random walks are on the log scale.
  X <- gp_data()
  base <- run_gp(X)
  scaled <- run_gp(X, time_points = list(2.5 * (0:5)))
  expect_equal(scaled$hypers[[1]]$length, 2.5 * base$hypers[[1]]$length, tolerance = 1e-6)
  expect_equal(scaled$hypers[[1]]$amplitude, base$hypers[[1]]$amplitude, tolerance = 1e-6)
  expect_equal(scaled$hypers[[1]]$noise, base$hypers[[1]]$noise, tolerance = 1e-6)
  expect_identical(scaled$allocations, base$allocations)
})

test_that("the length-scale floor and prior follow the positions", {
  X <- gp_data(P = 5)
  tp <- c(0, 0.5, 2, 3.5, 7)
  for (mult in c(1, 2)) {
    dp <- as.numeric(densityPrior(gp_min_length = mult))
    hp <- mdir:::densityHyperparameters(X, 1, 3, dp, tp)
    expect_equal(hp$min_length, mult * 0.5)
    extent <- 7
    ab <- c(hp$length_prior_shape, hp$length_prior_rate)
    expect_equal(pgamma(1 / hp$min_length, ab[1], rate = ab[2], lower.tail = FALSE), 0.01, tolerance = 1e-6)
    expect_equal(pgamma(1 / extent, ab[1], rate = ab[2]), 0.01, tolerance = 1e-6)
    expect_equal(hp$max_length, 10 * extent)
  }
  # the default positions give the old floor and extent
  hp0 <- mdir:::densityHyperparameters(X, 1, 3, numeric(0))
  expect_equal(hp0$min_length, 1)
  expect_equal(hp0$max_length, 10 * 4)
  # positions need not be sorted: the floor uses the sorted spacings
  hp1 <- mdir:::densityHyperparameters(X, 1, 3, numeric(0), rev(tp))
  expect_equal(hp1$min_length, 0.5)
})

test_that("the length scale never falls below the floor set by the closest pair", {
  skip_on_cran()
  X <- gp_data(N = 40, P = 6, seed = 3)
  tp <- c(0, 1, 1.5, 4, 8, 9)
  fit <- run_gp(X, time_points = list(tp))
  expect_gte(min(fit$hypers[[1]]$length), 0.5)
})

# The hyperparameters of a component whose items all have observed labels have a posterior that
# is computable on a grid: with the mean function integrated out, the cluster's sample mean is
# N(xi, K + noise / n) and its residual sum of squares is independent of the mean function.
# (xi is the column mean of all the data, so with a single cluster the sample mean carries no
# information on the kernel; two labelled clusters do.) The population of amplitudes and noises
# is fixed here (pool_gp = FALSE), so the clusters are independent.
gp_log_marginal_offset <- function(Xk, d, a, lambda, s2, positions, jitter = 1e-8) {
  n <- nrow(Xk); P <- ncol(Xk)
  SS <- sum(sweep(Xk, 2, colMeans(Xk))^2)
  D <- outer(positions, positions, function(i, j) -0.5 * (i - j)^2)
  ch <- chol(a * exp(D / lambda^2) + diag(jitter * a, P) + diag(s2 / n, P))
  z <- backsolve(ch, d, transpose = TRUE)
  -0.5 * (n - 1) * P * log(2 * pi * s2) - 0.5 * P * log(n) - SS / (2 * s2) -
    0.5 * (P * log(2 * pi) + 2 * sum(log(diag(ch))) + sum(z^2))
}

test_that("the hyperparameter sampler targets the exact marginal posterior at unequal positions", {
  skip_on_cran()
  set.seed(95)
  tp <- c(0, 0.5, 2, 3.5, 7, 8, 9.5, 12)
  P <- length(tp); n_per <- 20
  Kt <- 1 * exp(-0.5 * outer(tp, tp, "-")^2 / 3^2)
  f <- t(chol(Kt + diag(1e-8, P))) %*% matrix(rnorm(2 * P), P)
  X <- rbind(matrix(rnorm(n_per * P, rep(f[, 1], each = n_per), 0.3), n_per),
             matrix(rnorm(n_per * P, rep(f[, 2], each = n_per), 0.3), n_per))
  rownames(X) <- seq_len(2 * n_per)
  group <- rep(1:2, each = n_per)

  dp <- densityPrior(pool_gp = FALSE)
  hp <- mdir:::densityHyperparameters(X, 2, 3, as.numeric(dp), tp)
  c0 <- hp$log_variance_centre
  alpha <- hp$length_prior_shape; beta <- hp$length_prior_rate
  lo <- hp$min_length; hi <- hp$max_length
  Xk <- X[group == 1, ]
  d <- colMeans(Xk) - colMeans(X)
  s2hat <- sum(sweep(Xk, 2, colMeans(Xk))^2) / ((n_per - 1) * P)
  ya <- c0 + seq(-8, 4, length.out = 40)
  yn <- log(s2hat) + seq(-1.5, 1.5, length.out = 25)
  z <- seq(log(lo), log(hi), length.out = 30)
  eg <- expand.grid(ia = seq_along(ya), in_ = seq_along(yn), iz = seq_along(z))
  grid_posterior <- function(positions) {
    lp <- mapply(function(ia, in_, iz) {
      gp_log_marginal_offset(Xk, d, exp(ya[ia]), exp(z[iz]), exp(yn[in_]), positions) +
        dnorm(ya[ia], c0, 1, log = TRUE) + dnorm(yn[in_], c0, 1, log = TRUE) +
        (-alpha * z[iz] - beta / exp(z[iz]))
    }, eg$ia, eg$in_, eg$iz)
    list(amplitude = grid_moments(ya[eg$ia], lp), noise = grid_moments(yn[eg$in_], lp),
         length = grid_moments(z[eg$iz], lp))
  }
  gm <- grid_posterior(tp)

  # the comparison can tell the positions from equal spacing: the length scale's posterior
  # mean moves by more than a standard deviation if the positions are replaced by 0, 1, ...
  gm_index <- grid_posterior(seq_len(P) - 1)
  expect_gt(abs(gm$length["mean"] - gm_index$length["mean"]) / gm$length["sd"], 1)

  out <- callMDI(list(X), n_iter = 400000, thin = 10, types = "GP", K_max = 2, time_points = list(tp),
                 initial_labels = matrix(group, ncol = 1), fixed = matrix(1, nrow(X), 1),
                 density_prior = dp, check_prior = FALSE, control = mdirControl(auto_tune = FALSE))
  h_out <- out$hypers[[1]]
  keep <- -(1:200)
  draws <- list(amplitude = log(h_out$amplitude[keep, 1]), noise = log(h_out$noise[keep, 1]),
                length = log(h_out$length[keep, 1]))
  for (nm in names(draws)) {
    dr <- draws[[nm]]
    expect_lt(abs(mean(dr) - gm[[nm]]["mean"]), 5 * gm[[nm]]["sd"] / sqrt(ess_of(dr)) + 0.05 * gm[[nm]]["sd"],
              label = paste(nm, "mean"))
    expect_equal(sd(dr), unname(gm[[nm]]["sd"]), tolerance = 0.15, label = paste(nm, "sd"))
  }
})

test_that("the positions are checked", {
  X <- gp_data()
  Xm <- cbind(X[, 1:2], X[, 1:2]); rownames(Xm) <- rownames(X)
  fit <- function(tp, ...) callMDI(list(X), n_iter = 20, thin = 5, types = "GP", K_max = 2,
                                   check_prior = FALSE, time_points = tp, ...)
  expect_error(fit(list(1:5)), "one entry for each of the 6 columns")
  expect_error(fit(list(c(0, 1, 2, 3, 4, NA))), "finite")
  expect_error(fit(list(c(0, 1, 2, 3, 4, 4))), "must not repeat")
  expect_error(fit(list("a")), "finite numeric")
  expect_error(fit(list(1:6, 1:6)), "entry for each view")
  expect_error(
    callMDI(list(X, Xm), n_iter = 20, thin = 5, types = c("GP", "MVN"), K_max = c(2, 2),
            check_prior = FALSE, time_points = list(NULL, 1:4)),
    "applies to \"GP\" and \"TAGPM\" views only"
  )
  # a plain vector is accepted for one view, and equals the list
  set.seed(1); a <- fit(c(0, 1, 2, 3, 4, 5.5))
  set.seed(1); b <- fit(list(c(0, 1, 2, 3, 4, 5.5)))
  expect_identical(a$hypers, b$hypers)
  # and for a single mixture
  expect_no_error(callMixtureModel(X, n_iter = 20, thin = 5, type = "GP", K_max = 2, time_points = c(0, 1, 2, 3, 4, 5.5)))
  expect_error(callMixtureModel(X, n_iter = 20, thin = 5, type = "GP", K_max = 2, time_points = 1:3), "6 columns")
})

test_that("every entry point takes the positions, and the predictive functions use them", {
  skip_on_cran()
  X <- gp_data(N = 30, P = 6)
  tp <- list(c(0, 1, 2.5, 4, 8, 9))
  chains <- runMCMCChains(list(X), 2, n_iter = 120, thin = 4, types = "GP", K_max = 2, time_points = tp)
  expect_equal(chains[[2]]$time_points, list(tp[[1]]))
  fit <- fitMDI(list(X), 2, n_iter = 200, thin = 4, types = "GP", K_max = 2, time_points = tp, verbose = FALSE)
  expect_equal(fit[[1]]$time_points, list(tp[[1]]))
  sims <- simulatePosteriorPredictive(fit, list(X), n_draws = 5, burn = 50)
  expect_length(sims$replicates, 5)
  new <- predictMDI(fit, list(X), list(X[1:2, , drop = FALSE]), burn = 50)
  expect_true(all(is.finite(new$log_predictive)))
  prior_sims <- simulatePriorPredictive(list(X), "GP", K_max = 2, n_datasets = 3, time_points = tp)
  expect_length(prior_sims$replicates, 3)

  dir <- tempfile("gp_time_points")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  expect_no_error(callMDIWritingToFile(list(X), n_iter = 40, thin = 4, types = "GP", K_max = 2,
                                       dir_path = dir, time_points = tp))
  chain_dir <- file.path(dir, "chains")
  dir.create(chain_dir)
  expect_no_error(runMCMCChainsSavingToFile(list(X), 2, n_iter = 40, thin = 4, types = "GP", K_max = 2,
                                            dir_path = chain_dir, time_points = tp))
})

# Two curves observed at unequally spaced positions, partly missing: the mean curves are
# recovered at every position, including those where many items are unobserved.
test_that("curves are recovered at unequal positions, with missing measurements", {
  skip_on_cran()
  set.seed(31)
  tp <- c(0, 0.4, 1.5, 2, 3.7, 5, 5.3, 8, 9.5, 12)
  P <- length(tp)
  N <- 60
  truth <- rep(1:2, each = N / 2)
  curves <- rbind(2 * sin(tp / 2.5), -1.5 * cos(tp / 4) + 1)
  X <- curves[truth, ] + matrix(rnorm(N * P, sd = 0.4), N, P)
  rownames(X) <- seq_len(N)
  X[sample(length(X), 0.3 * length(X))] <- NA
  fit <- callMDI(list(X), n_iter = 6000, thin = 5, types = "GP", K_max = 3, time_points = list(tp),
                 control = mdirControl(n_burn = 3000))
  keep <- seq(floor(dim(fit$allocations)[1] / 2) + 1, dim(fit$allocations)[1])
  pred <- apply(fit$allocations[keep, , 1], 2, function(z) as.integer(names(which.max(table(z)))))
  tab <- table(pred, truth)
  expect_gte(sum(apply(tab, 2, max)) / N, 0.95)

  # posterior mean of each cluster's curve, from the draws in which it is the main cluster
  par <- fit$parameters[[1]][keep, , drop = FALSE]
  K <- 3
  f_draws <- lapply(seq_len(K), function(k) par[, (k - 1) * P + seq_len(P), drop = FALSE])
  size <- sapply(seq_len(K), function(k) mean(fit$N_k[k, 1, keep]))
  main <- order(size, decreasing = TRUE)[1:2]
  est <- t(sapply(main, function(k) colMeans(f_draws[[k]])))
  err <- sapply(1:2, function(g) {
    cand <- apply(est, 1, function(e) sqrt(mean((e - curves[g, ])^2)))
    min(cand)
  })
  expect_lt(max(err), 0.15)
})
