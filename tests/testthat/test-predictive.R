make_X <- function(N = 60, P = 2, seed = 41) {
  set.seed(seed)
  X <- matrix(rnorm(N * P, rep(c(0, 3), each = N / 2)), N, P)
  rownames(X) <- seq_len(N)
  X
}

test_that("prior predictive draws follow the prior (MVN, one view)", {
  X <- make_X()
  sims <- simulatePriorPredictive(X, "MVN", K = 3, n_datasets = 200)
  expect_s3_class(sims, "mdir_predictive")
  expect_length(sims$replicates, 200)
  expect_equal(dim(sims$replicates[[1]][[1]]), dim(X))
  # Data are drawn around the data-driven prior mean xi with large spread
  hp <- mdir:::densityHyperparameters(X, 3, 1, numeric(0))
  first_col_means <- vapply(sims$replicates, function(r) mean(r[[1]][, 1]), numeric(1))
  expect_lt(abs(mean(first_col_means) - hp$xi[1]), 4 * sd(first_col_means) / sqrt(200) + 0.5)
  expect_true(all(vapply(sims$replicates, function(r) all(is.finite(r[[1]])), logical(1))))
})

test_that("prior predictive agreement between views follows phi", {
  X <- list(make_X(), make_X(seed = 42))
  X[[2]] <- X[[2]][, 1, drop = FALSE]
  weak <- simulatePriorPredictive(X, c("MVN", "G"), K = c(3, 3), n_datasets = 300,
                                  prior = mdiPrior(phi_shape = 2, phi_rate = 200))
  strong <- simulatePriorPredictive(X, c("MVN", "G"), K = c(3, 3), n_datasets = 300,
                                    prior = mdiPrior(phi_shape = 20, phi_rate = 0.2))
  agree <- function(s) mean(vapply(s$parameters, function(p) mean(p$labels[, 1] == p$labels[, 2]), numeric(1)))
  expect_gt(agree(strong), agree(weak) + 0.2)
})

test_that("prior predictive works for every type and mimics missingness", {
  set.seed(43)
  N <- 30
  base <- matrix(rnorm(N * 3), N); rownames(base) <- seq_len(N)
  cat_x <- (base > 0) * 1; storage.mode(cat_x) <- "double"
  for (type in c("G", "MVN", "GP", "TAGM", "TAGPM", "C")) {
    X <- if (type == "C") cat_x else base
    X[2, 2] <- NA
    sims <- simulatePriorPredictive(X, type, K = 3, n_datasets = 3)
    expect_true(is.na(sims$replicates[[1]][[1]][2, 2]), info = type)
    expect_false(anyNA(sims$replicates[[1]][[1]][-2, ]), info = type)
    if (type == "C") expect_true(all(sims$replicates[[1]][[1]] %in% c(0, 1, NA)))
  }
})

test_that("posterior predictive checks detect a mis-specified variance", {
  skip_on_cran()
  X <- make_X(N = 80)
  # Data with heavy tails that a two-component Gaussian mixture cannot match at K = 1
  set.seed(44)
  Xt <- matrix(rt(80 * 2, 2) * 2, 80, 2); rownames(Xt) <- seq_len(80)
  fit <- callMDI(list(Xt), R = 800, thin = 4, types = "G", K = 1)
  post <- simulatePosteriorPredictive(fit, Xt, n_draws = 100, burn = 400)
  chk <- predictiveCheck(Xt, post, function(m) max(abs(m[, 1])))
  expect_gt(chk$observed, quantile(chk$replicated, 0.95))
  expect_lt(chk$p_value, 0.05)

  # A well specified model is not flagged
  X2 <- make_X(N = 100)
  fit2 <- callMDI(list(X2), R = 800, thin = 4, types = "MVN", K = 4)
  post2 <- simulatePosteriorPredictive(fit2, X2, n_draws = 100, burn = 400)
  chk2 <- predictiveCheck(X2, post2, function(m) sd(m[, 1]))
  expect_gt(chk2$p_value, 0.02)
  expect_lt(chk2$p_value, 0.98)
})

test_that("posterior predictive draws use saved parameters and respect missingness", {
  skip_on_cran()
  X <- make_X(N = 50); X[sample(length(X), 10)] <- NA
  fit <- suppressWarnings(callMDI(list(X), R = 200, thin = 10, types = "MVN", K = 3))
  post <- simulatePosteriorPredictive(fit, X, n_draws = 5, burn = 100)
  expect_length(post$replicates, 5)
  expect_equal(is.na(post$replicates[[1]][[1]]), is.na(X))
  fit_nopar <- suppressWarnings(callMDI(list(X), R = 100, thin = 10, types = "MVN", K = 3, save_parameters = FALSE))
  expect_error(simulatePosteriorPredictive(fit_nopar, X), "save_parameters")
})

test_that("plotPredictiveCheck returns ggplot objects", {
  X <- make_X()
  sims <- simulatePriorPredictive(X, "MVN", K = 3, n_datasets = 10)
  expect_s3_class(plotPredictiveCheck(X, sims, "density"), "ggplot")
  expect_s3_class(plotPredictiveCheck(X, sims, "statistic", statistic = sd, statistic_name = "SD"), "ggplot")
})
