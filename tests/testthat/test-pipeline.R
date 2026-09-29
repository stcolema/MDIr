two_view_data <- function(N = 50, seed = 51) {
  set.seed(seed)
  truth <- rep(1:2, each = N / 2)
  X <- lapply(1:2, function(v) {
    m <- matrix(rnorm(N * 2, ifelse(truth == 1, 0, 3)), N, 2); rownames(m) <- seq_len(N); m
  })
  list(X = X, truth = truth)
}

test_that("callMDI output can be processed and pooled across chains", {
  d <- two_view_data()
  chains <- runMCMCChains(d$X, 2, R = 200, thin = 5, types = c("MVN", "G"), K = c(4, 4))
  expect_length(chains, 2)
  proc <- processMCMCChain(chains[[1]], burn = 50)
  expect_equal(length(proc$complete_likelihood), floor(200 / 5) + 1 - (floor(50 / 5) + 1))
  expect_equal(nrow(proc$parameters[[1]]), length(proc$complete_likelihood))
  expect_equal(length(proc$pred), 2)
  pooled <- predictFromMultipleChains(chains, burn = 50)
  expect_equal(length(pooled$pred[[1]]), 50)
  expect_error(processMCMCChain(chains[[1]], burn = 200), "no saved samples")
  # a burn in smaller than the thinning interval must still work
  expect_silent(predictFromMultipleChains(chains, burn = 2))
})

test_that("two clearly separated clusters are recovered in both views", {
  skip_on_cran()
  d <- two_view_data(N = 60)
  chains <- runMCMCChains(d$X, 2, R = 1500, thin = 5, types = c("MVN", "MVN"), K = c(5, 5))
  pooled <- predictFromMultipleChains(chains, burn = 500)
  for (v in 1:2) {
    tab <- table(pooled$pred[[v]], d$truth)
    expect_gte(sum(apply(tab, 2, max)) / sum(tab), 0.95)
  }
})

test_that("semi-supervised views keep observed labels", {
  d <- two_view_data()
  fixed <- matrix(0, 50, 2)
  fixed[c(1:5, 30:35), 1] <- 1
  labels <- matrix(1, 50, 2)
  labels[1:5, 1] <- 1; labels[30:35, 1] <- 2
  out <- callMDI(d$X, R = 150, thin = 5, types = c("MVN", "MVN"), K = c(4, 4),
                 initial_labels = labels, fixed = fixed)
  expect_true(all(out$allocations[, c(1:5), 1] == 0))
  expect_true(all(out$allocations[, c(30:35), 1] == 1))
  # only two components are fixed; label swaps must never move them
  expect_true(all(out$allocations[, 1:5, 1] == 0))
})

test_that("the single view mixture wrapper works with missing data", {
  d <- two_view_data()
  X <- d$X[[1]]; X[3, 1] <- NA
  out <- callMixtureModel(X, R = 100, thin = 5, type = "MVN", K = 4)
  proc <- processMixtureModelChain(out, burn = 20)
  expect_equal(length(proc$pred[[1]]), 50)
})

test_that("results are reproducible with set.seed", {
  d <- two_view_data()
  set.seed(7); a <- callMDI(d$X, R = 60, thin = 6, types = c("MVN", "C_dummy")[c(1, 1)], K = c(3, 3))
  set.seed(7); b <- callMDI(d$X, R = 60, thin = 6, types = c("MVN", "MVN"), K = c(3, 3))
  expect_identical(a$allocations, b$allocations)
  expect_identical(a$phis, b$phis)
})

test_that("writing samples to file round-trips the saved quantities", {
  skip_on_cran()
  d <- two_view_data()
  dir <- tempfile(); dir.create(dir)
  details <- callMDIWritingToFile(d$X, R = 40, thin = 4, types = c("G", "MVN"), K = c(3, 4), dir_path = dir)
  out <- readInSavedSamples(details)
  expect_equal(dim(out$allocations), c(11, 50, 2))
  expect_equal(dim(out$weights), c(11, 4, 2))
  expect_true(all(out$weights[, 4, 1] == 0))
  expect_equal(dim(out$phis), c(11, 1))
  expect_true(all(out$weights[, 1:4, 2] > 0))
  expect_true(all(out$mass > 0))
})

test_that("input validation gives informative errors", {
  d <- two_view_data()
  expect_error(callMDI(d$X, R = 10, thin = 1, types = c("MVN", "MVN"), K = c(3, 3),
                       fixed = matrix(0, 10, 2)), "``fixed`` must have")
  Xc <- d$X; Xc[[1]][1, 1] <- Inf
  expect_error(callMDI(Xc, R = 10, thin = 1, types = c("MVN", "MVN"), K = c(3, 3)), "infinite")
  Xd <- d$X; Xd[[2]][, 1] <- c(0.5, rep(1, 49))
  expect_error(callMDI(Xd, R = 10, thin = 1, types = c("MVN", "C"), K = c(3, 3)), "non-negative integers")
  expect_error(mdiPrior(phi_rate = -1), "positive")
})

test_that("likelihood traces are extracted and plotted", {
  d <- two_view_data()
  chains <- runMCMCChains(d$X, 2, R = 100, thin = 5, types = c("MVN", "MVN"), K = c(3, 3))
  df <- getLikelihood(chains[[1]])
  expect_equal(sort(unique(df$type)), c("complete", "observed"))
  expect_equal(range(df$iteration), c(0, 100))
  expect_s3_class(plotLikelihoods(chains), "ggplot")
  proc <- processMCMCChain(chains[[1]], burn = 50)
  expect_equal(min(getLikelihood(proc)$iteration), 55)
})
