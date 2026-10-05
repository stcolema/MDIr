# Output post-processing and input helpers.

test_that("stick-breaking draws have the Beta(1, alpha) first stick", {
  set.seed(2)
  m <- mean(replicate(4000, stickBreakingPrior(10, 5)[1]))
  expect_equal(m, 1 / 11, tolerance = 0.1)
  m2 <- mean(replicate(4000, sampleStickBreakingPrior(10, 5)[1]))
  expect_equal(m2, 1 / 11, tolerance = 0.1)
})

test_that("consensus clustering keeps every chain's weights and accepts a short depth", {
  set.seed(1)
  X <- lapply(1:2, function(v) {
    m <- matrix(rnorm(40 * 2, rep(c(0, 3), each = 20)), 40, 2)
    rownames(m) <- 1:40
    m
  })
  ch <- runMCMCChains(X, 3, R = 100, thin = 5, types = c("G", "G"), K = c(4, 4))
  cc <- compileConsensusClustering(ch, D = 100, W = 3)
  expect_equal(dim(cc$weights[[1]]), c(3, 4))
  for (w in 1:3) expect_equal(unname(cc$weights[[1]][w, ]), unname(ch[[w]]$weights[21, 1:4, 1]))
  expect_equal(cc$normalising_constant[, 1], vapply(ch, function(x) x$evidence[21], numeric(1)))
  expect_silent(compileConsensusClustering(ch, D = 3, W = 3))
})

test_that("calcAllocProb drops the initial state and floor(burn / thin) samples", {
  set.seed(1)
  X <- lapply(1:2, function(v) {
    m <- matrix(rnorm(40 * 2, rep(c(0, 3), each = 20)), 40, 2)
    rownames(m) <- 1:40
    m
  })
  lab <- matrix(1, 40, 2)
  fixed <- matrix(0, 40, 2)
  fixed[c(1:5, 21:25), 1] <- 1
  lab[21:25, 1] <- 2
  f <- callMDI(X, R = 50, thin = 5, types = c("G", "G"), K = c(4, 4), initial_labels = lab,
               fixed = fixed, check_prior = FALSE)
  a <- f$allocation_probabilities[[1]]
  expect_equal(calcAllocProb(f, 1, burn = 0, method = "mean"), rowSums(a, dims = 2) / 11)
  expect_equal(calcAllocProb(f, 1, burn = 2, method = "mean"), rowSums(a[, , -1], dims = 2) / 10)
  expect_equal(calcAllocProb(f, 1, burn = 10, method = "mean"), rowSums(a[, , -(1:3)], dims = 2) / 8)
  expect_error(calcAllocProb(f, 1, burn = 50, method = "mean"), "no saved")
})

test_that("a single view round-trips through the saved-sample files", {
  set.seed(3)
  X <- list(matrix(rnorm(40), 20, 2))
  rownames(X[[1]]) <- 1:20
  d <- tempfile()
  dir.create(d)
  rd <- callMDIWritingToFile(X, R = 20, thin = 5, types = "G", K = 4, dir_path = d)
  out <- readInSavedSamples(rd)
  expect_equal(dim(out$allocations), c(5, 20, 1))
  expect_true(all(out$allocations %in% 0:3))
  expect_true(all(is.finite(out$complete_likelihood)))
  expect_true(all(out$mass > 0))
  expect_true(all(out$complete_likelihood[-1] < 0))
  expect_lt(abs(out$observed_likelihood[1]) , Inf)
  expect_true(all(out$observed_likelihood != 0))
})

test_that("likelihood traces are labelled with the iterations at which they were saved", {
  set.seed(3)
  X <- list(matrix(rnorm(40), 20, 2))
  rownames(X[[1]]) <- 1:20
  fit <- callMDI(X, R = 103, thin = 5, types = "G", K = 4, check_prior = FALSE)
  gl <- getLikelihood(fit)
  expect_equal(range(gl$iteration), c(0, 100))
  expect_equal(sum(gl$type == "complete"), 21)
  # the initial state has a recorded likelihood
  expect_true(all(gl$log_likelihood != 0))
})

test_that("pooled semi-supervised predictions are the most probable class", {
  set.seed(1)
  X <- lapply(1:2, function(v) {
    m <- matrix(rnorm(40 * 2, rep(c(0, 3), each = 20)), 40, 2)
    rownames(m) <- 1:40
    m
  })
  lab <- matrix(1, 40, 2)
  fixed <- matrix(0, 40, 2)
  fixed[c(1:5, 21:25), 1] <- 1
  lab[21:25, 1] <- 2
  ch <- runMCMCChains(X, 2, R = 60, thin = 5, types = c("G", "G"), K = c(4, 4), initial_labels = lab,
                      fixed = fixed)
  pooled <- predictFromMultipleChains(ch, burn = 20)
  expect_equal(unname(pooled$pred[[1]]), unname(apply(pooled$allocation_probability[[1]], 1, which.max)))
})

test_that("small helpers cope with one depth, one width and vector windows", {
  cms <- list(matrix(0, 3, 3), matrix(0.1, 3, 3))
  desc <- data.frame(Depth = c(10, 10), Width = c(1, 2))
  out <- makeCMComparisonSummaryDF(cms, desc)
  expect_equal(nrow(out), 1)
  expect_equal(out$Quantity_varied, "Width")
  expect_warning(w <- processProposalWindows(list(c(0.1, 0.2, 0.3)), "G"), "non-GP")
  expect_equal(w[[1]], 0)
})

test_that("the deprecated names forward to the current functions with a warning", {
  set.seed(1)
  X <- lapply(1:2, function(v) {
    m <- matrix(rnorm(40 * 2, rep(c(0, 3), each = 20)), 40, 2)
    rownames(m) <- 1:40
    m
  })
  fit <- callMDI(X, R = 20, thin = 5, types = c("G", "G"), K = c(3, 3), check_prior = FALSE)
  expect_warning(old <- calcFusionProbabiliy(fit, 1, 2), "deprecated")
  expect_equal(old, calcFusionProbability(fit, 1, 2))
  expect_warning(old_all <- calcFusionProbabiliyAllViews(fit), "deprecated")
  expect_equal(old_all, calcFusionProbabilityAllViews(fit))
  # `evidence` is the former name of `normalising_constant`
  expect_identical(fit$evidence, fit$normalising_constant)
  expect_equal(length(fit$normalising_constant), 5)
  expect_identical(processMCMCChain(fit, burn = 5)$evidence, processMCMCChain(fit, burn = 5)$normalising_constant)
})
