# Only semi-supervised views record allocation probabilities (N x K x draws)

test_that("allocation probabilities are recorded for semi-supervised views only", {
  set.seed(20)
  N <- 40; K <- c(5L, 6L)
  z <- sample(2, N, TRUE)
  X <- lapply(1:2, function(v) { m <- matrix(rnorm(N * 2, 3 * z), N, 2); rownames(m) <- seq_len(N); m })

  unsup <- callMDI(X, R = 30, thin = 10, types = c("G", "G"), K = K)
  expect_null(unsup$allocation_probabilities[[1]])
  expect_null(unsup$allocation_probabilities[[2]])
  expect_length(unsup$allocation_probabilities, 2)

  fixed <- matrix(0, N, 2); fixed[1:15, 1] <- 1
  labels <- matrix(1, N, 2); labels[1:15, 1] <- z[1:15]
  semi <- callMDI(X, R = 30, thin = 10, types = c("G", "G"), K = K, initial_labels = labels, fixed = fixed)
  expect_equal(dim(semi$allocation_probabilities[[1]]), c(N, K[1], 4))
  expect_null(semi$allocation_probabilities[[2]])
  # the recorded probabilities are probabilities
  expect_equal(unname(rowSums(semi$allocation_probabilities[[1]][, , 4])), rep(1, N), tolerance = 1e-10)

  # downstream: point estimates for both kinds of view, and a clear error otherwise
  proc <- processMCMCChain(semi, burn = 10)
  expect_equal(dim(proc$allocation_probability[[1]]), c(N, K[1]))
  expect_null(proc$allocation_probability[[2]])
  expect_length(proc$pred[[2]], N)
  expect_error(calcAllocProb(unsup, 1), "semi-supervised")
  expect_equal(dim(calcAllocProb(semi, 1)), c(N, K[1]))

  # nothing is stored in their place
  expect_lt(as.numeric(object.size(unsup$allocation_probabilities)), 1000)
})

test_that("a single-view mixture model has no allocation probabilities unless semi-supervised", {
  set.seed(21)
  N <- 30
  X <- matrix(rnorm(N * 2, rep(c(0, 3), each = N / 2)), N, 2); rownames(X) <- seq_len(N)
  m <- callMixtureModel(X, R = 20, thin = 10, type = "G", K = 4)
  expect_null(m$allocation_probabilities)
  fixed <- rep(0, N); fixed[1:10] <- 1
  labs <- rep(1, N); labs[1:10] <- rep(1:2, 5)
  s <- callMixtureModel(X, R = 20, thin = 10, type = "G", K = 4, initial_labels = labs, fixed = fixed)
  expect_equal(dim(s$allocation_probabilities), c(N, 4, 3))
})

test_that("joint likelihood and burn in: processed chains drop the burn in", {
  set.seed(22)
  N <- 20
  X <- lapply(1:2, function(v) { m <- matrix(rnorm(N * 2, rep(c(0, 3), each = N / 2)), N, 2); rownames(m) <- seq_len(N); m })
  fit <- callMDI(X, R = 50, thin = 5, types = c("G", "G"), K = c(3, 3), save_pointwise = TRUE)
  proc <- processMCMCChain(fit, burn = 20)
  expect_length(proc$joint_likelihood, 6)
  expect_equal(nrow(proc$pointwise_likelihood), 6)
  expect_equal(proc$joint_likelihood, fit$joint_likelihood[6:11])
  # a processed chain is accepted by pointwiseLogLik (its burn in is already removed)
  expect_equal(nrow(pointwiseLogLik(proc)), 6)
  expect_output(print(summary(fit)), "joint = ")
})
