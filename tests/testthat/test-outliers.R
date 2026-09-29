test_that("t-augmented model flags gross outliers and estimates their weight", {
  # K equals the number of clusters: with spare components a broad component can
  # absorb the outliers instead of the outlier distribution
  skip_on_cran()
  set.seed(71)
  N <- 120
  X <- rbind(
    matrix(rnorm(50 * 2, 0, 0.7), 50), matrix(rnorm(50 * 2, 5, 0.7), 50),
    matrix(runif(20 * 2, -15, 20), 20)
  )
  is_outlier <- rep(c(FALSE, TRUE), c(100, 20))
  rownames(X) <- seq_len(N)
  # Start at the true clusters: from random initial labels the sampler can settle
  # in a local mode where one component holds both clusters and the other absorbs
  # the outliers, which is what running several chains and checking Rhat is for.
  labels <- matrix(rep(1:2, c(50, 70)), N, 1)
  fit <- callMDI(list(X), R = 3000, thin = 5, types = "TAGM", K = 2,
                 initial_labels = labels, initial_labels_as_intended = TRUE)
  keep <- -(1:200)
  p_out <- colMeans(fit$outliers[keep, , 1])
  expect_gt(mean(p_out[is_outlier]), 0.7)
  expect_lt(mean(p_out[!is_outlier]), 0.1)
  w <- fit$outlier_weights[keep, 1]
  expect_true(abs(mean(w) - 0.17) < 0.1)     # 20 / 120 true outliers, Beta(2, 10) prior
  expect_true(all(w > 0 & w < 1))
})

test_that("outlier weight follows its Beta posterior (regression: hyperparameter was overwritten)", {
  # sampleOutlier() used to assign a random number to the Beta hyperparameter
  # `u`. When every item has an observed label none can be an outlier, so the
  # weight must be exactly Beta(2, N + 10) whatever the data are.
  skip_on_cran()
  set.seed(72)
  N <- 20
  X <- matrix(rnorm(N * 2), N); rownames(X) <- seq_len(N)
  labels <- matrix(rep(1:2, length.out = N), N, 1)
  fit <- callMDI(list(X), R = 12000, thin = 3, types = "TAGM", K = 2,
                 initial_labels = labels, fixed = matrix(1, N, 1))
  w <- fit$outlier_weights[-(1:50), 1]
  a <- 2; b <- N + 10
  expect_equal(mean(w), a / (a + b), tolerance = 0.03)
  expect_equal(var(w), a * b / ((a + b)^2 * (a + b + 1)), tolerance = 0.1)
})
