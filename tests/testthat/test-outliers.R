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

test_that("outlier weight follows its Beta posterior", {
  # When every item has an observed label none can be an outlier and none
  # carries a factor of the outlier weight, so the weight keeps its Beta(2, 10)
  # prior whatever the data are.
  skip_on_cran()
  set.seed(72)
  N <- 20
  X <- matrix(rnorm(N * 2), N); rownames(X) <- seq_len(N)
  labels <- matrix(rep(1:2, length.out = N), N, 1)
  fit <- callMDI(list(X), R = 12000, thin = 3, types = "TAGM", K = 2,
                 initial_labels = labels, fixed = matrix(1, N, 1))
  w <- fit$outlier_weights[-(1:50), 1]
  a <- 2; b <- 10
  expect_equal(mean(w), a / (a + b), tolerance = 0.03)
  expect_equal(var(w), a * b / ((a + b)^2 * (a + b + 1)), tolerance = 0.1)
})

test_that("the outlier weight update counts only items without an observed label", {
  skip_on_cran()
  set.seed(4)
  N <- 60
  X <- list(matrix(rnorm(N * 2, rep(c(-3, 3), each = N / 2)), N, 2))
  rownames(X[[1]]) <- seq_len(N)
  X[[1]][c(52, 56, 60), ] <- c(25, -25, 30, 30, -30, 20)
  lab <- matrix(rep(c(1, 2), each = N / 2), N, 1)
  fixed <- matrix(0, N, 1)
  fixed[c(1:15, 31:45), 1] <- 1
  fit <- callMDI(X, R = 4000, thin = 1, types = "TAGM", K = 2, initial_labels = lab,
                 fixed = fixed, check_prior = FALSE)
  n_free <- sum(fixed == 0)
  n_out <- rowSums(fit$outliers[, , 1])
  eps <- fit$outlier_weights[, 1]
  d <- data.frame(k = n_out[-length(n_out)], e = eps[-1])
  agg <- aggregate(e ~ k, d, function(x) c(mean = mean(x), n = length(x)))
  agg <- data.frame(k = agg$k, m = agg$e[, "mean"], n = agg$e[, "n"])
  agg <- agg[agg$n >= 150, ]
  expect_gt(nrow(agg), 1)
  prior <- c(a = 2, b = 10)            # Beta(2, 10) prior on the outlier weight
  # Beta(a + n_out, b + n_free - n_out): the mean given the previous count of outliers
  expected <- (prior["a"] + agg$k) / (prior["a"] + prior["b"] + n_free)
  expect_equal(agg$m, unname(expected), tolerance = 0.12)
  # and clearly not the mean when the observed items are counted as non-outliers
  wrong <- (prior["a"] + agg$k) / (prior["a"] + prior["b"] + N)
  expect_gt(min(agg$m - wrong), 0.01)
})
