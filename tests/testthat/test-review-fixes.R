# Regression tests for the defects found in the code review (see NEWS.md, "Corrections").

# Call the C++ sampler directly, bypassing the R-side recoding of the labels, so
# that the sampler's own handling of observed labels is exercised.
raw_run <- function(X, types, K, labels, fixed, R = 300, split_merge = 0L, betas = numeric(0)) {
  V <- length(X)
  mdir:::runMDI(
    R, 1L, X, as.integer(K), mdir:::translateTypes(types), mdir:::setupOutlierComponents(types),
    labels, fixed, rep(list(0), V), FALSE, FALSE, as.numeric(mdiPrior()), as.numeric(densityPrior()),
    rep(0L, V), FALSE, TRUE, betas, 0L, 1L, as.integer(split_merge)
  )
}

three_cluster <- function(N = 45, V = 2, seed = 7) {
  set.seed(seed)
  m <- rep(c(-4, 0, 4), each = N / 3)
  X <- lapply(seq_len(V), function(v) {
    x <- matrix(rnorm(N * 2, m), N, 2)
    rownames(x) <- seq_len(N)
    x
  })
  list(X = X, truth = rep(1:3, each = N / 3))
}

test_that("observed labels survive label swaps and split-merge whatever classes they use", {
  d <- three_cluster()
  N <- 45
  obs <- c(1:4, 31:34)          # classes 1 and 3 observed, class 2 never
  fixed <- cbind(rep(0, N), rep(0, N))
  fixed[obs, 1] <- 1
  fixed[c(1:3, 16:18), 2] <- 1

  # Observed classes {0, 2}: not the first components. The sampler must not move them.
  labels <- matrix(1L, N, 2)
  labels[1:4, 1] <- 0L
  labels[31:34, 1] <- 2L
  labels[c(1:3), 2] <- 4L
  labels[c(16:18), 2] <- 1L
  for (sm in c(0L, 2L)) {
    for (tp in list(c("G", "MVN"), c("TAGM", "G"), c("MVN", "TAGM"))) {
      out <- raw_run(d$X, tp, c(5, 5), labels, fixed, R = 200, split_merge = sm)
      a <- out$allocations
      for (v in 1:2) {
        idx <- which(fixed[, v] == 1)
        expect_true(all(apply(a[, idx, v], 1, function(r) all(r == labels[idx, v]))),
                    info = paste("view", v, "split_merge", sm, paste(tp, collapse = "/")))
      }
      # observed items are never flagged as outliers
      expect_true(all(out$outliers[, which(fixed[, 1] == 1), 1] == 0))
    }
  }
})

test_that("split-merge alone (one view, no swaps) never moves an observed item", {
  d <- three_cluster(V = 1)
  N <- 45
  fixed <- matrix(0, N, 1)
  fixed[c(1:4, 31:34), 1] <- 1
  labels <- matrix(1L, N, 1)
  labels[1:4, 1] <- 0L
  labels[31:34, 1] <- 3L
  for (tp in c("G", "MVN", "TAGM", "C")) {
    X <- d$X
    if (tp == "C") X <- list(matrix(sample(0:2, N * 3, TRUE), N, 3))
    out <- raw_run(X, tp, 5, labels, fixed, R = 300, split_merge = 3L)
    expect_gt(out$split_merge$attempts, 0)
    idx <- which(fixed[, 1] == 1)
    expect_true(all(apply(out$allocations[, idx, 1], 1, function(r) all(r == labels[idx, 1]))), info = tp)
    expect_true(all(out$outliers[, idx, 1] == 0), info = tp)
  }
})

test_that("a split-merge or swap run still moves the free items", {
  d <- three_cluster()
  N <- 45
  fixed <- matrix(0, N, 2)
  fixed[c(1:4, 31:34), 1] <- 1
  labels <- matrix(1L, N, 2)
  labels[1:4, 1] <- 0L
  labels[31:34, 1] <- 2L
  out <- raw_run(d$X, c("G", "G"), c(5, 5), labels, fixed, R = 200, split_merge = 2L)
  expect_gt(out$split_merge$accepts, 0)
  free <- setdiff(seq_len(N), which(fixed[, 1] == 1))
  expect_gt(length(unique(apply(out$allocations[, free, 1], 1, paste, collapse = ","))), 20)
})

test_that("callMDI recodes observed classes to contiguous labels and keeps them", {
  d <- three_cluster()
  N <- 45
  fixed <- matrix(0, N, 2)
  fixed[c(1:4, 31:34), 1] <- 1
  for (codes in list(c(0, 2), c(1, 3), c(2, 4))) {
    lab <- matrix(codes[1], N, 2)
    lab[1:4, 1] <- codes[1]
    lab[31:34, 1] <- codes[2]
    fit <- callMDI(d$X, R = 120, thin = 1, types = c("G", "G"), K = c(5, 5), initial_labels = lab,
                   fixed = fixed, check_prior = FALSE, split_merge = 1)
    a <- fit$allocations[, , 1]
    expect_true(all(a[, 1:4] == 0), info = paste(codes, collapse = ","))
    expect_true(all(a[, 31:34] == 1), info = paste(codes, collapse = ","))
    expect_equal(fit$allocation_probabilities[[1]][1, 1, ], rep(1, 121))
  }
  expect_equal(mdir:::.recodeSemiSupervised(c(7, 7, 3, 9, 1), c(1, 1, 1, 0, 0)), c(1, 1, 0, 3, 2))
})

test_that("a single observed class with any code is accepted", {
  d <- three_cluster()
  fixed <- matrix(0, 45, 2)
  fixed[1:5, 1] <- 1
  lab <- matrix(1, 45, 2)
  lab[1:5, 1] <- 2
  fit <- callMDI(d$X, R = 30, thin = 1, types = c("G", "G"), K = c(4, 4), initial_labels = lab,
                 fixed = fixed, check_prior = FALSE)
  expect_true(all(fit$allocations[, 1:5, 1] == 0))
})

test_that("initial semi-supervised labels follow the observed class proportions", {
  set.seed(1)
  lab <- c(rep(2, 9), 1, 1, 1)
  fixed <- c(rep(1, 10), 0, 0)
  draws <- replicate(3000, generateInitialSemiSupervisedLabels(lab, fixed)[11:12])
  expect_equal(mean(draws == 2), 0.9, tolerance = 0.02)
  expect_equal(generateInitialSemiSupervisedLabels(c(3, 3, 1, 1), c(1, 1, 0, 0)), c(3, 3, 3, 3))
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
  expect_equal(cc$evidence[, 1], vapply(ch, function(x) x$evidence[21], numeric(1)))
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
  # the likelihood of the initial state is recorded, not left at zero
  expect_true(all(gl$log_likelihood != 0))
})

test_that("Gaussian-process acceptance rates are proportions", {
  skip_on_cran()
  set.seed(5)
  N <- 30
  X <- list(matrix(rnorm(N * 6, rep(c(0, 2), each = 15)), N, 6))
  rownames(X[[1]]) <- 1:N
  fit <- callMDI(X, R = 200, thin = 5, types = "GP", K = 3, check_prior = FALSE)
  ac <- fit$acceptance_count[[1]]
  expect_true(all(ac >= 0 & ac <= 1))
  expect_gt(max(ac), 0)
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

test_that("the effective sample size agrees with the posterior package", {
  skip_if_not_installed("posterior")
  set.seed(1)
  for (rho in c(0, 0.6, 0.9, -0.5)) {
    x <- sapply(1:4, function(i) as.numeric(suppressWarnings(arima.sim(list(ar = rho), 800))))
    mine <- rankNormalizedRhat(x)
    suppressWarnings({
      expect_equal(mine$ess_bulk, posterior::ess_bulk(x), tolerance = 1e-3)
      expect_equal(mine$ess_tail, posterior::ess_tail(x), tolerance = 1e-3)
      expect_equal(mine$rhat, posterior::rhat(x), tolerance = 1e-3)
    })
  }
})
