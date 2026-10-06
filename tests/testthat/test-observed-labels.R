# Observed labels under label swaps and split-merge, and their initial coding.

# Call the C++ sampler directly, bypassing the R-side recoding of the labels, so
# that the sampler's own handling of observed labels is exercised.
raw_run <- function(X, types, K, labels, fixed, R = 300, split_merge = 0L, betas = numeric(0)) {
  V <- length(X)
  mdir:::runMDI(
    R, 1L, X, as.integer(K), mdir:::translateTypes(types), mdir:::setupOutlierComponents(types),
    labels, fixed, rep(list(0), V), FALSE, FALSE, as.numeric(mdiPrior()), as.numeric(densityPrior()),
    rep(0L, V), FALSE, TRUE, betas, 0L, 1L, as.integer(split_merge), 0L
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
