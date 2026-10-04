# Split-merge move (sequentially allocated proposal, collapsed marginal likelihoods).
# References are independent of the sampler: quadrature, a closed form, exact enumeration.
# See verification/tempering/ for the larger experiments and the proofs.

dp0 <- as.numeric(densityPrior(scale_pool_shape = 0))
marg <- mdir:::collapsedLogMarginalCpp

test_that("collapsed marginals match independent calculations", {
  set.seed(1)
  # G, one dimension: quadrature over (mu, sigma^2)
  X <- matrix(rnorm(8, 1, 1.3), 8, 1)
  hp <- mdir:::densityHyperparameters(X, 3, 0, dp0)
  rows <- c(0, 2, 3, 6)
  x <- X[rows + 1, 1]
  for (beta in c(1, 0.4)) {
    f <- function(mu, s2) {
      ll <- vapply(seq_along(mu), function(i) sum(dnorm(x, mu[i], sqrt(s2), log = TRUE)), 0)
      exp(beta * ll + dnorm(mu, hp$xi, sqrt(s2 / hp$kappa), log = TRUE) +
        (hp$nu / 2) * log(hp$scale / 2) - lgamma(hp$nu / 2) - (hp$nu / 2 + 1) * log(s2) - hp$scale / (2 * s2))
    }
    q <- integrate(function(s2) vapply(s2, function(v) suppressWarnings(integrate(function(mu) f(mu, v), -Inf, Inf, rel.tol = 1e-10)$value), 0),
      0, Inf, rel.tol = 1e-8)$value
    expect_equal(marg(X, 3, 0, dp0, rows, beta), log(q), tolerance = 1e-6)
  }
  # C, two categories: quadrature over the success probability
  Xc <- matrix(c(0, 1, 1, 0, 0, 1, 1, 1, 0, 0), ncol = 1)
  hpc <- mdir:::densityHyperparameters(Xc, 2, 2, dp0)
  al <- as.numeric(hpc$concentration[[1]])
  rows <- c(0, 1, 2, 4, 5)
  y <- Xc[rows + 1, 1]
  for (beta in c(1, 0.5)) {
    q <- integrate(function(th) exp(beta * (sum(y == 0) * log(th) + sum(y == 1) * log(1 - th))) * dbeta(th, al[1], al[2]), 0, 1,
      rel.tol = 1e-10)$value
    expect_equal(marg(Xc, 2, 2, dp0, rows, beta), log(q), tolerance = 1e-8)
  }
  # MVN, P = 2: normal-inverse-Wishart closed form (beta = 1)
  Xm <- matrix(rnorm(18), 9, 2)
  hm <- mdir:::densityHyperparameters(Xm, 3, 1, dp0)
  rows <- c(1, 2, 4, 7, 8)
  Xs <- Xm[rows + 1, , drop = FALSE]
  n <- nrow(Xs); P <- 2
  xb <- colMeans(Xs); S <- crossprod(sweep(Xs, 2, xb))
  kn <- hm$kappa + n; nun <- hm$nu + n
  Psin <- hm$scale + S + hm$kappa * n / kn * tcrossprod(xb - as.numeric(hm$xi))
  lmv <- function(a) P * (P - 1) / 4 * log(pi) + sum(lgamma(a + (1 - seq_len(P)) / 2))
  ld <- function(M) as.numeric(determinant(M, logarithm = TRUE)$modulus)
  ref <- -n * P / 2 * log(pi) + lmv(nun / 2) - lmv(hm$nu / 2) + hm$nu / 2 * ld(hm$scale) - nun / 2 * ld(Psin) +
    P / 2 * (log(hm$kappa) - log(kn))
  expect_equal(marg(Xm, 3, 1, dp0, rows, 1), ref, tolerance = 1e-9)
})

# exact distribution of the labels under the collapsed target for a tiny categorical problem
sm_exact <- function(X, K, w, beta, fixed, lab_fixed) {
  free <- which(fixed == 0)
  grid <- as.matrix(expand.grid(rep(list(0:(K - 1)), length(free))))
  lp <- apply(grid, 1, function(g) {
    lab <- lab_fixed; lab[free] <- g
    s <- sum(log(w[lab + 1]))
    for (k in 0:(K - 1)) {
      rows <- which(lab == k) - 1
      if (length(rows)) s <- s + marg(X, K, 2, dp0, rows, beta)
    }
    s
  })
  p <- exp(lp - max(lp))
  list(grid = grid, prob = p / sum(p), free = free)
}
sm_tv <- function(X, K, w, beta, fixed, lab_fixed, ref_beta = beta, n_chains = 4, n_iter = 100000) {
  ex <- sm_exact(X, K, w, ref_beta, fixed, lab_fixed)
  code <- function(m) as.vector(m %*% K^(seq_len(ncol(m)) - 1))  # expand.grid: first column fastest
  counts <- numeric(nrow(ex$grid))
  for (i in seq_len(n_chains)) {
    set.seed(300 + i)
    init <- sample(0:(K - 1), length(fixed), TRUE)
    init[fixed == 1] <- lab_fixed[fixed == 1]
    out <- mdir:::splitMergeOnlyCpp(X, K, 2, dp0, init, fixed, w, n_iter, beta)
    lab <- out$labels[-(1:200), ex$free, drop = FALSE]
    counts <- counts + tabulate(code(lab) + 1, nrow(ex$grid))
  }
  0.5 * sum(abs(counts / sum(counts) - ex$prob))
}

test_that("the move alone leaves the collapsed label posterior invariant (exact enumeration)", {
  set.seed(11)
  N <- 5; K <- 3
  X <- matrix(sample(0:2, N * 2, TRUE, c(.5, .3, .2)), N, 2)
  w <- c(0.5, 0.3, 0.2)
  expect_lt(sm_tv(X, K, w, 1, rep(0, N), rep(0, N)), 0.03)
  expect_lt(sm_tv(X, K, w, 0.5, rep(0, N), rep(0, N)), 0.03)
  # semi-supervised: items 1 and 2 observed in classes 0 and 1
  expect_lt(sm_tv(X, K, w, 1, c(1, 1, 0, 0, 0), c(0, 1, 0, 0, 0)), 0.03)
  # negative control: the wrong tempered target is rejected
  expect_gt(sm_tv(X, K, w, 0.5, rep(0, N), rep(0, N), ref_beta = 1), 0.05)
})

test_that("split-merge is accepted for G, MVN, C and combinations, and refused otherwise", {
  set.seed(2)
  X <- matrix(rnorm(60, rep(c(0, 4), each = 30)), 30, 2); rownames(X) <- 1:30
  Xc <- matrix(sample(0:2, 60, TRUE), 30, 2, dimnames = list(1:30, NULL))
  fit <- callMDI(list(X, Xc), R = 40, thin = 1, types = c("MVN", "C"), K = c(3, 3), split_merge = 2L, check_prior = FALSE)
  expect_equal(fit$split_merge$moves, 2L)
  expect_gt(fit$split_merge$attempts, 0)
  # semi-supervised: fixed labels are not moved
  lab <- matrix(c(rep(0:1, length.out = 6), rep(0, 24)), ncol = 1)
  fx <- matrix(c(rep(1, 6), rep(0, 24)), ncol = 1)
  fit2 <- callMDI(list(X), R = 60, thin = 1, types = "MVN", K = 3, initial_labels = lab, fixed = fx, split_merge = 1L,
    check_prior = FALSE)
  expect_true(all(apply(fit2$allocations[, 1:6, 1], 1, identical, as.numeric(lab[1:6, 1]))))
  expect_error(callMDI(list(X), R = 20, thin = 1, types = "TAGM", K = 3, split_merge = 1L, check_prior = FALSE), "split-merge")
  expect_error(callMDI(list(X), R = 20, thin = 1, types = "MVN", K = 3, split_merge = -1L, check_prior = FALSE), "split_merge")
})

test_that("split-merge changes the chain but not at zero moves", {
  set.seed(3)
  X <- matrix(rnorm(40, rep(c(0, 4), each = 20)), 20, 2); rownames(X) <- 1:20
  set.seed(5); a <- callMDI(list(X), R = 30, thin = 1, types = "MVN", K = 3, check_prior = FALSE)
  set.seed(5); b <- callMDI(list(X), R = 30, thin = 1, types = "MVN", K = 3, split_merge = 0L, check_prior = FALSE)
  expect_identical(a$allocations, b$allocations)
})
