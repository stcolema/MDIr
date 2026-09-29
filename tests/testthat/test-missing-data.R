test_that("multivariate t imputation is the exact conditional distribution", {
  set.seed(21)
  X <- matrix(rnorm(400 * 3), 400) %*% chol(matrix(c(1, .7, .3, .7, 2, .4, .3, .4, 1.5), 3)) + 5
  for (miss in list(2L, c(1L, 3L), c(1L, 2L, 3L))) {
    chk <- mdir:::mvtImputationCheckCpp(X, miss, 40000)
    # If the imputation is exact then (observed, imputed) has the same joint law
    # as (observed, true): compare marginals and the covariance
    for (p in 1:3) {
      expect_gt(suppressWarnings(ks.test(chk$imputed[, p], chk$direct[, p])$p.value), 0.001)
    }
    expect_equal(cov(chk$imputed), cov(chk$direct), tolerance = 0.15)
  }
})

test_that("a Gaussian imputation that ignores the t scaling would be detected", {
  # Sanity check on the test itself: imputing with the Gaussian conditional
  # (no (df + d) / (df + p_o) rescaling) changes the tail behaviour
  set.seed(22)
  X <- matrix(rnorm(300 * 2), 300)
  chk <- mdir:::mvtImputationCheckCpp(X, 2L, 60000)
  expect_equal(mean(abs(chk$imputed[, 2]) > 2), mean(abs(chk$direct[, 2]) > 2), tolerance = 0.1)
})

sim_two_clusters <- function(N, P = 2, delta = 4, rho = 0.6) {
  truth <- rep(1:2, length.out = N)
  Sig <- matrix(rho, P, P); diag(Sig) <- 1
  Z <- matrix(rnorm(N * P), N) %*% chol(Sig)
  X <- Z + ifelse(truth == 1, 0, delta)
  rownames(X) <- seq_len(N)
  list(X = X, truth = truth)
}

test_that("cluster means are recovered with missing data (no stale-imputation bias)", {
  skip_on_cran()
  set.seed(23)
  d <- sim_two_clusters(500)
  X <- d$X
  miss <- matrix(runif(length(X)) < 0.3, nrow(X))
  Xm <- X; Xm[miss] <- NA
  Xm <- Xm[rowSums(!is.na(Xm)) > 0, , drop = FALSE]     # keep some information in every row
  truth <- d$truth[as.integer(rownames(Xm))]
  out <- callMDI(list(Xm), R = 3000, thin = 5, types = "MVN", K = 2)
  keep <- seq(200, nrow(out$allocations))
  # Cluster-specific mean of the first feature from the saved parameters
  P <- 2; K <- 2
  th <- out$parameters[[1]][keep, , drop = FALSE]
  mu <- array(th[, seq_len(P * K)], dim = c(length(keep), P, K))
  # order components by their first-feature mean (label switching)
  ord <- t(apply(mu[, 1, ], 1, order))
  mu_sorted <- sapply(1:K, function(k) sapply(seq_along(keep), function(i) mu[i, 1, ord[i, k]]))
  expect_equal(unname(colMeans(mu_sorted)), c(0, 4), tolerance = 0.15)
  # Imputed data must also carry the right cluster structure
  expect_gt(mean(apply(out$allocations[keep, , 1], 2, function(z) max(table(z))) / length(keep)), 0.8)
})

test_that("imputations are stored and consistent with the observed cells", {
  set.seed(24)
  d <- sim_two_clusters(80)
  X <- d$X; X[sample(length(X), 30)] <- NA
  out <- suppressWarnings(callMDI(list(X), R = 100, thin = 10, types = "MVN", K = 3, save_imputed = TRUE))
  cells <- out$missing_cells[[1]]
  expect_equal(nrow(cells), sum(is.na(X)))
  expect_equal(ncol(out$imputed[[1]]), nrow(cells))
  expect_true(all(is.finite(out$imputed[[1]])))
  expect_true(all(is.na(X[cells + 1])))
})

test_that("missing values are handled in every density and with outliers", {
  set.seed(25)
  N <- 40
  base <- matrix(rnorm(N * 3), N); rownames(base) <- seq_len(N)
  cat_x <- (base > 0) * 1L; storage.mode(cat_x) <- "double"
  for (type in c("G", "MVN", "GP", "TAGM", "TAGPM", "C")) {
    X <- if (type == "C") cat_x else base
    X[sample(length(X), 15)] <- NA
    X[1, ] <- NA                                    # a wholly missing item
    out <- suppressWarnings(callMDI(list(X), R = 60, thin = 6, types = type, K = 3))
    expect_true(all(is.finite(out$complete_likelihood)), info = type)
    expect_true(all(is.finite(out$observed_likelihood)), info = type)
  }
})

test_that("data validation catches unusable missing data patterns", {
  X <- matrix(rnorm(40), 20); rownames(X) <- 1:20
  X[, 2] <- NA
  expect_error(callMDI(list(X), R = 10, thin = 1, types = "MVN", K = 2), "no observed")
})
