ar1_chains <- function(n, m, rho, shift = 0) {
  sapply(seq_len(m), function(j) {
    x <- if (rho == 0) stats::rnorm(n) else as.numeric(stats::arima.sim(list(ar = rho), n))
    x + shift * (j - 1)
  })
}

test_that("Rhat and ESS agree with the posterior package", {
  skip_if_not_installed("posterior")
  set.seed(31)
  cases <- list(
    iid = ar1_chains(1000, 4, 0),
    ar_moderate = ar1_chains(1000, 4, 0.6),
    ar_strong = ar1_chains(2000, 3, 0.9),
    offset = ar1_chains(500, 4, 0.3, shift = 0.5),
    heavy_tail = matrix(stats::rt(4000, 2), 1000)
  )
  for (nm in names(cases)) {
    x <- cases[[nm]]
    mine <- rankNormalizedRhat(x)
    expect_equal(mine$rhat, posterior::rhat(x), tolerance = 1e-6, info = nm)
    expect_equal(mine$ess_bulk, suppressWarnings(posterior::ess_bulk(x)), tolerance = 0.1, info = nm)
    expect_equal(mine$ess_tail, suppressWarnings(posterior::ess_tail(x)), tolerance = 0.15, info = nm)
  }
})

test_that("Rhat detects chains that have not mixed and chains that differ in scale", {
  set.seed(32)
  expect_lt(rankNormalizedRhat(matrix(rnorm(4000), ncol = 4))$rhat, 1.01)
  expect_gt(rankNormalizedRhat(cbind(rnorm(1000), rnorm(1000, 1), rnorm(1000, -1)))$rhat, 1.1)
  # Same location, different scale: only the folded (tail) Rhat sees this
  scale_diff <- cbind(rnorm(1000, sd = 1), rnorm(1000, sd = 4), rnorm(1000, sd = 1), rnorm(1000, sd = 4))
  out <- rankNormalizedRhat(scale_diff)
  expect_gt(out$rhat_tail, 1.05)
  expect_gt(out$rhat, out$rhat_bulk)
})

test_that("ESS of independent draws is close to the number of draws and shrinks with autocorrelation", {
  set.seed(33)
  iid <- rankNormalizedRhat(matrix(rnorm(8000), ncol = 4))
  expect_gt(iid$ess_bulk, 0.8 * 8000)
  ar <- rankNormalizedRhat(ar1_chains(2000, 4, 0.9))
  # Theoretical ESS for AR(1) with rho = 0.9 is n (1 - rho) / (1 + rho) = 8000 * 0.0526
  expect_equal(ar$ess_bulk, 8000 * 0.1 / 1.9, tolerance = 0.35)
})

test_that("non-finite input gives NA rather than an error", {
  x <- matrix(rnorm(400), ncol = 4); x[3, 2] <- NA
  expect_true(is.na(rankNormalizedRhat(x)$rhat))
})

test_that("assessConvergence monitors label-switching-invariant quantities", {
  skip_on_cran()
  set.seed(34)
  X <- lapply(1:2, function(v) {
    m <- matrix(rnorm(60 * 2, rep(c(0, 3), each = 30)), 60, 2); rownames(m) <- 1:60; m
  })
  chains <- runMCMCChains(X, 3, R = 600, thin = 3, types = c("MVN", "G"), K = c(4, 4))
  conv <- assessConvergence(chains, burn = 300)
  expect_s3_class(conv, "mdir_convergence")
  expect_true(all(c("complete_likelihood", "phi[1,2]", "mass[1]") %in% conv$quantity))
  expect_true(all(is.finite(conv$rhat)))
  expect_output(print(conv), "convergence diagnostics")
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
