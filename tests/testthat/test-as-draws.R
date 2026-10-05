draws_fit <- function() {
  set.seed(21)
  X <- lapply(1:2, function(v) {
    m <- matrix(rnorm(40 * 2, rep(c(0, 3), each = 20)), 40, 2)
    rownames(m) <- 1:40
    m
  })
  list(X = X, fit = runMCMCChains(X, 3, R = 200, thin = 5, types = c("G", "MVN"), K = c(4, 4)))
}

test_that("as_draws() gives a draws_array of the label-invariant quantities", {
  skip_if_not_installed("posterior")
  d <- draws_fit()
  dr <- posterior::as_draws(d$fit)
  expect_s3_class(dr, "draws_array")
  # 41 saved samples; the initial state and the first 20 (burn = R / 2 = 100) are dropped
  expect_equal(posterior::niterations(dr), 20)
  expect_equal(posterior::nchains(dr), 3)
  expect_equal(posterior::variables(dr),
    c("log_lik_complete", "log_lik_observed", "log_lik_joint", "log_normalising_constant",
      "mass[1]", "mass[2]", "phi[1,2]", "occupied[1]", "pooled[1,1]", "pooled[1,2]",
      "occupied[2]", "pooled[2,1]", "pooled[2,2]"))
  ch2 <- unclass(d$fit)[[2]]
  expect_equal(as.numeric(dr[, 2, "mass[1]"]), ch2$mass[22:41, 1])
  expect_equal(as.numeric(dr[, 2, "phi[1,2]"]), ch2$phis[22:41, 1])
  expect_equal(as.numeric(dr[, 2, "log_lik_joint"]), ch2$joint_likelihood[22:41])
  expect_equal(as.numeric(dr[, 2, "log_normalising_constant"]), log(ch2$normalising_constant[22:41]))
  expect_equal(as.numeric(dr[, 2, "occupied[1]"]),
               apply(ch2$allocations[22:41, , 1], 1, function(z) length(unique(z))))
})

test_that("the other formats and a single chain work, and burn is respected", {
  skip_if_not_installed("posterior")
  d <- draws_fit()
  expect_s3_class(posterior::as_draws_df(d$fit), "draws_df")
  expect_s3_class(posterior::as_draws_matrix(d$fit), "draws_matrix")
  expect_s3_class(posterior::as_draws_list(d$fit), "draws_list")
  expect_s3_class(posterior::as_draws_rvars(d$fit), "draws_rvars")
  one <- posterior::as_draws_df(unclass(d$fit)[[1]])
  expect_equal(nrow(one), 20)
  expect_equal(posterior::niterations(posterior::as_draws(d$fit, burn = 150)), 10)
  expect_error(posterior::as_draws(d$fit, burn = 200), "no saved samples")
})

test_that("processed and unprocessed chains give the same draws", {
  skip_if_not_installed("posterior")
  d <- draws_fit()
  raw <- posterior::as_draws(d$fit, burn = 100)
  pr <- posterior::as_draws(processMCMCChains(d$fit, burn = 100))
  expect_equal(unclass(raw), unclass(pr), ignore_attr = TRUE)
})

test_that("the posterior diagnostics agree with those of assessConvergence", {
  skip_if_not_installed("posterior")
  d <- draws_fit()
  s <- suppressWarnings(posterior::summarise_draws(posterior::as_draws(d$fit), "rhat", "ess_bulk", "ess_tail"))
  conv <- assessConvergence(d$fit, burn = 100)
  # the variables are named differently; compare the shared ones
  pairs <- c("phi[1,2]" = "phi[1,2]", "mass[1]" = "mass[1]", "occupied[1]" = "occupied_components[1]",
             "log_lik_complete" = "complete_likelihood", "log_lik_joint" = "joint_likelihood")
  for (v in names(pairs)) {
    a <- s[s$variable == v, ]
    b <- conv[conv$quantity == pairs[[v]], ]
    expect_equal(a$ess_bulk, b$ess_bulk, tolerance = 1e-6, info = v)
    # posterior gives NA when one of the two tail indicators is constant (a discrete variable)
    if (!is.na(a$ess_tail)) expect_equal(a$ess_tail, b$ess_tail, tolerance = 1e-6, info = v)
    expect_equal(a$rhat, b$rhat, tolerance = 1e-6, info = v)
  }
})

test_that("weights and allocations are added on request", {
  skip_if_not_installed("posterior")
  d <- draws_fit()
  dr <- posterior::as_draws(d$fit, weights = TRUE, allocations = TRUE)
  v <- posterior::variables(dr)
  expect_true(all(sprintf("sorted_weight[1,%d]", 1:4) %in% v))
  expect_equal(sum(grepl("^allocation", v)), 2 * 40)
  w <- posterior::subset_draws(dr, variable = sprintf("sorted_weight[1,%d]", 1:4))
  expect_true(all(abs(apply(w[, 1, ], 1, sum) - 1) < 1e-12))
  expect_true(all(apply(w[, 1, ], 1, function(r) !is.unsorted(rev(r)))))
  a <- as.numeric(dr[, 3, "allocation[2,7]"])
  expect_equal(a, unclass(d$fit)[[3]]$allocations[22:41, 7, 2] + 1)
})

test_that("an smcMDI ensemble becomes weighted draws", {
  skip_if_not_installed("posterior")
  set.seed(3)
  X <- list(matrix(rnorm(60, rep(c(0, 3), each = 15)), 30, 2), matrix(rnorm(60, rep(c(0, 3), each = 15)), 30, 2))
  rownames(X[[1]]) <- rownames(X[[2]]) <- 1:30
  s <- smcMDI(X, c("G", "G"), K = c(3, 3), n_particles = 20, final_sweeps = 4, thin = 2)
  dr <- posterior::as_draws(s)
  expect_s3_class(dr, "draws_matrix")
  expect_equal(posterior::ndraws(dr), 20 * 3)
  expect_equal(sum(stats::weights(dr)), 1)
  expect_equal(stats::weights(dr), smcWeights(s, "draw"))
  expect_equal(as.numeric(dr[1:3, "mass[1]"]), s$mass[, 1, 1])
  expect_true(all(c("log_lik_complete", "phi[1,2]", "occupied[2]") %in% posterior::variables(dr)))
})

test_that("pointwise log-likelihoods feed loo", {
  skip_if_not_installed("loo")
  set.seed(5)
  X <- lapply(1:2, function(v) {
    m <- matrix(rnorm(30 * 2, rep(c(0, 3), each = 15)), 30, 2)
    rownames(m) <- 1:30
    m
  })
  ch <- runMCMCChains(X, 2, R = 120, thin = 3, types = c("G", "G"), K = c(3, 3), save_pointwise = TRUE)
  ll <- pointwiseLogLik(ch, burn = 60)
  r_eff <- suppressWarnings(loo::relative_eff(exp(ll), chain_id = attr(ll, "chain_id")))
  fit <- suppressWarnings(loo::loo(ll, r_eff = r_eff))
  expect_s3_class(fit, "loo")
  expect_equal(nrow(fit$pointwise), 30)
  expect_true(is.finite(fit$estimates["elpd_loo", "Estimate"]))
})
