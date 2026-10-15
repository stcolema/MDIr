# Adaptation of the Gaussian process proposal windows during burn in: towards the
# acceptance rate of a one-dimensional random-walk Metropolis update (about 0.44), by a
# Robbins-Monro recursion that diminishes and is stopped after `n_burn` sweeps, so that
# the chain after burn in has fixed windows.

gp_views <- function(N = 30, P = 6) {
  X <- matrix(rnorm(N * P, rep(c(0, 2), each = N / 2)), N, P)
  rownames(X) <- seq_len(N)
  list(X)
}

# Mean acceptance rate of each of the amplitude, length and noise updates
gp_rates <- function(fit) {
  ac <- fit$acceptance_count[[1]]
  K <- length(ac) / 3
  # an empty component makes no proposals, and shows a rate of 0
  vapply(0:2, function(h) mean(ac[h * K + seq_len(K)][ac[h * K + seq_len(K)] > 0]), numeric(1))
}

test_that("badly scaled windows are moved towards the target, and a run started from them accepts at about that rate", {
  set.seed(1)
  X <- gp_views()
  for (start in list(rep(0.002, 3), rep(8, 3))) {
    fit <- callMDI(X, n_iter = 6000, thin = 10, types = "GP", K_max = 2, check_prior = FALSE,
                   control = mdirControl(proposal_windows = list(start), auto_tune = TRUE))
    final <- fit$final_proposal_windows[[1]]
    expect_length(final, 3)
    if (start[1] < 1) expect_true(all(final > 10 * start)) else expect_true(all(final < start / 3))

    # the adapted windows, kept fixed: acceptance close to the target
    set.seed(2)
    again <- callMDI(X, n_iter = 4000, thin = 10, types = "GP", K_max = 2, check_prior = FALSE,
                     control = mdirControl(proposal_windows = list(final), auto_tune = FALSE))
    rates <- gp_rates(again)
    expect_true(all(rates > 0.2 & rates < 0.7), label = paste("acceptance rates", paste(round(rates, 2), collapse = ", ")))
  }

  # without adaptation the same starting windows stay as they are, and accept at far from the target
  fixed <- callMDI(X, n_iter = 2000, thin = 10, types = "GP", K_max = 2, check_prior = FALSE,
                   control = mdirControl(proposal_windows = list(rep(0.002, 3)), auto_tune = FALSE))
  expect_equal(fixed$final_proposal_windows[[1]], rep(0.002, 3))
  expect_true(all(gp_rates(fixed) > 0.9))
})

test_that("auto_tune = FALSE, n_burn = 0 and views without hyperparameters leave the chain as it was", {
  X <- gp_views()
  run <- function(...) {
    set.seed(3)
    callMDI(X, n_iter = 300, thin = 5, types = "GP", K_max = 2, check_prior = FALSE, ...)
  }
  off <- run(control = mdirControl(auto_tune = FALSE))
  zero <- run(control = mdirControl(auto_tune = TRUE, n_burn = 0))
  expect_identical(off$hypers, zero$hypers)
  expect_identical(off$allocations, zero$allocations)
  expect_equal(zero$final_proposal_windows[[1]], c(0.125, 0.15, 0.10))
  on <- run(control = mdirControl(auto_tune = TRUE))
  expect_false(identical(off$hypers, on$hypers))

  # adaptation uses no random numbers: for views with no Gaussian process it changes nothing
  views <- lapply(1:2, function(v) {
    m <- matrix(rnorm(40 * 2, rep(c(0, 3), each = 20)), 40, 2)
    rownames(m) <- 1:40
    m
  })
  mvn <- function(auto_tune) {
    set.seed(4)
    callMDI(views, n_iter = 100, thin = 5, types = c("MVN", "G"), K_max = c(3, 3),
            control = mdirControl(auto_tune = auto_tune))
  }
  expect_identical(mvn(TRUE)$allocations, mvn(FALSE)$allocations)
  expect_identical(mvn(TRUE)$phis, mvn(FALSE)$phis)
})

test_that("the windows are fixed once the burn in is over", {
  X <- gp_views()
  run <- function(n_iter, n_burn) {
    set.seed(5)
    callMDI(X, n_iter = n_iter, thin = 5, types = "GP", K_max = 2, check_prior = FALSE,
            control = mdirControl(auto_tune = TRUE, n_burn = n_burn))
  }
  short <- run(n_iter = 500, n_burn = 500)
  long <- run(n_iter = 1500, n_burn = 500)
  # the first 500 sweeps use the same random numbers, and nothing adapts after them
  expect_equal(long$final_proposal_windows[[1]], short$final_proposal_windows[[1]])
  expect_false(isTRUE(all.equal(short$final_proposal_windows[[1]], c(0.125, 0.15, 0.10))))
  expect_equal(long$n_burn, 500)
  # a longer adaptation goes further
  longer <- run(n_iter = 1500, n_burn = 1500)
  expect_false(isTRUE(all.equal(longer$final_proposal_windows[[1]], long$final_proposal_windows[[1]])))
})

test_that("the arguments are checked and reach every entry point", {
  X <- gp_views()
  expect_error(mdirControl(auto_tune = NA), "auto_tune")
  expect_error(mdirControl(n_burn = -1), "n_burn")
  expect_error(mdirControl(n_burn = c(1, 2)), "n_burn")

  chains <- runMCMCChains(X, 2, n_iter = 100, thin = 5, types = "GP", K_max = 2,
                          control = mdirControl(auto_tune = FALSE, n_burn = 7))
  expect_false(chains[[1]]$auto_tune)
  expect_equal(chains[[2]]$n_burn, 7)
  fit <- fitMDI(X, 2, n_iter = 100, thin = 5, types = "GP", K_max = 2, control = mdirControl(n_burn = 30), verbose = FALSE)
  expect_equal(fit[[1]]$n_burn, 30)

  # writing the samples to file adapts too, and runs
  dir <- tempfile("mdir_autotune")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  expect_no_error(callMDIWritingToFile(X, n_iter = 60, thin = 5, types = "GP", K_max = 2,
                                       dir_path = dir, control = mdirControl(auto_tune = TRUE, n_burn = 30)))
})
