# print(), summary() and format() methods for the fit objects (see
# R/mdirFitMethods.R). The objects are plain lists underneath, so every
# existing `$` / `[[` access must keep working; that is checked alongside the
# console-facing behaviour.

make_views <- function(N = 40, V = 2) {
  lapply(seq_len(V), function(v) {
    m <- matrix(stats::rnorm(N * 2, rep(c(0, 3), each = N)), N, 2)
    rownames(m) <- seq_len(N)
    m
  })
}

# A synthetic mdir_convergence object, to test the wording of each verdict
make_convergence <- function(rhat, ess, n_chains = 4, threshold = 1.01) {
  out <- data.frame(
    quantity = paste0("q", seq_along(rhat)), rhat = rhat, rhat_bulk = rhat,
    rhat_tail = rhat, ess_bulk = ess, ess_tail = ess, stringsAsFactors = FALSE
  )
  out$converged <- rhat < threshold & ess >= 100 * n_chains
  attr(out, "threshold") <- threshold
  attr(out, "min_ess") <- 100 * n_chains
  attr(out, "burn") <- 100
  attr(out, "n_chains") <- n_chains
  class(out) <- c("mdir_convergence", "data.frame")
  out
}

test_that("callMDI() returns an mdir_fit that prints briefly and keeps `$` access", {
  set.seed(1)
  X <- make_views()
  fit <- callMDI(X, R = 100, thin = 5, types = c("MVN", "MVN"), K = c(3, 3))

  expect_s3_class(fit, "mdir_fit")
  expect_true(is.list(fit))
  expect_equal(fit$N, 40)
  expect_equal(dim(fit$allocations), c(21, 40, 2))

  out <- capture.output(print(fit))
  expect_lt(length(out), 12)
  expect_match(out[1], "<mdir fit> 2 views, 40 items", fixed = TRUE)
  expect_true(any(grepl("burn in not applied", out)))
  expect_false(any(grepl("allocations", out)))
  expect_identical(withVisible(capture.output(res <- print(fit)))$visible, TRUE)
  expect_identical(res, fit)
})

test_that("summary.mdir_fit() reports views, phi and, once processed, clustering tables", {
  set.seed(2)
  X <- make_views()
  fit <- callMDI(X, R = 200, thin = 5, types = c("MVN", "MVN"), K = c(3, 3))

  smry <- summary(fit)
  expect_s3_class(smry, "summary.mdir_fit")
  expect_equal(smry$burn, 100)
  expect_true(smry$burn_default)
  expect_equal(smry$n_samples, length(seq(floor(100 / 5) + 2, 41)))
  expect_equal(rownames(smry$phi), "phi[1,2]")
  expect_equal(colnames(smry$phi), c("mean", "sd", "2.5%", "50%", "97.5%"))
  expect_null(smry$clustering)

  out <- capture.output(print(smry))
  expect_true(any(grepl("MDI model fitted by MCMC", out)))
  expect_true(any(grepl("Occupied components", out)))
  expect_true(any(grepl("phi\\[1,2\\]", out)))
  expect_true(any(grepl("processMCMCChain", out)))

  # An explicit burn in replaces the default
  expect_equal(summary(fit, burn = 50)$burn, 50)
  expect_false(summary(fit, burn = 50)$burn_default)
  expect_error(summary(fit, burn = 1000), "leaves no saved samples")

  processed <- processMCMCChain(fit, burn = 100)
  expect_s3_class(processed, "mdir_fit")
  smry_p <- summary(processed)
  expect_false(smry_p$burn_default)
  expect_equal(smry_p$burn, 100)
  expect_length(smry_p$clustering, 2)
  expect_equal(sum(smry_p$clustering[[1]]), 40)
  out_p <- capture.output(print(smry_p))
  expect_true(any(grepl("Clustering table, view 1", out_p)))
  expect_true(any(grepl("burn = 100 applied", capture.output(print(processed)))))
})

test_that("a single view has no phi table", {
  set.seed(3)
  fit <- callMDI(make_views(V = 1), R = 60, thin = 5, types = "MVN", K = 3)
  expect_null(summary(fit)$phi)
  expect_no_error(capture.output(print(summary(fit))))
})

test_that("runMCMCChains() returns an mdir_fit_list, reports progress only when asked, and validates n_chains", {
  set.seed(4)
  X <- make_views()
  run <- function(...) {
    runMCMCChains(X, 2, R = 60, thin = 5, types = c("MVN", "MVN"), K = c(3, 3), ...)
  }

  expect_no_message(chains <- run())
  expect_s3_class(chains, "mdir_fit_list")
  expect_s3_class(chains[[1]], "mdir_fit")
  expect_equal(chains[[2]]$Chain, 2)

  msgs <- testthat::capture_messages(run(verbose = TRUE))
  expect_true(any(grepl("Chain 1/2: running 60 iterations", msgs)))
  expect_true(any(grepl("Chain 2/2: finished in", msgs)))

  expect_error(runMCMCChains(X, 0, R = 60, thin = 5, types = c("MVN", "MVN")), "n_chains")

  out <- capture.output(print(chains))
  expect_lt(length(out), 12)
  expect_match(out[1], "<mdir fit list> 2 chains", fixed = TRUE)
  expect_true(any(grepl("Convergence has not been assessed", out)))

  # `[` keeps the class, and drops diagnostics that refer to the full set
  attr(chains, "convergence") <- "stale"
  sub <- chains[1]
  expect_s3_class(sub, "mdir_fit_list")
  expect_null(attr(sub, "convergence"))
  expect_length(sub, 1)
})

test_that("fitMDI() attaches convergence, reports it when verbose, and stays quiet otherwise", {
  set.seed(5)
  X <- make_views()
  fit_args <- list(X, n_chains = 3, R = 300, thin = 5, types = c("MVN", "MVN"), K = c(3, 3), burn = 150)

  expect_no_message(fit <- do.call(fitMDI, c(fit_args, verbose = FALSE)))
  expect_s3_class(fit, "mdir_fit_list")
  conv <- attr(fit, "convergence")
  expect_s3_class(conv, "mdir_convergence")
  expect_equal(attr(conv, "burn"), 150)
  expect_length(attr(conv, "chain_loglik"), 3)

  msgs <- testthat::capture_messages(do.call(fitMDI, c(fit_args, verbose = TRUE)))
  expect_true(any(grepl("Fitting MDI to 40 items in 2 views (MVN, MVN): 3 chains of 300 iterations", msgs, fixed = TRUE)))
  expect_true(any(grepl("Chain 3/3: finished in", msgs)))
  expect_true(any(grepl("^Convergence:", msgs)))

  out <- capture.output(print(fit))
  expect_true(any(grepl("^Convergence:", out)))

  smry <- summary(fit)
  expect_s3_class(smry, "summary.mdir_fit_list")
  expect_length(smry$per_chain, 3)
  out_s <- capture.output(print(smry))
  expect_true(any(grepl("Mean complete log-lik", out_s)))
  expect_true(any(grepl("MDI convergence diagnostics: 3 chains", out_s)))

  # The convergence attribute survives processing
  processed <- processMCMCChains(fit, burn = 150)
  expect_s3_class(processed, "mdir_fit_list")
  expect_identical(attr(processed, "convergence"), conv)
  expect_true(all(vapply(processed, function(ch) !is.null(ch$pred), logical(1))))
})

test_that("fitMDI() warns, rather than failing silently, when diagnostics cannot be computed", {
  set.seed(6)
  X <- make_views()
  expect_warning(
    fit <- fitMDI(X, n_chains = 2, R = 20, thin = 5, types = c("MVN", "MVN"), K = c(3, 3), verbose = FALSE),
    "Convergence diagnostics could not be computed"
  )
  expect_s3_class(fit, "mdir_fit_list")
  expect_null(attr(fit, "convergence"))
})

test_that("the convergence verdict names what failed and what to do", {
  pass <- make_convergence(c(1.001, 1.004, 1.002), c(900, 800, 700))
  expect_match(format(pass), "^Convergence: All 3 monitored quantities have Rhat < 1.01 and ESS >= 400")
  expect_match(format(pass), "not proof")

  rhat_fail <- make_convergence(c(1.001, 1.30, 1.002), c(900, 800, 700))
  expect_match(format(rhat_fail), "1 of 3 quantities have Rhat >= 1.01 (worst: q2, 1.30)", fixed = TRUE)

  ess_fail <- make_convergence(c(1.001, 1.002, 1.002), c(900, 150, 700))
  expect_match(format(ess_fail), "1 further quantity has Rhat < 1.01 but ESS < 400", fixed = TRUE)
  expect_false(grepl("disagree", format(ess_fail)))

  single <- make_convergence(c(1.001, 1.002), c(900, 800), n_chains = 1)
  expect_match(format(single), "Only one chain was run")
})

test_that("print.mdir_convergence() lists quantities in aligned columns and truncates long tables", {
  conv <- make_convergence(c(1.001, 1.30, 1.002), c(900, 150, 700))
  out <- capture.output(print(conv))
  expect_true(any(grepl("^quantity +Rhat +ESS bulk +ESS tail", out)))
  expect_equal(sum(grepl("\\*$", out[grepl("^q[0-9]", out)])), 1)
  capture.output(res <- print(conv))
  expect_identical(res, conv)

  many <- make_convergence(seq(1.001, 1.5, length.out = 30), rep(1000, 30))
  out_many <- capture.output(print(many, max_rows = 10))
  expect_equal(sum(grepl("^q[0-9]", out_many)), 10)
  expect_true(any(grepl("10 of 30 quantities shown", out_many)))
  # the ten shown are the ten with the largest Rhat
  expect_true(all(paste0("q", 21:30) %in% sub(" .*", "", out_many[grepl("^q[0-9]", out_many)])))
})
