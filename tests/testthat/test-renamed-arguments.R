# `R` was renamed `n_iter` and `K` was renamed `K_max`; the old names still work
# with a warning. The exports are an explicit list rather than a pattern.

views <- function(N = 30, V = 2) {
  lapply(seq_len(V), function(v) {
    m <- matrix(stats::rnorm(N * 2, rep(c(0, 3), each = N)), N, 2)
    rownames(m) <- seq_len(N)
    m
  })
}

# Run `expr`, returning its value and the messages of every warning it raised
catch_warnings <- function(expr) {
  w <- character()
  value <- withCallingHandlers(
    expr,
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = w)
}

test_that("deprecated `R` and `K` give the same fit as `n_iter` and `K_max`, with warnings", {
  X <- views()
  set.seed(1)
  new <- callMDI(X, n_iter = 20, thin = 5, types = c("MVN", "MVN"), K_max = c(3, 3))
  set.seed(1)
  old <- catch_warnings(callMDI(X, R = 20, thin = 5, types = c("MVN", "MVN"), K = c(3, 3)))
  expect_length(old$warnings, 2)
  expect_match(old$warnings[1], "`R`.*deprecated.*`n_iter`")
  expect_match(old$warnings[2], "`K`.*deprecated.*`K_max`")
  expect_identical(old$value$allocations, new$allocations)
  expect_identical(old$value$K, new$K)
})

test_that("positional calls are unchanged and the fit records n_iter alongside R", {
  X <- views()
  expect_no_warning(fit <- callMDI(X, 20, 5, c("MVN", "MVN"), c(3, 3)))
  expect_identical(fit$n_iter, 20)
  expect_identical(fit$R, fit$n_iter)
})

test_that("a name and its deprecated alias together, or neither, are errors", {
  X <- views()
  expect_error(
    callMDI(X, n_iter = 20, R = 20, thin = 5, types = c("MVN", "MVN")),
    "both given"
  )
  expect_error(
    callMDI(X, n_iter = 20, thin = 5, types = c("MVN", "MVN"), K_max = 3, K = 3),
    "both given"
  )
  expect_error(callMDI(X, thin = 5, types = c("MVN", "MVN")), "n_iter")
})

test_that("every entry point accepts the new names without a warning", {
  X <- views()
  types <- c("MVN", "MVN")
  expect_no_warning(runMCMCChains(X, 2, n_iter = 20, thin = 5, types = types, K_max = c(3, 3)))
  # fitMDI() also warns that 20 iterations are too few for diagnostics; only the
  # deprecation matters here
  w <- catch_warnings(fitMDI(X, 2, n_iter = 20, thin = 5, types = types, K_max = c(3, 3), verbose = FALSE))$warnings
  expect_false(any(grepl("deprecated", w)))
  expect_no_warning(callMixtureModel(X[[1]], n_iter = 20, thin = 5, type = "MVN", K_max = 3))
  expect_no_warning(smcMDI(X, types, K_max = c(3, 3), n_particles = 10))
  expect_no_warning(simulatePriorPredictive(X, types, K_max = c(3, 3), n_datasets = 2))
})

test_that("the deprecated names are still accepted by the other entry points", {
  X <- views()
  types <- c("MVN", "MVN")
  expect_match(
    catch_warnings(runMCMCChains(X, 2, R = 20, thin = 5, types = types, K = c(3, 3)))$warnings[1],
    "deprecated"
  )
  expect_match(
    catch_warnings(callMixtureModel(X[[1]], R = 20, thin = 5, type = "MVN", K = 3))$warnings[1],
    "deprecated"
  )
  expect_match(
    catch_warnings(smcMDI(X, types, K = c(3, 3), n_particles = 10))$warnings[1],
    "`K`.*deprecated"
  )
  expect_match(
    catch_warnings(simulatePriorPredictive(X, types, K = c(3, 3), n_datasets = 2))$warnings[1],
    "`K`.*deprecated"
  )
})

test_that("only the documented interface is exported", {
  exports <- getNamespaceExports("mdir")
  expect_true(all(c("callMDI", "fitMDI", "runMCMCChains", "smcMDI", "predictMDI") %in% exports))
  # raw engines, validators and density helpers are internal
  expect_false(any(c("runMDI", "runMDISMC", "checkTypes", "gammaLogLikelihood", "wishartLogLikelihood") %in% exports))
  expect_false(any(grepl("Cpp$", exports)))
})
