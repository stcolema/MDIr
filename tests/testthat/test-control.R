# mdirControl(): the proposal windows and the auto-tuning schedule in one object, with
# the old individual `proposal_windows` argument deprecated.

gp_views <- function(N = 30, P = 6) {
  X <- matrix(rnorm(N * P, rep(c(0, 2), each = N / 2)), N, P)
  rownames(X) <- seq_len(N)
  list(X)
}

catch_warnings <- function(expr) {
  w <- character()
  value <- withCallingHandlers(expr, warning = function(cnd) {
    w <<- c(w, conditionMessage(cnd))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = w)
}

test_that("mdirControl() has the documented defaults and checks its arguments", {
  control <- mdirControl()
  expect_s3_class(control, "mdir_control")
  expect_null(control$proposal_windows)
  expect_true(control$auto_tune)
  expect_null(control$n_burn)

  ok <- mdirControl(proposal_windows = list(c(0.1, 0.2, 0.3), NULL), auto_tune = FALSE, n_burn = 50)
  expect_equal(ok$proposal_windows[[1]], c(0.1, 0.2, 0.3))
  expect_false(ok$auto_tune)

  expect_error(mdirControl(proposal_windows = c(0.1, 0.2, 0.3)), "proposal_windows")
  expect_error(mdirControl(proposal_windows = list("a")), "proposal_windows")
  expect_error(mdirControl(proposal_windows = list(c(0.1, -1, 0.3))), "proposal_windows")
  expect_error(mdirControl(auto_tune = "yes"), "auto_tune")
  expect_error(mdirControl(auto_tune = c(TRUE, FALSE)), "auto_tune")
  expect_error(mdirControl(n_burn = NA), "n_burn")
  expect_error(mdirControl(n_burn = "a"), "n_burn")
})

test_that("a control object prints", {
  expect_output(print(mdirControl()), "auto_tune = TRUE, n_burn = NULL")
  expect_output(print(mdirControl()), "defaults")
  expect_output(
    print(mdirControl(proposal_windows = list(c(0.05, 0.06, 0.07), NULL), n_burn = 10)),
    "view 1: 0.05, 0.06, 0.07"
  )
  expect_output(print(mdirControl(proposal_windows = list(NULL, c(1, 2, 3)))), "view 1: default")
})

test_that("control carries the windows into the fit, whichever entry point is used", {
  X <- gp_views()
  w <- c(0.2, 0.25, 0.3)
  control <- mdirControl(proposal_windows = list(w), auto_tune = FALSE)

  fit <- callMDI(X, n_iter = 60, thin = 5, types = "GP", K_max = 2, check_prior = FALSE, control = control)
  expect_equal(fit$final_proposal_windows[[1]], w)
  expect_equal(fit$proposal_windows[[1]], w)
  expect_false(fit$auto_tune)

  chains <- runMCMCChains(X, 2, n_iter = 60, thin = 5, types = "GP", K_max = 2, control = control)
  expect_equal(chains[[2]]$final_proposal_windows[[1]], w)
  fitted <- fitMDI(X, 2, n_iter = 200, thin = 2, types = "GP", K_max = 2, control = control, verbose = FALSE)
  expect_equal(fitted[[1]]$final_proposal_windows[[1]], w)
  mixture <- callMixtureModel(X[[1]], n_iter = 60, thin = 5, type = "GP", K_max = 2, control = control)
  expect_equal(mixture$final_proposal_windows[[1]], w)

  dir <- tempfile("mdir_control")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  expect_no_error(callMDIWritingToFile(X, n_iter = 60, thin = 5, types = "GP", K_max = 2, dir_path = dir,
                                       control = control))
  chain_dir <- file.path(dir, "chains")
  dir.create(chain_dir)
  expect_no_error(runMCMCChainsSavingToFile(X, 2, n_iter = 60, thin = 5, types = "GP", K_max = 2,
                                            dir_path = chain_dir, control = control))
})

test_that("a control object that is not an mdirControl() is refused", {
  X <- gp_views()
  expect_error(
    callMDI(X, n_iter = 20, thin = 5, types = "GP", K_max = 2, control = list(auto_tune = FALSE)),
    "mdirControl"
  )
})

test_that("the deprecated proposal_windows argument still works, with a warning, and control wins over it", {
  X <- gp_views()
  w <- c(0.05, 0.06, 0.07)
  run <- function(...) {
    set.seed(8)
    callMDI(X, n_iter = 200, thin = 5, types = "GP", K_max = 2, check_prior = FALSE, ...)
  }
  new <- run(control = mdirControl(proposal_windows = list(w)))
  old <- catch_warnings(run(proposal_windows = list(w)))
  expect_length(old$warnings, 1)
  expect_match(old$warnings, "`proposal_windows` argument of callMDI\\(\\) is deprecated")
  expect_identical(old$value$hypers, new$hypers)

  both <- catch_warnings(run(control = mdirControl(proposal_windows = list(c(0.3, 0.3, 0.3)), auto_tune = FALSE),
                             proposal_windows = list(w)))
  expect_match(both$warnings, "Both `control` and `proposal_windows`")
  expect_equal(both$value$final_proposal_windows[[1]], c(0.3, 0.3, 0.3))

  # the other entry points warn once for the call, not once per chain
  chains <- catch_warnings(runMCMCChains(X, 3, n_iter = 40, thin = 5, types = "GP", K_max = 2,
                                         proposal_windows = list(w)))
  expect_length(chains$warnings, 1)
  expect_match(chains$warnings, "runMCMCChains\\(\\)")
})

test_that("views without hyperparameters need no windows", {
  views <- lapply(1:2, function(v) {
    m <- matrix(rnorm(40 * 2, rep(c(0, 3), each = 20)), 40, 2)
    rownames(m) <- 1:40
    m
  })
  fit <- callMDI(views, n_iter = 40, thin = 5, types = c("MVN", "G"), K_max = c(3, 3),
                 control = mdirControl(auto_tune = TRUE, n_burn = 10))
  expect_equal(fit$final_proposal_windows, list(0, 0))
})
