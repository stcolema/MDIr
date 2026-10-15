# Parallel chains: opt-in, reproducible streams, capped under R CMD check.

small_views <- function(N = 30, L = 2) {
  lapply(seq_len(L), function(v) {
    m <- matrix(rnorm(N * 2, rep(c(0, 3), each = N / 2)), N, 2)
    rownames(m) <- seq_len(N)
    m
  })
}

with_env <- function(name, value, code) {
  old <- Sys.getenv(name, unset = NA)
  on.exit(if (is.na(old)) Sys.unsetenv(name) else do.call(Sys.setenv, setNames(list(old), name)))
  do.call(Sys.setenv, setNames(list(value), name))
  force(code)
}

test_that("the number of cores defaults to one and is capped by chains, cores and the check limit", {
  old <- options(mdir.cores = NULL)
  on.exit(options(old))
  with_env("_R_CHECK_LIMIT_CORES_", "", {
    expect_equal(mdir:::.mdirResolveCores(NULL, 4), 1L)
    expect_equal(mdir:::.mdirResolveCores(1, 4), 1L)
    # never more cores than chains
    expect_equal(mdir:::.mdirResolveCores(64, 1), 1L)
    expect_lte(mdir:::.mdirResolveCores(64, 8), 8L)
    options(mdir.cores = 2)
    expect_lte(mdir:::.mdirResolveCores(NULL, 4), 2L)
  })
  # as data.table does under R CMD check: at most two cores
  with_env("_R_CHECK_LIMIT_CORES_", "TRUE", {
    expect_lte(mdir:::.mdirResolveCores(64, 8), 2L)
  })
  for (bad in list(0, -1, NA, "a", c(1, 2))) {
    expect_error(mdir:::.mdirResolveCores(bad, 4), "n_cores")
  }
})

test_that("chain streams are reproducible from set.seed and leave the user's generator alone", {
  set.seed(10)
  a <- mdir:::.mdirStreamSeeds(3)
  set.seed(10)
  b <- mdir:::.mdirStreamSeeds(3)
  expect_identical(a, b)
  expect_length(unique(lapply(a, paste, collapse = " ")), 3)
  expect_true(all(vapply(a, function(s) s[1] %% 100 == 7, logical(1))))   # L'Ecuyer-CMRG kind code

  # The generator kind is unchanged, and its stream continues as if one draw had been taken
  set.seed(11)
  kind <- RNGkind()
  mdir:::.mdirStreamSeeds(3)
  expect_identical(RNGkind(), kind)
  after <- runif(1)
  set.seed(11)
  sample.int(.Machine$integer.max, 1L)
  expect_equal(after, runif(1))
})

test_that("serial runs (the default) are unchanged by the cores machinery", {
  X <- small_views()
  set.seed(12)
  a <- runMCMCChains(X, 2, R = 20, thin = 5, types = c("G", "G"), K = c(3, 3))
  set.seed(12)
  b <- lapply(1:2, function(i) callMDI(X, 20, 5, c("G", "G"), K = c(3, 3), check_prior = FALSE))
  for (i in 1:2) {
    expect_identical(a[[i]]$phis, b[[i]]$phis)
    expect_identical(a[[i]]$allocations, b[[i]]$allocations)
    expect_equal(a[[i]]$Chain, i)
  }
})

test_that("parallel chains are reproducible, differ from each other, and do not depend on the worker type", {
  skip_on_cran()
  X <- small_views()
  run <- function(type, n_cores = 2) {
    old <- options(mdir.parallel_type = type)
    on.exit(options(old))
    set.seed(13)
    runMCMCChains(X, 3, R = 20, thin = 5, types = c("G", "G"), K = c(3, 3), n_cores = n_cores)
  }
  if (.Platform$OS.type != "windows") {
    f1 <- run("FORK")
    f2 <- run("FORK")
    expect_s3_class(f1, "mdir_fit_list")
    expect_length(f1, 3)
    expect_equal(vapply(f1, function(ch) ch$Chain, numeric(1)), 1:3)
    for (i in 1:3) {
      expect_identical(f1[[i]]$phis, f2[[i]]$phis)
      expect_identical(f1[[i]]$allocations, f2[[i]]$allocations)
    }
    expect_false(identical(f1[[1]]$phis, f1[[2]]$phis))
  }

  # Sockets need the installed package in the workers
  skip_if(nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_", "")) && .Platform$OS.type == "windows")
  sock <- run("PSOCK")
  expect_length(sock, 3)
  if (.Platform$OS.type != "windows") {
    for (i in 1:3) expect_identical(sock[[i]]$phis, f1[[i]]$phis)
  }
})

test_that("a failing chain is reported", {
  skip_on_cran()
  skip_on_os("windows")
  X <- small_views()
  # fixed entries with labels beyond K make every chain fail inside callMDI
  expect_error(
    runMCMCChains(X, 2, R = 20, thin = 5, types = c("G", "G"), K = c(3, 3), n_cores = 2,
                  initial_labels = matrix(9, 30, 2), fixed = matrix(1, 30, 2),
                  initial_labels_as_intended = TRUE),
    "failed"
  )
})

test_that("fitMDI passes n_cores through", {
  skip_on_cran()
  skip_on_os("windows")
  X <- small_views()
  set.seed(14)
  fit <- fitMDI(X, 2, R = 60, thin = 2, types = c("G", "G"), K = c(3, 3), n_cores = 2, verbose = FALSE)
  expect_length(fit, 2)
  expect_false(is.null(attr(fit, "convergence")))
})
