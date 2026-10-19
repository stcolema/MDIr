# Guards against underflow of the component weights. With a mass prior that puts almost
# all of its probability near zero the prior weights Gamma(mass / K, rate) are astronomically
# small, the normalising constant is a few times the smallest double, the strategic latent
# variable overflows, and the first conditional used to fail with
# "randg(): incorrect distribution parameters" (about a third of the runs below).

underflow_data <- function(N = 20) {
  set.seed(99)
  matrix(c(rnorm(N / 2, 0, 1), rnorm(N / 2, 8, 1)), ncol = 1)
}

test_that("a mass prior concentrated near zero still gives a representable start", {
  X <- list(underflow_data())
  prior <- mdiPrior(mass_shape = 1, mass_rate = 1e3)
  for (seed in 1:20) {
    set.seed(seed)
    fit <- callMDI(X, n_iter = 10, thin = 1, types = "MVN", K_max = 5, prior = prior, check_prior = FALSE)
    expect_true(all(is.finite(fit$mass)))
    expect_true(all(fit$allocations >= 0 & fit$allocations < 5))
  }
})

test_that("smcMDI draws representable starts from a sparse mass prior", {
  X <- list(underflow_data())
  set.seed(3)
  fit <- smcMDI(X, "MVN", K_max = 5, n_particles = 20, prior = mdiPrior(mass_shape = 1, mass_rate = 1e3), check_prior = FALSE)
  expect_true(all(is.finite(fit$particle_weights)))
  expect_equal(sum(fit$particle_weights), 1)
})

test_that("a mass prior that cannot give a representable start stops with an explanation", {
  X <- list(underflow_data())
  set.seed(1)
  expect_error(
    callMDI(X, n_iter = 10, thin = 1, types = "MVN", K_max = 5,
            prior = mdiPrior(mass_shape = 1, mass_rate = 1e9), check_prior = FALSE),
    "numerically representable start"
  )
})

test_that("runs with the default and with a moderately sparse mass prior are unaffected", {
  X <- list(underflow_data())
  for (prior in list(mdiPrior(), mdiPrior(mass_shape = 1, mass_rate = 10))) {
    a <- { set.seed(5); callMDI(X, n_iter = 20, thin = 1, types = "MVN", K_max = 5, prior = prior, check_prior = FALSE) }
    b <- { set.seed(5); callMDI(X, n_iter = 20, thin = 1, types = "MVN", K_max = 5, prior = prior, check_prior = FALSE) }
    expect_identical(a$allocations, b$allocations)
  }
})
