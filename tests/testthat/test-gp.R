# The guarded, partially pooled GP priors (see ?densityPrior).

test_that("the inverse-gamma calibration puts the requested mass in each tail", {
  for (bounds in list(c(1, 2), c(1, 4), c(1, 19), c(2, 30))) {
    ab <- mdir:::calibrateInverseGammaCpp(bounds[1], bounds[2], 0.01)
    lower <- pgamma(1 / bounds[1], ab[1], rate = ab[2], lower.tail = FALSE)   # P(lambda < lower)
    upper <- pgamma(1 / bounds[2], ab[1], rate = ab[2])                        # P(lambda > upper)
    expect_equal(lower, 0.01, tolerance = 1e-6)
    expect_equal(upper, 0.01, tolerance = 1e-6)
  }
})

test_that("the population update matches the analytic posterior of (m, s)", {
  skip_on_cran()
  set.seed(91)
  y <- c(0.3, 1.1, -0.4, 0.8); m0 <- 0.2; t <- 2; h <- 1
  n <- length(y)
  chain <- mdir:::gpPopulationCheckCpp(y, m0, t, h, 60000)[-(1:200), ]
  # p(s | y) is proportional to the half-normal prior times the density of y with m integrated out
  logm <- function(s) {
    -(n - 1) * log(s) - 0.5 * log(s^2 + n * t^2) -
      0.5 * (sum((y - mean(y))^2) / s^2 + n * (mean(y) - m0)^2 / (s^2 + n * t^2)) - 0.5 * (s / h)^2
  }
  grid <- seq(1e-3, 8, length.out = 20000)
  w <- exp(logm(grid) - max(logm(grid))); w <- w / sum(w)
  post_m <- function(s) { prec <- 1 / t^2 + n / s^2; (m0 / t^2 + sum(y) / s^2) / prec }
  exp_s <- sum(w * grid); exp_m <- sum(w * post_m(grid)); sd_s <- sqrt(sum(w * grid^2) - exp_s^2)
  expect_lt(abs(mean(chain[, 2]) - exp_s), 5 * sd_s / sqrt(ess_of(chain[, 2])) + 0.01)
  expect_lt(abs(mean(chain[, 1]) - exp_m), 5 * sd(chain[, 1]) / sqrt(ess_of(chain[, 1])) + 0.01)
  expect_equal(sd(chain[, 2]), sd_s, tolerance = 0.08)
})

gp_log_marginal <- function(X, a, lambda, s2, jitter = 1e-8) {
  # X centred on its column means, so xbar = 0
  n <- nrow(X); P <- ncol(X)
  SS <- sum(sweep(X, 2, colMeans(X))^2)
  D <- outer(seq_len(P), seq_len(P), function(i, j) -0.5 * (i - j)^2)
  Kmat <- a * exp(D / lambda^2) + diag(jitter * a, P)
  ch <- chol(Kmat + diag(s2 / n, P))
  -0.5 * n * P * log(2 * pi * s2) + 0.5 * P * log(2 * pi * s2 / n) - SS / (2 * s2) -
    0.5 * (P * log(2 * pi) + 2 * sum(log(diag(ch))))
}

test_that("GP hyperparameter sampler targets the marginal posterior under the guarded, pooled priors", {
  skip_on_cran()
  set.seed(92)
  P <- 5; n <- 14
  Kt <- outer(1:P, 1:P, function(i, j) 0.8 * exp(-0.5 * (i - j)^2 / 2^2))
  mu <- as.numeric(t(chol(Kt)) %*% rnorm(P))
  X <- matrix(rnorm(n * P, mu, 0.5), n, P, byrow = TRUE); rownames(X) <- seq_len(n)

  hp <- mdir:::densityHyperparameters(X, 1, 3, numeric(0))
  c0 <- hp$log_variance_centre; t <- hp$center_sd; h <- hp$pool_sd_scale
  alpha <- hp$length_prior_shape; beta <- hp$length_prior_rate
  lo <- hp$min_length; hi <- hp$max_length

  # density of the log amplitude / log noise with the population (m, s) integrated out
  g_pop <- function(y) sapply(y, function(yy) {
    integrate(function(s) dnorm(yy, c0, sqrt(t^2 + s^2)) * 2 * dnorm(s, 0, h), 0, Inf, rel.tol = 1e-8)$value
  })
  # amplitude: wide grid (weakly identified); noise: fine grid around the data value
  s2hat <- sum(sweep(X, 2, colMeans(X))^2) / ((n - 1) * P)
  ya <- c0 + seq(-12, 5, length.out = 60)
  yn <- log(s2hat) + seq(-1.5, 1.5, length.out = 41)
  z <- seq(log(lo), log(hi), length.out = 30)
  ga <- log(g_pop(ya)); gn <- log(g_pop(yn))
  eg <- expand.grid(ia = seq_along(ya), in_ = seq_along(yn), iz = seq_along(z))
  lp <- mapply(function(ia, in_, iz) {
    gp_log_marginal(X, exp(ya[ia]), exp(z[iz]), exp(yn[in_])) + ga[ia] + gn[in_] +
      (-alpha * z[iz] - beta / exp(z[iz]))
  }, eg$ia, eg$in_, eg$iz)
  gm <- list(
    amplitude = grid_moments(ya[eg$ia], lp), noise = grid_moments(yn[eg$in_], lp),
    length = grid_moments(z[eg$iz], lp)
  )

  # one occupied component: the population (m, s) is informed by a single log amplitude, so the
  # amplitude mixes slowly and a long run is needed
  out <- callMDI(list(X), R = 400000, thin = 10, types = "GP", K = 1)
  h_out <- out$hypers[[1]]
  keep <- -(1:200)
  draws <- list(amplitude = log(h_out$amplitude[keep, 1]), noise = log(h_out$noise[keep, 1]),
                length = log(h_out$length[keep, 1]))
  for (nm in names(draws)) {
    d <- draws[[nm]]
    expect_lt(abs(mean(d) - gm[[nm]]["mean"]), 5 * gm[[nm]]["sd"] / sqrt(ess_of(d)) + 0.05 * gm[[nm]]["sd"],
              label = paste(nm, "mean"))
    expect_equal(sd(d), unname(gm[[nm]]["sd"]), tolerance = 0.15, label = paste(nm, "sd"))
  }
  # the guard holds
  expect_gte(min(h_out$length[keep, 1]), lo)
})

test_that("the length scale never falls below the floor, even for white noise data", {
  skip_on_cran()
  set.seed(93)
  X <- matrix(rnorm(60 * 8), 60); rownames(X) <- seq_len(60)
  for (floor_value in c(1, 2.5)) {
    fit <- callMDI(list(X), R = 3000, thin = 3, types = "GP", K = 3,
                   density_prior = densityPrior(gp_min_length = floor_value))
    len <- fit$hypers[[1]]$length
    expect_gte(min(len), floor_value)
    # and the posterior is not simply stuck at the floor
    expect_gt(median(len), 1.05 * floor_value)
  }
})

test_that("GP data need not be standardised and are centred at the column means", {
  skip_on_cran()
  set.seed(94)
  P <- 6; N <- 60
  truth <- rep(1:2, each = N / 2)
  f <- cbind(sin(seq_len(P) / 2), cos(seq_len(P) / 2))
  X <- 500 + 40 * (f[, truth] |> t()) + matrix(rnorm(N * P, sd = 5), N)
  rownames(X) <- seq_len(N)
  fit <- callMDI(list(X), R = 1500, thin = 3, types = "GP", K = 2)
  pred <- apply(fit$allocations[-(1:100), , 1], 2, function(z) as.integer(names(which.max(table(z)))))
  tab <- table(pred, truth)
  expect_gte(sum(apply(tab, 2, max)) / N, 0.95)
})
