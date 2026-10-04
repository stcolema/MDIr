# Parallel tempering. The statistical checks compare the sampler with exact
# references computed independently of it; see also verification/tempering/.

tp_data_G <- function(N = 30, L = 2, seed = 1) {
  set.seed(seed)
  lapply(seq_len(L), function(v) {
    m <- matrix(rnorm(N * 2, rep(c(0, 3), each = N / 2)), N, 2)
    rownames(m) <- seq_len(N)
    m
  })
}

test_that("the ladder helpers validate their input", {
  expect_equal(ptLadder(3, beta_min = 0.25, spacing = "linear"), c(0.25, 0.625, 1))
  g <- ptLadder(5, 0.01)
  expect_equal(g[c(1, 5)], c(0.01, 1))
  expect_equal(diff(log(g)), rep(diff(log(g))[1], 4))
  expect_error(ptLadder(1), "at least 2")
  expect_error(ptLadder(4, beta_min = 0), "geometric")
  expect_error(ptLadder(4, beta_min = 1), "beta_min")
  expect_error(mdir:::.mdirCheckLadder(c(0.5, 0.4, 1)), "strictly increasing")
  expect_error(mdir:::.mdirCheckLadder(c(0.2, 0.5)), "must be 1")
  expect_error(mdir:::.mdirCheckLadder(c(0.2, 1.2)), "\\[0, 1\\]")
  expect_equal(mdir:::.mdirCheckLadder(1), 1)
})

test_that("a ladder is refused where the tempered conditionals are not implemented", {
  X <- tp_data_G(20)
  expect_error(callMDI(X, R = 10, thin = 1, types = c("G", "G"), K = c(3, 3), betas = c(0.5, 1),
    check_prior = FALSE), NA)
  Xm <- X
  Xm[[1]][2, 1] <- NA
  expect_error(callMDI(Xm, R = 10, thin = 1, types = c("G", "G"), K = c(3, 3), betas = c(0.5, 1),
    check_prior = FALSE), "complete data")
  expect_error(callMDI(X, R = 10, thin = 1, types = c("TAGM", "G"), K = c(3, 3), betas = c(0.5, 1),
    check_prior = FALSE), "outlier")
  # but a missing value is fine without tempering
  expect_error(callMDI(Xm, R = 10, thin = 1, types = c("G", "G"), K = c(3, 3), check_prior = FALSE), NA)
})

test_that("betas = 1 is the untempered sampler, draw for draw", {
  X <- tp_data_G(30)
  set.seed(3)
  a <- callMDI(X, R = 60, thin = 2, types = c("G", "MVN"), K = c(3, 3), check_prior = FALSE)
  set.seed(3)
  b <- callMDI(X, R = 60, thin = 2, types = c("G", "MVN"), K = c(3, 3), check_prior = FALSE, betas = 1)
  for (nm in c("allocations", "phis", "weights", "mass", "joint_likelihood", "parameters")) {
    expect_identical(a[[nm]], b[[nm]])
  }
  expect_null(a$tempering)
})

test_that("a tempered component follows the tempered conjugate posterior", {
  skip_on_cran()
  one <- function(X, type, beta, R = 12000) {
    rownames(X) <- seq_len(nrow(X))
    callMDI(list(X), R = R, thin = 1, types = type, K = 1, initial_labels = matrix(1, nrow(X), 1),
      density_prior = densityPrior(scale_pool_shape = 0), betas = beta, check_prior = FALSE)$parameters[[1]][-(1:300), ]
  }
  set.seed(21)
  N <- 25
  X <- cbind(rnorm(N, 2, 1.5), rnorm(N, -1, 0.5))
  beta <- 0.35
  hp <- mdir:::densityHyperparameters(X, 1, 0, numeric(0))
  th <- one(X, "G", beta)
  n_eff <- beta * N; kn <- hp$kappa + n_eff; nun <- hp$nu + n_eff
  xb <- colMeans(X); mun <- (hp$kappa * hp$xi + n_eff * xb) / kn
  sc <- hp$scale + beta * colSums(sweep(X, 2, xb)^2) + n_eff * hp$kappa / kn * (xb - hp$xi)^2
  for (p in 1:2) {
    ev <- sc[p] / (nun - 2)
    expect_equal(mean(th[, 2 + p]), ev, tolerance = 0.06)
    expect_equal(var(th[, p]), ev / kn, tolerance = 0.12)
    expect_lt(abs(mean(th[, p]) - mun[p]), 0.1 * sqrt(ev / kn) * 3)
  }
  # untempered variance of mu would be ev / (kappa + N): the tempered one is larger
  expect_gt(var(th[, 1]), 1.5 * sc[1] / (nun - 2) / (hp$kappa + N))

  N <- 30; P <- 3
  Sig <- matrix(c(1, .5, .2, .5, 2, -.3, .2, -.3, .7), P)
  Xm <- matrix(rnorm(N * P), N) %*% chol(Sig) + rep(c(1, 0, -1), each = N)
  beta <- 0.4
  hp <- mdir:::densityHyperparameters(Xm, 1, 1, numeric(0))
  th <- one(Xm, "MVN", beta)
  n_eff <- beta * N; kn <- hp$kappa + n_eff; nun <- hp$nu + n_eff
  xb <- colMeans(Xm)
  sc <- hp$scale + beta * crossprod(sweep(Xm, 2, xb)) + n_eff * hp$kappa / kn * tcrossprod(xb - hp$xi)
  expect_equal(matrix(colMeans(th[, P + seq_len(P * P)]), P), sc / (nun - P - 1), tolerance = 0.06)

  Xc <- cbind(sample(0:2, 40, TRUE, c(.6, .3, .1)), sample(0:1, 40, TRUE))
  beta <- 0.3
  hp <- mdir:::densityHyperparameters(Xc, 1, 2, numeric(0))
  th <- one(Xc, "C", beta)
  start <- 0
  for (p in 1:2) {
    a <- as.numeric(hp$concentration[[p]]) + beta * tabulate(Xc[, p] + 1, hp$n_cat[p])
    cols <- start + seq_len(hp$n_cat[p])
    expect_equal(unname(colMeans(th[, cols, drop = FALSE])), a / sum(a), tolerance = 0.03)
    expect_equal(unname(apply(th[, cols, drop = FALSE], 2, sd)), sqrt(a / sum(a) * (1 - a / sum(a)) / (sum(a) + 1)), tolerance = 0.06)
    start <- start + hp$n_cat[p]
  }
})

test_that("the exchange log-likelihood equals an independent recomputation from the saved draws", {
  skip_on_cran()
  # The cold replica's saved allocation and parameters at each saved draw must
  # give the same data log-likelihood that the exchange step used.
  set.seed(5)
  N <- 24
  X <- tp_data_G(N, L = 2, seed = 5)
  fit <- callMDI(X, R = 30, thin = 3, types = c("G", "MVN"), K = c(3, 3), betas = c(0.3, 0.6, 1),
    check_prior = FALSE)
  ell_cpp <- fit$tempering$data_log_likelihood[, 3]
  ell_r <- vapply(seq_along(ell_cpp), function(s) {
    tot <- 0
    c1 <- fit$allocations[s, , 1] + 1
    th <- fit$parameters[[1]][s, ]
    mu <- matrix(th[1:6], 2, 3); vr <- matrix(th[7:12], 2, 3)
    for (n in 1:N) tot <- tot + sum(dnorm(X[[1]][n, ], mu[, c1[n]], sqrt(vr[, c1[n]]), log = TRUE))
    c2 <- fit$allocations[s, , 2] + 1
    th <- fit$parameters[[2]][s, ]
    mu <- matrix(th[1:6], 2, 3); cv <- array(th[6 + 1:12], c(2, 2, 3))
    for (n in 1:N) {
      S <- cv[, , c2[n]]
      d <- X[[2]][n, ] - mu[, c2[n]]
      tot <- tot - 0.5 * (2 * log(2 * pi) + log(det(S)) + drop(t(d) %*% solve(S, d)))
    }
    tot
  }, numeric(1))
  expect_equal(ell_cpp, ell_r, tolerance = 1e-8)

  set.seed(6)
  Xc <- list({ m <- matrix(sample(0:2, 24 * 2, TRUE), 24, 2); rownames(m) <- 1:24; m })
  fit <- callMDI(Xc, R = 20, thin = 2, types = "C", K = 3, betas = c(0.5, 1), check_prior = FALSE)
  ell_cpp <- fit$tempering$data_log_likelihood[, 2]
  ell_r <- vapply(seq_along(ell_cpp), function(s) {
    th <- fit$parameters[[1]][s, ]
    cl <- fit$allocations[s, , 1] + 1
    tot <- 0
    start <- 0
    for (p in 1:2) {
      ncat <- max(Xc[[1]][, p]) + 1
      pr <- matrix(th[start + seq_len(ncat * 3)], ncat, 3)
      tot <- tot + sum(log(pr[cbind(Xc[[1]][, p] + 1, cl)]))
      start <- start + ncat * 3
    }
    tot
  }, numeric(1))
  expect_equal(ell_cpp, ell_r, tolerance = 1e-8)
})

test_that("exchange bookkeeping is consistent", {
  X <- tp_data_G(30)
  betas <- c(0.2, 0.5, 0.8, 1)
  for (scheme in c("deo", "seo")) {
    fit <- callMDI(X, R = 200, thin = 5, types = c("G", "G"), K = c(3, 3), betas = betas,
      swap_scheme = scheme, check_prior = FALSE)
    tp <- fit$tempering
    expect_s3_class(tp, "mdir_tempering")
    expect_equal(tp$betas, betas)
    expect_length(tp$swap_attempts, 3)
    expect_true(all(tp$swap_accepts <= tp$swap_attempts))
    expect_true(all(tp$rejection_rate >= 0 & tp$rejection_rate <= 1))
    # rejection rate is one minus the mean acceptance probability
    expect_equal(tp$rejection_rate, 1 - tp$swap_accept_prob_sum / tp$swap_attempts)
    if (scheme == "deo") {
      # deterministic alternation: each round attempts the even pairs (1, 3) or the odd pair (2)
      expect_equal(tp$swap_attempts[1], tp$swap_attempts[3])
      expect_equal(tp$swap_attempts[1] + tp$swap_attempts[2], 200)
    } else {
      # random parity: a round attempts the pairs (1, 3) or the pair (2)
      expect_equal(tp$swap_attempts[1], tp$swap_attempts[3])
      expect_equal(tp$swap_attempts[1] + tp$swap_attempts[2], 200)
    }
    # every replica sits at exactly one temperature at each saved draw
    expect_true(all(apply(tp$replica, 1, function(r) all(sort(r) == 1:4))))
    d <- ptDiagnostics(fit)
    expect_s3_class(d, "mdir_pt_diagnostics")
    expect_equal(nrow(d$pairs), 3)
    expect_equal(d$barrier, sum(d$pairs$rejection_rate))
    expect_output(print(d), "Exchange diagnostics only")
  }
  # no ladder, no record
  expect_error(ptDiagnostics(callMDI(X, R = 20, thin = 2, types = c("G", "G"), K = c(3, 3), check_prior = FALSE)),
    "no tempering record")
})

test_that("the ladder update equalises the rejection rates (Gaussian path with known rates)", {
  skip_on_cran()
  set.seed(8)
  b0 <- ptLadder(7, 0.002)
  r0 <- vapply(seq_len(6), function(i) tp_gauss_rejection(b0[i], b0[i + 1]), numeric(1))
  expect_gt(max(r0) / min(r0), 3)             # the starting ladder is far from equal
  b1 <- tuneLadder(b0, r0)
  expect_equal(b1[c(1, 7)], c(0.002, 1))
  expect_true(all(diff(b1) > 0))
  r1 <- vapply(seq_len(6), function(i) tp_gauss_rejection(b1[i], b1[i + 1]), numeric(1))
  expect_lt(max(r1) / min(r1), 1.35)
  b2 <- tuneLadder(b1, r1)
  r2 <- vapply(seq_len(6), function(i) tp_gauss_rejection(b2[i], b2[i + 1]), numeric(1))
  expect_lt(max(r2) / min(r2), 1.2)
  # the sum of rejection rates (the barrier) is nearly ladder independent
  expect_equal(sum(r2), sum(r0), tolerance = 0.12)
  expect_error(tuneLadder(b0, r0[-1]), "neighbouring pair")
})

test_that("parallel tempering reproduces the exact label posterior of a tiny model", {
  skip_on_cran()
  set.seed(31)
  N <- 5; K <- 2
  X <- matrix(sample(0:2, N, TRUE, c(.5, .3, .2)), N, 1)
  rownames(X) <- seq_len(N)
  alpha <- mdir:::densityHyperparameters(X, K, 2, numeric(0))$concentration
  code <- function(a) as.numeric(a %*% K^(seq_len(N) - 1))
  run <- function(betas, seed, R = 30000) {
    set.seed(seed)
    a <- callMDI(list(X), R = R, thin = 1, types = "C", K = K, betas = betas, check_prior = FALSE,
      save_parameters = FALSE)$allocations[-(1:200), , 1]
    tabulate(code(a) + 1, K^N) / nrow(a)
  }
  ex1 <- tp_exact_L1(X, K, alpha, 1)
  freq_pt <- rowMeans(vapply(1:6, function(i) run(c(0.2, 0.5, 1), 100 + i), numeric(K^N)))
  tv <- 0.5 * sum(abs(freq_pt - ex1$prob))
  expect_lt(tv, 0.03)
  # the check has power: the same chain is far from the posterior at another temperature
  ex_wrong <- tp_exact_L1(X, K, alpha, 0.5)
  expect_gt(0.5 * sum(abs(freq_pt - ex_wrong$prob)), 0.06)
  # and a single tempered chain targets its own tempered posterior
  freq_t <- rowMeans(vapply(1:6, function(i) run(0.5, 200 + i), numeric(K^N)))
  expect_lt(0.5 * sum(abs(freq_t - ex_wrong$prob)), 0.03)
})
