# The phis are updated from their conditional with the strategic latent variable
# integrated out. Z is linear in phi(l, m), Z = A + B phi, so the conditional is
#   phi^(shape - 1) exp(-rate phi) (1 + phi)^N_lm (A + B phi)^(-N),
# and one slice-sampling update (Neal, 2003) draws from it without tuning.

phi_state <- function(L, K = rep(4L, L), seed = 1) {
  set.seed(seed)
  w <- matrix(0, max(K), L)
  for (l in seq_len(L)) w[seq_len(K[l]), l] <- rgamma(K[l], 1, 2)
  phi <- matrix(0, L, L)
  phi[upper.tri(phi)] <- rgamma(L * (L - 1) / 2, 2, 0.4)
  list(w = w, K = as.integer(K), phi = phi + t(phi))
}

test_that("Z is linear in a single phi, with the coefficient used by the update", {
  for (L in 3:5) {
    st <- phi_state(L, K = c(3L, 4L, 2L, 3L, 4L)[seq_len(L)], seed = L)
    pairs <- t(combn(L, 2))
    for (i in seq_len(nrow(pairs))) {
      l <- pairs[i, 1]; m <- pairs[i, 2]
      B <- mdir:::mdiPhiRateCpp(st$w, st$K, st$phi, l - 1, m - 1)
      zero <- st$phi; zero[l, m] <- zero[m, l] <- 0
      A <- mdir:::mdiNormalisingConstantCpp(st$w, st$K, zero)
      for (value in c(0.3, st$phi[l, m], 25)) {
        moved <- st$phi; moved[l, m] <- moved[m, l] <- value
        expect_equal(mdir:::mdiNormalisingConstantCpp(st$w, st$K, moved), A + B * value, tolerance = 1e-10)
      }
    }
  }
  # and against enumeration
  st <- phi_state(3, K = c(3L, 2L, 4L))
  expect_equal(mdir:::mdiNormalisingConstantCpp(st$w, st$K, st$phi), bf_Z(st$w, st$K, st$phi), tolerance = 1e-10)
})

test_that("the slice update draws from the exact collapsed conditional (L = 3 and 5, different scales)", {
  skip_on_cran()
  qs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  for (case in list(list(L = 3, N = 300, N_lm = 210, pair = c(1, 2)),
                    list(L = 5, N = 300, N_lm = 210, pair = c(2, 5)),
                    list(L = 4, N = 40, N_lm = 5, pair = c(1, 4)))) {
    st <- phi_state(case$L, seed = 100 + case$L)
    l <- case$pair[1]; m <- case$pair[2]
    draws <- mdir:::phiSliceChainCpp(st$w, st$K, st$phi, l - 1, m - 1, case$N, case$N_lm, 2, 0.2, 100000)
    grid <- seq(1e-5, max(draws) * 1.5, length.out = 200000)
    # The density written out in R from A and B, which are checked against
    # enumeration above (Z is linear in the phi), independently of the C++ code
    B <- mdir:::mdiPhiRateCpp(st$w, st$K, st$phi, l - 1, m - 1)
    zero <- st$phi; zero[l, m] <- zero[m, l] <- 0
    A <- mdir:::mdiNormalisingConstantCpp(st$w, st$K, zero)
    ld <- (2 - 1) * log(grid) - 0.2 * grid + case$N_lm * log1p(grid) - case$N * log(A + B * grid)
    expect_equal(ld[c(10, 1000, 100000)],
                 as.numeric(mdir:::phiConditionalLogDensityCpp(grid[c(10, 1000, 100000)], st$w, st$K, st$phi, l - 1, m - 1,
                                                               case$N, case$N_lm, 2, 0.2)),
                 tolerance = 1e-9)
    p <- exp(ld - max(ld)); p <- p / sum(p)
    exact_cdf <- cumsum(p)
    # the empirical CDF at the exact quantiles: draws are close to independent
    # (see below), so the standard error of each is about sqrt(q (1 - q) / 1e5)
    exact_q <- vapply(qs, function(q) grid[which(exact_cdf >= q)[1]], numeric(1))
    ecdf_at <- vapply(exact_q, function(x) mean(draws <= x), numeric(1))
    expect_lt(max(abs(ecdf_at - qs)), 0.01, label = paste("L =", case$L))
    expect_lt(abs(mean(draws) - sum(grid * p)), 5 * sd(draws) / sqrt(1e5 / 4))
  }
})

test_that("the slice update needs no tuning: draws are nearly independent at very different scales", {
  skip_on_cran()
  for (case in list(list(L = 3, N = 300, N_lm = 210), list(L = 5, N = 300, N_lm = 210), list(L = 3, N = 3000, N_lm = 2900))) {
    st <- phi_state(case$L, seed = 7)
    draws <- mdir:::phiSliceChainCpp(st$w, st$K, st$phi, 0, 1, case$N, case$N_lm, 2, 0.2, 20000)
    expect_lt(abs(acf(draws, plot = FALSE, lag.max = 1)$acf[2]), 0.1)
  }
})

test_that("the collapsed conditional is the marginal of the joint with v (Gibbs and slice agree for L = 2, fixed weights)", {
  skip_on_cran()
  set.seed(30)
  K <- c(5L, 5L); w <- matrix(rgamma(10, 1, 2), 5, 2)
  N <- 200; N_lm <- 150; shape <- 2; rate <- 0.2
  S1 <- sum(w[, 1]); S2 <- sum(w[, 2]); c1 <- sum(w[, 1] * w[, 2])
  # Gibbs with v, as the package did before: v | phi ~ Gamma(N, Z); phi | v is a mixture of Gammas
  phi <- 5; gibbs <- numeric(40000)
  r <- 0:N_lm
  for (s in seq_along(gibbs)) {
    v <- rgamma(1, N, S1 * S2 + phi * c1)
    lw <- lchoose(N_lm, r) + lgamma(r + shape) - (r + shape) * log(v * c1 + rate)
    pr <- exp(lw - max(lw))
    phi <- rgamma(1, shape + sample(r, 1, prob = pr), rate + v * c1)
    gibbs[s] <- phi
  }
  slice <- mdir:::phiSliceChainCpp(w, K, matrix(c(0, 5, 5, 0), 2), 0, 1, N, N_lm, shape, rate, 40000)
  qs <- c(0.1, 0.5, 0.9)
  expect_equal(unname(quantile(slice, qs)), unname(quantile(gibbs, qs)), tolerance = 0.03)
})

test_that("both phi updates run and record the same quantities; the choice is validated", {
  set.seed(31)
  N <- 30
  X <- lapply(1:3, function(v) { m <- matrix(rnorm(N * 2, rep(c(0, 3), each = N / 2)), N, 2); rownames(m) <- seq_len(N); m })
  for (upd in c("slice", "gibbs")) {
    fit <- callMDI(X, R = 40, thin = 10, types = rep("G", 3), K = rep(3, 3), phi_update = upd)
    expect_equal(fit$phi_update, upd)
    expect_equal(dim(fit$phis), c(5, 3))
    expect_true(all(is.finite(fit$phis)) && all(fit$phis > 0))
  }
  expect_error(callMDI(X, R = 40, thin = 10, types = rep("G", 3), K = rep(3, 3), phi_update = "mh"), "should be one of")
})

# Two samplers with one target: the Gibbs update given v and the collapsed slice
# update must give the same posterior of the phis. L = 3 views of 40 items, a
# small K, many independent chains, so that between-chain variation sets the
# Monte Carlo error.
test_that("slice and Gibbs updates give the same posterior of the phis (L = 3)", {
  skip_on_cran()
  set.seed(40)
  N <- 40; L <- 3; z <- sample(2, N, TRUE)
  X <- lapply(seq_len(L), function(v) {
    zz <- ifelse(runif(N) < 0.75, z, sample(2, N, TRUE))
    m <- matrix(rnorm(N * 2), N, 2) + 2 * zz; rownames(m) <- seq_len(N); m
  })
  summarise <- function(upd, seed) {
    set.seed(seed)
    out <- vapply(seq_len(16), function(i) {
      f <- callMDI(X, R = 12000, thin = 10, types = rep("G", L), K = rep(4, L), save_parameters = FALSE, phi_update = upd)
      k <- 101:1201
      c(colMeans(f$phis[k, ]), log_phi1 = mean(log(f$phis[k, 1])))
    }, numeric(4))
    list(mean = rowMeans(out), se = apply(out, 1, sd) / sqrt(ncol(out)))
  }
  a <- summarise("gibbs", 41)
  b <- summarise("slice", 42)
  z_score <- (a$mean - b$mean) / sqrt(a$se^2 + b$se^2)
  expect_lt(max(abs(z_score)), 4)
})

# With a constant likelihood the stationary distribution is the prior, whatever
# the number of views. Every phi must recover its Gamma(2, 0.2) prior (mean 10,
# sd 7.07), and so must the masses (Gamma(2, 0.1): mean 20, sd 14.1). Four and
# five views, which the sampler's partition recursion handles with the
# collapsed phi update at a higher cost per evaluation.
test_that("the sampler recovers the prior of every phi with four and five views", {
  skip_on_cran()
  ess_of_vec <- function(x) {
    n <- length(x) %/% 2
    rankNormalizedRhat(cbind(x[seq_len(n)], x[n + seq_len(n)]))$ess_bulk
  }
  for (K in list(c(2L, 3L, 2L, 2L), c(2L, 2L, 3L, 2L, 2L))) {
    L <- length(K)
    set.seed(60 + L)
    X <- lapply(seq_len(L), function(l) { m <- matrix(0, 12, 1); rownames(m) <- 1:12; m })
    out <- callMDI(X, R = 150000, thin = 10, types = rep("C", L), K = K)
    keep <- -(1:500)
    for (i in seq_len(ncol(out$phis))) {
      x <- out$phis[keep, i]
      e <- ess_of_vec(x)
      expect_lt(abs(mean(x) - 10), 5 * 7.07 / sqrt(e), label = paste("phi", i, "mean, L =", L))
      expect_lt(abs(sd(x) - 7.07), 0.15 * 7.07 + 5 * 7.07 / sqrt(e), label = paste("phi", i, "sd, L =", L))
    }
    for (l in seq_len(L)) {
      x <- out$mass[keep, l]
      e <- ess_of_vec(x)
      expect_lt(abs(mean(x) - 20), 5 * 14.14 / sqrt(e), label = paste("mass", l, "mean, L =", L))
    }
  }
})
