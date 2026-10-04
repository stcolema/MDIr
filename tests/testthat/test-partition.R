test_that("normalising constant and rates match brute-force enumeration", {
  set.seed(101)
  configs <- list(1L, c(4L, 4L), c(3L, 5L), c(5L, 3L), c(3L, 3L, 3L), c(2L, 4L, 3L),
                  c(3L, 2L, 4L, 3L), c(2L, 2L, 3L, 2L, 2L))
  for (K in configs) {
    for (rep in 1:3) {
      s <- random_mdi_state(K)
      L <- length(K)
      expect_equal(mdir:::mdiNormalisingConstantCpp(s$w, K, s$phi), bf_Z(s$w, K, s$phi),
                   tolerance = 1e-10, info = paste("Z, K =", paste(K, collapse = ",")))
      for (l in seq_len(L)) for (k in seq_len(K[l])) {
        expect_equal(mdir:::mdiWeightRateCpp(s$w, K, s$phi, l - 1, k - 1),
                     bf_weight_rate(s$w, K, s$phi, l, k), tolerance = 1e-10)
      }
      if (L > 1) for (l in 1:(L - 1)) for (m in (l + 1):L) {
        expect_equal(mdir:::mdiPhiRateCpp(s$w, K, s$phi, l - 1, m - 1),
                     bf_phi_rate(s$w, K, s$phi, l, m), tolerance = 1e-10,
                     info = paste("phi rate", l, m, "K =", paste(K, collapse = ",")))
      }
    }
  }
})

test_that("normalising constant is accurate when phi is very small or very large", {
  set.seed(102)
  K <- c(3L, 3L, 3L, 3L)
  for (scale in c(1e-8, 1e-4, 1e3)) {
    s <- random_mdi_state(K, phi_scale = scale)
    expect_equal(mdir:::mdiNormalisingConstantCpp(s$w, K, s$phi), bf_Z(s$w, K, s$phi), tolerance = 1e-9)
    expect_equal(mdir:::mdiPhiRateCpp(s$w, K, s$phi, 0, 3), bf_phi_rate(s$w, K, s$phi, 1, 4), tolerance = 1e-9)
  }
})

test_that("phi shape mixture weights are binomial-Gamma marginal likelihoods", {
  # sum_r C(N, r) Gamma(r + a) / (b + rate)^(r + a) is proportional to the
  # marginal of (1 + phi)^N under a Gamma(a, b) prior, up to exp(-rate phi)
  N_lm <- 7; a <- 2; b <- 0.2; rate <- 1.3
  direct <- integrate(function(phi) (1 + phi)^N_lm * dgamma(phi, a, b) * exp(-rate * phi), 0, Inf)$value
  mix <- sum(choose(N_lm, 0:N_lm) * gamma(0:N_lm + a) / (b + rate)^(0:N_lm + a)) * b^a / gamma(a)
  expect_equal(mix, direct, tolerance = 1e-6)
})

test_that("rates of all the weights of a view match brute-force enumeration", {
  set.seed(103)
  # Includes views with a single component, ragged K and up to eight views
  configs <- list(1L, c(1L, 3L), c(3L, 1L, 2L), c(4L, 4L), c(2L, 4L, 3L), c(3L, 2L, 4L, 3L),
                  c(2L, 2L, 2L, 2L, 2L, 2L), c(2L, 3L, 2L, 2L, 3L, 2L, 2L), rep(2L, 8L))
  for (K in configs) {
    s <- random_mdi_state(K)
    for (l in seq_along(K)) {
      got <- mdir:::mdiWeightRatesCpp(s$w, K, s$phi, l - 1)
      expect_length(got, K[l])
      expect_equal(got, vapply(seq_len(K[l]), function(k) bf_weight_rate(s$w, K, s$phi, l, k), numeric(1)),
                   tolerance = 1e-10, info = paste("view", l, "K =", paste(K, collapse = ",")))
    }
  }
})

test_that("weight rates agree with the single-weight path and with Euler's identity for many views", {
  set.seed(104)
  # Z is linear in the weights of each view, so sum_k w(k, l) dZ/dw(k, l) = Z
  for (K in list(rep(4L, 10L), c(3L, 5L, 2L, 4L, 3L, 4L, 2L, 3L, 5L, 4L, 3L, 2L))) {
    s <- random_mdi_state(K)
    Z <- mdir:::mdiNormalisingConstantCpp(s$w, K, s$phi)
    for (l in seq_along(K)) {
      rates <- mdir:::mdiWeightRatesCpp(s$w, K, s$phi, l - 1)
      expect_equal(sum(s$w[seq_len(K[l]), l] * rates), Z, tolerance = 1e-10)
    }
    if (length(K) == 10L) {
      for (l in c(1L, 6L, 10L)) {
        single <- vapply(seq_len(K[l]) - 1, function(k) mdir:::mdiWeightRateCpp(s$w, K, s$phi, l - 1, k), numeric(1))
        expect_equal(mdir:::mdiWeightRatesCpp(s$w, K, s$phi, l - 1), single, tolerance = 1e-10)
      }
    }
  }
})

test_that("weight rates are accurate when phi is very small or very large", {
  set.seed(105)
  K <- c(3L, 3L, 3L, 3L)
  for (scale in c(1e-8, 1e-4, 1e3)) {
    s <- random_mdi_state(K, phi_scale = scale)
    for (l in 1:4) {
      expect_equal(mdir:::mdiWeightRatesCpp(s$w, K, s$phi, l - 1),
                   vapply(1:3, function(k) bf_weight_rate(s$w, K, s$phi, l, k), numeric(1)), tolerance = 1e-9)
    }
  }
})

test_that("the MVN complete-data likelihood matches an independent calculation", {
  # complete_likelihood is the sum over items and views of the log density of
  # the item at its component's parameters; rebuild it from the saved draws
  log_mvn <- function(x, mu, S) {
    d <- x - mu
    -0.5 * (length(x) * log(2 * pi) + as.numeric(determinant(S, logarithm = TRUE)$modulus) +
              sum(d * solve(S, d)))
  }
  for (P in c(1L, 3L, 6L)) {
    set.seed(106 + P)
    N <- 30; K <- c(3L, 4L)
    X <- lapply(1:2, function(v) {
      m <- matrix(rnorm(N * P, rep(c(0, 3), each = N / 2)), N, P); rownames(m) <- seq_len(N); m
    })
    ch <- runMCMCChains(X, 1, R = 20, thin = 5, types = c("MVN", "MVN"), K = K)[[1]]
    for (s in 2:5) {
      total <- 0
      for (v in 1:2) {
        theta <- ch$parameters[[v]][s, ]
        mu <- matrix(theta[seq_len(P * K[v])], P, K[v])
        Sigma <- array(theta[-seq_len(P * K[v])], c(P, P, K[v]))
        lab <- ch$allocations[s, , v] + 1L  # raw chain labels are 0-based
        for (n in seq_len(N)) total <- total + log_mvn(X[[v]][n, ], mu[, lab[n]], matrix(Sigma[, , lab[n]], P, P))
      }
      expect_equal(ch$complete_likelihood[s], total, tolerance = 1e-9, info = paste("P =", P, "sample", s))
    }
  }
})

test_that("the normalising constant recorded by the sampler tracks the current phis and weights", {
  # With thin = 1, the evidence of a saved sample is Z at the start of the sweep
  # that follows the previous saved sample, so it must equal Z for the saved
  # weights and phis of that previous sample. This fails if Z (or the tables it
  # is built from) are not refreshed when the phis change.
  for (K in list(c(3L, 5L, 4L), c(4L, 4L, 4L, 3L))) {
    set.seed(107)
    L <- length(K); N <- 40
    X <- lapply(seq_len(L), function(v) {
      m <- matrix(rnorm(N * 2, rep(c(0, 3), each = N / 2)), N, 2); rownames(m) <- seq_len(N); m
    })
    ch <- runMCMCChains(X, 1, R = 25, thin = 1, types = rep("G", L), K = K)[[1]]
    pairs <- t(combn(L, 2))
    phi_at <- function(s) {
      phi <- matrix(0, L, L)
      for (i in seq_len(nrow(pairs))) phi[pairs[i, 1], pairs[i, 2]] <- phi[pairs[i, 2], pairs[i, 1]] <- ch$phis[s, i]
      phi
    }
    w_at <- function(s) {
      w <- matrix(0, max(K), L)
      for (l in seq_len(L)) w[, l] <- ch$weights[s, , l]
      w
    }
    expect_equal(ch$evidence[1], mdir:::mdiNormalisingConstantCpp(w_at(1), K, phi_at(1)), tolerance = 1e-10)
    for (s in 2:26) {
      expect_equal(ch$evidence[s], mdir:::mdiNormalisingConstantCpp(w_at(s - 1), K, phi_at(s - 1)),
                   tolerance = 1e-10, info = paste("K =", paste(K, collapse = ","), "sample", s))
    }
  }
})

test_that("weights beyond a view's number of components do not enter Z or its rates", {
  # The weight matrix is padded to K_max rows; after a label swap the padding of
  # a view with fewer components can hold non-zero values
  set.seed(108)
  for (K in list(c(2L, 4L, 3L), c(4L, 1L, 3L, 2L))) {
    s <- random_mdi_state(K)
    w_pad <- s$w
    for (l in seq_along(K)) if (K[l] < max(K)) w_pad[(K[l] + 1):max(K), l] <- runif(max(K) - K[l], 0.5, 2)
    expect_equal(mdir:::mdiNormalisingConstantCpp(w_pad, K, s$phi), bf_Z(s$w, K, s$phi), tolerance = 1e-10)
    for (l in seq_along(K)) {
      expect_equal(mdir:::mdiWeightRatesCpp(w_pad, K, s$phi, l - 1),
                   vapply(seq_len(K[l]), function(k) bf_weight_rate(s$w, K, s$phi, l, k), numeric(1)),
                   tolerance = 1e-10)
    }
  }
})

test_that("the label-swap ratio equals the exact ratio of the model's target", {
  # Target over (c, w) given the strategic latent variable v, from the MDI joint
  #   prod_n [prod_l w[c_nl, l] prod_{l<m}(1 + phi_lm 1[c_nl = c_nm])] exp(-v Z(w))
  # (Kirk et al. 2012; Coleman et al. 2025, S1 Text eq. 20). Exchanging
  # components k and k' of view l moves their labels and weights together; the
  # Gamma prior on the weights and the data likelihood (carried with the
  # parameters) are unchanged. Z is computed by enumeration.
  log_target <- function(cc, w, phi, K, v) {
    L <- length(K)
    lp <- -v * bf_Z(w, K, phi)
    for (n in seq_len(nrow(cc))) {
      lp <- lp + sum(log(w[cbind(cc[n, ], seq_len(L))]))
      if (L > 1) for (l in 1:(L - 1)) for (m in (l + 1):L) lp <- lp + log1p(phi[l, m] * (cc[n, l] == cc[n, m]))
    }
    lp
  }
  set.seed(109)
  for (K in list(c(3L, 3L), c(4L, 4L, 4L), c(2L, 4L, 3L), c(3L, 1L, 4L, 2L))) {
    L <- length(K); N <- 15
    s <- random_mdi_state(K, phi_scale = 4)
    cc <- sapply(seq_len(L), function(l) sample.int(K[l], N, TRUE))
    v <- rgamma(1, N, bf_Z(s$w, K, s$phi))
    for (l in seq_len(L)) if (K[l] > 1) for (rep in 1:4) {
      kk <- sample.int(K[l], 2); k <- kk[1]; kp <- kk[2]
      cs <- cc; cs[cc[, l] == k, l] <- kp; cs[cc[, l] == kp, l] <- k
      ws <- s$w; ws[c(k, kp), l] <- ws[c(kp, k), l]
      expected <- log_target(cs, ws, s$phi, K, v) - log_target(cc, s$w, s$phi, K, v)
      expect_equal(mdir:::mdiSwapLogRatioCpp(cc - 1L, s$phi, s$w, K, v, l - 1, k - 1, kp - 1), expected,
                   tolerance = 1e-10, info = paste("K =", paste(K, collapse = ","), "view", l, "swap", k, kp))
    }
  }
})

test_that("a label swap changes Z when the views are aligned differently", {
  # Z is unchanged by the same permutation of the components of every view, but
  # not by a permutation of one view's weights; the swap ratio must see the latter
  set.seed(110)
  K <- c(4L, 4L, 4L); s <- random_mdi_state(K, phi_scale = 4)
  w_all <- s$w; w_all[c(1, 2), ] <- w_all[c(2, 1), ]
  w_one <- s$w; w_one[c(1, 2), 1] <- w_one[c(2, 1), 1]
  Z <- bf_Z(s$w, K, s$phi)
  expect_equal(bf_Z(w_all, K, s$phi), Z, tolerance = 1e-12)
  expect_gt(abs(bf_Z(w_one, K, s$phi) - Z) / Z, 1e-3)
})
