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
