# With a constant likelihood (a categorical view with one category) the
# stationary distribution of the sampler is the prior. This exercises the
# weights, phi, mass, strategic latent variable and label-swap updates jointly
# and failed badly before the conditionals were corrected (phi mean 2.8 vs 10,
# mass mean 6.5 vs 20).

prior_recovery_run <- function(K, R, N = 12, seed = 1) {
  set.seed(seed)
  L <- length(K)
  X <- lapply(seq_len(L), function(l) { m <- matrix(0, N, 1); rownames(m) <- seq_len(N); m })
  callMDI(X, R = R, thin = 10, types = rep("C", L), K = K)
}

prior_simulator <- function(K, N = 12) {
  L <- length(K)
  mass <- rgamma(L, 2, 0.1)
  phi <- if (L > 1) rgamma(choose(L, 2), 2, 0.2) else 0
  w <- lapply(seq_len(L), function(l) rgamma(K[l], mass[l] / K[l], 2))
  grid <- as.matrix(expand.grid(lapply(K, seq_len)))
  lp <- rowSums(sapply(seq_len(L), function(l) log(w[[l]][grid[, l]])))
  idx <- 0
  if (L > 1) for (a in 1:(L - 1)) for (b in (a + 1):L) {
    idx <- idx + 1
    lp <- lp + log1p(phi[idx] * (grid[, a] == grid[, b]))
  }
  p <- exp(lp - max(lp)); p <- p / sum(p)
  g <- grid[sample.int(nrow(grid), N, TRUE, p), , drop = FALSE]
  c(agree = if (L > 1) mean(g[, 1] == g[, 2]) else 0, nocc = length(unique(g[, 1])))
}

test_that("sampler recovers the prior under a constant likelihood", {
  skip_on_cran()
  for (K in list(c(3L, 3L), c(2L, 4L, 3L))) {
    out <- prior_recovery_run(K, R = 60000)
    keep <- -(1:500)
    L <- length(K)
    # Exact prior moments: mass ~ Gamma(2, 0.1), phi ~ Gamma(2, 0.2)
    m <- out$mass[keep, 1]
    ess_m <- rankNormalizedRhat(cbind(m[1:(length(m) %/% 2 * 2)][c(TRUE, FALSE)], m[1:(length(m) %/% 2 * 2)][c(FALSE, TRUE)]))$ess_bulk
    expect_lt(abs(mean(m) - 20), 5 * 14.14 / sqrt(ess_m), label = "mass mean")
    phi <- out$phis[keep, 1]
    ess_p <- rankNormalizedRhat(cbind(phi[1:(length(phi) %/% 2 * 2)][c(TRUE, FALSE)], phi[1:(length(phi) %/% 2 * 2)][c(FALSE, TRUE)]))$ess_bulk
    expect_lt(abs(mean(phi) - 10), 5 * 7.07 / sqrt(ess_p), label = "phi mean")
    # Cross-view agreement and occupancy against direct simulation of the prior
    ref <- t(replicate(3000, prior_simulator(K)))
    al <- out$allocations[keep, , , drop = FALSE]
    agree <- rowMeans(al[, , 1] == al[, , 2])
    nocc <- apply(al[, , 1], 1, function(z) length(unique(z)))
    expect_lt(abs(mean(agree) - mean(ref[, "agree"])), 0.04)
    expect_lt(abs(mean(nocc) - mean(ref[, "nocc"])), 0.15)
  }
})
