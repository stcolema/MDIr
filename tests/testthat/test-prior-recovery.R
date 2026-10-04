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

# The label swap exchanges the labels, weights and parameters of two components
# of one view. If the weights of the other views are moved too, the acceptance
# ratio is wrong (Z is unchanged by the same permutation of every view) and
# items sit in the heavy components of their own view less often than the prior
# implies: the share of items in the largest-weight component was about 0.011
# too low with K = (3, 3), N = 10, and the mean normalised weight of the assigned
# component about 0.005 too low. The agreement and occupancy summaries above move
# by less than their tolerances, so they cannot see this.
test_that("sampler recovers the prior alignment of weights and labels", {
  skip_on_cran()
  K <- c(3L, 3L); N <- 10; L <- 2
  stats_of <- function(cc, w) {
    rho <- sweep(w, 2, colSums(w), "/")
    top <- apply(w, 2, which.max)
    c(own1 = mean(rho[cbind(cc[, 1], 1)]), own2 = mean(rho[cbind(cc[, 2], 2)]),
      top1 = mean(cc[, 1] == top[1]), top2 = mean(cc[, 2] == top[2]))
  }
  forward <- function() {
    mass <- rgamma(L, 2, 0.1); phi <- rgamma(1, 2, 0.2)
    w <- sapply(seq_len(L), function(l) rgamma(K[l], mass[l] / K[l], 2))
    grid <- as.matrix(expand.grid(lapply(K, seq_len)))
    lp <- rowSums(sapply(seq_len(L), function(l) log(w[grid[, l], l]))) + log1p(phi * (grid[, 1] == grid[, 2]))
    p <- exp(lp - max(lp)); p <- p / sum(p)
    stats_of(grid[sample.int(nrow(grid), N, TRUE, p), , drop = FALSE], w)
  }
  set.seed(2024)
  ref <- t(replicate(20000, forward()))
  out <- prior_recovery_run(K, R = 400000, N = N, seed = 2025)
  keep <- -(1:50)
  al <- out$allocations[keep, , , drop = FALSE] + 1L
  wt <- out$weights[keep, , , drop = FALSE]
  mc <- t(vapply(seq_len(dim(al)[1]), function(s) stats_of(al[s, , ], wt[s, , ]), numeric(4)))
  batch_se <- function(x, nb = 40) sd(vapply(split(x, cut(seq_along(x), nb, labels = FALSE)), mean, numeric(1))) / sqrt(nb)
  for (j in colnames(ref)) {
    se <- sqrt(batch_se(mc[, j])^2 + var(ref[, j]) / nrow(ref))
    expect_lt(abs(mean(mc[, j]) - mean(ref[, j])), 5 * se, label = paste("prior mean of", j))
  }
})
