# Exact (L = 1) and Monte Carlo (L = 2) reference posteriors over the labelled
# component assignments of tiny categorical MDI models, at inverse temperature
# beta. Used to check the tempered sampler and parallel tempering against
# something the sampler does not share code with.
#
# Target: pi_beta(c) proportional to prior(c) * prod_l prod_k B(alpha_l + beta n_lk) / B(alpha_l)
# where n_lk are category counts of the items assigned to component k of view l,
# alpha_l the Dirichlet concentration (data driven, as in the package) and
# prior(c) the marginal prior of the labels with weights, phis and masses
# integrated out.

lbeta_vec <- function(a) sum(lgamma(a)) - lgamma(sum(a))

# log marginal likelihood of the items of one component of one view under L^beta
log_marg_cat <- function(Xcol_items, alpha, beta) {
  # Xcol_items: matrix (n_k x P) of 0-based categories for the items in the component
  P <- length(alpha)
  out <- 0
  for (p in seq_len(P)) {
    cnt <- tabulate(Xcol_items[, p] + 1, length(alpha[[p]]))
    out <- out + lbeta_vec(alpha[[p]] + beta * cnt) - lbeta_vec(alpha[[p]])
  }
  out
}

# L = 1: prior(c) depends on the counts N_k only; Dirichlet-multinomial with
# concentration mass / K per component, integrated over mass ~ Gamma(shape, rate).
prior_counts_L1 <- function(Nk, K, mass_shape = 2, mass_rate = 0.1) {
  N <- sum(Nk)
  f <- function(m) {
    a <- m / K
    lp <- lgamma(m) - lgamma(m + N) + sum(lgamma(a + Nk) - lgamma(a)) +
      dgamma(m, mass_shape, mass_rate, log = TRUE)
    exp(lp)
  }
  stats::integrate(Vectorize(f), 0, Inf, rel.tol = 1e-10, subdivisions = 2000L)$value
}

exact_L1 <- function(X, K, alpha, beta) {
  N <- nrow(X)
  grid <- as.matrix(expand.grid(rep(list(0:(K - 1)), N)))
  cache <- new.env()
  lp <- numeric(nrow(grid))
  for (i in seq_len(nrow(grid))) {
    c_i <- grid[i, ]
    Nk <- tabulate(c_i + 1, K)
    key <- paste(Nk, collapse = "_")
    if (is.null(cache[[key]])) cache[[key]] <- log(prior_counts_L1(Nk, K))
    ll <- 0
    for (k in 0:(K - 1)) {
      idx <- which(c_i == k)
      if (length(idx)) ll <- ll + log_marg_cat(X[idx, , drop = FALSE], alpha, beta)
    }
    lp[i] <- cache[[key]] + ll
  }
  p <- exp(lp - max(lp)); p <- p / sum(p)
  # log of sum_c prior(c) prod marginal likelihoods: the evidence p(X) at this beta
  list(grid = grid, prob = p, log_evidence = max(lp) + log(sum(exp(lp - max(lp)))))
}

# L = 2, K_1 = K_2 = K: the prior of the joint labels is E[prod_n pi_{c_n}] over
# the prior of (masses, weights, phi), estimated by Monte Carlo (the items are
# independent given those). It depends on the counts of the K^2 joint cells only.
mc_prior_cells <- function(N, K, n_draws = 2e7, chunk = 1e6, w_rate = 2, mass_shape = 2, mass_rate = 0.1,
                           phi_shape = 2, phi_rate = 0.2) {
  cells <- as.matrix(expand.grid(rep(list(0:N), K * K)))
  cells <- cells[rowSums(cells) == N, , drop = FALSE]
  acc <- numeric(nrow(cells)); acc2 <- numeric(nrow(cells))
  n_chunks <- n_draws / chunk
  # Gamma(a, rate) draws on the log scale: G = Gamma(a + 1, rate) * U^(1 / a). A
  # prior with small mass / K puts weights below the smallest double, which would
  # otherwise give 0 / 0 in the cell probabilities.
  rlgamma <- function(n, a, rate) log(rgamma(n, a + 1, rate)) + log(runif(n)) / a
  for (b in seq_len(n_chunks)) {
    m1 <- rgamma(chunk, mass_shape, mass_rate); m2 <- rgamma(chunk, mass_shape, mass_rate)
    lW1 <- matrix(rlgamma(chunk * K, rep(m1 / K, K), w_rate), chunk)
    lW2 <- matrix(rlgamma(chunk * K, rep(m2 / K, K), w_rate), chunk)
    phi <- rgamma(chunk, phi_shape, phi_rate)
    lcell <- matrix(0, chunk, K * K)
    j <- 0
    for (k2 in 1:K) for (k1 in 1:K) {   # column-major cell order matches expand.grid below
      j <- j + 1
      lcell[, j] <- lW1[, k1] + lW2[, k2] + log1p(phi * (k1 == k2))
    }
    mx <- apply(lcell, 1, max)
    lc <- lcell - (mx + log(rowSums(exp(lcell - mx))))
    for (r in seq_len(nrow(cells))) {
      # only occupied cells enter; exp(-Inf) is fine, 0 * -Inf is avoided
      nz <- cells[r, ] > 0
      v <- exp(lc[, nz, drop = FALSE] %*% cells[r, nz])
      acc[r] <- acc[r] + sum(v); acc2[r] <- acc2[r] + sum(v^2)
    }
  }
  mean_p <- acc / n_draws
  se <- sqrt(pmax(acc2 / n_draws - mean_p^2, 0) / n_draws)
  list(cells = cells, prior = mean_p, se = se)
}

# Exact tempered posterior over all labelled (c1, c2) configurations, with the
# prior of the cell counts from mc_prior_cells(); returns the MC error too
exact_L2 <- function(Xlist, K, alphas, beta, pc) {
  N <- nrow(Xlist[[1]])
  grid1 <- as.matrix(expand.grid(rep(list(0:(K - 1)), N)))
  idx <- expand.grid(i1 = seq_len(nrow(grid1)), i2 = seq_len(nrow(grid1)))
  lik_view <- function(l) {
    vapply(seq_len(nrow(grid1)), function(i) {
      ll <- 0
      for (k in 0:(K - 1)) {
        it <- which(grid1[i, ] == k)
        if (length(it)) ll <- ll + log_marg_cat(Xlist[[l]][it, , drop = FALSE], alphas[[l]], beta)
      }
      ll
    }, numeric(1))
  }
  l1 <- lik_view(1); l2 <- lik_view(2)
  key <- apply(pc$cells, 1, paste, collapse = "_")
  lp <- numeric(nrow(idx)); lpse <- numeric(nrow(idx))
  for (r in seq_len(nrow(idx))) {
    c1 <- grid1[idx$i1[r], ]; c2 <- grid1[idx$i2[r], ]
    cnt <- tabulate(c1 + K * c2 + 1, K * K)      # cell (k1, k2) -> k1 + K k2 (column-major)
    row <- match(paste(cnt, collapse = "_"), key)
    lp[r] <- log(pc$prior[row]) + l1[idx$i1[r]] + l2[idx$i2[r]]
    lpse[r] <- pc$se[row] / pc$prior[row]
  }
  p <- exp(lp - max(lp)); p <- p / sum(p)
  list(c1 = grid1[idx$i1, , drop = FALSE], c2 = grid1[idx$i2, , drop = FALSE], prob = p, rel_se = lpse)
}
