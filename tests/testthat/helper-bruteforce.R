# Brute-force references for the MDI normalising constant and its derivatives.
# Enumerates every joint assignment (k_1, ..., k_L); only usable for small K^L.
bf_grid <- function(K) as.matrix(expand.grid(lapply(K, seq_len) ))

bf_terms <- function(w, K, phi, grid = bf_grid(K)) {
  L <- length(K)
  lw <- rowSums(matrix(vapply(seq_len(L), function(l) log(w[grid[, l], l]), numeric(nrow(grid))), nrow = nrow(grid)))
  up <- rep(1, nrow(grid))
  if (L > 1) for (l in 1:(L - 1)) for (m in (l + 1):L) up <- up * (1 + phi[l, m] * (grid[, l] == grid[, m]))
  list(grid = grid, term = exp(lw) * up)
}

bf_Z <- function(w, K, phi) sum(bf_terms(w, K, phi)$term)

# dZ/dw[k, l]: sum over assignments with k_l = k of Z-term / w[k, l]
bf_weight_rate <- function(w, K, phi, l, k) {
  tt <- bf_terms(w, K, phi)
  sel <- tt$grid[, l] == k
  sum(tt$term[sel]) / w[k, l]
}

# dZ/dphi[l, m]: terms with k_l = k_m, divided by (1 + phi[l, m])
bf_phi_rate <- function(w, K, phi, l, m) {
  tt <- bf_terms(w, K, phi)
  sel <- tt$grid[, l] == tt$grid[, m]
  sum(tt$term[sel]) / (1 + phi[l, m])
}

random_mdi_state <- function(K, phi_scale = 5, w_shape = 1) {
  L <- length(K)
  w <- matrix(0, max(K), L)
  for (l in seq_len(L)) w[seq_len(K[l]), l] <- rgamma(K[l], w_shape, 2)
  phi <- matrix(0, L, L)
  if (L > 1) for (l in 1:(L - 1)) for (m in (l + 1):L) phi[l, m] <- phi[m, l] <- rgamma(1, 2, 2 / phi_scale)
  list(w = w, phi = phi)
}

# Effective sample size of a single chain via the split-chain estimator
ess_of <- function(x) {
  n <- length(x) %/% 2
  rankNormalizedRhat(cbind(x[seq_len(n)], x[n + seq_len(n)]))$ess_bulk
}

# Posterior mean and sd of a density given on a grid
grid_moments <- function(x, logdens) {
  w <- exp(logdens - max(logdens)); w <- w / sum(w)
  m <- sum(w * x)
  c(mean = m, sd = sqrt(sum(w * x^2) - m^2))
}

# Brute-force marginal likelihood of one item: log of
#   sum_{k_1..k_L} prod_l w[k_l, l] g[k_l, l] prod_{l<m} (1 + phi[l, m] 1[k_l = k_m])
# minus log Z, with g the K_max x L matrix of per-view component likelihoods on
# the log scale (-Inf or NA beyond K[l] or for excluded components).
bf_log_marginal <- function(log_g, w, K, phi) {
  g <- exp(log_g); g[!is.finite(log_g)] <- 0
  list(
    log_marginal = log(sum(bf_terms(w * g, K, phi)$term)) - log(bf_Z(w, K, phi)),
    # p(c_l = k | x)
    class_prob = {
      tt <- bf_terms(w * g, K, phi)
      tt$term <- tt$term / sum(tt$term)
      sapply(seq_along(K), function(l) vapply(seq_len(max(K)), function(k) sum(tt$term[tt$grid[, l] == k]), numeric(1)))
    }
  )
}

# The state (weights, phi matrix) of saved draw s of a chain
bf_state <- function(ch, s) {
  L <- ch$V
  w <- matrix(0, max(ch$K), L)
  for (l in seq_len(L)) w[, l] <- ch$weights[s, , l]
  phi <- matrix(0, L, L)
  if (L > 1) {
    pairs <- t(utils::combn(L, 2))
    for (i in seq_len(nrow(pairs))) phi[pairs[i, 1], pairs[i, 2]] <- phi[pairs[i, 2], pairs[i, 1]] <- ch$phis[s, i]
  }
  list(w = w, phi = phi)
}

# Log-likelihood of the observed entries of item x (a vector, NA = missing) in
# each component of a diagonal Gaussian view ("G") at saved draw s
bf_gaussian_loglik <- function(ch, v, s, x) {
  P <- length(x); K <- ch$K[v]
  th <- ch$parameters[[v]][s, ]
  mu <- matrix(th[seq_len(P * K)], P, K)
  vr <- matrix(th[P * K + seq_len(P * K)], P, K)
  obs <- !is.na(x)
  vapply(seq_len(K), function(k) sum(stats::dnorm(x[obs], mu[obs, k], sqrt(vr[obs, k]), log = TRUE)), numeric(1))
}
