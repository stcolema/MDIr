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
