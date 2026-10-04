# Exact tempered label posterior for one view of a one-dimensional Gaussian
# model (MVN with P = 1 and the diagonal "G" density share the normal-inverse-
# gamma prior), variance scale held fixed (no pooling): the collapsed marginal
# likelihood of a block of items under L^beta is closed form (verified
# symbolically and by quadrature in sympy_checks.py):
#   m_beta = (2 pi)^(-beta n / 2) (kappa / (kappa + beta n))^(1/2) (psi/2)^(nu/2) / Gamma(nu/2)
#            * Gamma((nu + beta n)/2) / (psi_n / 2)^((nu + beta n)/2),
#   psi_n = psi + beta SS + kappa beta n / (kappa + beta n) (xbar - xi)^2.
args <- commandArgs(TRUE)
lib <- args[1]; n_chains <- as.integer(args[2]); R <- as.integer(args[3]); n_cores <- as.integer(args[4])
.libPaths(c(lib, .libPaths())); suppressMessages(library(mdir)); options(mdir.quiet = TRUE)
source("exact_reference.R"); source("compare.R")
set.seed(20260412)
N <- 6; K <- 2
X <- matrix(c(rnorm(3, -1.2, 0.8), rnorm(3, 1.5, 0.8)), N, 1); rownames(X) <- seq_len(N)
dp <- densityPrior(scale_pool_shape = 0)
cat("data:", round(X[, 1], 2), "\n")

exact_gauss <- function(type, beta) {
  mt <- if (type == "MVN") 1 else 0
  hp <- mdir:::densityHyperparameters(X, K, mt, as.numeric(dp))
  xi <- as.numeric(hp$xi); psi <- as.numeric(hp$scale); kappa <- hp$kappa; nu <- hp$nu
  lm <- function(x) {
    n <- length(x); if (n == 0) return(0)
    be <- beta * n; xb <- mean(x); SS <- sum((x - xb)^2)
    pn <- psi + beta * SS + kappa * be / (kappa + be) * (xb - xi)^2
    -be / 2 * log(2 * pi) + 0.5 * log(kappa / (kappa + be)) + nu / 2 * log(psi / 2) - lgamma(nu / 2) +
      lgamma((nu + be) / 2) - (nu + be) / 2 * log(pn / 2)
  }
  grid <- as.matrix(expand.grid(rep(list(0:(K - 1)), N)))
  cache <- new.env(); lp <- numeric(nrow(grid))
  for (i in seq_len(nrow(grid))) {
    c_i <- grid[i, ]; Nk <- tabulate(c_i + 1, K); key <- paste(Nk, collapse = "_")
    if (is.null(cache[[key]])) cache[[key]] <- log(prior_counts_L1(Nk, K))
    lp[i] <- cache[[key]] + sum(vapply(0:(K - 1), function(k) lm(X[c_i == k, 1]), numeric(1)))
  }
  p <- exp(lp - max(lp)); list(grid = grid, prob = p / sum(p))
}
run <- function(type, betas, seed) {
  set.seed(seed)
  a <- callMDI(list(X), R = R, thin = 1, types = type, K = K, betas = betas, check_prior = FALSE,
    save_parameters = FALSE, density_prior = dp)$allocations[-(1:200), , 1]
  config_code(a, K)
}
configs <- list(
  list(type = "MVN", name = "MVN (P = 1) single chain beta = 1", betas = 1, beta = 1),
  list(type = "MVN", name = "MVN (P = 1) single chain beta = 0.3", betas = 0.3, beta = 0.3),
  list(type = "MVN", name = "MVN (P = 1) PT cold chain (0.1, 0.4, 1)", betas = c(0.1, 0.4, 1), beta = 1),
  list(type = "G", name = "G single chain beta = 0.3", betas = 0.3, beta = 0.3),
  list(type = "G", name = "G PT cold chain (0.1, 0.4, 1)", betas = c(0.1, 0.4, 1), beta = 1),
  list(type = "MVN", name = "CONTROL: MVN single chain beta = 0.3 vs exact beta = 1", betas = 0.3, beta = 1)
)
for (cf in configs) {
  ex <- exact_gauss(cf$type, cf$beta)
  codes <- parallel::mclapply(seq_len(n_chains), function(i) run(cf$type, cf$betas, 3000 + i), mc.cores = n_cores)
  fr <- freq_by_chain(codes, nrow(ex$grid))
  res <- compare_to_exact(fr, ex$prob)
  cat(sprintf("%-58s TV %.4f (chain noise ~%.4f) | states %d | max|z| %.2f | mean z^2 %.2f (expect ~%.2f)\n",
    cf$name, res$tv, res$tv_chain_noise, res$n_states, res$max_abs_z, res$mean_z2, res$expected_mean_z2))
}
