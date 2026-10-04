# Semi-supervised: items 1 and 2 observed in classes 0 and 1. The target is the exact label
# posterior of the free items given the observed ones; checked for a plain chain, a PT cold
# chain and a chain with split-merge, plus a wrong-beta control.
args <- commandArgs(TRUE)
lib <- args[1]; n_chains <- as.integer(args[2]); R <- as.integer(args[3]); n_cores <- as.integer(args[4])
.libPaths(c(lib, .libPaths())); suppressMessages(library(mdir)); options(mdir.quiet = TRUE)
source("exact_reference.R"); source("compare.R")
set.seed(20260410)
N <- 7; K <- 3
X <- matrix(c(sample(0:2, N, TRUE, c(.5, .3, .2))), N, 1); rownames(X) <- seq_len(N)
hp <- mdir:::densityHyperparameters(X, K, 2, numeric(0)); alpha <- hp$concentration
fixed <- matrix(c(1, 1, rep(0, N - 2)), ncol = 1); lab <- matrix(c(0, 1, rep(0, N - 2)), ncol = 1)
run <- function(betas, sm, seed) {
  set.seed(seed)
  fit <- callMDI(list(X), R = R, thin = 1, types = "C", K = K, betas = betas, split_merge = sm, fixed = fixed,
    initial_labels = lab, check_prior = FALSE, save_parameters = FALSE)
  stopifnot(all(fit$allocations[, 1, 1] == 0), all(fit$allocations[, 2, 1] == 1))
  fit$allocations[-(1:200), 3:N, 1]
}
configs <- list(
  list("semi-supervised plain chain", 1, 0L, 1), list("semi-supervised PT (0.2, 0.5, 1)", c(0.2, 0.5, 1), 0L, 1),
  list("semi-supervised, split-merge", 1, 2L, 1), list("semi-supervised, tempered beta = 0.4", 0.4, 0L, 0.4),
  list("CONTROL: beta = 0.4 vs exact beta = 1", 0.4, 0L, 1))
for (cf in configs) {
  ex <- exact_L1(X, K, alpha, cf[[4]])
  keep <- ex$grid[, 1] == 0 & ex$grid[, 2] == 1
  grid <- ex$grid[keep, 3:N, drop = FALSE]; prob <- ex$prob[keep] / sum(ex$prob[keep])
  # grid rows are in expand.grid order (first column fastest) of the free items
  codes <- parallel::mclapply(seq_len(n_chains), function(i) config_code(run(cf[[2]], cf[[3]], 1000 + i), K), mc.cores = n_cores)
  fr <- freq_by_chain(codes, nrow(grid))
  stopifnot(all(config_code(grid, K) == seq_len(nrow(grid)) - 1))
  res <- compare_to_exact(fr, prob)
  cat(sprintf("%-42s TV %.4f (noise ~%.4f) | states %d | max|z| %.2f | mean z^2 %.2f (expect ~%.2f)\n", cf[[1]], res$tv,
    res$tv_chain_noise, res$n_states, res$max_abs_z, res$mean_z2, res$expected_mean_z2))
}
