args <- commandArgs(TRUE)
lib <- args[1]; n_chains <- as.integer(args[2]); R <- as.integer(args[3])
n_cores <- if (length(args) > 3) as.integer(args[4]) else 1L
.libPaths(c(lib, .libPaths())); suppressMessages(library(mdir)); options(mdir.quiet = TRUE)
source("exact_reference.R"); source("compare.R")
set.seed(20260410)
N <- 7; K <- 3
X <- matrix(c(sample(0:2, N, TRUE, c(.5, .3, .2))), N, 1); rownames(X) <- seq_len(N)
hp <- mdir:::densityHyperparameters(X, K, 2, numeric(0))
alpha <- hp$concentration
cat("data:", X[, 1], " alpha:", unlist(alpha), "\n")
run <- function(betas, seed) {
  set.seed(seed)
  fit <- callMDI(list(X), R = R, thin = 1, types = "C", K = K, betas = betas, check_prior = FALSE, save_parameters = FALSE)
  fit$allocations[-(1:200), , 1]
}
# name, ladder run, beta of the exact reference compared against; the last two
# are negative controls (a wrong reference), which must be rejected
configs <- list(
  list(name = "single chain beta = 1", betas = 1, beta = 1),
  list(name = "single chain beta = 0.4 (tempered target)", betas = 0.4, beta = 0.4),
  list(name = "PT cold chain, ladder (0.2, 0.5, 1)", betas = c(0.2, 0.5, 1), beta = 1),
  list(name = "PT cold chain, ladder (0.05, 0.15, 0.4, 1)", betas = c(0.05, 0.15, 0.4, 1), beta = 1),
  list(name = "CONTROL: single chain beta = 0.4 vs exact beta = 1", betas = 0.4, beta = 1),
  list(name = "CONTROL: single chain beta = 1 vs exact beta = 0.7", betas = 1, beta = 0.7)
)
for (cf in configs) {
  ex <- exact_L1(X, K, alpha, cf$beta)
  codes <- parallel::mclapply(seq_len(n_chains), function(i) config_code(run(cf$betas, 1000 + i), K), mc.cores = n_cores)
  fr <- freq_by_chain(codes, nrow(ex$grid))
  # exact$grid rows are in expand.grid order, whose code is the same mixed radix
  stopifnot(all(config_code(ex$grid, K) == seq_len(nrow(ex$grid)) - 1))
  res <- compare_to_exact(fr, ex$prob)
  cat(sprintf("%-62s TV %.4f (chain noise ~%.4f) | states %d | max|z| %.2f | mean z^2 %.2f (expect ~%.2f)\n",
    cf$name, res$tv, res$tv_chain_noise, res$n_states, res$max_abs_z, res$mean_z2, res$expected_mean_z2))
}
