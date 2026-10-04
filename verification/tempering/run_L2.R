args <- commandArgs(TRUE)
lib <- args[1]; n_chains <- as.integer(args[2]); R <- as.integer(args[3])
n_cores <- if (length(args) > 3) as.integer(args[4]) else 1L
.libPaths(c(lib, .libPaths())); suppressMessages(library(mdir)); options(mdir.quiet = TRUE)
source("exact_reference.R"); source("compare.R")
set.seed(20260411)
N <- 4; K <- 2
X1 <- matrix(sample(0:2, N, TRUE), N, 1); X2 <- matrix(sample(0:2, N, TRUE), N, 1)
rownames(X1) <- rownames(X2) <- seq_len(N)
cat("view 1:", X1[, 1], " view 2:", X2[, 1], "\n")
a1 <- mdir:::densityHyperparameters(X1, K, 2, numeric(0))$concentration
a2 <- mdir:::densityHyperparameters(X2, K, 2, numeric(0))$concentration
cat("alpha:", unlist(a1), "|", unlist(a2), "\n")
t0 <- Sys.time()
pc <- mc_prior_cells(N, K, n_draws = 4e7, chunk = 1e6)
cat(sprintf("prior MC done in %.0fs; max relative MC se of cell prior: %.4f\n", as.numeric(Sys.time() - t0, units = "secs"), max(pc$se / pc$prior)))
run <- function(betas, seed) {
  set.seed(seed)
  fit <- callMDI(list(X1, X2), R = R, thin = 1, types = c("C", "C"), K = c(K, K), betas = betas,
                 check_prior = FALSE, save_parameters = FALSE)
  a <- fit$allocations[-(1:200), , ]
  config_code(a[, , 1], K) + K^N * config_code(a[, , 2], K)
}
configs <- list(
  list(name = "single chain beta = 1", betas = 1, beta = 1),
  list(name = "single chain beta = 0.5 (tempered target)", betas = 0.5, beta = 0.5),
  list(name = "PT cold chain, ladder (0.2, 0.5, 1)", betas = c(0.2, 0.5, 1), beta = 1),
  list(name = "CONTROL: single chain beta = 1 vs exact beta = 0.5", betas = 1, beta = 0.5)
)
for (cf in configs) {
  ex <- exact_L2(list(X1, X2), K, list(a1, a2), cf$beta, pc)
  codes <- parallel::mclapply(seq_len(n_chains), function(i) run(cf$betas, 5000 + i), mc.cores = n_cores)
  fr <- freq_by_chain(codes, length(ex$prob))
  res <- compare_to_exact(fr, ex$prob)
  cat(sprintf("%-55s TV %.4f (chain noise ~%.4f) | states %d | max|z| %.2f | mean z^2 %.2f (expect ~%.2f)\n",
    cf$name, res$tv, res$tv_chain_noise, res$n_states, res$max_abs_z, res$mean_z2, res$expected_mean_z2))
  cat(sprintf("    exact-reference MC error: TV contribution <= %.4f\n", 0.5 * sum(ex$prob * ex$rel_se)))
}
