# Semi-supervised view: some items have observed labels (fixed). The exact tempered
# posterior over the free items' labels uses all items in the likelihood (the observed
# labels are data, their likelihood is tempered with the rest) and the Dirichlet-
# multinomial label prior over all items.
args <- commandArgs(TRUE)
lib <- args[1]; n_chains <- as.integer(args[2]); R <- as.integer(args[3]); n_cores <- as.integer(args[4])
.libPaths(c(lib, .libPaths())); suppressMessages(library(mdir)); options(mdir.quiet = TRUE)
source("exact_reference.R"); source("compare.R")
set.seed(20260413)
N <- 6; K <- 3
X <- matrix(sample(0:2, N, TRUE, c(.5, .3, .2)), N, 1); rownames(X) <- seq_len(N)
fixed <- c(1, 1, 0, 0, 0, 0); lab_fixed <- c(0, 1, 0, 0, 0, 0)   # items 1, 2 observed in classes 0, 1
alpha <- mdir:::densityHyperparameters(X, K, 2, numeric(0))$concentration
exact_fixed <- function(beta) {
  grid <- as.matrix(expand.grid(rep(list(0:(K - 1)), N - 2)))
  cache <- new.env(); lp <- numeric(nrow(grid))
  for (i in seq_len(nrow(grid))) {
    c_i <- c(lab_fixed[1:2], grid[i, ]); Nk <- tabulate(c_i + 1, K); key <- paste(Nk, collapse = "_")
    if (is.null(cache[[key]])) cache[[key]] <- log(prior_counts_L1(Nk, K))
    ll <- sum(vapply(0:(K - 1), function(k) if (any(c_i == k)) log_marg_cat(X[c_i == k, , drop = FALSE], alpha, beta) else 0, numeric(1)))
    lp[i] <- cache[[key]] + ll
  }
  p <- exp(lp - max(lp)); list(grid = grid, prob = p / sum(p))
}
init <- matrix(lab_fixed, N, 1)
run <- function(betas, seed) {
  set.seed(seed)
  a <- callMDI(list(X), R = R, thin = 1, types = "C", K = K, betas = betas, check_prior = FALSE, save_parameters = FALSE,
    initial_labels = init, fixed = matrix(fixed, N, 1), initial_labels_as_intended = FALSE)$allocations[-(1:200), , 1]
  stopifnot(all(a[, 1] == 0), all(a[, 2] == 1))     # observed labels kept
  config_code(a[, 3:N], K)
}
configs <- list(
  list(name = "semi-supervised, single chain beta = 1", betas = 1, beta = 1),
  list(name = "semi-supervised, single chain beta = 0.4", betas = 0.4, beta = 0.4),
  list(name = "semi-supervised, PT cold chain (0.2, 0.5, 1)", betas = c(0.2, 0.5, 1), beta = 1),
  list(name = "CONTROL: single chain beta = 0.4 vs exact beta = 1", betas = 0.4, beta = 1)
)
for (cf in configs) {
  ex <- exact_fixed(cf$beta)
  codes <- parallel::mclapply(seq_len(n_chains), function(i) run(cf$betas, 4000 + i), mc.cores = n_cores)
  fr <- freq_by_chain(codes, nrow(ex$grid))
  res <- compare_to_exact(fr, ex$prob)
  cat(sprintf("%-56s TV %.4f (chain noise ~%.4f) | states %d | max|z| %.2f | mean z^2 %.2f (expect ~%.2f)\n",
    cf$name, res$tv, res$tv_chain_noise, res$n_states, res$max_abs_z, res$mean_z2, res$expected_mean_z2))
}
