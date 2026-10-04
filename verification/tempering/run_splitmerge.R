# The split-merge move on its own (weights fixed, nothing else updated) has the collapsed
# target pi(c) proportional to prod_n w[c_n] prod_k m_beta(X_k) as its stationary law; compare
# long chains with exact enumeration. Observed labels (fixed items) stay put but count.
args <- commandArgs(TRUE)
lib <- args[1]; n_chains <- as.integer(args[2]); n_iter <- as.integer(args[3]); n_cores <- as.integer(args[4])
.libPaths(c(lib, .libPaths())); suppressMessages(library(mdir)); options(mdir.quiet = TRUE)
source("compare.R")
dp <- as.numeric(densityPrior(scale_pool_shape = 0))
cpp_marg <- mdir:::collapsedLogMarginalCpp; sm_only <- mdir:::splitMergeOnlyCpp
set.seed(20260414)
N <- 7; K <- 3; w <- c(0.5, 0.3, 0.2)

make_case <- function(type, P) {
  if (type == 2) X <- matrix(sample(0:2, N * P, TRUE, c(.5, .3, .2)), N, P)
  else X <- matrix(rnorm(N * P, rep(c(-1.5, 1.5), c(4, 3))), N, P)
  X
}
exact <- function(X, type, beta, fixed, lab_fixed) {
  free <- which(fixed == 0)
  grid <- as.matrix(expand.grid(rep(list(0:(K - 1)), length(free))))
  lp <- numeric(nrow(grid))
  for (i in seq_len(nrow(grid))) {
    lab <- lab_fixed; lab[free] <- grid[i, ]
    s <- sum(log(w[lab + 1]))
    for (k in 0:(K - 1)) { rows <- which(lab == k) - 1; if (length(rows)) s <- s + cpp_marg(X, K, type, dp, rows, beta) }
    lp[i] <- s
  }
  p <- exp(lp - max(lp)); list(grid = grid, prob = p / sum(p), free = free)
}
run_case <- function(name, type, P, beta, n_fixed = 0, ref_beta = beta) {
  X <- make_case(type, P)
  fixed <- c(rep(1, n_fixed), rep(0, N - n_fixed)); lab_fixed <- c(seq_len(n_fixed) - 1, rep(0, N - n_fixed))
  ex <- exact(X, type, ref_beta, fixed, lab_fixed)
  chains <- parallel::mclapply(seq_len(n_chains), function(i) {
    set.seed(8000 + i)
    init <- sample(0:(K - 1), N, TRUE); init[fixed == 1] <- lab_fixed[fixed == 1]
    out <- sm_only(X, K, type, dp, init, fixed, w, n_iter, beta, 0L, 0.1, integer(0))
    codes <- config_code(out$labels[-(1:500), ex$free, drop = FALSE], K)
    list(codes = codes, acc = out$acceptance)
  }, mc.cores = n_cores)
  fr <- freq_by_chain(lapply(chains, `[[`, "codes"), nrow(ex$grid))
  stopifnot(all(config_code(ex$grid, K) == seq_len(nrow(ex$grid)) - 1))
  res <- compare_to_exact(fr, ex$prob)
  cat(sprintf("%-52s TV %.4f (noise ~%.4f) | states %d | max|z| %.2f | mean z^2 %.2f (expect ~%.2f) | acceptance %.2f\n",
    name, res$tv, res$tv_chain_noise, res$n_states, res$max_abs_z, res$mean_z2, res$expected_mean_z2, mean(vapply(chains, `[[`, 1, "acc"))))
}
run_case("C (3 categories, 2 measurements), beta = 1", 2, 2, 1)
run_case("G (P = 2), beta = 1", 0, 2, 1)
run_case("MVN (P = 2), beta = 1", 1, 2, 1)
run_case("MVN (P = 2), beta = 0.5 (tempered)", 1, 2, 0.5)
run_case("C, beta = 1, items 1-2 observed in classes 0, 1", 2, 2, 1, n_fixed = 2)
run_case("MVN (P = 2), beta = 1, items 1-2 observed", 1, 2, 1, n_fixed = 2)
run_case("CONTROL: MVN beta = 0.5 chain vs exact beta = 1", 1, 2, 0.5, ref_beta = 1)
