# Split-merge with an outlier component (TAGM-type MVN view with a global multivariate t): the move
# alone, weights and outlier weight fixed, must leave the joint law of (label, outlier flag) of the free
# items, proportional to
#   prod_n w[c_n] * (o_n = 0: 1 - eps ; o_n = 1: eps * l_out(n)) * prod_k m(non-outlier members of k),
# invariant. Exact enumeration over (K * 2)^N_free states.
args <- commandArgs(TRUE)
lib <- args[1]; n_chains <- as.integer(args[2]); n_iter <- as.integer(args[3]); n_cores <- as.integer(args[4])
.libPaths(c(lib, .libPaths())); suppressMessages(library(mdir)); options(mdir.quiet = TRUE)
source("compare.R")
dp <- as.numeric(densityPrior(scale_pool_shape = 0))
cpp_marg <- mdir:::collapsedLogMarginalCpp; sm_only <- mdir:::splitMergeOnlyCpp
set.seed(20260416)
N <- 5; K <- 3; P <- 2; w <- c(0.5, 0.3, 0.2)
X <- rbind(matrix(rnorm(8, -1.5, 0.6), 4, 2), matrix(rnorm(6, 1.5, 0.6), 3, 2))[1:N, ]
X[5, ] <- c(6, -5)   # a gross outlier
X <- X[1:N, ]

run_case <- function(name, eps, ref_eps = eps, n_fixed = 0) {
  fixed <- c(rep(1, n_fixed), rep(0, N - n_fixed)); lab_fixed <- c(seq_len(n_fixed) - 1, rep(0, N - n_fixed))
  lo <- sm_only(X, K, 1, dp, lab_fixed, fixed, w, 0, 1, 1L, eps, integer(0))$outlier_loglik
  free <- which(fixed == 0)
  grid <- as.matrix(expand.grid(rep(list(0:(2 * K - 1)), length(free))))   # state = label + K * flag
  lp <- apply(grid, 1, function(g) {
    lab <- lab_fixed; flag <- rep(0, N)
    lab[free] <- g %% K; flag[free] <- g %/% K
    s <- sum(log(w[lab + 1])) + sum(ifelse(flag == 1, log(ref_eps) + lo, log(1 - ref_eps)))
    for (k in 0:(K - 1)) {
      rows <- which(lab == k & flag == 0) - 1
      if (length(rows)) s <- s + cpp_marg(X, K, 1, dp, rows, 1)
    }
    s
  })
  p <- exp(lp - max(lp)); p <- p / sum(p)
  chains <- parallel::mclapply(seq_len(n_chains), function(i) {
    set.seed(9000 + i)
    init <- sample(0:(K - 1), N, TRUE); init[fixed == 1] <- lab_fixed[fixed == 1]
    out <- sm_only(X, K, 1, dp, init, fixed, w, n_iter, 1, 1L, eps, sample(0:1, N, TRUE))
    st <- out$labels[-(1:500), free, drop = FALSE] + K * out$outliers[-(1:500), free, drop = FALSE]
    list(codes = as.vector(st %*% (2 * K)^(seq_along(free) - 1)), acc = out$acceptance)
  }, mc.cores = n_cores)
  fr <- freq_by_chain(lapply(chains, `[[`, "codes"), nrow(grid))
  res <- compare_to_exact(fr, p)
  cat(sprintf("%-50s TV %.4f (noise ~%.4f) | states %d | max|z| %.2f | mean z^2 %.2f (expect ~%.2f) | acc %.2f | outlier share %.3f\n",
    name, res$tv, res$tv_chain_noise, res$n_states, res$max_abs_z, res$mean_z2, res$expected_mean_z2,
    mean(vapply(chains, `[[`, 1, "acc")), sum(p * rowMeans(grid %/% K))))
}
run_case("MVN + t outlier, eps = 0.1", 0.1)
run_case("MVN + t outlier, eps = 0.4", 0.4)
run_case("MVN + t outlier, eps = 0.1, items 1-2 observed", 0.1, n_fixed = 2)
run_case("CONTROL: chain eps = 0.1 vs exact eps = 0.4", 0.1, ref_eps = 0.4)
