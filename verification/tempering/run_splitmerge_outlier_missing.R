# Outlier component with missing values. The state of the move includes the imputed values (the
# sampler draws them from the outlier law for flagged items and from the component otherwise); with the
# imputations held fixed the move must leave the joint law of (label, flag) given the imputed complete
# data invariant: proportional to
#   prod_n w[c_n] * (o_n = 0: 1 - eps ; o_n = 1: eps * t-density of the COMPLETE vector) * prod_k m(non-outlier members of k)
# (marginals from the hook's own imputed density). Exact enumeration.
args <- commandArgs(TRUE)
lib <- args[1]; n_chains <- as.integer(args[2]); n_iter <- as.integer(args[3]); n_cores <- as.integer(args[4])
.libPaths(c(lib, .libPaths())); suppressMessages(library(mdir)); options(mdir.quiet = TRUE)
source("compare.R")
dp <- as.numeric(densityPrior(scale_pool_shape = 0)); sm_only <- mdir:::splitMergeOnlyCpp
set.seed(20260418)
N <- 5; K <- 3; w <- c(0.5, 0.3, 0.2); eps <- 0.25
X <- rbind(matrix(rnorm(8, -1.5, 0.6), 4, 2), c(5, -4)); X <- rbind(X[1:3, ], matrix(rnorm(4, 1.5, 0.6), 2, 2)); X[5, ] <- c(5, -4)
X[2, 1] <- NA; X[4, 2] <- NA; X[5, 1] <- NA     # incl. the gross outlier
fixed <- rep(0, N)
set.seed(1)
info <- sm_only(X, K, 1, dp, rep(0, N), fixed, w, 0, 1, 1L, eps, integer(0), matrix(0, 0, 0))
lo <- info$outlier_loglik_complete; lm_sub <- info$subset_logmarg
cat("observed-only vs complete outlier log density:", round(info$outlier_loglik, 3), "|", round(lo, 3), "\n")
grid <- as.matrix(expand.grid(rep(list(0:(2 * K - 1)), N)))
lp <- apply(grid, 1, function(g) {
  lab <- g %% K; flag <- g %/% K
  s <- sum(log(w[lab + 1])) + sum(ifelse(flag == 1, log(eps) + lo, log(1 - eps)))
  for (k in 0:(K - 1)) {
    idx <- which(lab == k & flag == 0)
    if (length(idx)) s <- s + lm_sub[sum(2^(idx - 1)) + 1]
  }
  s
})
p <- exp(lp - max(lp)); p <- p / sum(p)
chains <- parallel::mclapply(seq_len(n_chains), function(i) {
  set.seed(9500 + i)
  init <- sample(0:(K - 1), N, TRUE); flags <- sample(0:1, N, TRUE)
  out <- sm_only(X, K, 1, dp, init, fixed, w, n_iter, 1, 1L, eps, flags, info$X_imputed)  # same imputed values in every chain
  stopifnot(isTRUE(all.equal(out$X_imputed, info$X_imputed)), isTRUE(all.equal(out$subset_logmarg, lm_sub)))
  st <- out$labels[-(1:500), , drop = FALSE] + K * out$outliers[-(1:500), , drop = FALSE]
  list(codes = as.vector(st %*% (2 * K)^(seq_len(N) - 1)), acc = out$acceptance)
}, mc.cores = n_cores)
fr <- freq_by_chain(lapply(chains, `[[`, "codes"), nrow(grid))
res <- compare_to_exact(fr, p)
cat(sprintf("outlier + missing values: TV %.4f (noise ~%.4f) | states %d | max|z| %.2f | mean z^2 %.2f (expect ~%.2f) | acc %.2f | outlier share %.3f\n",
  res$tv, res$tv_chain_noise, res$n_states, res$max_abs_z, res$mean_z2, res$expected_mean_z2, mean(vapply(chains, `[[`, 1, "acc")),
  sum(p * rowMeans(grid %/% K))))
