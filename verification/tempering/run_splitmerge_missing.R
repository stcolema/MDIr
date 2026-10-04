# Missing data (data augmentation): the split-merge move conditions on the current imputation.
# No closed-form reference exists, so compare pairwise co-clustering probabilities of a plain
# chain with those of a chain that also runs split-merge; both must target the same posterior.
args <- commandArgs(TRUE)
lib <- args[1]; n_chains <- as.integer(args[2]); R <- as.integer(args[3]); n_cores <- as.integer(args[4])
.libPaths(c(lib, .libPaths())); suppressMessages(library(mdir)); options(mdir.quiet = TRUE)
set.seed(20260415)
N <- 8; P <- 2
X <- matrix(rnorm(N * P, rep(c(-1.2, 1.2), each = 4)), N, P); rownames(X) <- seq_len(N)
X[2, 1] <- NA; X[5, 2] <- NA; X[7, 1] <- NA
pairs <- t(combn(N, 2))
run <- function(sm, seed) {
  set.seed(seed)
  f <- callMDI(list(X), R = R, thin = 1, types = "MVN", K = 3, split_merge = sm, check_prior = FALSE, save_parameters = FALSE)
  a <- f$allocations[-(1:500), , 1]
  colMeans(apply(pairs, 1, function(p) a[, p[1]] == a[, p[2]]))
}
est <- function(sm, base) do.call(cbind, parallel::mclapply(seq_len(n_chains), function(i) run(sm, base + i), mc.cores = n_cores))
A <- est(0L, 100); B <- est(2L, 200)
d <- rowMeans(A) - rowMeans(B); se <- sqrt(apply(A, 1, var) / n_chains + apply(B, 1, var) / n_chains)
z <- d / se
cat(sprintf("%d chains x %d sweeps each; %d pairwise co-clustering probabilities\n", n_chains, R, nrow(pairs)))
cat(sprintf("max |difference| %.4f; max |z| %.2f; mean z^2 %.2f (expect about 1 if the targets agree)\n", max(abs(d)), max(abs(z)), mean(z^2)))
cat(sprintf("co-clustering range (plain): %.3f to %.3f\n", min(rowMeans(A)), max(rowMeans(A))))
