# TAGM (MVN + global t outlier component), full sampler: pairwise co-clustering and per-item outlier
# probabilities of plain chains vs chains with split-merge; both must target the same posterior.
# Two settings: one TAGM view, and a TAGM view coupled to an MVN view (MDI).
args <- commandArgs(TRUE)
lib <- args[1]; n_chains <- as.integer(args[2]); R <- as.integer(args[3]); n_cores <- as.integer(args[4])
.libPaths(c(lib, .libPaths())); suppressMessages(library(mdir)); options(mdir.quiet = TRUE)
set.seed(20260417)
N <- 9; P <- 2
X <- matrix(rnorm(N * P, rep(c(-1.5, 1.5), c(5, 4)), 0.7), N, P); rownames(X) <- seq_len(N)
X[9, ] <- c(5, -4); X[4, ] <- c(-4.5, 4)
X2 <- matrix(rnorm(N * P, rep(c(-1.5, 1.5), c(5, 4)), 0.7), N, P); rownames(X2) <- seq_len(N)
pairs <- t(combn(N, 2))
summ <- function(f, view) {
  a <- f$allocations[-(1:500), , view]
  co <- colMeans(apply(pairs, 1, function(p) a[, p[1]] == a[, p[2]]))
  o <- colMeans(f$outliers[-(1:500), , view])
  c(co, o)
}
run <- function(Xs, types, sm, seed, view) {
  set.seed(seed)
  f <- callMDI(Xs, R = R, thin = 1, types = types, K = rep(3, length(Xs)), split_merge = sm, check_prior = FALSE, save_parameters = FALSE)
  summ(f, view)
}
est <- function(Xs, types, sm, base) do.call(cbind, parallel::mclapply(seq_len(n_chains), function(i) run(Xs, types, sm, base + i, 1), mc.cores = n_cores))
for (setting in list(list("one TAGM view", list(X), "TAGM"), list("TAGM + MVN views (view 1 reported)", list(X, X2), c("TAGM", "MVN")))) {
  A <- est(setting[[2]], setting[[3]], 0L, 100); B <- est(setting[[2]], setting[[3]], 2L, 200)
  d <- rowMeans(A) - rowMeans(B); se <- sqrt(apply(A, 1, var) / n_chains + apply(B, 1, var) / n_chains)
  z <- d / pmax(se, 1e-8)
  cat(sprintf("%-36s %d chains x %d sweeps; %d co-clustering + %d outlier probabilities: max |diff| %.4f; max |z| %.2f; mean z^2 %.2f; outlier prob range %.2f to %.2f\n",
    setting[[1]], n_chains, R, nrow(pairs), N, max(abs(d)), max(abs(z)), mean(z^2), min(rowMeans(A)[-seq_len(nrow(pairs))]), max(rowMeans(A)[-seq_len(nrow(pairs))])))
}
