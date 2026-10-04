# Posterior mass of each merge pattern in a genuinely multimodal problem, exact
# against plain chains and parallel tempering.
#
# Four well separated Gaussian clusters, K = 3 components, so one merge is forced.
# With the variance scale held fixed (no pooling) the collapsed marginal likelihood
# of a block of items is the closed-form normal-inverse-Wishart one, so the
# posterior mass of each partition of the four true groups into at most three
# blocks can be computed exactly, treating each true group as a unit. That is
# exact up to the (negligible, checked) probability that an item leaves its
# true group, so it is an independent reference for the mass of every mode.
args <- commandArgs(TRUE)
lib <- args[1]; sep <- as.numeric(args[2]); n_runs <- as.integer(args[3]); R <- as.integer(args[4])
n_cores <- as.integer(args[5]); out_rds <- args[6]
.libPaths(c(lib, .libPaths())); suppressMessages(library(mdir)); options(mdir.quiet = TRUE)
source("bimodal_setup.R"); source("exact_reference.R")
d <- make_data(sep = sep); X <- d$X; N <- nrow(X); K <- 3; P <- 2
dp <- densityPrior(scale_pool_shape = 0)

hp <- mdir:::densityHyperparameters(X, K, 1, as.numeric(dp))
# closed-form log marginal likelihood of a set of items under the NIW prior
log_marg_niw <- function(Xs) {
  n <- nrow(Xs); if (n == 0) return(0)
  xb <- colMeans(Xs); S <- crossprod(sweep(Xs, 2, xb))
  kn <- hp$kappa + n; nun <- hp$nu + n
  Psin <- hp$scale + S + hp$kappa * n / kn * tcrossprod(xb - hp$xi)
  lmvgamma <- function(a) P * (P - 1) / 4 * log(pi) + sum(lgamma(a + (1 - seq_len(P)) / 2))
  -n * P / 2 * log(pi) + lmvgamma(nun / 2) - lmvgamma(hp$nu / 2) +
    hp$nu / 2 * determinant(hp$scale, logarithm = TRUE)$modulus - nun / 2 * determinant(Psin, logarithm = TRUE)$modulus +
    P / 2 * (log(hp$kappa) - log(kn))
}
# all set partitions of the four groups with at most K blocks
parts <- list(c(1,1,1,1), c(1,1,1,2), c(1,1,2,1), c(1,2,1,1), c(1,2,2,2), c(1,1,2,2), c(1,2,1,2), c(1,2,2,1),
  c(1,1,2,3), c(1,2,1,3), c(1,2,3,1), c(1,2,2,3), c(1,2,3,2), c(1,2,3,3))
pat_name <- function(p) paste(vapply(split(1:4, p), paste, character(1), collapse = ""), collapse = "|")
exact <- vapply(parts, function(p) {
  # sum over the labelled assignments of blocks to components (K! / (K - b)!)
  b <- max(p)
  blocks <- lapply(1:b, function(j) which(truth_rows <- d$truth %in% which(p == j)))
  lm <- sum(vapply(blocks, function(ix) log_marg_niw(X[ix, , drop = FALSE]), numeric(1)))
  counts <- vapply(blocks, length, integer(1))
  perms <- if (b == 1) matrix(0:(K - 1), ncol = 1) else as.matrix(expand.grid(rep(list(0:(K - 1)), b)))
  perms <- perms[apply(perms, 1, function(r) !anyDuplicated(r)), , drop = FALSE]
  lp <- vapply(seq_len(nrow(perms)), function(i) {
    Nk <- integer(K); Nk[perms[i, ] + 1] <- counts
    log(prior_counts_L1(Nk, K))
  }, numeric(1))
  matrixStats_logsum <- function(x) max(x) + log(sum(exp(x - max(x))))
  matrixStats_logsum(lp) + lm
}, numeric(1))
exact <- exp(exact - max(exact)); exact <- exact / sum(exact)
names(exact) <- vapply(parts, pat_name, "")
cat("exact pattern masses (sep =", sep, "):\n"); print(round(sort(exact, decreasing = TRUE), 4))

pattern_freq <- function(f) {
  a <- f$allocations[-(1:50), , 1]
  pat <- apply(a, 1, group_partition, truth = d$truth)
  tab <- table(factor(pat, levels = names(exact)))
  as.numeric(tab) / sum(tab)
}
settings <- readRDS(paste0(out_rds, ".ladder.rds"))
run_plain <- function(i) { set.seed(900 + i); f <- callMDI(list(X), R = R, thin = 10, types = "MVN", K = K, check_prior = FALSE,
  save_parameters = FALSE, density_prior = dp); pattern_freq(f) }
run_pt <- function(i) { set.seed(1900 + i); f <- callMDI(list(X), R = R, thin = 10, types = "MVN", K = K, check_prior = FALSE,
  save_parameters = FALSE, density_prior = dp, betas = settings$betas)
  list(freq = pattern_freq(f), diag = mdir:::ptDiagnostics(f)) }
t0 <- Sys.time()
plain <- parallel::mclapply(seq_len(n_runs), run_plain, mc.cores = n_cores)
cat("plain runs done:", format(Sys.time() - t0), "\n")
pt <- parallel::mclapply(seq_len(n_runs), run_pt, mc.cores = n_cores)
cat("PT runs done:", format(Sys.time() - t0), "\n")
saveRDS(list(exact = exact, plain = plain, pt = pt, betas = settings$betas, sep = sep, R = R), out_rds)

summarise <- function(fl, name) {
  F <- do.call(cbind, fl)
  m <- rowMeans(F); se <- apply(F, 1, sd) / sqrt(ncol(F))
  cat(sprintf("\n%s: TV to exact %.3f; between-run sd of the total-variation distance %.3f\n", name, 0.5 * sum(abs(m - exact)),
    sd(apply(F, 2, function(f) 0.5 * sum(abs(f - exact))))))
  tab <- data.frame(pattern = names(exact), exact = round(exact, 4), estimate = round(m, 4), se = round(se, 4),
    z = round((m - exact) / pmax(se, 1e-6), 1))
  print(tab[order(-tab$exact), ][1:8, ], row.names = FALSE)
  cat("runs with a single dominant pattern (>95% of draws):", sum(apply(F, 2, max) > 0.95), "of", ncol(F), "\n")
}
summarise(plain, "Plain chains")
summarise(lapply(pt, `[[`, "freq"), "Parallel tempering (cold chain)")
cat("\nPT diagnostics (first run):\n"); print(pt[[1]]$diag)
