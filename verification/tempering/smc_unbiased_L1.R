# Unbiasedness of the SMC / AIS estimators on a model whose posterior and evidence are exact.
#   E[Zhat] = Z                       (evidence, up to the normalisation of the prior)
#   E[Zhat * fhat(c)] = Z * pi(c)     (unnormalised posterior measure)
# for a fixed schedule; the ratio estimator is only consistent.
args <- commandArgs(TRUE)
lib <- args[1]; n_runs <- as.integer(args[2]); n_part <- as.integer(args[3]); n_cores <- as.integer(args[4])
.libPaths(c(lib, .libPaths())); suppressMessages(library(mdir)); options(mdir.quiet = TRUE)
source("exact_reference.R"); source("compare.R")
set.seed(20260410)
N <- 7; K <- 3
X <- matrix(c(sample(0:2, N, TRUE, c(.5, .3, .2))), N, 1); rownames(X) <- seq_len(N)
alpha <- mdir:::densityHyperparameters(X, K, 2, numeric(0))$concentration
ex <- exact_L1(X, K, alpha, 1)
n_cfg <- nrow(ex$grid)
betas <- exp(seq(log(0.002), 0, length.out = 25))
one <- function(seed, thr, scheme, schedule) {
  set.seed(seed)
  fit <- if (schedule == "fixed") {
    smcMDI(list(X), "C", K = K, n_particles = n_part, schedule = "fixed", betas = betas, resample_threshold = thr,
      resample = scheme, check_prior = FALSE)
  } else {
    smcMDI(list(X), "C", K = K, n_particles = n_part, schedule = "adaptive", cess = 0.9, resample_threshold = thr,
      resample = scheme, check_prior = FALSE)
  }
  codes <- config_code(t(fit$allocations[1, , 1, ]), K)
  f <- as.numeric(tapply(fit$particle_weights, factor(codes + 1, levels = seq_len(n_cfg)), sum))
  f[is.na(f)] <- 0
  list(logZ = fit$log_evidence, f = f, steps = nrow(fit$trace))
}
settings <- list(
  list("AIS (never resample), fixed", 0, "systematic", "fixed"),
  list("SMC systematic, fixed", 0.5, "systematic", "fixed"),
  list("SMC multinomial, fixed", 0.5, "multinomial", "fixed"),
  list("SMC systematic, every step, fixed", 1, "systematic", "fixed"),
  list("SMC systematic, ADAPTIVE schedule", 0.5, "systematic", "adaptive")
)
top <- order(-ex$prob)[1:12]
for (st in settings) {
  res <- parallel::mclapply(seq_len(n_runs), function(i) one(7000 + i, st[[2]], st[[3]], st[[4]]), mc.cores = n_cores)
  Zr <- vapply(res, function(r) exp(r$logZ - ex$log_evidence), numeric(1))
  Fz <- vapply(res, function(r) r$f * exp(r$logZ - ex$log_evidence), numeric(n_cfg))
  m <- rowMeans(Fz); se <- apply(Fz, 1, sd) / sqrt(n_runs)
  z <- (m[top] - ex$prob[top]) / se[top]
  cat(sprintf("%-36s E[Zhat]/Z = %.4f (se %.4f; z = %.2f) | E[Zhat fhat(c)]/Z vs pi(c) over top 12 configs: max|z| %.2f, mean z^2 %.2f | mean steps %.0f\n",
    st[[1]], mean(Zr), sd(Zr) / sqrt(n_runs), (mean(Zr) - 1) / (sd(Zr) / sqrt(n_runs)), max(abs(z)), mean(z^2),
    mean(vapply(res, `[[`, 1, "steps"))))
  # ratio estimator pooled by evidence
  a <- vapply(res, function(r) exp(r$logZ), numeric(1)); a <- a / sum(a)
  pooled <- as.numeric(vapply(res, `[[`, numeric(n_cfg), "f") %*% a)
  cat(sprintf("    pooled (evidence-weighted) estimate: TV to exact posterior %.4f\n", 0.5 * sum(abs(pooled - ex$prob))))
}
