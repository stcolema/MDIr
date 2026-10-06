# Can weights from a sequential Monte Carlo ensemble recover the posterior mass of modes?
# Study 1: a problem whose mode masses are known exactly.
#
# Four well-separated Gaussian clusters fitted with K = 3 components (one merge is forced, and
# which pair merges is a mode of the partition posterior that Gibbs chains rarely leave; exact
# masses of the merge patterns from verification/tempering/mode_mass.R). Compared, for the
# same data and replicated runs:
#   * ensembles of independent plain chains with equal weights, as in consensus clustering
#     (the weight of a mode is its basin of attraction), by depth and width;
#   * weighted ensembles from smcMDI() on a schedule fixed by a pilot run, with more sweeps per
#     temperature step ("deeper" moves), with and without split-merge moves;
#   * pooling independent SMC runs by their evidence estimates (the consistent estimator).
#   Rscript study1_exact.R <lib> <n_rep> <cores> <out.rds>
args <- commandArgs(TRUE)
lib <- args[1]; n_rep <- as.integer(args[2]); n_cores <- as.integer(args[3]); out_rds <- args[4]
.libPaths(c(lib, .libPaths())); suppressMessages(library(mdir)); options(mdir.quiet = TRUE)
source("../tempering/bimodal_setup.R")
exact <- readRDS("../tempering/mode_mass_sep2.5_results.rds")$exact
d <- make_data(sep = 2.5); dp <- densityPrior(scale_pool_shape = 0)
K <- 3
top4 <- names(sort(exact, decreasing = TRUE))[1:4]
tv <- function(est) 0.5 * sum(abs(est[names(exact)] - exact))
as_vec <- function(x) { v <- setNames(rep(0, length(exact)), names(exact)); v[names(x)] <- x; v }

# the schedule from a pilot run, shared by all SMC runs
set.seed(5)
pilot <- smcMDI(list(d$X), "MVN", K = K, n_particles = 300, cess = 0.9, density_prior = dp, check_prior = FALSE)
betas <- pilot$trace$beta
cat("pilot schedule:", length(betas), "temperatures\n")

chain_patterns <- function(f, last_only) {
  a <- f$allocations[, , 1]
  a <- if (last_only) a[nrow(a), , drop = FALSE] else a[-(1:20), , drop = FALSE]
  apply(a, 1, group_partition, truth = d$truth)
}
ensemble_run <- function(job) {
  set.seed(job$seed)
  pats <- lapply(seq_len(job$W), function(i) {
    f <- callMDI(list(d$X), R = job$D, thin = max(1, job$D %/% 100), types = "MVN", K = K, check_prior = FALSE,
                 save_parameters = FALSE, density_prior = dp, split_merge = job$sm)
    list(last = chain_patterns(f, TRUE), pooled = chain_patterns(f, FALSE))
  })
  fr <- function(x) as_vec(table(factor(unlist(x), levels = names(exact))) / length(unlist(x)))
  data.frame(method = "ensemble", sm = job$sm, D = job$D, W = job$W, sweeps = job$W * job$D, seed = job$seed,
             tv_last = tv(fr(lapply(pats, `[[`, "last"))), tv_pooled = tv(fr(lapply(pats, `[[`, "pooled"))))
}
res_ens <- NULL
jobs <- do.call(c, lapply(c(0L, 2L), function(sm) lapply(c(250L, 1000L, 4000L), function(D)
  lapply(seq_len(n_rep), function(r) list(sm = sm, D = D, W = 100, seed = 7000 + r)))))
jobs <- unlist(jobs, recursive = FALSE)
res_ens <- do.call(rbind, parallel::mclapply(jobs, ensemble_run, mc.cores = n_cores))
cat("\nEQUAL-WEIGHT ENSEMBLES of W = 100 plain chains (basin weights), TV to the exact masses\n")
print(aggregate(cbind(tv_last, tv_pooled) ~ sm + D + sweeps, res_ens, function(x) c(mean = mean(x), sd = sd(x))), digits = 3)

smc_est <- function(fit) {
  pm <- smcPosterior(fit, function(l) group_partition(l, d$truth), view = 1)
  as_vec(pm$estimate)
}
smc_block <- function(sm, m, P, n_runs) {
  t0 <- Sys.time()
  rr <- smcReplicates(list(d$X), "MVN", n_runs = n_runs, n_particles = P, betas = betas, K = K,
                      n_cores = n_cores, density_prior = dp, check_prior = FALSE, split_merge = sm,
                      sweeps_per_step = m, resample_threshold = 0.5)
  est <- t(vapply(rr$runs, smc_est, numeric(length(exact))))
  pooled <- as_vec(smcPosterior(rr$combined, function(l) group_partition(l, d$truth), view = 1)$estimate)
  ess <- vapply(rr$runs, function(r) 1 / sum(r$particle_weights^2), numeric(1))
  logz <- vapply(rr$runs, `[[`, numeric(1), "log_evidence")
  list(row = data.frame(method = "SMC", sm = sm, sweeps_per_step = m, P = P, runs = n_runs,
                        sweeps_per_run = P * length(betas) * m,
                        tv_run_mean = mean(apply(est, 1, tv)), tv_run_sd = sd(apply(est, 1, tv)),
                        tv_pooled = tv(pooled), median_final_ess = median(ess), sd_log_evidence = sd(logz),
                        seconds = as.numeric(Sys.time() - t0, units = "secs")),
       top = rbind(exact = exact[top4], run_mean = colMeans(est)[top4], run_sd = apply(est, 2, sd)[top4],
                   pooled = pooled[top4]))
}
blocks <- list()
for (sm in c(0L, 2L)) for (m in c(1L, 4L, 16L)) {
  P <- if (m == 16L) 500L else 2000L
  b <- smc_block(sm, m, P, n_rep)
  blocks[[length(blocks) + 1]] <- b
  cat(sprintf("\nSMC: split_merge = %d, sweeps per step = %d, particles = %d (%d runs; %.0f s)\n", sm, m, P, n_rep, b$row$seconds))
  print(round(b$top, 3)); print(b$row[, -1], row.names = FALSE, digits = 3)
}
res_smc <- do.call(rbind, lapply(blocks, `[[`, "row"))
cat("\nSUMMARY OF SMC BLOCKS\n"); print(res_smc, row.names = FALSE, digits = 3)
saveRDS(list(ensemble = res_ens, smc = res_smc, top = lapply(blocks, `[[`, "top"), exact = exact, betas = betas), out_rds)
