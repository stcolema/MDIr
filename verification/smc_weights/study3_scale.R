# Study 3: do SMC weights remain usable with many clusters? 30 true clusters shared by three
# views, K = 45 components per view. There is no exact reference; the question is whether the
# weights are usable at all: how many temperatures the path needs, how far the effective sample
# size falls along it, how much the evidence estimate varies between independent runs (the
# variance of the weights), and what the weighted ensemble finds against the truth.
#   Rscript study3_scale.R <lib> <n_particles> <n_runs> <cores> <out.rds>
args <- commandArgs(TRUE)
lib <- args[1]; P <- as.integer(args[2]); n_runs <- as.integer(args[3]); n_cores <- as.integer(args[4]); out_rds <- args[5]
.libPaths(c(lib, .libPaths())); suppressMessages(library(mdir)); options(mdir.quiet = TRUE)
source("../joint_allocation/setup.R")
d <- make_large(); V <- 3
t0 <- Sys.time()
set.seed(21)
pilot <- smcMDI(d$X, d$types, K = d$K, n_particles = P, cess = 0.9, prior = d$prior, check_prior = FALSE,
                joint_allocation = V, split_merge = 2L, sweeps_per_step = 2L)
cat(sprintf("pilot: %d temperatures, %.0f s, min ESS %.1f of %d, final ESS %.1f, forced final step: %s\n",
            nrow(pilot$trace), as.numeric(Sys.time() - t0, units = "secs"), min(pilot$trace$ess), P,
            1 / sum(pilot$particle_weights^2), pilot$forced_final_step))
print(pilot$trace[round(seq(1, nrow(pilot$trace), length.out = 10)), c("step", "beta", "ess", "cess", "resampled")], row.names = FALSE)
rr <- smcReplicates(d$X, d$types, n_runs = n_runs, n_particles = P, betas = pilot$trace$beta, K = d$K,
                    n_cores = n_cores, prior = d$prior, check_prior = FALSE, joint_allocation = V, split_merge = 2L,
                    sweeps_per_step = 2L, resample_threshold = 0.5)
logz <- vapply(rr$runs, `[[`, numeric(1), "log_evidence")
ess <- vapply(rr$runs, function(r) 1 / sum(r$particle_weights^2), numeric(1))
cat("\nlog evidence per run:", round(logz, 1), "\n")
cat("sd of log evidence:", round(sd(logz), 2), "; final ESS per run:", round(ess, 1), "of", P, "\n")
cat("evidence weights of the runs:", round(exp(logz - max(logz)) / sum(exp(logz - max(logz))), 3), "\n")
occ <- function(fit, v = 1) vapply(seq_len(fit$n_particles), function(p) length(unique(fit$allocations[1, , v, p])), numeric(1))
for (i in seq_along(rr$runs)) {
  r <- rr$runs[[i]]
  w <- r$particle_weights
  a1 <- vapply(seq_len(r$n_particles), function(p) ari(r$allocations[1, , 1, p], d$truth[[1]]), numeric(1))
  cat(sprintf("run %d: weighted ARI view 1 %.3f (unweighted %.3f), weighted occupied %.1f (unweighted %.1f), distinct ancestors %d\n",
              i, sum(w * a1), mean(a1), sum(w * occ(r)), mean(occ(r)), length(unique(r$root))))
}
saveRDS(list(pilot = pilot$trace, logz = logz, ess = ess), out_rds)
cat(sprintf("total %.0f s\n", as.numeric(Sys.time() - t0, units = "secs")))
