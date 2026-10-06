# Study 2: SMC weights against equal-weight ensembles in a multi-view problem with aligned modes.
#
# The three-view data of the multimodal vignette. Plain chains sometimes sit in a mode in which a
# true cluster is split in two in both associated views (a "split mode"): a higher likelihood and more
# components but, on well-separated Gaussian clusters, almost no posterior mass. An equal-weight
# ensemble gives such a mode the weight of its basin of attraction. There is no exact reference here,
# so the posterior mass of the split mode is judged by agreement between methods that do not
# share the failure (tempered chains with joint allocation and split-merge moves, SMC), and by the
# agreement of independent replicates with each other.
#   Rscript study2_mdi.R <lib> <n_rep> <cores> <out.rds>
args <- commandArgs(TRUE)
lib <- args[1]; n_rep <- as.integer(args[2]); n_cores <- as.integer(args[3]); out_rds <- args[4]
.libPaths(c(lib, .libPaths())); suppressMessages(library(mdir)); options(mdir.quiet = TRUE)
source("../joint_allocation/setup.R")
d <- make_small(); V <- 3; X <- d$X; types <- d$types; K <- d$K; prior <- d$prior

# A view is in a split mode if it has more than three components that hold at least 8 items
n_big <- function(z) sum(table(z) >= 8)
split_mode <- function(lab) as.numeric(n_big(lab[, 1]) >= 4 || n_big(lab[, 2]) >= 4)
chain_last <- function(f) {
  a <- f$allocations[nrow(f$allocations), , ]
  split_mode(a)
}

# equal-weight ensembles of chains: the final state of each chain
ensemble <- function(job) {
  set.seed(job$seed)
  flags <- vapply(seq_len(job$W), function(i) {
    f <- callMDI(X, R = job$D, thin = max(1, job$D %/% 20), types = types, K = K, prior = prior,
                 check_prior = FALSE, save_parameters = FALSE, joint_allocation = job$joint)
    chain_last(f)
  }, numeric(1))
  data.frame(method = "ensemble, equal weights", joint = job$joint, D = job$D, W = job$W, seed = job$seed,
             p_split = mean(flags), sweeps = job$W * job$D)
}
jobs <- unlist(unlist(lapply(c(0L, V), function(j) lapply(c(2000L, 6000L), function(D)
  lapply(seq_len(n_rep), function(r) list(joint = j, D = D, W = 24, seed = 9000 + r)))), recursive = FALSE), recursive = FALSE)
res_ens <- parallel::mclapply(jobs, function(j) tryCatch(ensemble(j), error = function(e) conditionMessage(e)), mc.cores = n_cores)
if (any(vapply(res_ens, is.character, TRUE))) stop(unique(unlist(res_ens[vapply(res_ens, is.character, TRUE)])))
res_ens <- do.call(rbind, res_ens)
cat("\nEQUAL-WEIGHT ENSEMBLES: proportion of chains whose final state is in a split mode\n")
print(aggregate(p_split ~ joint + D + sweeps, res_ens, function(x) c(mean = mean(x), sd = sd(x))), digits = 3)

# the reference: tempered chains with joint allocation and split-merge moves
pt <- function(job) {
  set.seed(job$seed)
  f <- callMDI(X, R = 3000, thin = 10, types = types, K = K, prior = prior, check_prior = FALSE,
               save_parameters = FALSE, betas = ptLadder(4, 0.1), joint_allocation = V, split_merge = 2L)
  a <- f$allocations[-(1:151), , , drop = FALSE]
  data.frame(method = "PT + joint + split-merge", seed = job$seed, p_split_draws = mean(apply(a, 1, function(m) split_mode(matrix(m, nrow = dim(a)[2])))),
             final = chain_last(f))
}
res_pt <- do.call(rbind, parallel::mclapply(lapply(seq_len(max(8, n_rep)), function(r) list(seed = 9500 + r)), pt,
                                            mc.cores = n_cores))
cat("\nTEMPERED CHAINS (reference): proportion of post-burn draws in a split mode, per chain\n")
print(res_pt, row.names = FALSE, digits = 3)

# SMC on a schedule from a pilot
set.seed(11)
pilot <- smcMDI(X, types, K = K, n_particles = 100, cess = 0.9, prior = prior, check_prior = FALSE,
                joint_allocation = V)
betas <- pilot$trace$beta
cat("\npilot schedule:", length(betas), "temperatures\n")
smc_block <- function(joint, sm, m, P) {
  t0 <- Sys.time()
  rr <- smcReplicates(X, types, n_runs = n_rep, n_particles = P, betas = betas, K = K, n_cores = n_cores,
                      prior = prior, check_prior = FALSE, joint_allocation = joint, split_merge = sm,
                      sweeps_per_step = m, resample_threshold = 0.5)
  est <- vapply(rr$runs, function(r) smcPosterior(r, split_mode)$estimate["1"], numeric(1))
  est[is.na(est)] <- 0
  pooled <- smcPosterior(rr$combined, split_mode)$estimate["1"]
  if (is.na(pooled)) pooled <- 0
  ess <- vapply(rr$runs, function(r) 1 / sum(r$particle_weights^2), numeric(1))
  logz <- vapply(rr$runs, `[[`, numeric(1), "log_evidence")
  data.frame(method = "SMC", joint = joint, split_merge = sm, sweeps_per_step = m, P = P, runs = n_rep,
             p_split_run_mean = mean(est), p_split_run_sd = sd(est), p_split_pooled = unname(pooled),
             median_final_ess = median(ess), sd_log_evidence = sd(logz),
             seconds = as.numeric(Sys.time() - t0, units = "secs"))
}
res_smc <- do.call(rbind, list(smc_block(0L, 0L, 2L, 300L), smc_block(V, 0L, 2L, 300L), smc_block(V, 2L, 2L, 300L),
                               smc_block(V, 2L, 8L, 100L)))
cat("\nSMC (fixed pilot schedule), proportion of weighted particles in a split mode\n")
print(res_smc, row.names = FALSE, digits = 3)
saveRDS(list(ensemble = res_ens, pt = res_pt, smc = res_smc, betas = betas), out_rds)
