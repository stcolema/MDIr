# Does the joint allocation reduce the chains' disagreement? Repeats a multi-chain run for several
# seeds with and without the joint allocation, and reports, per run, the largest Rhat among the
# quantities that describe the partitions and their association, and the spread across chains of
# the mean number of occupied components (a chain in a different mode of the partition shows up
# as a spread).
#   Rscript mixing.R <lib> <scenario: small|large> <n_seeds> <R> <cores>
args <- commandArgs(TRUE)
lib <- args[1]; scenario <- args[2]; n_seeds <- as.integer(args[3]); R <- as.integer(args[4])
n_cores <- as.integer(args[5])
.libPaths(c(lib, .libPaths())); suppressMessages(library(mdir)); options(mdir.quiet = TRUE)
source("setup.R")
d <- if (scenario == "small") make_small() else make_large()
V <- length(d$X)
configs <- list(
  plain = list(joint_allocation = 0L),
  joint_all = list(joint_allocation = V),
  split_merge = list(joint_allocation = 0L, split_merge = 2L),
  joint_all_split_merge = list(joint_allocation = V, split_merge = 2L)
)
one <- function(job) {
  cfg <- configs[[job$cfg]]
  set.seed(job$seed)
  t0 <- Sys.time()
  f <- do.call(fitMDI, c(list(X = d$X, n_chains = 4, R = R, thin = max(1, R %/% 600), types = d$types,
                              K = d$K, prior = d$prior, burn = R / 2, verbose = FALSE), cfg))
  cv <- as.data.frame(unclass(attr(f, "convergence")))
  # the two associated views, whose shared structure is the object of interest, and the rest
  informative <- cv[cv$quantity %in% c("occupied_components[1]", "occupied_components[2]", "phi[1,2]",
                                       "agreement[1,2]"), ]
  keep <- seq(floor((R / 2) / f[[1]]$thin) + 2, nrow(f[[1]]$mass))
  occ <- vapply(1:V, function(v) vapply(unclass(f)[seq_along(f)], function(ch)
    mean(apply(ch$allocations[keep, , v], 1, function(z) length(unique(z)))), numeric(1)), numeric(4))
  data.frame(cfg = job$cfg, seed = job$seed, rhat_views12 = max(informative$rhat),
             rhat_view3_occupied = cv$rhat[cv$quantity == "occupied_components[3]"],
             rhat_all = max(cv$rhat),
             chains_with_extra_cluster_v1 = sum(occ[, 1] > 3.8), chains_with_extra_cluster_v2 = sum(occ[, 2] > 3.6),
             seconds = as.numeric(Sys.time() - t0, units = "secs"))
}
jobs <- expand.grid(cfg = names(configs), seed = 1000 + seq_len(n_seeds), stringsAsFactors = FALSE)
jobs <- split(jobs, seq_len(nrow(jobs)))
res <- do.call(rbind, parallel::mclapply(jobs, one, mc.cores = n_cores))
print(res, row.names = FALSE)
cat("\nBy configuration: runs whose views-1-2 quantities have Rhat > 1.05; runs in which some chain sits in\n",
    "a mode with an extra cluster in view 1 or 2; runs whose view-3 occupancy has Rhat > 1.05; median seconds\n", sep = "")
print(do.call(rbind, lapply(split(res, res$cfg), function(r)
  data.frame(cfg = r$cfg[1], runs = nrow(r), rhat12_fail = sum(r$rhat_views12 > 1.05),
             any_chain_extra_cluster = sum(r$chains_with_extra_cluster_v1 + r$chains_with_extra_cluster_v2 > 0),
             chains_extra_cluster = sum(r$chains_with_extra_cluster_v1 + r$chains_with_extra_cluster_v2),
             view3_fail = sum(r$rhat_view3_occupied > 1.05), median_seconds = median(r$seconds)))), row.names = FALSE)
saveRDS(res, paste0("mixing_", scenario, ".rds"))
