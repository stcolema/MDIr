# The same comparison with many clusters: 30 true clusters shared by three views, K = 45 components
# per view. Reports, for each configuration and seed, the number of occupied components and the
# adjusted Rand index to the truth of the last draw of each chain (are the chains finding the
# structure, and do they agree?), and the Rhat of the quantities describing the associated views.
#   Rscript mixing_large.R <lib> <n_seeds> <R> <cores>
args <- commandArgs(TRUE)
lib <- args[1]; n_seeds <- as.integer(args[2]); R <- as.integer(args[3]); n_cores <- as.integer(args[4])
.libPaths(c(lib, .libPaths())); suppressMessages(library(mdir)); options(mdir.quiet = TRUE)
source("setup.R")
d <- make_large(); V <- 3
configs <- list(
  plain = list(joint_allocation = 0L),
  joint_all = list(joint_allocation = V),
  split_merge = list(joint_allocation = 0L, split_merge = 2L),
  joint_all_split_merge = list(joint_allocation = V, split_merge = 2L)
)
one <- function(job) {
  set.seed(job$seed)
  t0 <- Sys.time()
  f <- do.call(fitMDI, c(list(X = d$X, n_chains = 4, R = R, thin = max(1, R %/% 400), types = d$types,
                              K = d$K, prior = d$prior, burn = R / 2, verbose = FALSE), configs[[job$cfg]]))
  cv <- as.data.frame(unclass(attr(f, "convergence")))
  inf <- cv[cv$quantity %in% c("occupied_components[1]", "occupied_components[2]", "phi[1,2]", "agreement[1,2]",
                               "joint_likelihood"), ]
  ch <- unclass(f)[seq_along(f)]
  last <- nrow(ch[[1]]$mass)
  data.frame(cfg = job$cfg, seed = job$seed, rhat_key = max(inf$rhat),
             occupied_v1 = mean(vapply(ch, function(x) length(unique(x$allocations[last, , 1])), numeric(1))),
             ari_v1 = mean(vapply(ch, function(x) ari(x$allocations[last, , 1], d$truth[[1]]), numeric(1))),
             ari_v2 = mean(vapply(ch, function(x) ari(x$allocations[last, , 2], d$truth[[2]]), numeric(1))),
             joint_ll = mean(vapply(ch, function(x) x$joint_likelihood[last], numeric(1))),
             joint_ll_sd = sd(vapply(ch, function(x) x$joint_likelihood[last], numeric(1))),
             seconds = as.numeric(Sys.time() - t0, units = "secs"))
}
jobs <- expand.grid(cfg = names(configs), seed = 2000 + seq_len(n_seeds), stringsAsFactors = FALSE)
res <- do.call(rbind, parallel::mclapply(split(jobs, seq_len(nrow(jobs))), one, mc.cores = n_cores))
print(res, row.names = FALSE, digits = 3)
cat("\nMeans by configuration\n")
print(aggregate(cbind(rhat_key, occupied_v1, ari_v1, ari_v2, joint_ll, joint_ll_sd, seconds) ~ cfg, res, mean), digits = 4)
saveRDS(res, "mixing_large.rds")
