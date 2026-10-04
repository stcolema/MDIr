# Pool many independent plain chains (equal weight per chain, as in consensus
# clustering of chain outputs) and compare the pooled pattern frequencies with
# the exact posterior masses of mode_mass.R.
args <- commandArgs(TRUE)
lib <- args[1]; rds <- args[2]; n_chains <- as.integer(args[3]); R <- as.integer(args[4]); n_cores <- as.integer(args[5])
.libPaths(c(lib, .libPaths())); suppressMessages(library(mdir)); options(mdir.quiet = TRUE)
source("bimodal_setup.R")
mm <- readRDS(rds); exact <- mm$exact
d <- make_data(sep = mm$sep)
dp <- densityPrior(scale_pool_shape = 0)
one <- function(i) {
  set.seed(7000 + i)
  f <- callMDI(list(d$X), R = R, thin = 10, types = "MVN", K = 3, split_merge = 2L, check_prior = FALSE, save_parameters = FALSE, density_prior = dp)
  pat <- apply(f$allocations[-(1:50), , 1], 1, group_partition, truth = d$truth)
  as.numeric(table(factor(pat, levels = names(exact)))) / length(pat)
}
F <- do.call(cbind, parallel::mclapply(seq_len(n_chains), one, mc.cores = n_cores))
m <- rowMeans(F)
cat(sprintf("%d independent plain chains of %d sweeps, pooled with equal weights\n", n_chains, R))
cat(sprintf("TV to exact posterior masses: %.3f\n", 0.5 * sum(abs(m - exact))))
print(data.frame(pattern = names(exact), exact = round(exact, 3), pooled = round(m, 3),
  se = round(apply(F, 1, sd) / sqrt(ncol(F)), 3))[order(-exact), ][1:5, ], row.names = FALSE)
cat("share of chains whose dominant pattern is the exact-modal one:", mean(apply(F, 2, which.max) == which.max(exact)), "\n")
