args <- commandArgs(TRUE); lib <- args[1]; sep <- as.numeric(args[2]); out <- args[3]
.libPaths(c(lib, .libPaths())); suppressMessages(library(mdir)); options(mdir.quiet = TRUE)
source("bimodal_setup.R")
d <- make_data(sep = sep)
dp <- densityPrior(scale_pool_shape = 0)
# the transition between the one-component and the three-component phase sits
# near beta = 0.6 - 0.85 for these separations, so start with a ladder that is
# dense there
b0 <- c(0.05, 0.3, 0.55, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1)
set.seed(42)
ad <- adaptLadder(list(d$X), "MVN", b0, rounds = 5, R_pilot = 2000, K = 3, density_prior = dp)
for (i in seq_along(ad$history)) cat(sprintf("round %d rejection: %s\n", i, paste(round(ad$history[[i]]$rejection_rate, 2), collapse = " ")))
cat("ladder:", round(ad$betas, 3), "\n")
saveRDS(ad, out)
