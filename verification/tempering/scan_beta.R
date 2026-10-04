.libPaths(c("/tmp/claude-0/-home-user-MDIr/7d714999-05c1-566d-84c6-cb4710ea967b/scratchpad/lib_new", .libPaths())); suppressMessages(library(mdir)); options(mdir.quiet = TRUE)
source("bimodal_setup.R")
args <- commandArgs(TRUE); sep <- as.numeric(args[1])
d <- make_data(sep = sep)
bs <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
res <- parallel::mclapply(bs, function(b) {
  set.seed(500 + round(b * 100))
  f <- callMDI(list(d$X), R = 4000, thin = 1, types = "MVN", K = 3, betas = b, check_prior = FALSE, save_parameters = FALSE)
  ll <- f$complete_likelihood[-(1:1000)]
  occ <- mean(apply(f$allocations[-(1:1000), , 1][seq(1, 3000, 30), ], 1, function(l) length(unique(l))))
  c(mean(ll), sd(ll), occ)
}, mc.cores = 3)
cat("sep", sep, "\n"); print(round(cbind(beta = bs, ll_mean = sapply(res, `[`, 1), ll_sd = sapply(res, `[`, 2), occupied = sapply(res, `[`, 3)), 2), row.names = FALSE)
