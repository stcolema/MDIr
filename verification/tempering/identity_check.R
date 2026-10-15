# Compare the sampler built before and after the tempering code at betas = 1, one replica.
# Install the two versions into different libraries, then
#   Rscript identity_check.R <lib_before> before.rds; Rscript identity_check.R <lib_after> after.rds
#   Rscript -e 'a <- readRDS("before.rds"); b <- readRDS("after.rds"); print(sapply(seq_along(a), function(i) identical(a[[i]], b[[i]])))'
# (all TRUE when run for this change: G/G, MVN/MVN, G, MVN/G/MVN with missing values, TAGM/G, C/C)
args <- commandArgs(TRUE)
lib <- args[1]; out <- args[2]
.libPaths(c(lib, .libPaths()))
library(mdir)
options(mdir.quiet = TRUE)
res <- list()
mk <- function(seed, N = 40, L = 2) {
  set.seed(seed)
  z <- sample(1:3, N, TRUE)
  X <- lapply(seq_len(L), function(l) { m <- matrix(rnorm(N * 2, z * 2), N, 2); rownames(m) <- seq_len(N); m })
  X
}
cfg <- list(
  list(types = c("G", "G"), L = 2, K = c(4, 4)),
  list(types = c("MVN", "MVN"), L = 2, K = c(3, 3)),
  list(types = c("G"), L = 1, K = 4),
  list(types = c("MVN", "G", "MVN"), L = 3, K = c(3, 4, 3)),
  list(types = c("TAGM", "G"), L = 2, K = c(3, 3))
)
for (i in seq_along(cfg)) {
  cf <- cfg[[i]]
  X <- mk(100 + i, L = cf$L)
  if (i == 4) X[[1]][c(3, 9), 1] <- NA
  set.seed(7 + i)
  fit <- callMDI(X, R = 60, thin = 2, types = cf$types, K = cf$K, save_pointwise = TRUE)
  res[[i]] <- fit[c("allocations", "phis", "weights", "mass", "complete_likelihood", "observed_likelihood", "joint_likelihood", "parameters", "pooled_hyperparameters", "imputed", "outliers", "evidence")]
}
# categorical
set.seed(5)
Xc <- list({ m <- matrix(sample(0:2, 60 * 2, TRUE), 60, 2); rownames(m) <- 1:60; m }, { m <- matrix(sample(0:1, 60 * 3, TRUE), 60, 3); rownames(m) <- 1:60; m })
set.seed(3)
fit <- callMDI(Xc, R = 60, thin = 2, types = c("C", "C"), K = c(3, 3))
res[[6]] <- fit[c("allocations", "phis", "weights", "mass", "complete_likelihood", "parameters")]
saveRDS(res, out)
