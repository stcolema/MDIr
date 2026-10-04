args <- commandArgs(TRUE); lib <- args[1]
.libPaths(c(lib, .libPaths())); suppressMessages(library(mdir)); options(mdir.quiet = TRUE)
source("bimodal_setup.R")
mm <- readRDS("mode_mass_sep2.5_results.rds"); exact <- mm$exact
top4 <- names(sort(exact, decreasing = TRUE))[1:4]
d <- make_data(sep = 2.5); dp <- densityPrior(scale_pool_shape = 0)
cat("exact:", paste(top4, round(exact[top4], 3), collapse = "  "), "\n")
one <- function(label, bstar, P, start_sweeps, thr = 0.5, seed = 1) {
  set.seed(seed); t0 <- Sys.time()
  fit <- smcMDI(list(d$X), "MVN", K = 3, n_particles = P, cess = 0.9, resample_threshold = thr, density_prior = dp,
    check_prior = FALSE, beta_start = bstar, start_sweeps = start_sweeps)
  pm <- smcPosterior(fit, function(l) group_partition(l, d$truth), view = 1)
  est <- setNames(rep(0, 4), top4); k <- intersect(names(pm$estimate), top4); est[k] <- pm$estimate[k]
  # the pattern frequencies of the starting population are not recorded; report the final ones
  cat(sprintf("%-30s P=%5d steps=%3d finalESS=%7.1f  %s  (%.0fs)\n", label, P, nrow(fit$trace), 1 / sum(fit$particle_weights^2),
    paste(sprintf("%.3f", est), collapse = " "), as.numeric(Sys.time() - t0, units = "secs")))
  invisible(fit)
}
for (b in c(0.7, 0.8, 0.9)) {
  one(sprintf("start beta* = %.1f, 300 sweeps", b), b, 1000, 300)
}
one("start beta* = 0.8, 2000 sweeps", 0.8, 1000, 2000)
one("start beta* = 1.0-eps (plain chains)", 0.999, 1000, 300)
