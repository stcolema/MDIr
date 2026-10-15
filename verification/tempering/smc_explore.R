args <- commandArgs(TRUE); lib <- args[1]
.libPaths(c(lib, .libPaths())); suppressMessages(library(mdir)); options(mdir.quiet = TRUE)
source("bimodal_setup.R")
mm <- readRDS("mode_mass_sep2.5_results.rds"); exact <- mm$exact
top4 <- names(sort(exact, decreasing = TRUE))[1:4]
d <- make_data(sep = 2.5); dp <- densityPrior(scale_pool_shape = 0)
pilot <- { set.seed(5); smcMDI(list(d$X), "MVN", K = 3, n_particles = 300, cess = 0.9, density_prior = dp, check_prior = FALSE) }
betas <- pilot$trace$beta
cat("exact:", paste(top4, round(exact[top4], 3), collapse = "  "), "\n")
run <- function(label, P, thr, scheme = "systematic", sched = "fixed", cess = 0.9, seed = 1) {
  set.seed(seed)
  t0 <- Sys.time()
  fit <- if (sched == "fixed") smcMDI(list(d$X), "MVN", K = 3, n_particles = P, schedule = "fixed", betas = betas,
      resample_threshold = thr, resample = scheme, density_prior = dp, check_prior = FALSE)
    else smcMDI(list(d$X), "MVN", K = 3, n_particles = P, cess = cess, resample_threshold = thr, density_prior = dp, check_prior = FALSE)
  pm <- smcPosterior(fit, function(l) group_partition(l, d$truth), view = 1)
  est <- setNames(rep(0, 4), top4); est[intersect(names(pm$estimate), top4)] <- pm$estimate[intersect(names(pm$estimate), top4)]
  cat(sprintf("%-34s P=%5d steps=%3d finalESS=%7.1f  %s  (%.0fs)\n", label, P, nrow(fit$trace), 1 / sum(fit$particle_weights^2),
    paste(sprintf("%.3f", est), collapse = " "), as.numeric(Sys.time() - t0, units = "secs")))
}
run("SMC adaptive cess .9, thr .5", 300, 0.5, sched = "adaptive")
run("SMC adaptive cess .9, thr .5", 3000, 0.5, sched = "adaptive")
run("SMC adaptive cess .99, thr .5", 3000, 0.5, sched = "adaptive", cess = 0.99)
run("SMC pilot schedule, thr .5", 3000, 0.5)
run("SMC pilot schedule, resample always", 3000, 1)
run("AIS pilot schedule (no resampling)", 3000, 0)
run("AIS pilot schedule (no resampling)", 20000, 0, seed = 2)
