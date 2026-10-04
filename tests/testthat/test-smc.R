# Weighted ensembles (annealed importance sampling / sequential Monte Carlo). The
# statistical checks compare with exact references computed independently of the
# sampler (helper-tempering.R); see also verification/tempering/.

smc_small <- function(N = 40, seed = 1) {
  set.seed(seed)
  m <- matrix(rnorm(N * 2, rep(c(0, 4), each = N / 2)), N, 2)
  rownames(m) <- seq_len(N)
  list(m)
}

test_that("input checks", {
  X <- smc_small()
  expect_error(smcMDI(X, "MVN", K = 3, n_particles = 1), "n_particles")
  expect_error(smcMDI(X, "MVN", K = 3, n_particles = 10, cess = 1), "cess")
  expect_error(smcMDI(X, "MVN", K = 3, n_particles = 10, schedule = "fixed"), "needs `betas`")
  expect_error(smcMDI(X, "MVN", K = 3, n_particles = 10, schedule = "fixed", betas = c(0.5, 0.4, 1)), "increasing")
  expect_error(smcMDI(X, "MVN", K = 3, n_particles = 10, schedule = "fixed", betas = c(0.2, 0.5)), "must be 1")
  expect_error(smcMDI(X, "MVN", K = 3, n_particles = 10, schedule = "fixed", betas = c(0.6, 1), beta_start = 0.5, start_sweeps = 5),
    NA)
  expect_error(smcMDI(X, "MVN", K = 3, n_particles = 10, schedule = "fixed", betas = c(0.4, 1), beta_start = 0.5, start_sweeps = 5),
    "above `beta_start`")
  expect_error(smcMDI(X, "MVN", K = 3, n_particles = 10, beta_start = 1), "beta_start")
  expect_warning(smcMDI(X, "MVN", K = 3, n_particles = 10, beta_start = 0.5), "not draws from the tempered")
  expect_error(smcReplicates(X, "MVN", n_runs = 1, K = 3), "n_runs")
})

test_that("models without tempered conditionals are refused", {
  X <- smc_small()
  Xm <- X
  Xm[[1]][3, 1] <- NA
  expect_error(smcMDI(Xm, "MVN", K = 3, n_particles = 10, check_prior = FALSE), "complete data")
  expect_error(smcMDI(X, "TAGM", K = 3, n_particles = 10, check_prior = FALSE), "outlier")
})

test_that("a run has the documented structure", {
  X <- smc_small()
  set.seed(2)
  fit <- smcMDI(X, "MVN", K = 3, n_particles = 30, final_sweeps = 6, thin = 3, check_prior = FALSE)
  expect_s3_class(fit, "mdir_smc")
  expect_equal(dim(fit$allocations), c(3, 40, 1, 30))
  expect_equal(sum(fit$particle_weights), 1)
  expect_true(all(fit$allocations >= 0 & fit$allocations < 3))
  expect_equal(fit$schedule[1], 0)
  expect_equal(utils::tail(fit$schedule, 1), 1)
  expect_true(all(diff(fit$schedule) > 0))
  expect_true(all(fit$trace$ess <= 30 + 1e-9 & fit$trace$ess >= 1))
  expect_true(all(fit$trace$cess >= 0.9 - 1e-6 | fit$trace$beta == 1))
  expect_equal(dim(fit$data_log_likelihood), c(3, 30))
  expect_output(print(fit), "particles")
  expect_s3_class(smcDiagnostics(fit), "mdir_smc_diagnostics")
  expect_output(print(smcDiagnostics(fit)), "log evidence")
  expect_identical(summary(fit)$n_particles, 30L)
  expect_equal(smcWeights(fit, "draw"), rep(fit$particle_weights, each = 3) / 3)
  # resampling thresholds: never resample keeps unequal weights, always keeps equal
  never <- smcMDI(X, "MVN", K = 3, n_particles = 30, resample_threshold = 0, check_prior = FALSE)
  expect_false(any(never$trace$resampled))
  expect_gt(stats::sd(never$particle_weights), 0)
  always <- smcMDI(X, "MVN", K = 3, n_particles = 30, resample_threshold = 1.01, check_prior = FALSE)
  expect_true(all(always$trace$resampled))
  expect_equal(always$particle_weights, rep(1 / 30, 30))
})

test_that("a fixed schedule is honoured and the multi-view model runs", {
  set.seed(3)
  X <- lapply(1:2, function(v) { m <- matrix(rnorm(30 * 2, rep(c(0, 4), each = 15)), 30, 2); rownames(m) <- 1:30; m })
  betas <- ptLadder(8, 0.01)
  fit <- smcMDI(X, c("G", "MVN"), K = c(3, 3), n_particles = 25, schedule = "fixed", betas = betas, check_prior = FALSE)
  expect_equal(fit$trace$beta, betas)
  expect_equal(dim(fit$allocations), c(1, 30, 2, 25))
  expect_equal(dim(fit$phis), c(1, 1, 25))
  expect_true(all(fit$phis > 0))
})

test_that("evidence and unnormalised posterior are unbiased against the exact values (fixed schedule)", {
  skip_on_cran()
  set.seed(31)
  N <- 5; K <- 2
  X <- matrix(sample(0:2, N, TRUE, c(.5, .3, .2)), N, 1)
  rownames(X) <- seq_len(N)
  alpha <- mdir:::densityHyperparameters(X, K, 2, numeric(0))$concentration
  ex <- tp_exact_L1(X, K, alpha, 1)
  # exact log evidence: log sum_c prior(c) prod marginal likelihood
  lp <- vapply(seq_len(nrow(ex$grid)), function(i) {
    c_i <- ex$grid[i, ]
    log(tp_prior_counts(tabulate(c_i + 1, K), K)) +
      sum(vapply(0:(K - 1), function(k) if (any(c_i == k)) tp_log_marg_cat(X[c_i == k, , drop = FALSE], alpha, 1) else 0, numeric(1)))
  }, numeric(1))
  log_Z <- max(lp) + log(sum(exp(lp - max(lp))))
  code <- function(a) as.numeric(a %*% K^(seq_len(N) - 1))
  betas <- ptLadder(12, 0.005)
  runs <- lapply(1:150, function(i) {
    set.seed(500 + i)
    f <- smcMDI(list(X), "C", K = K, n_particles = 40, schedule = "fixed", betas = betas, check_prior = FALSE)
    codes <- code(t(f$allocations[1, , 1, ]))
    w <- tapply(f$particle_weights, factor(codes + 1, levels = seq_len(2^N)), sum)
    w[is.na(w)] <- 0
    list(Z = exp(f$log_evidence - log_Z), fz = as.numeric(w) * exp(f$log_evidence - log_Z))
  })
  Zr <- vapply(runs, `[[`, 1, "Z")
  expect_lt(abs(mean(Zr) - 1) / (stats::sd(Zr) / sqrt(length(Zr))), 4)
  Fz <- vapply(runs, `[[`, numeric(2^N), "fz")
  top <- order(-ex$prob)[1:6]
  z <- (rowMeans(Fz)[top] - ex$prob[top]) / (apply(Fz, 1, stats::sd)[top] / sqrt(ncol(Fz)))
  expect_lt(max(abs(z)), 4)
})

test_that("weighted helpers agree with direct calculations", {
  X <- smc_small()
  set.seed(4)
  fit <- smcMDI(X, "MVN", K = 3, n_particles = 25, final_sweeps = 4, thin = 2, check_prior = FALSE,
    resample_threshold = 0)
  # weighted PSM equals the direct weighted sum over draws
  psm <- weightedPSM(fit)
  expect_equal(dim(psm), c(40, 40))
  expect_equal(unname(diag(psm)), rep(1, 40))
  expect_true(isSymmetric(psm))
  w <- smcWeights(fit, "draw")
  lab <- matrix(aperm(fit$allocations[, , 1, ], c(1, 3, 2)), nrow = 3 * 25)
  direct <- matrix(0, 40, 40)
  for (r in seq_len(nrow(lab))) direct <- direct + w[r] * outer(lab[r, ], lab[r, ], "==")
  expect_equal(unname(psm), unname(direct))
  # equal weights give the unweighted matrix of the draws
  expect_equal(unname(weightedPSM(lab[1:5, ], weights = rep(1, 5))), unname(createSimilarityMat(lab[1:5, ])))
  expect_error(weightedPSM(lab, weights = 1:2), "one value")

  # functional: weighted frequencies sum to one and agree with a manual sum
  f2 <- function(l) paste(sort(table(l)), collapse = "-")
  pm <- smcPosterior(fit, f2, view = 1)
  expect_equal(sum(pm$estimate), 1)
  manual <- tapply(w, factor(apply(lab, 1, f2)), sum)
  expect_equal(unname(pm$estimate[names(manual)]), as.numeric(manual))
  num <- smcPosterior(fit, function(l) length(unique(l)), view = 1)
  expect_equal(unname(num$estimate), sum(w * apply(lab, 1, function(r) length(unique(r)))))
  expect_output(print(num), "mean")
  expect_error(smcPosterior(fit, function(l) c(1, 2), view = 1), "single value")
  expect_error(weightedPSM(1:3 + 0), NA)

  # resampling to equal weights follows the weights
  idx <- resampleSMC(fit, n_draws = 4000)
  freq <- tabulate(idx$particle, 25) / 4000
  expect_lt(max(abs(freq - fit$particle_weights)), 0.03)

  wc <- weightedConsensus(fit, n_draws = 200)
  expect_length(wc$clustering, 40)
  expect_equal(wc$n_clusters, length(unique(wc$clustering)))
  expect_equal(wc$psm, psm)
})

test_that("a weighted ensemble feeds the existing chain workflow", {
  X <- smc_small()
  set.seed(5)
  fit <- smcMDI(X, "MVN", K = 3, n_particles = 40, final_sweeps = 5, thin = 1, check_prior = FALSE)
  ch <- smcAsChain(fit, n_draws = 60)
  expect_s3_class(ch, "mdir_fit")
  expect_equal(dim(ch$allocations), c(61, 40, 1))
  pr <- processMCMCChain(ch, burn = 0, construct_psm = TRUE)
  expect_equal(dim(pr$allocations)[1], 60)
  expect_length(pr$pred[[1]], 40)
  expect_equal(dim(pr$psms[[1]]), c(40, 40))
})

test_that("independent runs combine by their evidence estimates", {
  X <- smc_small()
  betas <- ptLadder(10, 0.01)
  runs <- lapply(1:5, function(i) { set.seed(60 + i); smcMDI(X, "MVN", K = 3, n_particles = 15, schedule = "fixed",
    betas = betas, check_prior = FALSE) })
  cb <- combineSMC(runs)
  expect_equal(sum(cb$particle_weights), 1)
  expect_equal(cb$n_particles, 75)
  logZ <- vapply(runs, `[[`, 1, "log_evidence")
  a <- exp(logZ - max(logZ)); a <- a / sum(a)
  expect_equal(cb$run_weights, a)
  expect_equal(cb$log_evidence, log(mean(exp(logZ))))
  expect_equal(cb$particle_weights[1:15], runs[[1]]$particle_weights * a[1])
  eq <- combineSMC(runs, "equal")
  expect_equal(eq$run_weights, rep(1 / 5, 5))
  expect_error(combineSMC(list(runs[[1]], smcMDI(smc_small(60), "MVN", K = 3, n_particles = 5, check_prior = FALSE))), "same numbers")
  expect_error(combineSMC(list(1, 2)), "mdir_smc")

  # jackknife standard errors need at least four runs and give a pooled estimate in [0, 1]
  expect_error(smcSE(runs[1:3], function(l) length(unique(l)), view = 1, discrete = TRUE), "four")
  se <- smcSE(runs, function(l) paste0("k", length(unique(l))), view = 1)
  expect_equal(sum(se$estimate), 1)
  expect_true(all(se$se >= 0))
  expect_equal(sum(se$equal_pool), 1)
})

test_that("smcReplicates fixes the schedule before the runs and agrees with its parts", {
  skip_on_cran()
  X <- smc_small(30)
  set.seed(8)
  rr <- smcReplicates(X, "MVN", n_runs = 4, n_particles = 12, K = 3, check_prior = FALSE)
  expect_s3_class(rr, "mdir_smc_runs")
  expect_length(rr$runs, 4)
  expect_true(all(vapply(rr$runs, function(r) identical(r$trace$beta, rr$pilot$trace$beta), logical(1))))
  expect_equal(rr$combined$n_particles, 48)
  expect_output(print(rr), "independent runs")
  ag <- compareRuns(rr)
  expect_s3_class(ag, "mdir_run_agreement")
  expect_gte(ag$max_range, 0)
  expect_output(print(ag), "Agreement")
})

test_that("compareRuns detects disagreement and accepts chains", {
  set.seed(9)
  X <- smc_small(30)
  a <- callMDI(X, R = 60, thin = 2, types = "MVN", K = 3, check_prior = FALSE)
  expect_equal(compareRuns(list(a, a))$max_range, 0)
  # two fits that cluster the items differently disagree
  b <- a
  b$allocations <- 2 - a$allocations
  b$allocations[, 1:15, ] <- (b$allocations[, 1:15, ] + 1) %% 3
  d <- compareRuns(list(a, b))
  expect_gt(d$max_range, 0.5)
  expect_gt(d$prop_disagree, 0)
  expect_error(compareRuns(list(a)), "two fits")
  expect_error(compareRuns(list(a, list(allocations = 1))), "mdir_fit or mdir_smc")
})

test_that("a start above the prior moves the schedule and keeps weights equal at the start", {
  X <- smc_small()
  set.seed(10)
  fit <- smcMDI(X, "MVN", K = 3, n_particles = 20, beta_start = 0.5, start_sweeps = 10, check_prior = FALSE)
  expect_equal(fit$schedule[1], 0.5)
  expect_true(all(fit$trace$beta > 0.5))
  expect_equal(fit$beta_start, 0.5)
})

test_that("the weight tail diagnostic separates light and heavy tails and refuses resampled runs", {
  skip_if_not_installed("loo")
  fake <- function(w, resampled = FALSE) {
    structure(list(particle_weights = w / sum(w), trace = data.frame(resampled = resampled), n_particles = length(w)),
      class = c("mdir_smc", "list"))
  }
  set.seed(12)
  light <- exp(rnorm(4000, 0, 0.5))            # all moments finite
  heavy <- (1 - runif(4000))^(-1 / 0.9)         # Pareto tail index 0.9: infinite variance (and nearly infinite mean)
  dl <- smcWeightDiagnostic(fake(light))
  dh <- smcWeightDiagnostic(fake(heavy))
  expect_lt(dl$pareto_k, 0.5)
  expect_true(dl$reliable)
  expect_gt(dh$pareto_k, 0.7)
  expect_false(dh$reliable)
  expect_equal(dl$threshold, 0.7)
  expect_error(smcWeightDiagnostic(fake(light, resampled = TRUE)), "resampled")
  # an annealed-importance-sampling run goes through
  X <- smc_small()
  set.seed(13)
  fit <- smcMDI(X, "MVN", K = 3, n_particles = 60, resample_threshold = 0, check_prior = FALSE)
  d <- smcWeightDiagnostic(fit)
  expect_true(is.finite(d$pareto_k))
})
