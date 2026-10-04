#' @title Weighted ensembles of tempered chains (annealed importance sampling and
#' sequential Monte Carlo)
#' @description Runs many chains from the prior to the posterior along the
#' likelihood-tempered path \eqn{\pi_\beta \propto L^\beta P} (\eqn{\beta} from 0
#' to 1; see \code{\link{ptLadder}}) and gives every chain an importance weight.
#' Giving every chain a weight is the principled form of pooling independent
#' chains; see the guarantees and the limits below.
#'
#' \strong{Why the weights.} Pooling independent chains with equal weights
#' converges to \eqn{\sum_m b_m \pi_m}, where \eqn{\pi_m} is the posterior
#' restricted to mode \eqn{m} and \eqn{b_m} is the probability that a chain
#' lands in it. That is the posterior only if \eqn{b_m} equals the posterior
#' mass of \eqn{m}, which is not true in general, and more chains do not change
#' it. Here each particle starts at an exact draw from the prior, is moved by
#' the tempered Gibbs sweep at each \eqn{\beta} (a kernel that leaves
#' \eqn{\pi_\beta} invariant), and is weighted by the product of the incremental
#' ratios \eqn{\exp\{(\beta_k - \beta_{k-1})\ell\}}, \eqn{\ell} the data
#' log-likelihood (Neal, 2001; Del Moral, Doucet and Jasra, 2006). The weights
#' carry the information about the mass of each mode that the chains
#' cannot find from their own moves.
#'
#' \strong{What is guaranteed, and under what assumptions} (proofs and checks in
#' \code{verification/tempering/}; each item says which kind of statement it is).
#' \itemize{
#'   \item \emph{Exactly, for any number of particles (finite sample).} Suppose the
#'     initial draw is exactly from \eqn{\pi_0} (the prior), every move leaves its
#'     \eqn{\pi_\beta} invariant, and the schedule is fixed in advance (or chosen
#'     independently of the particles) and the resampling rule depends only on the
#'     weights. Then the average of the weights is an unbiased estimate of
#'     \eqn{Z_1 / Z_0}, the marginal likelihood (up to the normalisation of the
#'     prior), and, for any function \eqn{f} of the state, the weighted sum of
#'     \eqn{f} is unbiased for \eqn{Z \, E_\pi f}. Without resampling this is proved
#'     in Lean 4 (\code{Tempering/AIS.lean}); with resampling it is Del Moral (2004,
#'     Theorem 7.4.2) and I confirmed it in exact rational arithmetic for tiny
#'     systems (systematic and multinomial resampling, adaptive resampling).
#'     This does \emph{not} require the Gibbs sweeps to mix.
#'   \item \emph{Asymptotically.} The self-normalised estimator (the weighted
#'     frequency of a mode, the weighted similarity matrix) is a ratio of two
#'     unbiased estimates: it is biased at finite size but consistent as the number
#'     of particles, or of independent runs, grows, with a central limit theorem
#'     when the weights have finite variance (Del Moral, Doucet and Jasra, 2006;
#'     Chopin, 2004).
#'   \item \emph{Adaptive schedule.} With \code{schedule = "adaptive"} the temperatures
#'     depend on the particles. The estimates are consistent as the number of
#'     particles grows (Beskos et al., 2016) but are \strong{not} unbiased: the exact
#'     enumeration shows a non-zero bias. Use a fixed schedule (for example one taken
#'     from a pilot run, as \code{\link{smcReplicates}} does) for unbiased estimates.
#'   \item \emph{Starting above the prior} (\code{beta_start > 0}) is a
#'     \strong{conditional} procedure. If the particles at \code{beta_start} were exactly
#'     draws from \eqn{\pi_{\beta_{start}}} the statements above hold with
#'     \eqn{Z_{\beta_{start}}} in place of \eqn{Z_0}. If they were drawn from some other
#'     law \eqn{\mu_0}, the bias of the weighted sum of a bounded \eqn{f} is at most
#'     \eqn{\|f\|_\infty \prod_k R_k \sum_x |\mu_0(x) - \pi_{\beta_{start}}(x)|}, where
#'     \eqn{R_k} bounds the incremental weights (proved in Lean 4,
#'     \code{Tempering/AIS.lean}, \code{bias_bound}). The total-variation distance
#'     of the start from \eqn{\pi_{\beta_{start}}} is exactly what cannot be
#'     checked from the output when the start is made by running the sampler, so
#'     there is \strong{no guarantee} for \code{start_sweeps} sweeps unless that sampler
#'     is known to have mixed at \code{beta_start}. In the example in
#'     \code{verification/tempering/} the masses moved toward the exact values as
#'     \code{start_sweeps} grew, which is an illustration and not a guarantee.
#' }
#'
#' \strong{What is not guaranteed.} Nothing bounds the variance for a finite
#' number of particles, and the variance can be astronomically large. When the
#' tempering path has a sharp transition (see \code{\link{ptDiagnostics}}) the
#' weights degenerate there; in the example in \code{verification/tempering/}
#' (four clusters, three components), with 20,000 particles and no resampling the
#' effective sample size was 1.6, and with resampling the effective sample size was
#' high but the mode masses were wrong and differed between runs: all the
#' particles inherited the pattern chosen by the few lineages that crossed the
#' transition first. The effective sample size therefore does not detect this
#' failure. \code{\link{smcDiagnostics}} and \code{\link{compareRuns}} are
#' diagnostics with no guarantee. Memory grows with \code{n_particles} times the
#' data size and the particles run in sequence.
#' @inheritParams callMDI
#' @param n_particles Number of particles (chains).
#' @param schedule \code{"adaptive"} chooses each next inverse temperature so that
#' the conditional effective sample size of the incremental weights is
#' \code{cess} times \code{n_particles} (Zhou, Johansen and Aston, 2016, eq. 3.16);
#' \code{"fixed"} uses \code{betas}.
#' @param betas For \code{schedule = "fixed"}, a strictly increasing vector in
#' (0, 1] ending at 1. The run starts at 0 (the prior).
#' @param cess Target conditional ESS fraction for the adaptive schedule, in
#' (0, 1). Larger values take smaller steps.
#' @param resample_threshold Resample when the effective sample size of the
#' weights falls below this fraction of \code{n_particles}. \code{0} never
#' resamples (annealed importance sampling, with the largest weight variance),
#' \code{1} resamples at every step.
#' @param resample \code{"systematic"} (default) or \code{"multinomial"}.
#' @param sweeps_per_step Gibbs sweeps of every particle after each change of
#' temperature.
#' @param final_sweeps Further sweeps at \eqn{\beta = 1} after the last
#' temperature. Each particle's draws share its weight. Every extra sweep is
#' valid for the weighted sample, but a draw is highly dependent on the particle
#' it comes from, so these add little against more particles.
#' @param thin Thinning of the extra sweeps.
#' @param split_merge Split-merge style moves in every sweep of every particle; see
#' \code{\link{callMDI}}. The moves leave each tempered target invariant, so the
#' guarantees below are unchanged; they can only change the variance.
#' @param fixed Not supported; any non-\code{NULL} value is an error (the prior start needs unsupervised models).
#' @param beta_start,start_sweeps Start the particles at inverse temperature
#' \code{beta_start} after \code{start_sweeps} sweeps there, instead of at the prior. Conditional
#' procedure: read the corresponding item under \emph{What is guaranteed}.
#' @param max_steps Upper bound on the number of temperatures. If it is reached
#' the next temperature is set to 1 (flagged in \code{\link{smcDiagnostics}}).
#' @return An object of class \code{mdir_smc}: a list with
#' \code{allocations} (an array of draws x items x views x particles, labels from
#' 0), \code{particle_weights} (normalised), \code{log_evidence}, \code{schedule}
#' (the temperatures from 0), \code{trace} (a data frame of the run: temperature,
#' ESS, conditional ESS, whether it resampled, distinct prior ancestors),
#' \code{phis}, \code{mass}, \code{mdi_weights}, \code{data_log_likelihood} and
#' the settings.
#' Supported: complete data, \code{"G"}, \code{"MVN"} and \code{"C"} views, no
#' outlier component and no observed labels (the prior of the labels cannot be
#' sampled exactly given observed labels).
#' @references Neal, R. M. (2001). Annealed importance sampling. \emph{Statistics
#' and Computing}, 11, 125-139.
#'
#' Del Moral, P., Doucet, A. and Jasra, A. (2006). Sequential Monte Carlo
#' samplers. \emph{Journal of the Royal Statistical Society B}, 68, 411-436.
#'
#' Del Moral, P. (2004). \emph{Feynman-Kac Formulae}. Springer.
#'
#' Zhou, Y., Johansen, A. M. and Aston, J. A. D. (2016). Toward automatic model
#' comparison: an adaptive sequential Monte Carlo approach. \emph{Journal of
#' Computational and Graphical Statistics}, 25, 701-726.
#'
#' Beskos, A., Jasra, A., Kantas, N. and Thiery, A. (2016). On the convergence of
#' adaptive sequential Monte Carlo methods. \emph{Annals of Applied Probability},
#' 26, 1111-1146.
#'
#' Chopin, N. (2004). Central limit theorem for sequential Monte Carlo methods and
#' its application to Bayesian inference. \emph{Annals of Statistics}, 32,
#' 2385-2411.
#' @seealso \code{\link{smcDiagnostics}}, \code{\link{smcReplicates}},
#' \code{\link{weightedConsensus}}, \code{\link{smcPosterior}}
#' @examples
#' set.seed(1)
#' X <- list(matrix(rnorm(60 * 2, rep(c(0, 4), each = 30)), 60, 2))
#' rownames(X[[1]]) <- 1:60
#' fit <- smcMDI(X, "MVN", K = 3, n_particles = 20)
#' smcDiagnostics(fit)
#' @export
smcMDI <- function(X,
                   types,
                   K = NULL,
                   n_particles = 200,
                   schedule = c("adaptive", "fixed"),
                   betas = NULL,
                   cess = 0.9,
                   resample_threshold = 0.5,
                   resample = c("systematic", "multinomial"),
                   sweeps_per_step = 1L,
                   final_sweeps = 0L,
                   thin = 1L,
                   prior = mdiPrior(),
                   density_prior = densityPrior(),
                   phi_update = c("slice", "gibbs"),
                   max_steps = 2000L,
                   check_prior = TRUE,
                   beta_start = 0,
                   start_sweeps = 0L,
                   split_merge = 0L,
                   fixed = NULL) {
  if (!is.null(fixed)) {
    stop("smcMDI is unsupervised only: the particles start from exact draws of the prior, which are not available ",
      "given observed labels. Use callMDI(..., betas = ) (parallel tempering) for semi-supervised models.")
  }
  schedule <- match.arg(schedule)
  resample <- match.arg(resample)
  phi_update <- match.arg(phi_update)
  checkDataCorrectInput(X, types)
  .mdirCheckCount(n_particles, "n_particles", 2)
  .mdirCheckCount(sweeps_per_step, "sweeps_per_step", 1)
  .mdirCheckCount(thin, "thin", 1)
  .mdirCheckCount(max_steps, "max_steps", 1)
  .mdirCheckCount(split_merge, "split_merge", 0)
  if (!is.numeric(final_sweeps) || length(final_sweeps) != 1 || is.na(final_sweeps) || final_sweeps < 0) {
    stop("`final_sweeps` must be a non-negative integer.", call. = FALSE)
  }
  if (!is.numeric(cess) || length(cess) != 1 || is.na(cess) || cess <= 0 || cess >= 1) {
    stop("`cess` must lie in (0, 1).", call. = FALSE)
  }
  if (!is.numeric(resample_threshold) || length(resample_threshold) != 1 || is.na(resample_threshold) ||
    resample_threshold < 0) {
    stop("`resample_threshold` must be a non-negative number.", call. = FALSE)
  }
  if (!is.numeric(beta_start) || length(beta_start) != 1 || is.na(beta_start) || beta_start < 0 || beta_start >= 1) {
    stop("`beta_start` must lie in [0, 1).", call. = FALSE)
  }
  if (beta_start > 0 && start_sweeps >= 1 && !isTRUE(getOption("mdir.quiet"))) {
    message("beta_start > 0: the weights are exact only if the particles follow the tempered distribution at ",
      "beta_start, which cannot be checked from the output (see ?smcMDI).")
  }
  if (beta_start > 0 && start_sweeps < 1) {
    warning("`beta_start` > 0 without `start_sweeps`: the particles are prior draws, not draws from the tempered ",
      "distribution at `beta_start`, so the weights do not correct anything below `beta_start`.", call. = FALSE)
  }
  if (schedule == "fixed") {
    if (is.null(betas)) {
      stop("`schedule = \"fixed\"` needs `betas`.", call. = FALSE)
    }
    betas <- .mdirCheckLadder(betas)
    if (betas[1] <= beta_start) {
      stop("A fixed schedule must start above `beta_start` (0 unless set).", call. = FALSE)
    }
  } else {
    betas <- numeric(0)
  }

  N <- nrow(X[[1]])
  V <- length(X)
  P <- lapply(X, ncol)
  if (is.null(K)) {
    K <- rep(floor(N / 2), V)
  }
  if (check_prior) {
    for (m in .checkSparsity(X, types, K, prior)) message(m)
  }
  density_types <- translateTypes(types)
  outlier_types <- setupOutlierComponents(types)
  fixed <- matrix(0, N, V)

  t_0 <- Sys.time()
  raw <- runMDISMC(
    n_particles = as.integer(n_particles),
    Y = X,
    K = K,
    mixture_types = density_types,
    outlier_types = outlier_types,
    fixed = fixed,
    prior = as.numeric(prior),
    density_prior = as.numeric(density_prior),
    phi_slice = (phi_update == "slice"),
    betas = betas,
    adaptive = (schedule == "adaptive"),
    cess_target = cess,
    resample_threshold = resample_threshold,
    resample_scheme = as.integer(resample == "multinomial"),
    sweeps_per_step = as.integer(sweeps_per_step),
    max_steps = as.integer(max_steps),
    final_sweeps = as.integer(final_sweeps),
    final_thin = as.integer(thin),
    beta_start = beta_start,
    start_sweeps = as.integer(start_sweeps),
    split_merge = as.integer(split_merge)
  )
  time_taken <- Sys.time() - t_0

  P_n <- as.integer(n_particles)
  n_draws <- nrow(raw$data_log_likelihood)
  K_max <- max(K)
  LC2 <- if (V > 1) choose(V, 2) else 1L
  stack <- function(lst, dims) array(unlist(lapply(lst, as.numeric)), dim = c(dims, P_n))

  trace <- data.frame(
    step = seq_along(raw$beta),
    beta = as.numeric(raw$beta),
    ess = as.numeric(raw$ess),
    cess = as.numeric(raw$cess),
    resampled = as.logical(raw$resampled),
    unique_ancestors = as.integer(raw$unique_roots),
    log_evidence_increment = as.numeric(raw$evidence_increment)
  )
  out <- list(
    allocations = stack(raw$allocations, c(n_draws, N, V)),
    particle_weights = as.numeric(raw$particle_weights),
    log_evidence = as.numeric(raw$log_evidence),
    schedule = c(beta_start, trace$beta),
    trace = trace,
    phis = stack(raw$phis, c(n_draws, LC2)),
    mass = stack(raw$mass, c(n_draws, V)),
    mdi_weights = stack(raw$weights_record, c(n_draws, K_max, V)),
    data_log_likelihood = unname(raw$data_log_likelihood),
    root = as.integer(raw$root),
    forced_final_step = isTRUE(raw$forced_final_step),
    N = N, V = V, P = P, K = K, types = types, n_particles = P_n,
    n_draws = n_draws, schedule_type = schedule, cess = cess,
    beta_start = beta_start, start_sweeps = as.integer(start_sweeps),
    resample_threshold = resample_threshold, resample = resample,
    sweeps_per_step = as.integer(sweeps_per_step),
    final_sweeps = as.integer(final_sweeps), thin = as.integer(thin),
    sample_ids = row.names(X[[1]]), Time = time_taken
  )
  class(out) <- c("mdir_smc", "list")
  if (out$forced_final_step) {
    warning("`max_steps` was reached and the last temperature was set to 1; ",
      "the weights are probably degenerate.", call. = FALSE)
  }
  out
}

.mdirCheckCount <- function(x, name, minimum) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || x < minimum || x != round(x)) {
    stop(sprintf("`%s` must be an integer of at least %d.", name, minimum), call. = FALSE)
  }
  invisible(x)
}

# Internal: log(sum(exp(x)))
.mdirLogSumExp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}

# Internal: systematic resampling of n indices from normalised weights w
.mdirSystematic <- function(w, n = length(w)) {
  cum <- cumsum(w / sum(w))
  cum[length(cum)] <- 1
  pos <- (stats::runif(1) + seq_len(n) - 1) / n
  pmin(findInterval(pos, c(0, cum), left.open = TRUE, rightmost.closed = TRUE), length(w))
}

# Internal: the draws of one view as a matrix with one row per (particle, draw),
# the draw varying fastest, with the matching normalised weights
.smcFlatten <- function(fit, view = 1) {
  D <- dim(fit$allocations)[1]
  P <- dim(fit$allocations)[4]
  N <- dim(fit$allocations)[2]
  sub <- fit$allocations[, , view, , drop = FALSE]
  dim(sub) <- c(D, N, P)
  labels <- matrix(aperm(sub, c(1, 3, 2)), nrow = D * P, ncol = N)
  w <- rep(fit$particle_weights, each = D) / D
  list(labels = labels, weights = w / sum(w), D = D, P = P)
}

.smcCheck <- function(x) {
  if (!inherits(x, "mdir_smc")) {
    stop("`fit` must be the result of smcMDI() (class mdir_smc).", call. = FALSE)
  }
  invisible(x)
}

#' @title Print and summarise a weighted ensemble
#' @description \code{print} gives a short report; \code{summary} is
#' \code{\link{smcDiagnostics}}.
#' @param x,object An \code{mdir_smc} object from \code{\link{smcMDI}}.
#' @param ... Unused.
#' @return \code{x} invisibly (\code{print}); the diagnostics (\code{summary}).
#' @export
print.mdir_smc <- function(x, ...) {
  cat(sprintf(
    "<mdir smc> %d particles, %d items, %d view%s (%s)\n",
    x$n_particles, x$N, x$V, if (x$V == 1) "" else "s", paste(x$types, collapse = ", ")
  ))
  cat(sprintf(
    "  %d temperatures (%s schedule), %d resampling steps, log evidence %.2f\n",
    nrow(x$trace), x$schedule_type, sum(x$trace$resampled), x$log_evidence
  ))
  cat(sprintf(
    "  final effective sample size %.1f of %d; %d distinct prior ancestors\n",
    1 / sum(x$particle_weights^2), x$n_particles, length(unique(x$root))
  ))
  cat(sprintf("  run time %s. See smcDiagnostics().\n", format(x$Time, digits = 3)))
  invisible(x)
}

#' @rdname print.mdir_smc
#' @export
summary.mdir_smc <- function(object, ...) {
  smcDiagnostics(object)
}

#' @title Normalised particle weights
#' @description The weights of the particles of a weighted ensemble.
#' @param fit An \code{mdir_smc} object.
#' @param per \code{"particle"} (one weight per particle) or \code{"draw"} (one
#' per recorded draw, a particle's weight divided over its draws).
#' @return A numeric vector summing to one.
#' @export
smcWeights <- function(fit, per = c("particle", "draw")) {
  .smcCheck(fit)
  per <- match.arg(per)
  w <- fit$particle_weights / sum(fit$particle_weights)
  if (per == "draw") {
    w <- rep(w, each = fit$n_draws) / fit$n_draws
  }
  w
}

#' @title Diagnostics of a weighted ensemble
#' @description Summarises how reliable the weights are: the effective sample
#' size of the final weights, the lowest effective sample size and conditional ESS
#' along the schedule, the number of resampling steps, the number of temperatures
#' (a sharp transition in the path shows as very many small steps in an adaptive
#' schedule), the number of distinct prior draws the final particles descend from
#' (reported for information: after resampling this is usually small, and it is
#' harmless if the moves mix at the hot end, so it is not an alarm), and advisories.
#'
#' \strong{Reading it.} A small final ESS means the estimate rests on a handful of
#' particles and can be badly wrong with no sign in the estimate itself,
#' particularly if modes are separated by a sharp transition in the tempering
#' path. A large ESS does not prove the weights are adequate
#' either (a mode that no particle reached cannot be weighed), so run independent
#' replicates and compare them (\code{\link{smcReplicates}},
#' \code{\link{compareRuns}}).
#' @param fit An \code{mdir_smc} object.
#' @return An object of class \code{mdir_smc_diagnostics}.
#' @export
smcDiagnostics <- function(fit) {
  .smcCheck(fit)
  W <- fit$particle_weights / sum(fit$particle_weights)
  P <- fit$n_particles
  tr <- fit$trace
  out <- list(
    n_particles = P,
    n_temperatures = nrow(tr),
    schedule = fit$schedule,
    log_evidence = fit$log_evidence,
    final_ess = 1 / sum(W^2),
    min_ess = min(tr$ess),
    min_cess_fraction = min(tr$cess),
    n_resampled = sum(tr$resampled),
    distinct_ancestors = length(unique(fit$root)),
    forced_final_step = fit$forced_final_step
  )
  adv <- character(0)
  if (out$final_ess < 0.1 * P) {
    adv <- c(adv, sprintf("Final effective sample size %.1f is under 10%% of the particles.", out$final_ess))
  }
  if (out$forced_final_step) {
    adv <- c(adv, "The step limit was reached and the last temperature was set to 1.")
  }
  if (fit$schedule_type == "adaptive") {
    adv <- c(adv, "Adaptive schedule: estimates are consistent as the number of particles grows but not exactly unbiased; use a fixed schedule for unbiased estimates.")
  }
  out$advisories <- adv
  class(out) <- "mdir_smc_diagnostics"
  out
}

#' @rdname smcDiagnostics
#' @param x An object returned by the function this page documents.
#' @param ... Unused.
#' @export
print.mdir_smc_diagnostics <- function(x, ...) {
  cat(sprintf("Weighted ensemble: %d particles, %d temperatures\n", x$n_particles, x$n_temperatures))
  cat(sprintf("  log evidence estimate: %.3f\n", x$log_evidence))
  cat(sprintf("  final ESS %.1f; lowest ESS along the path %.1f; lowest conditional ESS fraction %.2f\n",
    x$final_ess, x$min_ess, x$min_cess_fraction))
  cat(sprintf("  resampling steps: %d; distinct prior ancestors of the final particles: %d\n",
    x$n_resampled, x$distinct_ancestors))
  if (length(x$advisories) > 0) {
    cat("Advisories:\n")
    for (a in x$advisories) cat("  -", a, "\n")
  }
  invisible(x)
}

#' @title Tail diagnostic for the importance weights of an ensemble
#' @description The effective sample size of weights can look healthy while the
#' estimate is unreliable, and can look poor for estimates that are fine, so it is
#' not a guarantee in either direction. The tail index of the weights is a more
#' informative diagnostic: the Pareto \eqn{\hat k} of Vehtari et al. (2024) fits a
#' generalised Pareto distribution to the largest weights. \eqn{\hat k < 0.5}
#' indicates a finite variance of the weights (so a central limit theorem for the
#' weighted estimates); the weighted estimates are reliable for
#' \eqn{\hat k < \min(1 - 1/\log_{10} S, 0.7)}, with \eqn{S} the number of weights;
#' above that, the estimate converges very slowly, or not at all in practice.
#'
#' It applies to \emph{independent} importance weights, that is to an ensemble run
#' without resampling (\code{resample_threshold = 0}: annealed importance
#' sampling). After resampling the particles are dependent and equally weighted
#' and the diagnostic is meaningless, so it refuses. It also has no guarantee: it
#' estimates the tail from the weights that were drawn and cannot see a region
#' that no particle reached.
#' @param fit An \code{mdir_smc} object run with \code{resample_threshold = 0}.
#' @return A list: \code{pareto_k}, the \code{threshold} above which the estimates
#' are reported as unreliable, \code{finite_variance_plausible} (\eqn{\hat k <
#' 0.5}), \code{reliable} (below the threshold) and \code{ess}.
#' @references Vehtari, A., Simpson, D., Gelman, A., Yao, Y. and Gabry, J. (2024).
#' Pareto smoothed importance sampling. \emph{Journal of Machine Learning
#' Research}, 25(72), 1-58.
#' @export
smcWeightDiagnostic <- function(fit) {
  .smcCheck(fit)
  if (any(fit$trace$resampled)) {
    stop("The ensemble was resampled, so its weights are not independent importance weights; ",
      "rerun with `resample_threshold = 0` (annealed importance sampling) to use this diagnostic.", call. = FALSE)
  }
  if (!requireNamespace("loo", quietly = TRUE)) {
    stop("This diagnostic needs the 'loo' package.", call. = FALSE)
  }
  W <- fit$particle_weights / sum(fit$particle_weights)
  S <- length(W)
  psis <- suppressWarnings(loo::psis(log(W), r_eff = 1))
  k <- as.numeric(loo::pareto_k_values(psis))
  thr <- min(1 - 1 / log10(S), 0.7)
  list(pareto_k = k, threshold = thr, finite_variance_plausible = k < 0.5, reliable = k < thr,
    ess = 1 / sum(W^2), n_weights = S)
}

#' @title Resample a weighted ensemble to equal weights
#' @description Draws \code{n_draws} (particle, draw) pairs with probability
#' proportional to their weights, so that an unweighted summary of the result
#' estimates the weighted one. This adds Monte Carlo variance; prefer the weighted
#' functions (\code{\link{weightedPSM}}, \code{\link{smcPosterior}}) for estimates
#' and use this to feed methods that need unweighted draws.
#' @param fit An \code{mdir_smc} object.
#' @param n_draws Number of draws.
#' @param method \code{"systematic"} or \code{"multinomial"}.
#' @return A data frame with columns \code{particle} and \code{draw}.
#' @export
resampleSMC <- function(fit, n_draws = 1000, method = c("systematic", "multinomial")) {
  .smcCheck(fit)
  method <- match.arg(method)
  .mdirCheckCount(n_draws, "n_draws", 1)
  w <- smcWeights(fit, "draw")
  idx <- if (method == "systematic") {
    .mdirSystematic(w, n_draws)
  } else {
    sample.int(length(w), n_draws, replace = TRUE, prob = w)
  }
  data.frame(
    particle = (idx - 1) %/% fit$n_draws + 1,
    draw = (idx - 1) %% fit$n_draws + 1
  )
}

#' @title Weighted posterior similarity matrix
#' @description The posterior similarity matrix (the probability that two items
#' share a component) of one view, with every draw carrying its weight. Unlike
#' the matrix from pooled independent chains, this estimates the posterior
#' co-clustering probabilities (under the assumptions in \code{\link{smcMDI}}).
#' @param x An \code{mdir_smc} object, or a matrix of labels with one row per
#' draw and one column per item.
#' @param view The view (for an \code{mdir_smc} object).
#' @param weights Weights of the rows when \code{x} is a matrix (default equal).
#' @return A symmetric items x items matrix.
#' @export
weightedPSM <- function(x, view = 1, weights = NULL) {
  if (inherits(x, "mdir_smc")) {
    fl <- .smcFlatten(x, view)
    labels <- fl$labels
    w <- fl$weights
    ids <- x$sample_ids
  } else {
    labels <- as.matrix(x)
    w <- if (is.null(weights)) rep(1 / nrow(labels), nrow(labels)) else weights / sum(weights)
    ids <- colnames(labels)
    if (length(w) != nrow(labels)) {
      stop("`weights` needs one value for each row of `x`.", call. = FALSE)
    }
  }
  N <- ncol(labels)
  psm <- matrix(0, N, N)
  for (k in sort(unique(as.vector(labels)))) {
    ind <- (labels == k) * sqrt(w)
    psm <- psm + crossprod(ind)
  }
  if (!is.null(ids)) {
    dimnames(psm) <- list(ids, ids)
  }
  psm
}

#' @title Point estimate clustering from a weighted ensemble
#' @description Resamples the ensemble to equal weights (see
#' \code{\link{resampleSMC}}) and finds a point estimate of the partition of one
#' view with \code{salso::salso} (the same method \code{\link{processMCMCChain}}
#' uses), together with the weighted similarity matrix.
#' @param fit An \code{mdir_smc} object.
#' @param view The view.
#' @param n_draws Number of resampled draws used for the point estimate.
#' @param ... Passed to \code{salso::salso} (for example \code{loss}).
#' @return A list with \code{clustering} (a named vector of labels),
#' \code{psm} (the weighted similarity matrix), \code{n_clusters} and
#' \code{ess} (the effective sample size of the particle weights).
#' @export
weightedConsensus <- function(fit, view = 1, n_draws = 1000, ...) {
  .smcCheck(fit)
  fl <- .smcFlatten(fit, view)
  idx <- .mdirSystematic(fl$weights, n_draws)
  est <- suppressWarnings(salso::salso(fl$labels[idx, , drop = FALSE], ...))
  names(est) <- fit$sample_ids
  list(
    clustering = est,
    psm = weightedPSM(fit, view),
    n_clusters = length(unique(est)),
    ess = 1 / sum(fit$particle_weights^2)
  )
}

#' @title Weighted posterior summary of any function of the partition
#' @description Applies \code{fun} to the labels of every recorded draw and
#' returns the weighted posterior distribution of the result: weighted
#' frequencies if it is discrete (a character, factor or integer-valued
#' label such as "which clusters merged"), the weighted mean and standard
#' deviation if it is numeric. This is the tool for the mass of a mode, defined by
#' whatever feature distinguishes the modes.
#' @param fit An \code{mdir_smc} object.
#' @param fun A function of the labels of one draw: an items x views matrix, or a
#' vector of items if \code{view} is given. It must return a single value.
#' @param view The view passed to \code{fun}, or \code{NULL} for all views.
#' @param discrete Treat the result as categories (default: if it is not numeric).
#' @return An object of class \code{mdir_smc_functional}: \code{estimate} (named
#' weighted frequencies, or the mean), \code{sd} (numeric case), \code{ess}
#' (effective sample size of the particle weights) and the per-draw
#' \code{values} and \code{weights}.
#' @export
smcPosterior <- function(fit, fun, view = NULL, discrete = NULL) {
  .smcCheck(fit)
  D <- fit$n_draws
  P <- fit$n_particles
  w <- smcWeights(fit, "draw")
  vals <- vector("list", D * P)
  k <- 0
  for (p in seq_len(P)) {
    for (d in seq_len(D)) {
      k <- k + 1
      lab <- fit$allocations[d, , , p]
      lab <- matrix(lab, nrow = fit$N, ncol = fit$V)
      vals[[k]] <- if (is.null(view)) fun(lab) else fun(lab[, view])
    }
  }
  if (any(lengths(vals) != 1)) {
    stop("`fun` must return a single value for each draw.", call. = FALSE)
  }
  vals <- unlist(vals, use.names = FALSE)
  if (is.null(discrete)) discrete <- !is.numeric(vals) || is.factor(vals)
  ess <- 1 / sum(fit$particle_weights^2)
  if (discrete) {
    vals <- as.character(vals)
    est <- tapply(w, factor(vals, levels = unique(vals)), sum)
    est <- stats::setNames(as.numeric(est), names(est))
    est <- sort(est[!is.na(est)], decreasing = TRUE)
    out <- list(estimate = est, sd = NULL, ess = ess, values = vals, weights = w, discrete = TRUE)
  } else {
    m <- sum(w * vals)
    out <- list(
      estimate = c(mean = m), sd = sqrt(sum(w * (vals - m)^2)), ess = ess,
      values = vals, weights = w, discrete = FALSE
    )
  }
  class(out) <- "mdir_smc_functional"
  out
}

#' @rdname smcPosterior
#' @param x An object returned by the function this page documents.
#' @param digits Number of digits shown.
#' @param ... Unused.
#' @export
print.mdir_smc_functional <- function(x, digits = 3, ...) {
  cat(sprintf("Weighted posterior of the function (effective sample size of the weights %.1f)\n", x$ess))
  if (x$discrete) {
    print(round(x$estimate, digits))
  } else {
    cat(sprintf("  mean %.*f, sd %.*f\n", digits, x$estimate, digits, x$sd))
  }
  invisible(x)
}

#' @title Combine independent weighted ensembles
#' @description Pools the particles of independent runs of \code{\link{smcMDI}}.
#' With \code{weight_by = "evidence"} each run is weighted by its estimate of the
#' marginal likelihood, which is what makes the pooled estimate a ratio of
#' unbiased numerator and denominator (for a fixed schedule), so it estimates the
#' posterior whatever each run's own mixture of modes (the guarantees are those
#' of \code{\link{smcMDI}}). \code{"equal"} gives runs equal weight: that is the
#' usual pooling of chains, whose limit is a mixture of the modes weighted by the
#' runs' own (basin) weights, which is not the posterior in general; it is
#' offered for comparison only.
#' @param fits A list of \code{mdir_smc} objects with the same data and settings.
#' @param weight_by \code{"evidence"} (default) or \code{"equal"}.
#' @return An \code{mdir_smc} object with all the particles, \code{log_evidence}
#' the log of the average of the runs' evidence estimates, and \code{run_weights}.
#' @export
combineSMC <- function(fits, weight_by = c("evidence", "equal")) {
  weight_by <- match.arg(weight_by)
  if (length(fits) < 1 || !all(vapply(fits, inherits, logical(1), "mdir_smc"))) {
    stop("`fits` must be a list of mdir_smc objects.", call. = FALSE)
  }
  d1 <- dim(fits[[1]]$allocations)[1:3]
  same <- vapply(fits, function(f) identical(dim(f$allocations)[1:3], d1), logical(1))
  if (!all(same)) {
    stop("All runs must have the same numbers of draws, items and views.", call. = FALSE)
  }
  J <- length(fits)
  logZ <- vapply(fits, `[[`, numeric(1), "log_evidence")
  a <- if (weight_by == "evidence") exp(logZ - .mdirLogSumExp(logZ)) else rep(1 / J, J)
  a <- a / sum(a)
  cat_dim4 <- function(nm, dims) {
    array(unlist(lapply(fits, function(f) as.numeric(f[[nm]]))), dim = c(dims, sum(vapply(fits, `[[`, 1L, "n_particles"))))
  }
  f1 <- fits[[1]]
  out <- f1
  out$allocations <- cat_dim4("allocations", d1)
  out$phis <- cat_dim4("phis", dim(f1$phis)[1:2])
  out$mass <- cat_dim4("mass", dim(f1$mass)[1:2])
  out$mdi_weights <- cat_dim4("mdi_weights", dim(f1$mdi_weights)[1:3])
  out$data_log_likelihood <- do.call(cbind, lapply(fits, `[[`, "data_log_likelihood"))
  out$particle_weights <- unlist(Map(function(f, aj) f$particle_weights / sum(f$particle_weights) * aj, fits, a))
  out$n_particles <- length(out$particle_weights)
  out$root <- unlist(Map(function(f, j) f$root + (j - 1) * 1e6L, fits, seq_len(J)))
  out$log_evidence <- .mdirLogSumExp(logZ) - log(J)
  out$run_log_evidence <- logZ
  out$run_weights <- a
  out$runs <- J
  out$weight_by <- weight_by
  out$trace <- f1$trace
  out$Time <- sum(vapply(fits, function(f) as.numeric(f$Time, units = "secs"), numeric(1)))
  class(out) <- c("mdir_smc", "list")
  out
}

#' @title Independent weighted ensembles with a common schedule
#' @description Runs \code{\link{smcMDI}} \code{n_runs} times independently on
#' the same temperature schedule, combines them by evidence
#' (\code{\link{combineSMC}}), and keeps the runs so that
#' the agreement between them and the standard errors of any summary can be
#' assessed (\code{\link{compareRuns}}, \code{\link{smcSE}}).
#'
#' The common schedule is fixed before the runs: either given (\code{betas}) or
#' taken from an adaptive pilot run. With a fixed schedule independent of the runs,
#' each run's evidence estimate is unbiased and the pooled estimate is a ratio of
#' unbiased numerator and denominator, consistent as the number of runs or of
#' particles grows. Independent runs also give an honest measure of Monte Carlo
#' error, which a single weighted ensemble cannot (the weights are dependent).
#' @inheritParams smcMDI
#' @param n_runs Number of independent runs (at least 2).
#' @param betas A fixed schedule, or \code{NULL} to take it from a pilot run.
#' @param n_cores Cores for the runs (see \code{\link{runMCMCChains}}).
#' @param verbose Report the pilot and run progress as messages.
#' @param check_prior Report prior warnings once (for the pilot run).
#' @param ... Further arguments passed to \code{\link{smcMDI}} (for example
#' \code{resample_threshold}).
#' @return An object of class \code{mdir_smc_runs}: \code{runs} (a list of
#' \code{mdir_smc}), \code{combined}, \code{schedule}, \code{pilot} (or \code{NULL}).
#' @export
smcReplicates <- function(X, types, n_runs = 8, n_particles = 200, betas = NULL,
                          K = NULL, n_cores = NULL, verbose = FALSE, check_prior = TRUE, ...) {
  .mdirCheckCount(n_runs, "n_runs", 2)
  pilot <- NULL
  if (is.null(betas)) {
    if (verbose) message("Pilot run to choose the schedule...")
    pilot <- smcMDI(X, types, K = K, n_particles = n_particles, schedule = "adaptive",
      check_prior = check_prior, ...)
    betas <- pilot$trace$beta
  } else {
    betas <- .mdirCheckLadder(betas)
  }
  n_cores <- .mdirResolveCores(n_cores, n_runs)
  fit_one <- function() {
    smcMDI(X, types, K = K, n_particles = n_particles, schedule = "fixed", betas = betas,
      check_prior = FALSE, ...)
  }
  if (verbose) message(sprintf("%d runs of %d particles on %d temperatures...", n_runs, n_particles, length(betas)))
  runs <- if (n_cores > 1) {
    .mdirRunChainsParallel(n_runs, n_cores, fit_one)
  } else {
    lapply(seq_len(n_runs), function(i) fit_one())
  }
  out <- list(runs = runs, combined = combineSMC(runs), schedule = c(0, betas), pilot = pilot)
  class(out) <- c("mdir_smc_runs", "list")
  out
}

#' @rdname smcReplicates
#' @param x An object returned by the function this page documents.
#' @param ... Unused.
#' @export
print.mdir_smc_runs <- function(x, ...) {
  ev <- vapply(x$runs, `[[`, numeric(1), "log_evidence")
  ess <- vapply(x$runs, function(r) 1 / sum(r$particle_weights^2), numeric(1))
  cat(sprintf("<mdir smc runs> %d independent runs of %d particles on a common schedule of %d temperatures\n",
    length(x$runs), x$runs[[1]]$n_particles, length(x$schedule) - 1))
  cat(sprintf("  log evidence per run: mean %.2f, sd %.2f (pooled log mean evidence %.2f)\n",
    mean(ev), stats::sd(ev), x$combined$log_evidence))
  cat(sprintf("  final ESS per run: median %.1f (range %.1f to %.1f)\n", stats::median(ess), min(ess), max(ess)))
  cat("  See smcSE() for standard errors and compareRuns() for agreement between runs.\n")
  invisible(x)
}

#' @title Standard errors of a summary from independent weighted ensembles
#' @description Applies a function of the partition to every run
#' (\code{\link{smcPosterior}}), pools the runs by their evidence estimates and
#' gives a jackknife standard error over runs, together with the range across
#' runs. The jackknife is appropriate for the ratio estimator
#' \eqn{\sum_j \hat Z_j \hat f_j / \sum_j \hat Z_j}; it assumes the runs are
#' independent and have finite variance, and it is justified asymptotically as the
#' number of runs grows (the estimator is a smooth function of means of
#' independent terms). There is no finite-sample guarantee for a small number of
#' runs; at least 4 are required and 8 or more advisable.
#' @param runs An \code{mdir_smc_runs} object from \code{\link{smcReplicates}}, or a
#' list of \code{mdir_smc} objects.
#' @param fun,view,discrete As in \code{\link{smcPosterior}}.
#' @return A data frame with one row per category (or one row, \code{mean}, for
#' a numeric function): \code{estimate} (pooled), \code{se} (jackknife),
#' \code{run_min}, \code{run_max} (over the runs' own estimates) and
#' \code{equal_pool} (the estimate when the runs are pooled with equal weight,
#' for comparison).
#' @export
smcSE <- function(runs, fun, view = NULL, discrete = NULL) {
  fits <- if (inherits(runs, "mdir_smc_runs")) runs$runs else runs
  J <- length(fits)
  if (J < 4) {
    stop("At least four independent runs are needed for a standard error.", call. = FALSE)
  }
  res <- lapply(fits, smcPosterior, fun = fun, view = view, discrete = discrete)
  cats <- unique(unlist(lapply(res, function(r) names(r$estimate))))
  F <- vapply(res, function(r) {
    v <- stats::setNames(numeric(length(cats)), cats)
    v[names(r$estimate)] <- r$estimate
    v
  }, numeric(length(cats)))
  F <- matrix(F, nrow = length(cats), dimnames = list(cats, NULL))
  logZ <- vapply(fits, `[[`, numeric(1), "log_evidence")
  pool <- function(idx) {
    a <- exp(logZ[idx] - .mdirLogSumExp(logZ[idx]))
    drop(F[, idx, drop = FALSE] %*% (a / sum(a)))
  }
  est <- pool(seq_len(J))
  jk <- vapply(seq_len(J), function(j) pool(seq_len(J)[-j]), numeric(length(cats)))
  jk <- matrix(jk, nrow = length(cats))
  se <- sqrt((J - 1) / J * rowSums((jk - rowMeans(jk))^2))
  out <- data.frame(
    category = cats, estimate = est, se = se,
    run_min = apply(F, 1, min), run_max = apply(F, 1, max),
    equal_pool = rowMeans(F), row.names = NULL
  )
  out[order(-out$estimate), ]
}

#' @title Agreement between independent runs
#' @description Compares the posterior similarity matrices of independent runs
#' of the same analysis: parallel tempering chains (\code{\link{callMDI}} with
#' \code{betas}), plain chains or weighted ensembles. Disagreement between
#' runs on co-clustering probabilities is the practical sign that they have not
#' all found the same posterior, in particular when modes have been weighted
#' differently. Agreement is necessary and not sufficient: runs can agree and
#' all miss the same mode.
#' @param runs A list of fits (\code{mdir_fit}, \code{mdir_smc}), or the result
#' of \code{\link{runMCMCChains}} or \code{\link{smcReplicates}}.
#' @param view The view.
#' @param burn Burn in (iterations) removed from \code{mdir_fit} runs; the first
#' saved draw is always dropped.
#' @param tolerance Co-clustering probabilities that differ by more than this
#' between runs count as disagreeing.
#' @return An object of class \code{mdir_run_agreement}: \code{max_range} (the
#' largest range across runs of a co-clustering probability), \code{mean_sd},
#' \code{prop_disagree} (the proportion of item pairs whose range exceeds
#' \code{tolerance}), \code{pairwise} (the mean absolute difference between each
#' pair of runs) and the matrices \code{psms}.
#' @export
compareRuns <- function(runs, view = 1, burn = 0, tolerance = 0.1) {
  fits <- if (inherits(runs, "mdir_smc_runs")) runs$runs else unclass(runs)
  if (!is.list(fits) || length(fits) < 2) {
    stop("`runs` must hold at least two fits.", call. = FALSE)
  }
  psms <- lapply(fits, function(f) {
    if (inherits(f, "mdir_smc")) {
      weightedPSM(f, view)
    } else if (!is.null(f$allocations) && length(dim(f$allocations)) == 3) {
      thin <- if (is.null(f$thin)) 1 else f$thin
      drop_n <- floor(burn / thin) + 1
      n_saved <- dim(f$allocations)[1]
      if (drop_n >= n_saved) stop("The burn in leaves no saved draws.", call. = FALSE)
      weightedPSM(f$allocations[-seq_len(drop_n), , view])
    } else {
      stop("Each run must be an mdir_fit or mdir_smc object.", call. = FALSE)
    }
  })
  N <- nrow(psms[[1]])
  if (!all(vapply(psms, function(m) nrow(m) == N, logical(1)))) {
    stop("The runs must describe the same items.", call. = FALSE)
  }
  arr <- array(unlist(psms), dim = c(N, N, length(psms)))
  ut <- upper.tri(matrix(0, N, N))
  vals <- apply(arr, 3, function(m) m[ut])
  vals <- matrix(vals, ncol = length(psms))
  rng <- apply(vals, 1, function(v) max(v) - min(v))
  J <- length(psms)
  pw <- matrix(NA_real_, J, J)
  for (a in seq_len(J)) for (b in seq_len(J)) pw[a, b] <- mean(abs(vals[, a] - vals[, b]))
  out <- list(
    max_range = max(rng),
    mean_sd = mean(apply(vals, 1, stats::sd)),
    prop_disagree = mean(rng > tolerance),
    tolerance = tolerance,
    pairwise = pw,
    n_runs = J,
    psms = psms
  )
  class(out) <- "mdir_run_agreement"
  out
}

#' @rdname compareRuns
#' @param x An object returned by the function this page documents.
#' @param ... Unused.
#' @export
print.mdir_run_agreement <- function(x, ...) {
  cat(sprintf("Agreement of %d runs on co-clustering probabilities\n", x$n_runs))
  cat(sprintf("  largest range across runs: %.3f\n", x$max_range))
  cat(sprintf("  mean between-run sd: %.4f\n", x$mean_sd))
  cat(sprintf("  item pairs whose range exceeds %.2f: %.1f%%\n", x$tolerance, 100 * x$prop_disagree))
  cat("Runs that disagree have not all found the posterior; agreement does not show that they have.\n")
  invisible(x)
}

#' @title Present a weighted ensemble as an MCMC chain
#' @description Resamples an ensemble to equal weights and arranges the draws as
#' a chain in the format of \code{\link{callMDI}}, so that the existing
#' post-processing (\code{\link{processMCMCChain}}, similarity matrices,
#' plotting) can be applied. Because the draws are an unweighted resample the
#' result carries the extra Monte Carlo variance of resampling and the draws
#' are not a Markov chain (convergence diagnostics for chains do not apply);
#' the first row is a copy of the second, since \code{processMCMCChain} always
#' drops the first saved sample.
#' @param fit An \code{mdir_smc} object.
#' @param n_draws Number of resampled draws.
#' @return An object of class \code{mdir_fit} with \code{R = n_draws},
#' \code{thin = 1}. Quantities not recorded by the ensemble (likelihood traces
#' other than the data log-likelihood, outliers) are filled with zeros or \code{NA}.
#' @export
smcAsChain <- function(fit, n_draws = 1000) {
  .smcCheck(fit)
  idx <- resampleSMC(fit, n_draws)
  idx <- idx[c(1, seq_len(n_draws)), ]            # the first row copies the second
  N <- fit$N; V <- fit$V; K_max <- max(fit$K)
  alloc <- array(0L, dim = c(nrow(idx), N, V))
  for (r in seq_len(nrow(idx))) {
    alloc[r, , ] <- fit$allocations[idx$draw[r], , , idx$particle[r]]
  }
  mass <- t(vapply(seq_len(nrow(idx)), function(r) fit$mass[idx$draw[r], , idx$particle[r]], numeric(V)))
  mass <- matrix(mass, nrow = nrow(idx), ncol = V)
  weights <- array(0, dim = c(nrow(idx), K_max, V))
  for (r in seq_len(nrow(idx))) {
    weights[r, , ] <- fit$mdi_weights[idx$draw[r], , , idx$particle[r]]
  }
  chain <- list(
    allocations = alloc,
    outliers = array(0L, dim = dim(alloc)),
    mass = mass,
    weights = weights,
    complete_likelihood = vapply(seq_len(nrow(idx)), function(r) fit$data_log_likelihood[idx$draw[r], idx$particle[r]], numeric(1)),
    evidence = rep(NA_real_, nrow(idx)),
    N_k = array(vapply(seq_len(nrow(idx) * V), function(i) {
      r <- (i - 1) %% nrow(idx) + 1; v <- (i - 1) %/% nrow(idx) + 1
      tabulate(alloc[r, , v] + 1, K_max)
    }, numeric(K_max)), dim = c(K_max, nrow(idx), V)),
    N = N, P = fit$P, K = fit$K, V = V, types = fit$types,
    R = n_draws, thin = 1, burn = 0, Semisupervised = rep(FALSE, V),
    Overfitted = rep(TRUE, V), sample_ids = fit$sample_ids,
    missing_proportion = rep(0, V), Time = fit$Time
  )
  chain$N_k <- aperm(chain$N_k, c(1, 3, 2))   # K_max x V x draws
  if (V > 1) {
    phis <- matrix(0, nrow(idx), dim(fit$phis)[2])
    for (r in seq_len(nrow(idx))) phis[r, ] <- fit$phis[idx$draw[r], , idx$particle[r]]
    chain$phis <- phis
  }
  chain$allocation_probabilities <- vector("list", V)
  class(chain) <- c("mdir_fit", "list")
  chain
}
