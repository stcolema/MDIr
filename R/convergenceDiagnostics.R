#' @title Rank-normalised split-Rhat and effective sample sizes
#' @description Computes the rank-normalised, folded, split
#' \eqn{\hat{R}} and the bulk and tail effective sample sizes (ESS) for a
#' scalar quantity sampled by several MCMC chains, following Vehtari, Gelman,
#' Simpson, Carpenter and Burkner (2021). The classical Gelman-Rubin
#' \eqn{\hat{R}} can miss non-convergence in the tails and differences in
#' scale between chains; rank-normalising, folding and splitting address that.
#' Vehtari et al. recommend treating \eqn{\hat{R} > 1.01} as a sign of
#' non-convergence and requiring a bulk and tail ESS of at least 100 per chain.
#'
#' The implementation follows the definitions in the paper and the reference
#' implementation in Stan (the same as the \pkg{posterior} package): draws are
#' split in half, transformed to normal scores using average ranks with the
#' fractional offset of 3/8, the between/within chain variance ratio gives
#' \eqn{\hat{R}}, and the ESS uses Geyer's initial monotone positive sequence
#' estimator of the autocorrelation time from the multi-chain autocorrelation
#' estimate. The reported \eqn{\hat{R}} is the maximum of the bulk (rank
#' normalised) and tail (folded) values; the tail ESS is the minimum of the
#' ESS for the 5\% and 95\% quantile indicators.
#' @param chains A numeric matrix with iterations in rows and chains in
#' columns, or a list of numeric vectors (one per chain, truncated to the
#' shortest).
#' @return A named list: \code{rhat} (the maximum of \code{rhat_bulk} and
#' \code{rhat_tail}), \code{rhat_bulk}, \code{rhat_tail}, \code{ess_bulk} and
#' \code{ess_tail}.
#' @references Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and
#' Burkner, P.-C. (2021). Rank-normalization, folding, and localization: an
#' improved Rhat for assessing convergence of MCMC. \emph{Bayesian Analysis},
#' 16(2), 667-718. \doi{10.1214/20-BA1221}
#'
#' Geyer, C. J. (1992). Practical Markov chain Monte Carlo. \emph{Statistical
#' Science}, 7(4), 473-483.
#' @export
#' @examples
#' set.seed(1)
#' # Well mixed chains
#' rankNormalizedRhat(matrix(rnorm(4000), ncol = 4))
#'
#' # Chains that have not mixed
#' rankNormalizedRhat(cbind(rnorm(1000, 0), rnorm(1000, 1), rnorm(1000, -1)))
rankNormalizedRhat <- function(chains) {
  if (is.list(chains)) {
    min_len <- min(vapply(chains, length, integer(1)))
    chains <- vapply(chains, function(x) utils::tail(x, min_len), numeric(min_len))
  }
  if (!is.matrix(chains) || !is.numeric(chains)) {
    stop("`chains` must be a numeric matrix (iterations x chains) or a list of numeric vectors.")
  }
  if (anyNA(chains) || any(!is.finite(chains))) {
    return(list(rhat = NA_real_, rhat_bulk = NA_real_, rhat_tail = NA_real_,
                ess_bulk = NA_real_, ess_tail = NA_real_))
  }
  n <- nrow(chains)
  if (n < 4) {
    stop("Need at least 4 iterations per chain to split the chains.")
  }

  split <- .splitChains(chains)

  bulk <- .rhatAndEss(.rankNormalize(split))
  folded <- .rhatAndEss(.rankNormalize(abs(split - stats::median(split))))

  # Tail ESS: minimum over the 5% and 95% quantile indicator functions
  q <- stats::quantile(split, c(0.05, 0.95), names = FALSE, type = 7)
  ess_q <- vapply(q, function(cut) {
    ind <- (split <= cut) * 1
    if (all(ind == ind[1])) {
      return(NA_real_)
    }
    .rhatAndEss(ind, normalise = FALSE)$ess
  }, numeric(1))

  list(
    rhat = max(bulk$rhat, folded$rhat),
    rhat_bulk = bulk$rhat,
    rhat_tail = folded$rhat,
    ess_bulk = bulk$ess,
    ess_tail = suppressWarnings(min(ess_q, na.rm = TRUE))
  )
}

# Split every chain into two halves, dropping the middle draw if n is odd
.splitChains <- function(x) {
  n <- nrow(x)
  half <- n %/% 2
  cbind(x[seq_len(half), , drop = FALSE], x[n - half + seq_len(half), , drop = FALSE])
}

# Rank-normalise jointly over all draws: z = qnorm((r - 3/8) / (S + 1/4))
.rankNormalize <- function(x) {
  r <- rank(x, ties.method = "average")
  z <- stats::qnorm((r - 3 / 8) / (length(r) + 1 / 4))
  matrix(z, nrow = nrow(x), ncol = ncol(x))
}

# Rhat and ESS for (already transformed, split) draws, iterations x chains
.rhatAndEss <- function(x, normalise = TRUE) {
  n <- nrow(x)
  m <- ncol(x)
  chain_means <- colMeans(x)
  chain_vars <- apply(x, 2, stats::var)

  W <- mean(chain_vars)
  B <- n * stats::var(chain_means)
  var_plus <- ((n - 1) / n) * W + B / n

  rhat <- if (W > 0) sqrt(var_plus / W) else NA_real_
  list(rhat = rhat, ess = .effectiveSampleSize(x, W, var_plus))
}

# Multi-chain ESS: rho_t = 1 - (W - mean_j acov_j(t)) / var_plus with Geyer's
# initial monotone positive sequence (Vehtari et al. 2021, Section 3.2; Stan
# reference manual, "Effective sample size").
.effectiveSampleSize <- function(x, W, var_plus) {
  n <- nrow(x)
  m <- ncol(x)
  total <- m * n
  if (!is.finite(var_plus) || var_plus <= 0) {
    return(NA_real_)
  }

  # Autocovariance with the 1 / n normalisation, via the FFT-free acf()
  acov <- vapply(seq_len(m), function(j) {
    stats::acf(x[, j], lag.max = n - 1, plot = FALSE, type = "covariance", demean = TRUE)$acf[, 1, 1]
  }, numeric(n))
  # acf() uses the 1/n divisor; Stan rescales so the lag-0 term is the unbiased variance
  mean_acov <- rowMeans(acov) * n / (n - 1)

  rho <- 1 - (W - mean_acov) / var_plus
  rho[1] <- 1

  # Geyer: sums of adjacent pairs P_t = rho_{2t} + rho_{2t+1} (t = 0, 1, ...)
  # while positive, made monotone non-increasing
  n_pairs <- (n - 1) %/% 2
  if (n_pairs < 1) {
    return(total)
  }
  pairs <- rho[2 * (seq_len(n_pairs) - 1) + 1] + rho[2 * (seq_len(n_pairs) - 1) + 2]
  positive <- which(pairs <= 0)
  last <- if (length(positive) == 0) n_pairs else positive[1] - 1
  if (last < 1) {
    # P_0 <= 0: the chain is anti-correlated at lag one; tau = 1 by the estimator
    tau <- 1
  } else {
    pairs <- cummin(pairs[seq_len(last)])
    tau <- -1 + 2 * sum(pairs)
  }
  # Stan also caps tau below by 1 / log10(S) to avoid absurdly large ESS
  tau <- max(tau, 1 / log10(total))
  total / tau
}

#' @title Assess MCMC convergence of an MDI run
#' @description Convergence diagnostics for several chains from
#' \code{\link{runMCMCChains}}. The rank-normalised split \eqn{\hat{R}} and bulk
#' and tail ESS (\code{\link{rankNormalizedRhat}}) are computed for
#' quantities that are invariant to the label switching that mixture models
#' exhibit (cluster labels are only identified up to permutation, so
#' diagnostics on raw component parameters can be misleading): the
#' complete-data and observed-data log-likelihoods, the MDI dataset
#' association parameters \eqn{\phi}, the concentration masses, and the number
#' of occupied components in each view. For the pairwise agreement of the
#' clusterings it also reports \eqn{\hat{R}} for the fusion probabilities.
#' @param mcmc_chains Output of \code{\link{runMCMCChains}} (two or more chains
#' are needed for \eqn{\hat{R}}; with one chain split-\eqn{\hat{R}} still uses
#' the two halves of the chain).
#' @param burn Number of iterations to discard as burn in. Defaults to half of
#' the iterations run.
#' @param threshold The \eqn{\hat{R}} above which a quantity is flagged
#' (1.01 following Vehtari et al., 2021).
#' @param min_ess_per_chain Minimum bulk and tail ESS per chain for a quantity to pass
#' (100, following Vehtari et al., 2021).
#' @return A data frame with a row per monitored quantity and columns
#' \code{quantity}, \code{rhat}, \code{rhat_bulk}, \code{rhat_tail},
#' \code{ess_bulk}, \code{ess_tail} and \code{converged}, with class
#' \code{mdir_convergence} and the chosen \code{threshold} and \code{burn} as
#' attributes.
#' @references Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and
#' Burkner, P.-C. (2021). Rank-normalization, folding, and localization: an
#' improved Rhat for assessing convergence of MCMC. \emph{Bayesian Analysis},
#' 16(2), 667-718.
#'
#' Stephens, M. (2000). Dealing with label switching in mixture models.
#' \emph{Journal of the Royal Statistical Society B}, 62(4), 795-809.
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' X <- lapply(1:2, function(v) {
#'   m <- matrix(rnorm(60 * 2, rep(c(0, 3), each = 30)), 60, 2)
#'   rownames(m) <- 1:60
#'   m
#' })
#' chains <- runMCMCChains(X, 3, R = 1000, thin = 5, types = c("MVN", "MVN"), K = c(4, 4))
#' assessConvergence(chains, burn = 500)
#' }
assessConvergence <- function(mcmc_chains,
                              burn = NULL,
                              threshold = 1.01,
                              min_ess_per_chain = 100) {
  n_chains <- length(mcmc_chains)
  if (n_chains < 1) {
    stop("No chains supplied.")
  }
  first <- mcmc_chains[[1]]
  thin <- first$thin
  R <- first$R
  if (is.null(burn)) {
    burn <- floor(R / 2)
  }
  # Saved sample 1 is the initial state; sample s corresponds to iteration (s - 1) * thin
  first_kept <- floor(burn / thin) + 2
  n_saved <- length(first$complete_likelihood)
  if (first_kept > n_saved - 7) {
    stop("Too few saved iterations remain after the burn in to compute diagnostics.")
  }
  keep <- seq(first_kept, n_saved)

  V <- first$V
  traces <- list()
  add_trace <- function(name, extract) {
    traces[[name]] <<- vapply(mcmc_chains, function(ch) extract(ch)[keep], numeric(length(keep)))
  }
  add_trace("complete_likelihood", function(ch) ch$complete_likelihood)
  if (!is.null(first$observed_likelihood)) {
    add_trace("observed_likelihood", function(ch) ch$observed_likelihood)
  }
  for (v in seq_len(V)) {
    add_trace(paste0("mass[", v, "]"), function(ch) ch$mass[, v])
    add_trace(paste0("occupied_components[", v, "]"), function(ch) {
      apply(ch$allocations[, , v, drop = TRUE], 1, function(z) length(unique(z)))
    })
  }
  if (V > 1) {
    pairs <- utils::combn(V, 2)
    for (i in seq_len(ncol(pairs))) {
      add_trace(paste0("phi[", pairs[1, i], ",", pairs[2, i], "]"), function(ch) ch$phis[, i])
      add_trace(paste0("agreement[", pairs[1, i], ",", pairs[2, i], "]"), function(ch) {
        rowMeans(ch$allocations[, , pairs[1, i]] == ch$allocations[, , pairs[2, i]])
      })
    }
  }

  # A single chain is split into two halves by rankNormalizedRhat
  rows <- lapply(names(traces), function(nm) {
    d <- rankNormalizedRhat(traces[[nm]])
    data.frame(quantity = nm, rhat = d$rhat, rhat_bulk = d$rhat_bulk, rhat_tail = d$rhat_tail,
               ess_bulk = d$ess_bulk, ess_tail = d$ess_tail, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  n_effective_chains <- n_chains
  out$converged <- with(out, !is.na(rhat) & rhat < threshold &
    ess_bulk >= min_ess_per_chain * n_effective_chains & ess_tail >= min_ess_per_chain * n_effective_chains)
  attr(out, "threshold") <- threshold
  attr(out, "burn") <- burn
  attr(out, "n_chains") <- n_chains
  class(out) <- c("mdir_convergence", "data.frame")
  out
}

#' @export
print.mdir_convergence <- function(x, digits = 3, ...) {
  cat("MDI convergence diagnostics (", attr(x, "n_chains"), " chain(s), burn = ", attr(x, "burn"),
      ", Rhat threshold ", attr(x, "threshold"), ")\n", sep = "")
  df <- as.data.frame(x)
  num <- vapply(df, is.numeric, logical(1))
  df[num] <- lapply(df[num], signif, digits)
  print(df, row.names = FALSE)
  if (!all(x$converged)) {
    cat("\nQuantities not meeting the Rhat / ESS criteria: ",
        paste(x$quantity[!x$converged], collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}
