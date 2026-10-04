#' @title Parallel tempering ladder
#' @description Inverse temperatures for the \code{betas} argument of
#' \code{\link{callMDI}}, \code{\link{runMCMCChains}} and \code{\link{fitMDI}}.
#'
#' The replica at inverse temperature \eqn{\beta} samples
#' \eqn{\pi_\beta \propto L^\beta P}, where \eqn{L} is the likelihood of the data
#' given the component assignments and component parameters and \eqn{P} is
#' the rest of the model (the priors and the coupling of the views). \eqn{\beta =
#' 1} is the posterior and \eqn{\beta = 0} the prior.
#' @param n_temperatures Number of replicas (at least 2).
#' @param beta_min The smallest inverse temperature. \code{0} is the prior, which
#' mixes fastest, but a geometric ladder needs a positive value.
#' @param spacing \code{"geometric"} (the default) or \code{"linear"}. A
#' geometric ladder is a reasonable first guess when the log-likelihood scales
#' with the number of items; \code{\link{tuneLadder}} and
#' \code{\link{adaptLadder}} refine it from the observed exchange rates.
#' @return A strictly increasing numeric vector ending at 1.
#' @examples
#' ptLadder(5)
#' ptLadder(5, beta_min = 0, spacing = "linear")
#' @export
ptLadder <- function(n_temperatures, beta_min = 0.01, spacing = c("geometric", "linear")) {
  spacing <- match.arg(spacing)
  if (!is.numeric(n_temperatures) || length(n_temperatures) != 1 || is.na(n_temperatures) ||
    n_temperatures < 2 || n_temperatures != round(n_temperatures)) {
    stop("`n_temperatures` must be an integer of at least 2.", call. = FALSE)
  }
  if (!is.numeric(beta_min) || length(beta_min) != 1 || is.na(beta_min) || beta_min < 0 || beta_min >= 1) {
    stop("`beta_min` must lie in [0, 1).", call. = FALSE)
  }
  if (spacing == "geometric") {
    if (beta_min <= 0) {
      stop("A geometric ladder needs `beta_min` > 0.", call. = FALSE)
    }
    return(exp(seq(log(beta_min), 0, length.out = n_temperatures)))
  }
  seq(beta_min, 1, length.out = n_temperatures)
}

.mdirCheckLadder <- function(betas) {
  if (!is.numeric(betas) || length(betas) < 1 || anyNA(betas) || any(!is.finite(betas))) {
    stop("`betas` must be a numeric vector of inverse temperatures.", call. = FALSE)
  }
  if (any(betas < 0 | betas > 1)) {
    stop("`betas` must lie in [0, 1].", call. = FALSE)
  }
  if (length(betas) > 1) {
    if (any(diff(betas) <= 0)) {
      stop("`betas` must be strictly increasing.", call. = FALSE)
    }
    if (betas[length(betas)] != 1) {
      stop("The last element of `betas` must be 1 (the posterior).", call. = FALSE)
    }
  }
  as.numeric(betas)
}

.mdirTidyTempering <- function(tempering, betas, swap_scheme) {
  if (length(betas) < 2 || length(tempering) == 0) {
    return(NULL)
  }
  out <- list(
    betas = as.numeric(tempering$betas),
    swap_scheme = swap_scheme,
    swap_every = as.integer(tempering$swap_every),
    swap_attempts = as.numeric(tempering$swap_attempts),
    swap_accepts = as.numeric(tempering$swap_accepts),
    swap_accept_prob_sum = as.numeric(tempering$swap_accept_prob_sum),
    rejection_rate = as.numeric(tempering$rejection_rate),
    round_trips = as.numeric(tempering$round_trips),
    data_log_likelihood = unname(tempering$data_log_likelihood),
    replica = unname(tempering$replica)
  )
  class(out) <- "mdir_tempering"
  out
}

#' @title Update a parallel tempering ladder from exchange rejection rates
#' @description Places the inverse temperatures so that the estimated rejection
#' rates of all neighbouring exchanges are equal (Syed et al., 2022, Algorithm
#' 2). The cumulative communication barrier \eqn{\Lambda(\beta)} is estimated by
#' the running sum of the rejection rates at the current ladder, interpolated
#' monotonically (Fritsch-Carlson), and the new ladder has equal increments of
#' \eqn{\Lambda}. The end points are kept.
#'
#' The rates come from the exchange attempts of a finished run
#' (\code{\link{ptDiagnostics}}); each is estimated as one minus the mean
#' acceptance probability over the attempts, not the fraction accepted.
#' Estimates from a short run that has not reached equilibrium are noisy and can
#' be biased; \code{\link{adaptLadder}} repeats the update.
#' @param betas The ladder the rejection rates were measured at.
#' @param rejection_rate Estimated rejection rate of each neighbouring pair
#' (\code{length(betas) - 1} values in [0, 1]).
#' @param n_temperatures Number of temperatures of the new ladder (defaults to
#' the current number).
#' @return The new ladder.
#' @references Syed, S., Bouchard-Cote, A., Deligiannidis, G. and Doucet, A.
#' (2022). Non-reversible parallel tempering: a scalable highly parallel MCMC
#' scheme. \emph{Journal of the Royal Statistical Society Series B}, 84(2),
#' 321-350.
#' @export
tuneLadder <- function(betas, rejection_rate, n_temperatures = length(betas)) {
  betas <- .mdirCheckLadder(betas)
  if (length(betas) < 2) {
    stop("A ladder needs at least two temperatures.", call. = FALSE)
  }
  if (length(rejection_rate) != length(betas) - 1 || anyNA(rejection_rate) ||
    any(rejection_rate < 0 | rejection_rate > 1)) {
    stop("`rejection_rate` needs one value in [0, 1] for each neighbouring pair.", call. = FALSE)
  }
  if (n_temperatures < 2) {
    stop("`n_temperatures` must be at least 2.", call. = FALSE)
  }
  # The estimated cumulative barrier at each current temperature. A small floor
  # keeps it strictly increasing so that it can be inverted.
  barrier <- c(0, cumsum(pmax(rejection_rate, 1e-8)))
  total <- barrier[length(barrier)]
  interp <- stats::splinefun(betas, barrier, method = "monoH.FC")
  target <- seq(0, total, length.out = n_temperatures)
  new_betas <- numeric(n_temperatures)
  new_betas[1] <- betas[1]
  new_betas[n_temperatures] <- 1
  if (n_temperatures > 2) {
    for (k in 2:(n_temperatures - 1)) {
      new_betas[k] <- stats::uniroot(
        function(b) interp(b) - target[k], c(betas[1], 1), tol = 1e-10
      )$root
    }
  }
  new_betas
}

#' @title Diagnostics of the exchange process in parallel tempering
#' @description Summarises how well replicas exchange: for each neighbouring
#' pair the acceptance and (Rao-Blackwellised) rejection rate of the exchange
#' proposals, the cumulative communication barrier, the round trips that
#' replicas made between the hottest and coldest temperature, and the round-trip
#' rate that Syed et al. (2022) predict from the rejection rates.
#'
#' \strong{What this does and does not show.} A round trip means a state
#' travelled from the prior end of the ladder to the posterior and back, so a
#' high rate is evidence that states are carried between temperatures. It does
#' not show that the cold chain has converged: the guarantee is only
#' asymptotic (the cold chain targets the posterior whatever the ladder, once
#' the whole ladder has mixed), and slow mixing within a temperature, a phase
#' transition in the path, or a ladder too short to bridge the prior and the
#' posterior all defeat tempering without being visible here. Few or no round
#' trips, or one pair with a rejection rate near 1, are reliable signs of
#' trouble. The predicted rates assume that every replica is in
#' equilibrium at its temperature when it is proposed an exchange, and
#' that the prior end of the ladder is sampled exactly; neither holds exactly
#' for this sampler.
#' @param fit A fit from \code{\link{callMDI}} run with a ladder in \code{betas},
#' or a list of such fits (the exchange counts are pooled).
#' @return An object of class \code{mdir_pt_diagnostics}: a list with
#' \code{pairs} (a data frame with one row per neighbouring pair),
#' \code{barrier} (the estimated global communication barrier, the sum of the
#' rejection rates), \code{round_trips}, \code{round_trip_rate_observed},
#' \code{round_trip_rate_predicted} (finite ladder) and
#' \code{round_trip_rate_limit} (many temperatures).
#' @references Syed, S., Bouchard-Cote, A., Deligiannidis, G. and Doucet, A.
#' (2022). Non-reversible parallel tempering: a scalable highly parallel MCMC
#' scheme. \emph{Journal of the Royal Statistical Society Series B}, 84(2),
#' 321-350. Rates: Corollary 1 and Theorem 3.
#' @export
ptDiagnostics <- function(fit) {
  fits <- if (inherits(fit, "mdir_fit")) list(fit) else unclass(fit)
  temp <- lapply(fits, function(f) f$tempering)
  if (any(vapply(temp, is.null, logical(1)))) {
    stop("The fit has no tempering record; run it with `betas` of length two or more.", call. = FALSE)
  }
  betas <- temp[[1]]$betas
  for (t in temp) {
    if (!isTRUE(all.equal(t$betas, betas))) {
      stop("All fits must use the same ladder.", call. = FALSE)
    }
  }
  attempts <- Reduce(`+`, lapply(temp, `[[`, "swap_attempts"))
  accepts <- Reduce(`+`, lapply(temp, `[[`, "swap_accepts"))
  prob_sum <- Reduce(`+`, lapply(temp, `[[`, "swap_accept_prob_sum"))
  trips <- sum(vapply(temp, function(t) sum(t$round_trips), numeric(1)))
  rounds <- sum(vapply(fits, function(f) floor(f$R / temp[[1]]$swap_every), numeric(1)))

  rej <- ifelse(attempts > 0, 1 - prob_sum / attempts, NA_real_)
  pairs <- data.frame(
    beta_hot = betas[-length(betas)],
    beta_cold = betas[-1],
    attempts = attempts,
    accepted = accepts,
    acceptance_rate = ifelse(attempts > 0, accepts / attempts, NA_real_),
    rejection_rate = rej
  )
  barrier <- sum(rej)
  finite_ineff <- sum(rej / (1 - rej))
  out <- list(
    betas = betas,
    pairs = pairs,
    barrier = barrier,
    round_trips = trips,
    swap_rounds = rounds,
    round_trip_rate_observed = trips / rounds,
    round_trip_rate_predicted = 1 / (2 + 2 * finite_ineff),
    round_trip_rate_limit = 1 / (2 + 2 * barrier)
  )
  class(out) <- "mdir_pt_diagnostics"
  out
}

#' @rdname ptDiagnostics
#' @param x An object returned by \code{ptDiagnostics}.
#' @param digits Number of significant digits shown.
#' @param ... Unused.
#' @export
print.mdir_pt_diagnostics <- function(x, digits = 3, ...) {
  cat(sprintf(
    "Parallel tempering: %d temperatures, beta from %s to 1\n",
    length(x$betas), format(signif(min(x$betas), 3))
  ))
  p <- x$pairs
  p$beta_hot <- signif(p$beta_hot, digits)
  p$beta_cold <- signif(p$beta_cold, digits)
  p$acceptance_rate <- round(p$acceptance_rate, digits)
  p$rejection_rate <- round(p$rejection_rate, digits)
  print(p, row.names = FALSE)
  cat(sprintf(
    "\nEstimated communication barrier (sum of rejection rates): %.2f\n", x$barrier
  ))
  cat(sprintf(
    "Round trips (all replicas): %d in %d exchange rounds\n",
    as.integer(x$round_trips), as.integer(x$swap_rounds)
  ))
  cat(sprintf(
    "Round-trip rate per round: observed %.4f; predicted %.4f (this ladder), %.4f (limit)\n",
    x$round_trip_rate_observed, x$round_trip_rate_predicted, x$round_trip_rate_limit
  ))
  cat(
    "Exchange diagnostics only: they do not show that the cold chain has converged.\n"
  )
  invisible(x)
}

#' @title Adapt a parallel tempering ladder with pilot runs
#' @description Runs short pilot fits, estimates the rejection rate of each
#' neighbouring exchange and moves the ladder with \code{\link{tuneLadder}}, so
#' that all rates become equal. The pilot runs start from the supplied initial
#' labels (or the default) and are not used for inference.
#'
#' Rejection rates from a pilot that has not reached equilibrium are noisy and
#' can be biased, and the equal-rate ladder is only guaranteed to maximise the
#' round-trip rate under the assumptions listed in \code{\link{ptDiagnostics}}.
#' Check the final ladder with \code{\link{ptDiagnostics}} on the real run.
#' @inheritParams callMDI
#' @param betas The starting ladder, for example from \code{\link{ptLadder}}.
#' @param rounds Number of pilot runs.
#' @param R_pilot Iterations of each pilot run.
#' @param ... Further arguments passed to \code{\link{callMDI}} (for example
#' \code{K}, \code{initial_labels} or \code{prior}).
#' @return A list with \code{betas} (the adapted ladder) and \code{history} (a
#' list with the ladder and the estimated rejection rates of each round).
#' @export
adaptLadder <- function(X, types, betas, rounds = 3, R_pilot = 500, ...) {
  betas <- .mdirCheckLadder(betas)
  if (length(betas) < 2) {
    stop("A ladder needs at least two temperatures.", call. = FALSE)
  }
  history <- vector("list", rounds)
  for (i in seq_len(rounds)) {
    pilot <- callMDI(X, R = R_pilot, thin = max(1L, R_pilot %/% 2L), types = types,
      betas = betas, save_parameters = FALSE, check_prior = FALSE, ...)
    d <- ptDiagnostics(pilot)
    rates <- d$pairs$rejection_rate
    history[[i]] <- list(betas = betas, rejection_rate = rates)
    if (anyNA(rates)) {
      stop("A pilot run made no exchange attempts for some pair; increase `R_pilot`.", call. = FALSE)
    }
    betas <- tuneLadder(betas, rates)
  }
  list(betas = betas, history = history)
}
