# Conversion of fits to the draws formats of the posterior package, so that the
# summaries, diagnostics and plots of the Stan workflow (posterior, bayesplot,
# loo) apply to mdir output. The methods are registered when posterior is loaded.

# Internal: the scalar quantities that do not depend on how clusters are
# labelled (or, optionally, labelled-dependent ones the user asks for), as a draws
# x variables matrix for the stored samples `index` of one chain.
.mdirDrawsMatrix <- function(x, index, weights = FALSE, allocations = FALSE) {
  V <- x$V
  cols <- list()
  cols[["log_lik_complete"]] <- as.numeric(x$complete_likelihood)[index]
  if (!is.null(x$observed_likelihood)) {
    cols[["log_lik_observed"]] <- as.numeric(x$observed_likelihood)[index]
  }
  if (!is.null(x$joint_likelihood)) {
    cols[["log_lik_joint"]] <- as.numeric(x$joint_likelihood)[index]
  }
  for (v in seq_len(V)) {
    cols[[sprintf("mass[%d]", v)]] <- x$mass[index, v]
  }
  if (V > 1 && !is.null(x$phis)) {
    pairs <- utils::combn(V, 2)
    for (i in seq_len(ncol(pairs))) {
      cols[[sprintf("phi[%d,%d]", pairs[1, i], pairs[2, i])]] <- x$phis[index, i]
    }
  }
  for (v in seq_len(V)) {
    cols[[sprintf("occupied[%d]", v)]] <- .mdirOccupied(x, v, index)
    if (!is.null(x$outlier_weights) && any(x$outlier_weights[, v] != 0)) {
      cols[[sprintf("outlier_weight[%d]", v)]] <- x$outlier_weights[index, v]
    }
    pooled <- x$pooled_hyperparameters[[v]]
    if (!is.null(pooled) && length(pooled) > 0) {
      for (j in seq_len(ncol(pooled))) {
        cols[[sprintf("pooled[%d,%d]", v, j)]] <- pooled[index, j]
      }
    }
  }
  if (weights) {
    # Normalised weights sorted within each draw, so that they do not depend on the labelling
    for (v in seq_len(V)) {
      w <- matrix(x$weights[index, seq_len(x$K[v]), v], nrow = length(index))
      w <- t(apply(w / rowSums(w), 1, sort, decreasing = TRUE))
      w <- matrix(w, nrow = length(index))
      for (k in seq_len(ncol(w))) {
        cols[[sprintf("sorted_weight[%d,%d]", v, k)]] <- w[, k]
      }
    }
  }
  if (allocations) {
    # Labelled-dependent unless the view is semi-supervised; numbered from 1
    for (v in seq_len(V)) {
      z <- matrix(x$allocations[index, , v], nrow = length(index)) + 1
      for (n in seq_len(ncol(z))) {
        cols[[sprintf("allocation[%d,%d]", v, n)]] <- z[, n]
      }
    }
  }
  do.call(cbind, cols)
}

.mdirDrawsArray <- function(chains, burn = NULL, weights = FALSE, allocations = FALSE) {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("This needs the 'posterior' package.", call. = FALSE)
  }
  mats <- lapply(chains, function(ch) {
    idx <- .mdirRetained(ch, burn)$index
    # The initial state is not a draw from the chain
    if (!.mdirIsProcessed(ch)) idx <- idx[idx > 1]
    .mdirDrawsMatrix(ch, idx, weights, allocations)
  })
  n_iter <- vapply(mats, nrow, integer(1))
  if (length(unique(n_iter)) != 1) {
    stop("The chains retain different numbers of draws.", call. = FALSE)
  }
  arr <- array(
    unlist(mats, use.names = FALSE),
    dim = c(n_iter[1], ncol(mats[[1]]), length(mats)),
    dimnames = list(iteration = NULL, variable = colnames(mats[[1]]), chain = NULL)
  )
  posterior::as_draws_array(aperm(arr, c(1, 3, 2)))
}

.mdirChains <- function(x) {
  if (inherits(x, "mdir_fit_list")) unclass(x)[seq_along(x)] else list(x)
}

#' @title Convert an mdir fit to the draws formats of the posterior package
#' @description Methods of \code{posterior::as_draws()} and its relatives for
#' the output of \code{\link{callMDI}}, \code{\link{runMCMCChains}},
#' \code{\link{fitMDI}} and \code{\link{processMCMCChains}}, and for weighted
#' ensembles from \code{\link{smcMDI}}. With the result, the Stan workflow carries
#' over: \code{posterior::summarise_draws()} (mean, quantiles, Rhat, bulk and
#' tail ESS), \code{posterior::subset_draws()}, \code{posterior::resample_draws()},
#' and the \pkg{bayesplot} functions that take draws. Per-item log-likelihoods
#' for \pkg{loo} come from \code{\link{pointwiseLogLik}}.
#'
#' A mixture model has no parameters that mean the same thing in every draw
#' (components can swap labels), so the variables are restricted to quantities
#' that do not depend on the labelling:
#' \describe{
#'   \item{\code{log_lik_complete}, \code{log_lik_observed}, \code{log_lik_joint}}{
#'   The log-likelihoods described in \code{\link{callMDI}}. These are not
#'   \code{lp__}: they omit the priors.}
#'   \item{\code{mass[v]}}{Concentration of view \code{v}.}
#'   \item{\code{phi[l,m]}}{Association of views \code{l} and \code{m}.}
#'   \item{\code{occupied[v]}}{Number of occupied components in view \code{v}.}
#'   \item{\code{outlier_weight[v]}}{Outlier weight of a view with an outlier component.}
#'   \item{\code{pooled[v,j]}}{The pooled hyperparameters of view \code{v} (variance
#'   scales, or the Gaussian-process population parameters; see \code{\link{densityPrior}}).}
#'   \item{\code{sorted_weight[v,k]}}{With \code{weights = TRUE}: the normalised
#'   component weights of a view sorted from largest to smallest.}
#'   \item{\code{allocation[v,n]}}{With \code{allocations = TRUE}: the component
#'   (from 1) of item \code{n} in view \code{v}. These depend on the labelling, so
#'   their Rhat and ESS are meaningful only for semi-supervised views, whose
#'   observed classes fix the labels.}
#' }
#' The normalising constant \eqn{Z} is left out: the overall scale of each view's
#' unnormalised weights does not affect the model, so \eqn{Z} wanders with the scale and
#' its Rhat and ESS say nothing about the posterior.
#'
#' The initial state of an unprocessed chain is not a draw and is dropped, as
#' is the burn in.
#'
#' For an \code{\link{smcMDI}} ensemble the draws of all particles are stacked
#' (particle by particle) and carry the particle weights, which
#' \pkg{posterior} then uses in its summaries; the quantities are the data
#' log-likelihood, \code{mass}, \code{phi} and \code{occupied}. The
#' draws of one particle are not independent, so \pkg{posterior}'s effective sample sizes do not apply.
#' See \code{\link{smcDiagnostics}}.
#' @param x Output of \code{\link{callMDI}} (\code{mdir_fit}), of
#' \code{\link{runMCMCChains}} / \code{\link{fitMDI}} /
#' \code{\link{processMCMCChains}} (\code{mdir_fit_list}), or of
#' \code{\link{smcMDI}} (\code{mdir_smc}).
#' @param burn Number of iterations treated as burn in. The default is half of
#' \code{R}, as in \code{\link{summary.mdir_fit}}. Ignored for chains that have
#' been through \code{\link{processMCMCChain}}.
#' @param weights Include the sorted normalised component weights.
#' @param allocations Include the item allocations (large: one variable per item
#' and view).
#' @param ... Unused.
#' @return A draws object of the requested class: \code{as_draws()} and
#' \code{as_draws_array()} give a \code{draws_array} (iterations x chains x
#' variables), the others the corresponding format.
#' @seealso \code{\link{plot.mdir_fit_list}}, \code{\link{assessConvergence}},
#' \code{\link{pointwiseLogLik}}
#' @name as_draws.mdir_fit
#' @examples
#' if (requireNamespace("posterior", quietly = TRUE)) {
#'   set.seed(1)
#'   X <- lapply(1:2, function(v) {
#'     m <- matrix(rnorm(40 * 2, rep(c(0, 3), each = 20)), 40, 2)
#'     rownames(m) <- 1:40
#'     m
#'   })
#'   fit <- runMCMCChains(X, n_chains = 2, R = 400, thin = 5, types = c("G", "G"), K = c(4, 4))
#'   draws <- posterior::as_draws(fit)
#'   posterior::summarise_draws(draws)
#' }
NULL

#' @rdname as_draws.mdir_fit
#' @exportS3Method posterior::as_draws
as_draws.mdir_fit <- function(x, burn = NULL, weights = FALSE, allocations = FALSE, ...) {
  .mdirDrawsArray(.mdirChains(x), burn, weights, allocations)
}

#' @rdname as_draws.mdir_fit
#' @exportS3Method posterior::as_draws
as_draws.mdir_fit_list <- as_draws.mdir_fit

#' @rdname as_draws.mdir_fit
#' @exportS3Method posterior::as_draws_array
as_draws_array.mdir_fit <- as_draws.mdir_fit

#' @rdname as_draws.mdir_fit
#' @exportS3Method posterior::as_draws_array
as_draws_array.mdir_fit_list <- as_draws.mdir_fit

#' @rdname as_draws.mdir_fit
#' @exportS3Method posterior::as_draws_df
as_draws_df.mdir_fit <- function(x, ...) posterior::as_draws_df(as_draws.mdir_fit(x, ...))

#' @rdname as_draws.mdir_fit
#' @exportS3Method posterior::as_draws_df
as_draws_df.mdir_fit_list <- as_draws_df.mdir_fit

#' @rdname as_draws.mdir_fit
#' @exportS3Method posterior::as_draws_matrix
as_draws_matrix.mdir_fit <- function(x, ...) posterior::as_draws_matrix(as_draws.mdir_fit(x, ...))

#' @rdname as_draws.mdir_fit
#' @exportS3Method posterior::as_draws_matrix
as_draws_matrix.mdir_fit_list <- as_draws_matrix.mdir_fit

#' @rdname as_draws.mdir_fit
#' @exportS3Method posterior::as_draws_list
as_draws_list.mdir_fit <- function(x, ...) posterior::as_draws_list(as_draws.mdir_fit(x, ...))

#' @rdname as_draws.mdir_fit
#' @exportS3Method posterior::as_draws_list
as_draws_list.mdir_fit_list <- as_draws_list.mdir_fit

#' @rdname as_draws.mdir_fit
#' @exportS3Method posterior::as_draws_rvars
as_draws_rvars.mdir_fit <- function(x, ...) posterior::as_draws_rvars(as_draws.mdir_fit(x, ...))

#' @rdname as_draws.mdir_fit
#' @exportS3Method posterior::as_draws_rvars
as_draws_rvars.mdir_fit_list <- as_draws_rvars.mdir_fit

# Weighted draws of an smcMDI ensemble: one row per (particle, draw), particle
# varying slowest, with the particle weights
.mdirSmcDraws <- function(x) {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("This needs the 'posterior' package.", call. = FALSE)
  }
  D <- x$n_draws
  P <- x$n_particles
  V <- x$V
  per <- function(arr3) matrix(aperm(arr3, c(1, 3, 2)), nrow = D * P)
  cols <- list(log_lik_complete = as.numeric(x$data_log_likelihood))
  mass <- per(x$mass)
  for (v in seq_len(V)) cols[[sprintf("mass[%d]", v)]] <- mass[, v]
  if (V > 1) {
    pairs <- utils::combn(V, 2)
    phi <- per(x$phis)
    for (i in seq_len(ncol(pairs))) {
      cols[[sprintf("phi[%d,%d]", pairs[1, i], pairs[2, i])]] <- phi[, i]
    }
  }
  for (v in seq_len(V)) {
    lab <- x$allocations[, , v, , drop = FALSE]
    dim(lab) <- c(D, x$N, P)
    lab <- matrix(aperm(lab, c(1, 3, 2)), nrow = D * P, ncol = x$N)
    cols[[sprintf("occupied[%d]", v)]] <- apply(lab, 1, function(r) length(unique(r)))
  }
  m <- posterior::as_draws_matrix(do.call(cbind, cols))
  posterior::weight_draws(m, weights = smcWeights(x, "draw"))
}

#' @rdname as_draws.mdir_fit
#' @exportS3Method posterior::as_draws
as_draws.mdir_smc <- function(x, ...) .mdirSmcDraws(x)

#' @rdname as_draws.mdir_fit
#' @exportS3Method posterior::as_draws_matrix
as_draws_matrix.mdir_smc <- function(x, ...) .mdirSmcDraws(x)

#' @rdname as_draws.mdir_fit
#' @exportS3Method posterior::as_draws_df
as_draws_df.mdir_smc <- function(x, ...) posterior::as_draws_df(.mdirSmcDraws(x))
