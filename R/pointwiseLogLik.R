#' @title Pointwise log-likelihood of the fitted items
#' @description The log-likelihood of each item at each saved draw, with the
#' item's component assignments in every view summed out. The result is the
#' draws x items matrix that \code{loo::waic()} and \code{loo::loo()} take
#' (Vehtari, Gelman and Gabry, 2017); for several chains it carries the chain
#' of each row in \code{attr(., "chain_id")}, for \code{loo::relative_eff()}.
#'
#' Because the assignments are summed out, the likelihood of an item depends
#' on the parameters only (the weights, \eqn{\phi} and the component
#' parameters), not on any item-specific latent variable. This is the form that
#' information criteria and leave-one-out cross-validation need. It is the
#' marginal likelihood of the whole model, including the coupling between views:
#' with weights \eqn{w}, \eqn{\phi} and per-view component likelihoods
#' \eqn{g_l(k) = p(x_{nl} \mid \theta_{lk})}, it is \eqn{Z(w g) / Z(w)}, where
#' \eqn{Z} is the MDI normalising constant and \eqn{w g} the elementwise product
#' (outliers are marginalised within \eqn{g}). It is computed exactly and
#' reuses the recursion that gives \eqn{Z}.
#'
#' Items with an observed label in a view (semi-supervised) contribute the
#' joint density of their data and of that label.
#'
#' Pointwise values are only recorded when the model is fitted with
#' \code{save_pointwise = TRUE}.
#' @param object Output of \code{\link{callMDI}} (one chain) or
#' \code{\link{runMCMCChains}} / \code{\link{fitMDI}} (several chains), before
#' \code{\link{processMCMCChain}}.
#' @param burn Number of iterations to discard as burn in (the initial state is
#' always dropped).
#' @return A matrix with a row per saved draw (all chains stacked) and a column
#' per item.
#' @references Vehtari, A., Gelman, A. and Gabry, J. (2017). Practical Bayesian
#' model evaluation using leave-one-out cross-validation and WAIC.
#' \emph{Statistics and Computing}, 27, 1413-1432.
#' @seealso \code{\link{predictMDI}}
#' @export
#' @examples
#' set.seed(1)
#' X <- lapply(1:2, function(v) {
#'   m <- matrix(rnorm(40 * 2, rep(c(0, 3), each = 20)), 40, 2)
#'   rownames(m) <- 1:40
#'   m
#' })
#' fit <- callMDI(X, R = 100, thin = 5, types = c("MVN", "MVN"), K = c(4, 4),
#'                save_pointwise = TRUE)
#' ll <- pointwiseLogLik(fit, burn = 20)
#' dim(ll)
#'
#' # With the loo package: loo::waic(ll)
pointwiseLogLik <- function(object, burn = 0) {
  chains <- .mdirAsChains(object)
  out <- lapply(seq_along(chains), function(i) {
    ch <- chains[[i]]
    if (is.null(ch$pointwise_likelihood)) {
      stop("Pointwise log-likelihoods were not saved. Re-run with `save_pointwise = TRUE`.",
           call. = FALSE)
    }
    draws <- .mdirDraws(ch, burn)
    ch$pointwise_likelihood[draws, , drop = FALSE]
  })
  chain_id <- rep(seq_along(out), vapply(out, nrow, integer(1)))
  ll <- do.call(rbind, out)
  colnames(ll) <- chains[[1]]$sample_ids
  attr(ll, "chain_id") <- chain_id
  ll
}

# Internal: a single chain or a list of chains as a list of chains
.mdirAsChains <- function(object) {
  if (!is.null(object$allocations)) {
    return(list(object))
  }
  if (is.list(object) && length(object) > 0 && !is.null(object[[1]]$allocations)) {
    return(unclass(object)[seq_along(object)])
  }
  stop("`object` must be the output of callMDI(), runMCMCChains() or fitMDI().", call. = FALSE)
}

# Internal: indices of the saved draws kept after a burn in. The initial state
# (saved draw 1) is always dropped. A chain that has been through
# processMCMCChain() has already dropped its burn in.
.mdirDraws <- function(chain, burn) {
  n_saved <- length(chain$joint_likelihood)
  if (!is.null(chain$pred)) {
    return(seq_len(n_saved))
  }
  first_kept <- floor(burn / chain$thin) + 2
  if (first_kept > n_saved) {
    stop("No saved iterations remain after the burn in of ", burn, ".", call. = FALSE)
  }
  seq(first_kept, n_saved)
}
