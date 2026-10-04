#' @title Predict new items from a fitted MDI model
#' @description For items that were not in the fitted data, computes at each
#' saved draw the posterior probability that each belongs to each component of
#' each view, and the log-likelihood of its data. The new items are treated as
#' further exchangeable draws from the model at that draw's weights,
#' \eqn{\phi} and component parameters, so the calculation is exact given the
#' draw (it uses the same recursion as the normalising constant) and averages
#' over the posterior by averaging over draws. Nothing is refitted.
#'
#' A new item may be missing entries or whole views; the likelihood then uses its
#' observed entries only, and an item observed in a single view is classified
#' using the other views' information through \eqn{\phi} only to the extent that
#' the coupling allows (with no observed data in a view its probabilities there
#' come from the weights and the coupling alone).
#'
#' Component labels of a view are only identified through its observed items
#' when it is semi-supervised. The averaged \code{class_probability} is
#' therefore returned for semi-supervised views only; for other views it is
#' \code{NULL}, because averaging a label across draws (and chains) that permute
#' the labels would be meaningless. Use \code{coclustering} instead, the
#' probability that a new item shares a component with each fitted item, which
#' does not depend on labelling. In an over-fitted semi-supervised view only the
#' columns of the observed classes are identified.
#' @param object Output of \code{\link{callMDI}}, \code{\link{runMCMCChains}} or
#' \code{\link{fitMDI}}, fitted with \code{save_parameters = TRUE}, before
#' \code{\link{processMCMCChain}}.
#' @param X The data the model was fitted to (list of matrices). The chains do not
#' store it; it is needed to rebuild the densities and their data-driven
#' hyperparameters.
#' @param newdata The new items: a list with a matrix per view (same columns as
#' \code{X}, any number of rows common to all views), or a matrix for a
#' single view.
#' @param burn Number of iterations to discard as burn in (the initial state is
#' always dropped).
#' @param n_draws If given, use this many of the kept draws, evenly spaced
#' across chains. Defaults to all of them.
#' @param coclustering Also return the co-clustering probabilities with the
#' fitted items. Costs memory (new items x fitted items, per view).
#' @return A list with
#' \describe{
#'   \item{\code{class_probability}}{a list with, for each semi-supervised view, a
#'   matrix of the probability of each component (columns) for each new item
#'   (rows); \code{NULL} for other views.}
#'   \item{\code{coclustering}}{(if requested) a list with, for each view, the
#'   matrix of the probability that each new item (rows) is in the same component
#'   as each fitted item (columns).}
#'   \item{\code{log_predictive}}{the log of the posterior predictive density of
#'   each new item, the log of the mean of the likelihood over draws.}
#'   \item{\code{log_likelihood}}{the draws x new items matrix of log-likelihoods
#'   behind it.}
#' }
#' @seealso \code{\link{pointwiseLogLik}}
#' @export
#' @examples
#' set.seed(2)
#' make_view <- function(n) {
#'   m <- matrix(rnorm(n * 2, rep(c(0, 4), each = n / 2)), n, 2)
#'   rownames(m) <- seq_len(n)
#'   m
#' }
#' X <- list(make_view(40), make_view(40))
#' fit <- callMDI(X, R = 200, thin = 5, types = c("MVN", "MVN"), K = c(4, 4))
#'
#' new <- list(matrix(c(0.2, -0.1, 4.1, 3.8), 2, 2, byrow = TRUE),
#'             matrix(c(0.1, 0.3, 4.2, 4.0), 2, 2, byrow = TRUE))
#' pred <- predictMDI(fit, X, new, burn = 100, coclustering = TRUE)
#' round(pred$coclustering[[1]][, 1:5], 2)
#' pred$log_predictive
predictMDI <- function(object,
                       X,
                       newdata,
                       burn = 0,
                       n_draws = NULL,
                       coclustering = FALSE) {
  chains <- .mdirAsChains(object)
  first <- chains[[1]]
  if (!is.null(first$pred)) {
    stop("`object` has been through processMCMCChain(); pass the chains as fitted.", call. = FALSE)
  }
  V <- first$V
  X <- .asViewList(X)
  if (is.matrix(newdata) || is.data.frame(newdata)) {
    newdata <- list(as.matrix(newdata))
  }
  newdata <- .asViewList(newdata)
  if (length(X) != V || length(newdata) != V) {
    stop("`X` and `newdata` must each have one matrix for each of the ", V, " views.", call. = FALSE)
  }
  if (nrow(X[[1]]) != first$N) {
    stop("`X` does not match the data the chains were fitted to.", call. = FALSE)
  }
  for (v in seq_len(V)) {
    if (ncol(newdata[[v]]) != ncol(X[[v]])) {
      stop("View ", v, " of `newdata` has ", ncol(newdata[[v]]), " columns; the model was fitted to ",
           ncol(X[[v]]), ".", call. = FALSE)
    }
    if (nrow(newdata[[v]]) != nrow(newdata[[1]])) {
      stop("All views of `newdata` must have the same number of items.", call. = FALSE)
    }
  }
  if (any(vapply(newdata, function(m) any(is.infinite(m)), logical(1)))) {
    stop("`newdata` contains infinite values. Use NA to mark missing entries.", call. = FALSE)
  }
  for (ch in chains) {
    if (is.null(ch$parameters) || any(vapply(ch$parameters, function(p) is.null(p) || length(p) == 0, logical(1)))) {
      stop("Component parameters were not saved. Re-run with `save_parameters = TRUE`.", call. = FALSE)
    }
  }

  # Pool the kept draws of all chains
  pool <- do.call(rbind, lapply(seq_along(chains), function(i) {
    cbind(chain = i, draw = .mdirDraws(chains[[i]], burn))
  }))
  if (!is.null(n_draws)) {
    if (!is.numeric(n_draws) || length(n_draws) != 1 || n_draws < 1) {
      stop("`n_draws` must be a single positive integer.", call. = FALSE)
    }
    if (n_draws < nrow(pool)) {
      pool <- pool[unique(round(seq(1, nrow(pool), length.out = n_draws))), , drop = FALSE]
    }
  }
  S <- nrow(pool)

  K <- first$K
  K_max <- max(K)
  weights <- array(0, c(S, K_max, V))
  for (s in seq_len(S)) {
    ch <- chains[[pool[s, "chain"]]]
    weights[s, , ] <- ch$weights[pool[s, "draw"], , ]
  }
  allocations <- array(0, c(S, first$N, V))
  for (s in seq_len(S)) {
    ch <- chains[[pool[s, "chain"]]]
    allocations[s, , ] <- ch$allocations[pool[s, "draw"], , ]
  }
  phis <- do.call(rbind, lapply(seq_len(S), function(s) {
    chains[[pool[s, "chain"]]]$phis[pool[s, "draw"], ]
  }))
  outlier_weights <- do.call(rbind, lapply(seq_len(S), function(s) {
    chains[[pool[s, "chain"]]]$outlier_weights[pool[s, "draw"], ]
  }))
  parameters <- lapply(seq_len(V), function(v) {
    do.call(rbind, lapply(seq_len(S), function(s) {
      chains[[pool[s, "chain"]]]$parameters[[v]][pool[s, "draw"], ]
    }))
  })
  phis <- matrix(phis, nrow = S)
  outlier_weights <- matrix(outlier_weights, nrow = S)

  prior <- if (is.null(first$prior)) mdiPrior() else first$prior
  density_prior <- if (is.null(first$density_prior)) densityPrior() else first$density_prior
  codes <- .typeCodes(first$types)
  raw <- predictNewItemsCpp(
    X, newdata, as.integer(K), codes$density, codes$outlier, parameters,
    weights, phis, outlier_weights, allocations, isTRUE(coclustering),
    as.numeric(prior), as.numeric(density_prior)
  )

  new_ids <- rownames(newdata[[1]])
  if (is.null(new_ids)) new_ids <- seq_len(nrow(newdata[[1]]))
  fitted_ids <- first$sample_ids
  if (is.null(fitted_ids)) fitted_ids <- seq_len(first$N)

  log_lik <- raw$log_likelihood
  colnames(log_lik) <- new_ids
  log_predictive <- apply(log_lik, 2, function(x) {
    m <- max(x)
    if (!is.finite(m)) return(m)
    m + log(mean(exp(x - m)))
  })

  semisupervised <- first$Semisupervised
  class_probability <- lapply(seq_len(V), function(v) {
    if (!isTRUE(semisupervised[v])) return(NULL)
    pr <- t(raw$class_probability[[v]][seq_len(K[v]), , drop = FALSE])
    dimnames(pr) <- list(new_ids, NULL)
    pr
  })

  out <- list(class_probability = class_probability)
  if (isTRUE(coclustering)) {
    out$coclustering <- lapply(raw$coclustering, function(m) {
      dimnames(m) <- list(new_ids, fitted_ids)
      m
    })
  }
  out$log_predictive <- log_predictive
  out$log_likelihood <- log_lik
  out
}
