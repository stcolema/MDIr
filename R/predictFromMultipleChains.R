#' @title Predict from multiple MCMC chains
#' @description Applies a burn in to and finds a point estimate by combining
#' multiple chains of ``callMDI``.
#' @param mcmc_outputs Output from ``runMCMCChains``
#' @param burn The number of MCMC samples to drop as part of a burn in.
#' @param point_estimate_method Summary statistic used to define the point
#' estimate. Must be ``'mean'`` or ``'median'``. ``'median'`` is the default.
#' @returns A named list of quantities related to prediction/clustering:
#'
#'  * ``allocation_probability``: List with an $(N x K)$ matrix for each
#'    semi-supervised view. The point estimate of the allocation probabilities for
#'    each data point to each class.
#'
#'  * ``prob``: List with an $N$ vector for each semi-supervised view. The point
#'    estimate of the probability of being allocated to the class with the
#'    highest probability.
#'
#'  * ``pred``: List of $N$ vectorsfor each semi-supervised view. The predicted
#'    class for each sample.
#'
#'  * ``allocations``: List of sampled allocations for each view. Columns
#'    correspond to items being clustered, rows to MCMC samples.
#'
#' @export
predictFromMultipleChains <- function(mcmc_outputs,
                                      burn,
                                      point_estimate_method = "median",
                                      chains_already_processed = FALSE) {

  
  
  if(chains_already_processed ) {
    processed_chains <- mcmc_outputs
  } else {
    # Process the chains, making point estimates and applying a burn-in
    processed_chains <- processMCMCChains(mcmc_outputs, burn, point_estimate_method)
  }

  # The number of chains
  n_chains <- length(processed_chains)
  chain_indices <- seq(1, n_chains)

  # This is used to derive some characteristics of the model run
  first_chain <- processed_chains[[1]]

  # Dimensions of the dataset
  V <- first_chain$V
  N <- first_chain$N
  K <- first_chain$K

  # MCMC call
  R <- first_chain$R
  thin <- first_chain$thin

  # The type of mixture model used
  types <- first_chain$types

  # Indices for views
  view_inds <- seq(1, V)

  # Is the output semisupervised
  is_semisupervised <- first_chain$Semisupervised

  # What summary statistic is used to define our point estimates
  use_median <- point_estimate_method == "median"
  use_mean <- point_estimate_method == "mean"
  wrong_method <- !(use_median | use_mean)
  if (wrong_method) {
    stop("Wrong point estimate method given. Must be one of 'mean' or 'median'")
  }

  # We burn the floor of burn / thin of these
  eff_burn <- floor(burn / thin)

  # We record only the floor of R / thin samples
  eff_R <- floor(R / thin) - eff_burn

  # The indices dropped as part of the burn in
  dropped_indices <- seq(1, eff_burn)

  # Setup the output list
  merged_outputs <- list()
  merged_outputs$allocations <- vector("list", V)
  merged_outputs$allocation_probability <- vector("list", V)
  merged_outputs$prob <- vector("list", V)
  merged_outputs$pred <- vector("list", V)

  merged_outputs$R <- R
  merged_outputs$thin <- thin
  merged_outputs$burn <- burn
  merged_outputs$n_chains <- n_chains

  merged_outputs$Point_estimate <- point_estimate_method

  merged_outputs$N <- N
  merged_outputs$V <- V
  merged_outputs$K <- K

  first_chain <- TRUE
  for (v in view_inds) {
    current_view_is_semi_supervised <- is_semisupervised[v]
    
    merged_outputs$allocation_probability[[v]] <- .alloc_prob <- matrix(
      0,
      N,
      K[v]
    )
    for (ii in chain_indices) {
      .curr_chain <- processed_chains[[ii]]
      in_first_chain <- (ii == 1)

      if (in_first_chain) {
        merged_outputs$allocations[[v]] <- .curr_chain$allocations[, , v, drop = TRUE]
      } else {
        .prev <- merged_outputs$allocations[[v]]
        .current <- .curr_chain$allocations[, , v, drop = TRUE]
        merged_outputs$allocations[[v]] <- rbind(.prev, .current)
      }

      if (current_view_is_semi_supervised) {
        .prev <- .alloc_prob
        .curr <- .curr_chain$allocation_probability[[v]]

        .alloc_prob <- .prev + .curr
      }
    }

    if (current_view_is_semi_supervised) {

      # Normalise the probabilities
      .alloc_prob <- .alloc_prob / n_chains

      merged_outputs$allocation_probability[[v]] <- .alloc_prob

      merged_outputs$prob[[v]] <- .prob <- apply(.alloc_prob, 1, max)
      merged_outputs$pred[[v]] <- apply(.alloc_prob, 1, which.max)
    }
  }

  merged_outputs
}
