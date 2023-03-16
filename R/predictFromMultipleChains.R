#' @title Predict from multiple MCMC chains
#' @description Applies a burn in to and finds a point estimate by combining
#' multiple chains of ``callMDI``.
#' @param mcmc_outputs Output from ``runMCMCChains``
#' @param burn The number of MCMC samples to drop as part of a burn in.
#' @param point_estimate_method Summary statistic used to define the point
#' estimate. Must be ``'mean'`` or ``'median'``. ``'median'`` is the default.
#' @param construct_psm Logical indicating if PSMs be constructed in the
#' unsupervised views. Defaults to FALSE. If TRUE the PSM is constructed and
#' this is used to infer the point estimate rather than the sampled partitions.
#' @param chains_already_processed Logical indicating if the the chains have
#' been processed already.
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
#' @importFrom salso salso
#' @examples
#' \donttest{
#' N <- 100
#' X <- matrix(c(rnorm(N, 0, 1), rnorm(N, 3, 1)), ncol = 2, byrow = TRUE)
#' Y <- matrix(c(rnorm(N, 0, 1), rnorm(N, 3, 1)), ncol = 2, byrow = TRUE)
#'
#' truth <- c(rep(1, N / 2), rep(2, N / 2))
#' data_modelled <- list(X, Y)
#'
#' V <- length(data_modelled)
#'
#' # MCMC parameters
#' R <- 5000
#' thin <- 50
#' burn <- 1000
#'
#' K_max <- 10
#' K <- rep(K_max, V)
#' types <- rep("G", V)
#'
#' n_chains <- 3
#' mcmc_out <- runMCMCChains(data_modelled, n_chains, R, thin, types, K = K)
#' predictFromMultipleChains(mcmc_out, burn)
#' }
#' @export
predictFromMultipleChains <- function(mcmc_outputs,
                                      burn,
                                      point_estimate_method = "median",
                                      construct_psm = FALSE,
                                      chains_already_processed = FALSE) {
  if (chains_already_processed) {
    processed_chains <- mcmc_outputs
  } else {
    # Process the chains, making point estimates and applying a burn-in
    processed_chains <- processMCMCChains(mcmc_outputs, burn, point_estimate_method, construct_psm)
  }

  # The number of chains
  n_chains <- length(processed_chains)
  chain_indices <- seq(1, n_chains)

  # This is used to derive some characteristics of the model run
  first_chain <- processed_chains[[1]]

  # Dimensions of the dataset
  V <- first_chain$V
  N <- first_chain$N
  P <- first_chain$P
  K <- first_chain$K

  # Flag for mixture model vs MDI
  multiple_views <- V > 1

  # MCMC call
  R <- first_chain$R
  thin <- first_chain$thin

  # The type of mixture model used
  types <- first_chain$types

  # Flag indicating if model contains GP
  gp_used <- types %in% c("GP", "TAGPM")

  # The prior on the concentration
  alpha <- first_chain$alpha

  # Indices for views
  view_inds <- seq(1, V)

  # Is the output semisupervised and overfitted
  is_semisupervised <- first_chain$Semisupervised
  is_overfitted <- first_chain$Overfitted

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
  merged_outputs$Semisupervised <- is_semisupervised
  merged_outputs$Overfitted <- is_overfitted

  # I use lists rather than arrays for many of the outputs. This change is not
  # ideal (should build more consistently), but I prefer lists.
  merged_outputs$allocations <- vector("list", V)
  merged_outputs$allocation_probability <- vector("list", V)
  merged_outputs$prob <- vector("list", V)
  merged_outputs$pred <- vector("list", V)
  merged_outputs$outliers <- vector("list", V)
  merged_outputs$N_k <- vector("list", V)

  if (construct_psm) {
    merged_outputs$psm <- vector("list", V)
  }

  # Various parameters of model runs and processing
  merged_outputs$R <- R
  merged_outputs$thin <- thin
  merged_outputs$burn <- burn
  merged_outputs$n_chains <- n_chains

  merged_outputs$Point_estimate <- point_estimate_method

  # Key statistics of the data modelled
  merged_outputs$N <- N
  merged_outputs$P <- P
  merged_outputs$V <- V
  merged_outputs$K <- K

  # Some modelling decisions
  merged_outputs$types <- types
  merged_outputs$alpha <- alpha

  # Hyperparameters for GPs
  merged_outputs$hypers <- vector("list", V)

  # The mass hyper parameter of the component weights and the phi parameters for
  # MDI
  merged_outputs$phis <- do.call(rbind, lapply(processed_chains, function(x) x$phis))
  merged_outputs$mass <- do.call(rbind, lapply(processed_chains, function(x) x$mass))

  # Some model fit measures (not really working for MDI)
  merged_outputs$complete_likelihood <- as.matrix(do.call(
    c,
    lapply(processed_chains, function(x) x$complete_likelihood)
  ))

  merged_outputs$evidence <- as.matrix(do.call(
    c,
    lapply(processed_chains, function(x) x$evidence)
  ))

  # The total time for running each chain (this assumes chains are run serially)
  for (ii in chain_indices) {
    if (ii == 1) {
      merged_outputs$Time <- processed_chains[[ii]]$Time
    } else {
      merged_outputs$Time <- merged_outputs$Time + processed_chains[[ii]]$Time
    }
  }

  merged_outputs$weights <- list()
  for (v in view_inds) {
    merged_outputs$weights[[v]] <- do.call(rbind, lapply(processed_chains, function(x) x$weights[, , v, drop = TRUE]))
  }

  first_chain <- TRUE
  for (v in view_inds) {
    current_view_is_semi_supervised <- is_semisupervised[v]
    current_view_is_overfitted <- is_overfitted[v]

    merged_outputs$allocation_probability[[v]] <- .alloc_prob <- matrix(
      0,
      N,
      K[v]
    )
    for (ii in chain_indices) {
      .curr_chain <- processed_chains[[ii]]
      in_first_chain <- (ii == 1)

      if (in_first_chain) {
        merged_outputs$allocations[[v]] <- .alloc <- .curr_chain$allocations[, , v, drop = TRUE]
        merged_outputs$outliers[[v]] <- .outliers <- .curr_chain$outliers[, , v, drop = TRUE]
        merged_outputs$N_k[[v]] <- .n_k <- t(.curr_chain$N_k[, v, , drop = TRUE])
      } else {
        .prev <- .alloc
        .current <- .curr_chain$allocations[, , v, drop = TRUE]
        merged_outputs$allocations[[v]] <- .alloc <- rbind(.prev, .current)

        .out_prev <- .outliers
        .out_current <- .curr_chain$outliers[, , v, drop = TRUE]
        merged_outputs$allocations[[v]] <- .outliers <- rbind(.out_prev, .out_current)

        .n_k_prev <- .n_k
        .n_k_current <- t(.curr_chain$N_k[, v, , drop = TRUE])
        merged_outputs$N_k[[v]] <- .n_k <- rbind(.n_k_prev, .n_k_current)
      }

      if (current_view_is_semi_supervised) {
        .prev <- .alloc_prob
        .curr <- .curr_chain$allocation_probability[[v]]

        .alloc_prob <- .prev + .curr
      }
    }

    if (construct_psm) {
      merged_outputs$psm[[v]] <- .psm <- createSimilarityMat(.alloc)
    }

    if (current_view_is_semi_supervised) {

      # Normalise the probabilities
      .alloc_prob <- .alloc_prob / n_chains

      merged_outputs$allocation_probability[[v]] <- .alloc_prob
      merged_outputs$prob[[v]] <- .prob <- apply(.alloc_prob, 1, max)
      merged_outputs$pred[[v]] <- apply(.alloc_prob, 1, which.max)

      if (current_view_is_overfitted) {
        merged_outputs$pred[[v]] <- suppressWarnings(salso::salso(.alloc))
      }
    } else {
      merged_outputs$pred[[v]] <- suppressWarnings(salso::salso(.alloc))
    }

    if (gp_used[v]) {
      merged_outputs$hypers[[v]]$amplitude <- do.call(rbind, lapply(processed_chains, function(x) x$hypers[[v]]$amplitude))
      merged_outputs$hypers[[v]]$length <- do.call(rbind, lapply(processed_chains, function(x) x$hypers[[v]]$length))
      merged_outputs$hypers[[v]]$noise <- do.call(rbind, lapply(processed_chains, function(x) x$hypers[[v]]$noise))
    } else {
      # merged_outputs$hypers[[v]] <- NULL
    }
  }

  # Find the fusion probabilities across views
  # If using a mixture model (i.e. V = 1), this returns a NULL
  merged_outputs$fusion_probabilities <- NULL
  if (multiple_views) {
    VC2 <- choose(V, 2)
    merged_outputs$fusion_probabilities <- vector("list", VC2)
    entry <- 0
    names <- c()
    for (v in seq(1, V - 1)) {
      for (w in seq(v + 1, V)) {
        name <- paste0("fused_items_", v, w)
        names <- c(names, name)
        entry <- entry + 1
        merged_outputs$fusion_probabilities[[entry]] <- rep(0, N)
        for (ii in seq(1, n_chains)) {
          merged_outputs$fusion_probabilities[[entry]] <- merged_outputs$fusion_probabilities[[entry]] +
            calcFusionProbabiliy(processed_chains[[ii]], v, w)
        }
        merged_outputs$fusion_probabilities[[entry]] <- merged_outputs$fusion_probabilities[[entry]] / n_chains
      }
    }
    names(merged_outputs$fusion_probabilities) <- names
  }

  merged_outputs
}
