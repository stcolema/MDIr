#' @title Compile Consensus Clustering
#' @description Performs consensus clustering on a list of MCMC chains from
#' ``runMCMCChains``.
#' @param mcmc_lst Output of ``runMCMCChains``.
#' @param D The iteration to use from within the chains. Defaults to the largest
#' possible value, i.e., the length of the chains.
#' @param W The number of chains to use in compiling the consensus clustring.
#' Defaults to the length of ``mcmc_lst`` but can be smaller.
#' @param construct_cm Logical indicating if the consensus matrix should be
#' constructed.
#' @param point_estimate_method Point estimate method for allocation
#' probabilities in semi-supervised views. One of "median" or "mean".
#' @param ... Arguments passed to ``salso::salso`` for inferring a point
#' estimate in unsupervised views.
#' @returns A named list similar to the output of ``runMCMCChains``.
#' @examples
#' \donttest{
#' N <- 500
#' K <- 4
#' P <- 10
#'
#' X <- generateSimulationDataset(K, N, P)
#' Y <- generateSimulationDataset(K, N, P)
#' Z <- generateSimulationDataset(K, N, P)
#'
#' row.names(Z$data) <- row.names(Y$data) <- row.names(X$data)
#'
#' data_modelled <- list(X$data, Y$data, Z$data)
#' V <- length(data_modelled)
#'
#' # This R is much too low for real applications
#' R <- 1000
#' thin <- 100
#'
#' K_max <- 50
#' K <- rep(K_max, V)
#'
#' types <- rep("G", V)
#'
#' n_chains <- 500
#' mcmc_lst <- runMCMCChains(data_modelled, n_chains, R, thin, types, K = K)
#' cc <- compileConsensusClustering(mcmc_lst)
#' }
#' @export
compileConsensusClustering <- function(mcmc_lst, D = NULL, W = NULL,
                                       point_estimate_method = "mean",
                                       construct_cm = TRUE,
                                       ...) {
  first_chain <- mcmc_lst[[1]]
  D_ub <- first_chain$R
  thin <- first_chain$thin

  N <- first_chain$N
  P <- first_chain$P
  V <- first_chain$V
  VC2 <- choose(V, 2)
  view_inds <- seq(1, V)

  K <- first_chain$K

  types <- first_chain$types
  gp_used <- types %in% c("GP", "TAGPM")

  is_semisupervised <- first_chain$Semisupervised
  multiple_views <- (V > 1)

  if (is.null(D)) {
    D <- D_ub
  }
  if (D > D_ub) {
    stop("D cannot be greater than ", D_ub, " for this list of chains.")
  }

  W_ub <- length(mcmc_lst)
  if (is.null(W)) {
    W <- W_ub
  }
  if (W > W_ub) {
    stop("W cannot be greater than ", W_ub, " for this list of chains.")
  }

  sample_used <- floor(D / thin) + 1
  width_inds <- seq(1, W)

  cc_lst <- list()
  cc_lst <- first_chain

  cc_lst$allocations <- list()
  cc_lst$weights <- list()
  cc_lst$outliers <- list()
  cc_lst$allocation_probabilities <- list()
  cc_lst$N_k <- list()
  cc_lst$hypers <- list()
  cc_lst$acceptance_count <- list()

  cc_lst$n_chains <- W

  cc_lst$mass <- matrix(0, W, V)
  cc_lst$phis <- matrix(0, W, VC2)
  cc_lst$Time <- cc_lst$evidence <- cc_lst$complete_likelihood <- matrix(0, W, 1)

  for (v in view_inds) {
    cc_lst$outliers[[v]] <- cc_lst$allocations[[v]] <- matrix(0, W, N)
    cc_lst$N_k[[v]] <- cc_lst$weights[[v]] <- matrix(0, W, K[v])
    cc_lst$allocation_probabilities[[v]] <- array(0, c(N, K[v], W))
    cc_lst$hypers[[v]] <- cc_lst$acceptance_count[[v]] <- NA
    if (gp_used[v]) {
      cc_lst$acceptance_count[[v]] <- matrix(0, W, 3 * K[v])
      cc_lst$hypers[[v]] <- vector("list", 3)
      names(cc_lst$hypers[[v]]) <- names(first_chain$hypers[[v]])
      for (ii in seq(1, 3)) {
        cc_lst$hypers[[v]][[ii]] <- matrix(0, W, K[v])
      }
    }
  }

  for (w in width_inds) {
    .mcmc <- mcmc_lst[[w]]

    cc_lst$Time[w] <- .mcmc$Time
    cc_lst$evidence[w] <- .mcmc$evidence[sample_used - 1]
    cc_lst$complete_likelihood[w] <- .mcmc$complete_likelihood[sample_used]
    cc_lst$mass[w, ] <- .mcmc$mass[sample_used, ]
    cc_lst$phis[w, ] <- .mcmc$phis[sample_used, ]

    for (v in view_inds) {
      cc_lst$allocations[[v]][w, ] <- .mcmc$allocations[sample_used, , v]
      cc_lst$outliers[[v]][w, ] <- .mcmc$outliers[sample_used, , v]
      cc_lst$N_k[[v]][w, ] <- .mcmc$N_k[seq(1, K[v]), v, sample_used]
      cc_lst$weights[[v]] <- .mcmc$weights[sample_used, , v]
      cc_lst$allocation_probabilities[[v]][, , w] <- .mcmc$allocation_probabilities[[v]][, , sample_used]

      if (gp_used[v]) {
        cc_lst$acceptance_count[[v]][w, ] <- .mcmc$acceptance_count[[v]]
        for (ii in seq(1, 3)) {
          cc_lst$hypers[[v]][[ii]][w, ] <- .mcmc$hypers[[v]][[ii]][sample_used, ]
        }
      }
    }
  }

  cc_lst$cm <- list()
  for (v in view_inds) {
    if (is_semisupervised[v]) {

      # The estimate of the allocation probability matrix, the probability of the
      # most probable class and the predicted class
      cc_lst$allocation_probability[[v]] <- .alloc_prob <-
        calcAllocProb(cc_lst, v,
          method = point_estimate_method
        )

      cc_lst$prob[[v]] <- apply(.alloc_prob, 1, max)
      cc_lst$pred[[v]] <- apply(.alloc_prob, 1, which.max)
    } else {
      cc_lst$pred[[v]] <- suppressWarnings(salso::salso(cc_lst$allocations[[v]], ...))
    }
    if (construct_cm) {
      cc_lst$cm[[v]] <- createSimilarityMat(cc_lst$allocations[[v]])
    }
  }
  if (multiple_views) {
    cc_lst$fusion_probabilities <- calcFusionProbabiliyAllViews(cc_lst, processed = TRUE)
  }
  cc_lst
}
