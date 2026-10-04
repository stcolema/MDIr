#' @title Call Multiple Dataset Integration
#' @description Runs a MCMC chain of the integrative clustering method,
#' Multiple Dataset Integration (MDI), to V datasets.
#' @param X Data to cluster. List of matrices with the N items to cluster held
#' in rows.
#' @param R The number of iterations in the sampler.
#' @param thin The factor by which the samples generated are thinned, e.g. if
#' ``thin=50`` only every 50th sample is kept.
#' @param types Character vector indicating density types to use. 'G' (Gaussian
#' with diagonal covariance matrix) 'MVN' (multivariate normal), 'TAGM'
#' (t-adjust Gaussian mixture), 'GP' (MVN with Gaussian process prior on the
#' mean), 'TAGPM' (TAGM with GP prior on the mean), 'C' (categorical).
#' @param K Vector indicating the number of components to include (the upper
#' bound on the number of clusters in each dataset).
#' @param initial_labels Initial clustering. $N x V$ matrix.
#' @param fixed Which items are fixed in their initial label. $N x V$ matrix.
#' @param alpha The concentration parameter for the stick-breaking prior and the
#' weights in the model.
#' @param initial_labels_as_intended Logical indicating if the passed initial
#' labels are as intended or should ``generateInitialLabels`` be called.
#' @param proposal_windows List of the proposal windows for the Metropolis-Hastings
#' sampling of Gaussian process hyperparameters. Each entry corresponds to a
#' view. For views modelled using a Gaussian process, the first entry is the
#' proposal window for the ampltiude, the second is for the length-scale and the
#' third is for the noise. These are not used in other mixture types.
#' @param save_parameters Logical. Record the component parameters (means,
#' covariances, category probabilities, ...) at every saved iteration. These are
#' needed for posterior predictive checks (see ``simulatePosteriorPredictive``).
#' Set to ``FALSE`` to save memory for large models.
#' @param prior MDI-level prior hyperparameters, created by
#' \code{\link{mdiPrior}}. See its documentation for the defaults and how to
#' check them with \code{\link{simulatePriorPredictive}}.
#' @param density_prior Density-level (hierarchical) prior options, created by
#' \code{\link{densityPrior}}: pooling of the variance scale across components and
#' the priors of Gaussian process views, including a floor on the GP length scale.
#' @param check_prior Logical. Report a message if the prior on the weights is on
#' the side of the Rousseau and Mengersen (2011) threshold that duplicates rather
#' than empties superfluous components (see \code{\link{mdiPrior}}). Set
#' \code{options(mdir.quiet = TRUE)} to silence it globally.
#' @param save_imputed Logical. Record the imputed value of every missing entry
#' at every saved iteration (``FALSE`` by default).
#' @param save_pointwise Logical. Record the log-likelihood of every item at every
#' saved iteration (an iterations x N matrix, ``pointwise_likelihood``), for
#' \code{\link{pointwiseLogLik}} and information criteria. ``FALSE`` by default;
#' the total over items (``joint_likelihood``) is always recorded.
#' @param phi_update How the \eqn{\phi} parameters are updated. ``"slice"`` (the
#' default) draws each \eqn{\phi_{lm}} from its conditional with the strategic
#' latent variable integrated out, using one slice-sampling update (Neal, 2003)
#' that needs no tuning. ``"gibbs"`` is the update conditional on the
#' strategic latent variable used before version 0.11. Both have the same target;
#' the slice update mixes faster (see ``NEWS.md``).
#' @param betas Inverse temperatures for parallel tempering: a strictly
#' increasing vector in (0, 1] whose last element is 1 (the posterior), for
#' example from \code{\link{ptLadder}}. One replica of the sampler runs at
#' each temperature and neighbouring replicas exchange states, which lets the
#' chain at \code{beta = 1} cross barriers between modes that a single chain
#' crosses rarely. The default, \code{1}, is the ordinary sampler. A single
#' value other than 1 samples the likelihood-tempered target and is only useful
#' for testing. Tempering requires complete data, \code{"G"}, \code{"MVN"} or
#' \code{"C"} views and no outlier component. It multiplies the run time by
#' \code{length(betas)} and guarantees nothing about finite-time mixing; see
#' \code{\link{ptDiagnostics}}.
#' @param swap_scheme How replicas are paired for exchange. \code{"deo"}
#' (default) alternates deterministically between even and odd neighbouring
#' pairs, the non-reversible scheme of Syed et al. (2022); \code{"seo"} chooses
#' the parity at random (reversible). Both leave the tempered targets invariant.
#' @param swap_every Attempt exchanges after every \code{swap_every} sweeps.
#' @return An object of class \code{mdir_fit}: a named list containing the
#' sampled partitions, component weights, phi and mass parameters, model fit
#' measures and some details on the model call. It prints as a short report
#' (see \code{\link{print.mdir_fit}}); use \code{summary()} for posterior
#' summaries. Missing data (``NA`` entries in ``X``) are treated as missing at random and
#' imputed within the sampler.
#'
#' The per-item allocation probabilities (``allocation_probabilities``, an
#' \eqn{N \times K \times} draws array per view) are recorded only for
#' semi-supervised views, where they give the class probabilities. For other
#' views the entry is ``NULL``; the allocations are in ``allocations``.
#'
#' ``joint_likelihood`` is the log-likelihood of the data under the MDI model
#' with the component assignments of every item summed out, at each saved
#' draw. It differs from ``observed_likelihood``, which sums over each view's
#' components separately with that view's own normalised weights and so ignores
#' the coupling between views (the two coincide when every \eqn{\phi} is zero).
#' In a semi-supervised view the observed labels are treated as data: an item
#' with an observed label contributes the joint density of its data and label.
#' @references Neal, R. M. (2003). Slice sampling. \emph{Annals of Statistics},
#' 31(3), 705-767.
#' @examples
#'
#' N <- 100
#' X <- matrix(c(rnorm(N, 0, 1), rnorm(N, 3, 1)), ncol = 2, byrow = TRUE)
#' Y <- matrix(c(rnorm(N, 0, 1), rnorm(N, 3, 1)), ncol = 2, byrow = TRUE)
#'
#' truth <- c(rep(1, N / 2), rep(2, N / 2))
#' data_modelled <- list(X, Y)
#'
#' V <- length(data_modelled)
#'
#' # This R is much too low for real applications
#' R <- 100
#' thin <- 5
#' burn <- 10
#'
#' K_max <- 10
#' K <- rep(K_max, V)
#' types <- rep("G", V)
#'
#' mcmc_out <- callMDI(data_modelled, R, thin, types, K = K)
#'
#' @export
callMDI <- function(X,
                    R,
                    thin,
                    types,
                    K = NULL,
                    initial_labels = NULL,
                    fixed = NULL,
                    alpha = NULL,
                    initial_labels_as_intended = FALSE,
                    proposal_windows = NULL,
                    save_parameters = TRUE,
                    save_imputed = FALSE,
                    prior = mdiPrior(),
                    density_prior = densityPrior(),
                    check_prior = TRUE,
                    save_pointwise = FALSE,
                    phi_update = c("slice", "gibbs"),
                    betas = 1,
                    swap_scheme = c("deo", "seo"),
                    swap_every = 1L) {

  phi_update <- match.arg(phi_update)
  swap_scheme <- match.arg(swap_scheme)
  betas <- .mdirCheckLadder(betas)
  if (!is.numeric(swap_every) || length(swap_every) != 1 || is.na(swap_every) || swap_every < 1) {
    stop("`swap_every` must be a positive integer.", call. = FALSE)
  }
  if (!is.logical(save_pointwise) || length(save_pointwise) != 1 || is.na(save_pointwise)) {
    stop("`save_pointwise` must be TRUE or FALSE.", call. = FALSE)
  }

  # Check that the R > thin
  checkNumberOfSamples(R, thin)

  # Check inputs and translate to C++ inputs
  checkDataCorrectInput(X, types)

  # The number of items modelled
  N <- nrow(X[[1]])

  # The number of views modelled
  V <- length(X)

  # The number of measurements in each view
  P <- lapply(X, ncol)

  # If no upper bound on K is passed set K to half of N
  if (is.null(K)) {
    K <- rep(floor(N / 2), V)
  }

  no_labels_passed <- is.null(initial_labels)
  if (no_labels_passed) {
    initial_labels <- matrix(1, N, V)
    if (initial_labels_as_intended) {
      stop("Cannot fail to pass initial labels and pass ``initial_labels_as_intended=TRUE``.")
    }
  }

  if (is.null(fixed)) {
    fixed <- matrix(0, N, V)
  }

  # Check that the matrix indicating observed labels is correctly formatted.
  checkFixedInput(fixed, N, V)

  if (check_prior) {
    for (m in .checkSparsity(X, types, K, prior)) message(m)
  }

  # Translate user input into appropriate types for C++ function
  density_types <- translateTypes(types)
  outlier_types <- setupOutlierComponents(types)
  gp_used <- types %in% c("GP", "TAGPM")

  if (is.null(alpha)) {
    alpha <- rep(1, V)
  }

  # Generate initial labels. Uses the stick-breaking prior if unsupervised,
  # proportions of observed classes is semi-supervised.
  initial_labels <- generateInitialLabels(initial_labels, fixed, K, alpha,
    labels_as_intended = initial_labels_as_intended
  )
  # for(v in seq(V))
  #   checkLabels(initial_labels[, v], K[v])

  proposal_windows <- processProposalWindows(proposal_windows, types)

  # The allocation probabilities (N x K x draws) are large and only used for
  # semi-supervised views
  is_semisupervised <- apply(fixed, 2, function(x) any(x == 1))

  t_0 <- Sys.time()

  # Pull samples from the MDI model
  mcmc_output <- runMDI(
    R,
    thin,
    X,
    K,
    density_types,
    outlier_types,
    initial_labels,
    fixed,
    proposal_windows,
    save_parameters,
    save_imputed,
    as.numeric(prior),
    as.numeric(density_prior),
    save_allocation_probabilities = as.integer(is_semisupervised),
    save_pointwise = save_pointwise,
    phi_slice = (phi_update == "slice"),
    betas = betas,
    swap_scheme = as.integer(swap_scheme == "seo"),
    swap_every = as.integer(swap_every)
  )
  
  # Traces are returned as one-column matrices; use plain vectors
  for (nm in c("complete_likelihood", "observed_likelihood", "joint_likelihood", "evidence", "mass_acceptance_rate")) {
    mcmc_output[[nm]] <- as.numeric(mcmc_output[[nm]])
  }
  if (!save_pointwise) {
    mcmc_output["pointwise_likelihood"] <- list(NULL)
  }
  mcmc_output$allocation_probabilities[!is_semisupervised] <- list(NULL)
  mcmc_output$sample_ids <- row.names(X[[1]])

  t_1 <- Sys.time()
  time_taken <- t_1 - t_0

  # Record details of model run to output
  # MCMC details
  mcmc_output$thin <- thin
  mcmc_output$R <- R
  mcmc_output$burn <- 0

  # Density choice
  mcmc_output$types <- types

  # Proportion of missing entries in each view
  mcmc_output$missing_proportion <- vapply(X, function(x) mean(is.na(x)), numeric(1))

  # Dimensions of data
  mcmc_output$P <- P
  mcmc_output$N <- N
  mcmc_output$V <- V

  # Number of components modelled
  mcmc_output$K <- K

  # Record hyperparameter choice
  mcmc_output$alpha <- alpha
  mcmc_output$prior <- prior
  mcmc_output$density_prior <- density_prior
  mcmc_output$phi_update <- phi_update
  mcmc_output$betas <- betas
  mcmc_output$tempering <- .mdirTidyTempering(mcmc_output$tempering, betas, swap_scheme)

  # Indicate if the model was semi-supervised or unsupervised
  mcmc_output$Semisupervised <- is_semisupervised
  mcmc_output$Overfitted <- rep(TRUE, V)

  # Proposal windows if any used
  mcmc_output$proposal_windows <- proposal_windows

  for (v in seq(1, V)) {
    if (is_semisupervised[v]) {
      known_labels <- which(fixed[, v] == 1)
      K_fix <- length(unique(initial_labels[known_labels, v]))
      is_overfitted <- (K[v] > K_fix)
      mcmc_output$Overfitted[v] <- is_overfitted
    }
    if (gp_used[v]) {
      hypers <- vector("list", 3)
      names(hypers) <- c("amplitude", "length", "noise")
      hypers$amplitude <- mcmc_output$hypers[[v]][, seq(1, K[v]), drop = FALSE]
      hypers$length <- mcmc_output$hypers[[v]][, seq(K[v] + 1, 2 * K[v]), drop = FALSE]
      hypers$noise <- mcmc_output$hypers[[v]][, seq(2 * K[v] + 1, 3 * K[v]), drop = FALSE]
      mcmc_output$hypers[[v]] <- hypers
    } else {
      mcmc_output$hypers[[v]] <- NA
    }
  }
  
  # Record how long the algorithm took
  mcmc_output$Time <- time_taken

  # A classed list (every `$` access still works): print() and summary() show
  # a report instead of the sampled arrays, see R/mdirFitMethods.R
  class(mcmc_output) <- c("mdir_fit", class(mcmc_output))

  mcmc_output
}
