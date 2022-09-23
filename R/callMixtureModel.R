#' @title Call mixture model
#' @description Runs a MCMC chain of a Bayesian mixture model. Essentially a
#' wrapper to allow more intuitive inputs for the single dataset case of MDI.
#' @param X Data to cluster. Matrix with the N items to cluster held
#' in rows.
#' @param R The number of iterations in the sampler.
#' @param thin The factor by which the samples generated are thinned, e.g. if
#' ``thin=50`` only every 50th sample is kept.
#' @param type Character vector indicating density type to use. 'G' (Gaussian
#' with diagonal covariance matrix) 'MVN' (multivariate normal), 'TAGM'
#' (t-adjust Gaussian mixture), 'GP' (MVN with Gaussian process prior on the
#' mean), 'TAGPM' (TAGM with GP prior on the mean), 'C' (categorical).
#' @param K Integer indicating the number of components to include (the upper
#' bound on the number of clusters).
#' @param initial_labels Initial clustering. $N$-vector.
#' @param fixed Which items are fixed in their initial label. $N$-vector.
#' @param alpha The concentration parameter for the stick-breaking prior and the
#' weights in the model.
#' @param initial_labels_as_intended Logical indicating if the passed initial
#' labels are as intended or should ``generateInitialLabels`` be called.
#' @param proposal_windows List of the proposal windows for the Metropolis-Hastings
#' sampling of Gaussian process hyperparameters. Each entry corresponds to a
#' view. For views modelled using a Gaussian process, the first entry is the
#' proposal window for the ampltiude, the second is for the length-scale and the
#' third is for the noise. These are not used in other mixture types.
#' @return A named list containing the sampled partitions, component weights,
#' and mass parameters, model fit measures and some details on the model call.
#' @examples
#' N <- 100
#' X <- matrix(c(rnorm(N, 0, 1), rnorm(N, 3, 1)), ncol = 2, byrow = TRUE)
#'
#' # This R is much too low for real applications
#' R <- 100
#' thin <- 5
#'
#' alpha <- 1
#' K <- 10
#' type <- "MVN"
#'
#' mcmc_out <- callMixtureModel(X, R, thin, type, K = K)
#'
#' @export
callMixtureModel <- function(X,
                             R,
                             thin,
                             type,
                             K = 75,
                             initial_labels = NULL,
                             fixed = NULL,
                             alpha = NULL,
                             initial_labels_as_intended = FALSE,
                             proposal_windows = NULL) {

  # Check that the R > thin
  checkNumberOfSamples(R, thin)

  # Check inputs and translate to C++ inputs
  inputDataValid <- is.matrix(X)
  if (!inputDataValid) {
    stop("Data must be a matrix.")
  }

  # The number of items modelled
  N <- nrow(X)

  # The number of measurements
  P <- ncol(X)

  # Convert to a list for passing to C++
  X <- list(X)

  # Translate user input into appropriate types for C++ function
  density_type <- translateTypes(type)
  outlier_type <- setupOutlierComponents(type)
  gp_used <- type %in% c("GP", "TAGPM")
  
  if (is.null(alpha)) {
    alpha <- 1
  }

  if (is.null(initial_labels)) {
    initial_labels <- rep(1, N)
  }

  if (is.null(fixed)) {
    fixed <- rep(0, N)
  }

  # Generate initial labels. Uses the stick-breaking prior if unsupervised,
  # proportions of observed classes is semi-supervised.
  initial_labels <- matrix(initial_labels, ncol = 1)
  fixed <- matrix(fixed, ncol = 1)
  initial_labels <- generateInitialLabels(initial_labels,
    fixed,
    K,
    alpha,
    labels_as_intended = initial_labels_as_intended
  )
  # for(v in seq(V))
  checkLabels(initial_labels, K)

  proposal_windows <- processProposalWindows(proposal_windows, type)
  
  t_0 <- Sys.time()

  # Pull samples from the MDI model
  mcmc_output <- runMDI(
    R,
    thin,
    X,
    K,
    density_type,
    outlier_type,
    initial_labels,
    fixed,
    proposal_windows
  )

  t_1 <- Sys.time()
  time_taken <- t_1 - t_0

  # Put the outputs into a more intuitive object for a single dataset and remove
  # the phis
  mcmc_output$allocations <- mcmc_output$allocations[, , 1]
  mcmc_output$phis <- NULL
  mcmc_output$weights <- mcmc_output$weights[, , 1]
  mcmc_output$outliers <- mcmc_output$outliers[, , 1]
  mcmc_output$allocation_probabilities <- mcmc_output$allocation_probabilities[[1]]
  mcmc_output$N_k <- mcmc_output$N_k[, 1, ]
  mcmc_output$complete_likelihood <- mcmc_output$complete_likelihood[, 1]


  # Record details of model run to output
  # MCMC details
  mcmc_output$thin <- thin
  mcmc_output$R <- R
  mcmc_output$burn <- 0

  # Density choice
  mcmc_output$type <- type

  # Dimensions of data
  mcmc_output$P <- P
  mcmc_output$N <- N

  # Number of components modelled
  mcmc_output$K <- K

  # Record hyperparameter choice
  mcmc_output$alpha <- alpha

  # Indicate if the model was semi-supervised or unsupervised
  mcmc_output$Semisupervised <- is_semisupervised <- any(fixed == 1)

  mcmc_output$Overfitted <- TRUE
  if (is_semisupervised) {
    known_labels <- which(fixed == 1)
    K_fix <- length(unique(initial_labels[known_labels]))
    is_overfitted <- (K > K_fix)
    mcmc_output$Overfitted <- is_overfitted
  }

  if (gp_used) {
    hypers <- vector("list", 3)
    names(hypers) <- c("amplitude", "length", "noise")
    hypers$amplitude <- mcmc_output$hypers[[1]][, seq(1, K), drop = FALSE]
    hypers$length <- mcmc_output$hypers[[1]][, seq(K + 1, 2 * K), drop = FALSE]
    hypers$noise <- mcmc_output$hypers[[1]][, seq(2 * K + 1, 3 * K), drop = FALSE]
    mcmc_output$hypers <- hypers
  } else {
    mcmc_output$hypers <- NA
  }
  
  # Record how long the algorithm took
  mcmc_output$Time <- time_taken

  mcmc_output
}
