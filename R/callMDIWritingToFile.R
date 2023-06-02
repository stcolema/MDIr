#' @title Call Multiple Dataset Integration Writing to File
#' @description Runs a MCMC chain of the integrative clustering method,
#' Multiple Dataset Integration (MDI), to V datasets and writes samples to file.
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
#' @param dir_path A string for the directory to save samples to. This cannot be
#' a shorthand containing the `~` symbol.
#' @return A named list containing the sampled partitions, component weights,
#' phi and mass parameters, model fit measures and some details on the model call.
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
#' mcmc_out <- callMDIWritingToFile(data_modelled, R, thin, types, K = K)
#' @export
callMDIWritingToFile <- function(X,
                    R,
                    thin,
                    types,
                    K = NULL,
                    initial_labels = NULL,
                    fixed = NULL,
                    alpha = NULL,
                    initial_labels_as_intended = FALSE,
                    proposal_windows = NULL,
                    dir_path = tempdir()) {
  
  # Check that the R > thin
  checkNumberOfSamples(R, thin)
  
  # Check inputs and translate to C++ inputs
  checkDataCorrectInput(X)
  
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
  
  t_0 <- Sys.time()
  
  # Pull samples from the MDI model
  runMDIWriteToFile(
    R,
    thin,
    X,
    K,
    density_types,
    outlier_types,
    initial_labels,
    fixed,
    proposal_windows,
    dir_path
  )
  
  run_details <- list()
  
  run_details$n_samples <- floor(R / thin) + 1
  run_details$n_param <- N * V + sum(K) + V + choose(V, 2);
  
  t_1 <- Sys.time()
  time_taken <- t_1 - t_0
  
  # Record details of model run to output
  # MCMC details
  run_details$thin <- thin
  run_details$R <- R
  run_details$burn <- 0
  
  # Density choice
  run_details$types <- types
  
  # Dimensions of data
  run_details$P <- P
  run_details$N <- N
  run_details$V <- V
  
  # Number of components modelled
  run_details$K <- K
  
  # Record hyperparameter choice
  run_details$alpha <- alpha
  
  # Indicate if the model was semi-supervised or unsupervised
  run_details$Semisupervised <- is_semisupervised <- apply(fixed, 2, function(x) any(x == 1))
  run_details$Overfitted <- rep(TRUE, V)
  
  # Proposal windows if any used
  run_details$proposal_windows <- proposal_windows
  
  for (v in seq(1, V)) {
    if (is_semisupervised[v]) {
      known_labels <- which(fixed[, v] == 1)
      K_fix <- length(unique(initial_labels[known_labels, v]))
      is_overfitted <- (K[v] > K_fix)
      run_details$Overfitted[v] <- is_overfitted
    }
  }
  
  # Record how long the algorithm took
  run_details$Time <- time_taken
  
  # Where are the samples saved to
  run_details$Save_dir <- dir_path
  
  run_details
}
