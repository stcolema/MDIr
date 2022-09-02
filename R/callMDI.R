#' @title Call MDI
#' @description Runs a MCMC chain of the integrative clustering method,
#' Multiple Dataset Integration (MDI), to L datasets. It is recommended that
#' L < 5.
#' @param X Data to cluster. List of matrices with the N items to cluster held
#' in rows.
#' @param initial_labels Initial clustering. $N x L$ matrix.
#' @param fixed Which items are fixed in their initial label. $N x L$ matrix.
#' @param R The number of iterations in the sampler.
#' @param thin The factor by which the samples generated are thinned, e.g. if
#' ``thin=50`` only every 50th sample is kept.
#' @param type Character vector indicating density type to use. 'MVN'
#' (multivariate normal), 'TAGM' (t-adjust Gaussian mixture) or 'C' (categorical).
#' @param K_max Vector indicating the number of components to include (the upper
#' bound on the number of clusters in each dataset).
#' @param alpha The concentration parameter for the stick-breaking prior and the
#' weights in the model.
#' @param initial_labels_as_intended Logical indicating if the passed initial 
#' labels are as intended or should ``generateInitialLabels`` be called.
#' @return A named list containing the sampled partitions, component weights,
#' phi and mass parameters, model fit measures and some details on the model call.
#' @examples 
#' 
#' N <- 100
#' X <- matrix(c(rnorm(100, 0, 1), rnorm(100, 3, 1)), ncol = 2)
#' Y <- matrix(c(rnorm(100, 0, 1), rnorm(100, 3, 1)), ncol = 2)
#' data_modelled <- list(X, Y)
#' 
#' R <- 100
#' thin <- 5
#' 
#' alpha <- c(1, 1)
#' K <- c(10, 15)
#' labels <- matrix(0, nrow = N, ncol = 2)
#' labels[, 1] <- generateInitialUnsupervisedLabels(N, alpha[1], K[1])
#' labels[, 2] <- generateInitialUnsupervisedLabels(N, alpha[1], K[2])
#' 
#' fixed <- matrix(0, nrow = N, ncol = 2)
#' types <- c("MVN", "G")
#' 
#' mcmc_out <- callMDI(data_modelled, R, thin, labels, fixed, types, K, alpha)
#' 
#' @export
callMDI <- function(X,
                    R,
                    thin,
                    initial_labels,
                    fixed,
                    types,
                    K,
                    alpha = NULL,
                    initial_labels_as_intended = FALSE) {

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

  # Check that the matrix indicating observed labels is correctly formatted.
  checkFixedInput(fixed, N, V)

  # Translate user input into appropriate types for C++ function
  density_types <- translateTypes(types)
  outlier_types <- setupOutlierComponents(types)

  if(is.null(alpha))
    alpha <- rep(1, V)
  
  # Generate initial labels. Uses the stick-breaking prior if unsupervised,
  # proportions of observed classes is semi-supervised.
  initial_labels <- generateInitialLabels(initial_labels, fixed, K, alpha,
    labels_as_intended = initial_labels_as_intended
  )
  # for(v in seq(V))
  #   checkLabels(initial_labels[, v], K[v])
  
  t_0 <- Sys.time()
  
  # Pull samples from the MDI model
  mcmc_output <- runAltMDI(
    R,
    thin,
    X,
    K,
    density_types,
    outlier_types,
    initial_labels,
    fixed
  )

  t_1 <- Sys.time()
  time_taken <- t_1 - t_0
  
  # Record details of model run to output
  # MCMC details
  mcmc_output$thin <- thin
  mcmc_output$R <- R
  mcmc_output$burn <- 0

  # Density choice
  mcmc_output$types <- types

  # Dimensions of data
  mcmc_output$P <- P
  mcmc_output$N <- N
  mcmc_output$V <- V

  # Number of components modelled
  mcmc_output$K <- K

  # Record hyperparameter choice
  mcmc_output$alpha <- alpha

  # Indicate if the model was semi-supervised or unsupervised
  mcmc_output$Semisupervised <- is_semisupervised <- apply(fixed, 2, function(x) any(x == 1))
  mcmc_output$Overfitted <- rep(TRUE, V)
  
  for(v in seq(1, V)) {
    if(is_semisupervised[v]) {
      known_labels <- which(fixed[, v] == 1)
      K_fix <- length(unique(initial_labels[known_labels, v]))
      is_overfitted <- (K[v] > K_fix)
      mcmc_output$Overfitted[v] <- is_overfitted
    }
  }

  # Record how long the algorithm took
  mcmc_output$Time <- time_taken
  
  mcmc_output
}
