#' @title Process Mixture Model Chain
#' @description Process MCMC output from callMixtureModel, applying burn-in and 
#' computing point estimates. Wrapper around processMCMCChain that ensures 
#' correct handling of single-view outputs.
#' @param mcmc_output Output from `callMixtureModel`
#' @param burn Number of MCMC samples to discard as burn-in
#' @param point_estimate_method Summary statistic for point estimates: 
#' `'mean'` or `'median'` (default)
#' @param construct_psm Logical. If TRUE, construct posterior similarity matrix 
#' and use it for point estimates rather than sampled partitions
#' @returns List with processed MCMC output including point estimates
#' @export
#' @examples
#' N <- 100
#' X <- matrix(rnorm(N * 2, c(0, 3)), ncol = 2)
#' mcmc_out <- callMixtureModel(X, R = 100, thin = 5, type = "MVN", K = 10)
#' processed <- processMixtureModelChain(mcmc_out, burn = 10)
processMixtureModelChain <- function(mcmc_output,
                                     burn,
                                     point_estimate_method = "median",
                                     construct_psm = FALSE) {
  
  # Validate input is from callMixtureModel
  is_mixture_output <- is.null(mcmc_output$V) || mcmc_output$V == 1
  
  if (!is_mixture_output) {
    stop(paste(
      "Input appears to be from callMDI (V > 1).",
      "Use processMCMCChain() instead."
    ))
  }
  
  # Reshape arrays to match processMCMCChain expectations
  # callMixtureModel flattens arrays, we need to add view dimension back
  
  # Allocations: should be [iterations, N, V]
  if (length(dim(mcmc_output$allocations)) == 2) {
    mcmc_output$allocations <- array(
      mcmc_output$allocations,
      dim = c(nrow(mcmc_output$allocations), 
              ncol(mcmc_output$allocations), 
              1)
    )
  }
  
  # Weights: should be [iterations, K, V]
  if (length(dim(mcmc_output$weights)) == 2) {
    mcmc_output$weights <- array(
      mcmc_output$weights,
      dim = c(nrow(mcmc_output$weights), 
              ncol(mcmc_output$weights), 
              1)
    )
  }
  
  # Outliers: should be [iterations, N, V]
  if (length(dim(mcmc_output$outliers)) == 2) {
    mcmc_output$outliers <- array(
      mcmc_output$outliers,
      dim = c(nrow(mcmc_output$outliers), 
              ncol(mcmc_output$outliers), 
              1)
    )
  }
  
  # N_k: should be [iterations, V, K]
  if (length(dim(mcmc_output$N_k)) == 2) {
    # Current shape is [K, iterations], need [iterations, V=1, K]
    mcmc_output$N_k <- array(
      t(mcmc_output$N_k),  # Transpose to [iterations, K]
      dim = c(nrow(mcmc_output$weights), 1, ncol(mcmc_output$weights))
    )
  }
  
  # Complete likelihood: should be [iterations, V]
  if (length(mcmc_output$complete_likelihood) == nrow(mcmc_output$weights)) {
    mcmc_output$complete_likelihood <- matrix(
      mcmc_output$complete_likelihood,
      ncol = 1
    )
  }
  
  # Allocation probabilities: should be list of length V
  if (!is.list(mcmc_output$allocation_probabilities)) {
    mcmc_output$allocation_probabilities <- list(mcmc_output$allocation_probabilities)
  }
  
  # Add V if missing
  if (is.null(mcmc_output$V)) {
    mcmc_output$V <- 1
  }
  
  # Make types a vector
  if (!is.null(mcmc_output$type) && length(mcmc_output$type) == 1) {
    mcmc_output$types <- mcmc_output$type
  }
  
  # Process using existing function
  processed <- processMCMCChain(
    mcmc_output,
    burn,
    point_estimate_method,
    construct_psm
  )
  
  # Flatten PSM back to matrix if constructed
  if (construct_psm && is.list(processed$psm)) {
    processed$psm <- processed$psm[[1]]
  }
  
  processed
}

#' @title Process Mixture Model Chains
#' @description Process multiple MCMC chains from mixture model runs
#' @param mcmc_lst List of outputs from `callMixtureModel` or 
#' `runMCMCChains` with V=1
#' @param burn Number of MCMC samples to discard as burn-in
#' @param point_estimate_method Summary statistic: `'mean'` or `'median'`
#' @param construct_psm Logical. Construct PSMs for point estimates
#' @returns List of processed chains
#' @export
processMixtureModelChains <- function(mcmc_lst,
                                      burn,
                                      point_estimate_method = "median",
                                      construct_psm = FALSE) {
  
  lapply(
    mcmc_lst,
    processMixtureModelChain,
    burn,
    point_estimate_method,
    construct_psm
  )
}