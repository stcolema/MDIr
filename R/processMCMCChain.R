#' @title Process MCMC chain
#' @description Applies a burn in to and finds a point estimate for the output
#' of ``batchSemiSupervisedMixtureModel``. 
#' @param mcmc_output Output from ``batchSemiSupervisedMixtureModel``
#' @param burn The number of MCMC samples to drop as part of a burn in.
#' @param point_estimate_method Summary statistic used to define the point 
#' estimate. Must be ``'mean'`` or ``'median'``. ``'median'`` is the default.
#' @returns A named list similar to the output of 
#' ``batchSemiSupervisedMixtureModel`` with some additional entries:
#'  * ``allocation_probability``: $(N x K)$ matrix. The point estimate of 
#'  the allocation probabilities for each data point to each class.
#'  
#'  * ``prob``: $N$ vector. The point estimate of the probability of being 
#'  allocated to the class with the highest probability.
#'  
#'  * ``pred``: $N$ vector. The predicted class for each sample.
#'  
#' @export
#' @importFrom stats median
processMCMCChain <- function(mcmc_output, burn, point_estimate_method = "median") {
  
  # Dimensions of the dataset
  N <- mcmc_output$N
  P <- mcmc_output$P
  K <- mcmc_output$K
  V <- mcmc_output$V
  
  multiple_views <- V > 1
  
  # The type of mixture model used
  types <- mcmc_output$types
  
  # Indices for views
  view_inds <- seq(1, V)
  
  # MCMC iterations and thinning
  R <- mcmc_output$R
  thin <- mcmc_output$thin
  
  # Is the output semisupervised
  is_semisupervised <- mcmc_output$Semisupervised
  
  # What summary statistic is used to define our point estimates
  use_median <- point_estimate_method == "median"
  use_mean <- point_estimate_method == "mean"
  wrong_method <- ! (use_median | use_mean)
  if(wrong_method)
    stop("Wrong point estimate method given. Must be one of 'mean' or 'median'")
  
  # We burn the floor of burn / thin of these 
  eff_burn <- floor(burn / thin) + 1
  
  # We record only the floor of R / thin samples
  eff_R <- floor(R / thin) - eff_burn
  
  # The indices dropped as part of the burn in
  dropped_indices <- seq(1, eff_burn)
  
  new_output <- mcmc_output

  new_output$mass <- mcmc_output$mass[-dropped_indices, , drop = F]
  if(multiple_views) {
    # The information sharing parameters
    new_output$phis <- mcmc_output$phis[-dropped_indices, , drop = F]
  }
    
  # The model fit
  new_output$complete_likelihood <- mcmc_output$complete_likelihood[-dropped_indices] # , , drop = F]
  new_output$evidence <- mcmc_output$evidence[-dropped_indices[-eff_burn]]
  
  # The allocations and allocation probabilities
  new_output$allocations <- mcmc_output$allocations[-dropped_indices, , , drop = F]
  
  new_output$allocation_probability <- vector("list", V)
  new_output$allocation_probabilities <- vector("list", V)
  new_output$prob <- vector("list", V)
  new_output$pred <- vector("list", V)
  
  for(v in view_inds) {
    if(is_semisupervised[v]) {
      
      # Drop the allocation probabilities from the warm up period
      new_output$allocation_probabilities[[v]] <- mcmc_output$allocation_probabilities[[v]][, , -dropped_indices, drop = F]
      
      # The estimate of the allocation probability matrix, the probability of the 
      # most probable class and the predicted class
      new_output$allocation_probability[[v]] <- .alloc_prob <- 
        calcAllocProb(new_output, v, 
          method = point_estimate_method
        )
      
      new_output$prob[[v]] <- apply(.alloc_prob, 1, max)
      new_output$pred[[v]] <- apply(.alloc_prob, 1, which.max)
    }
  }
  
  # Record the applied burn in
  new_output$burn <- burn 
  
  # Return the MCMC object with burn in applied and point estimates found
  new_output
}
