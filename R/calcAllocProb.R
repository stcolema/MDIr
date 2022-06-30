#' @title Calculate allocation probabilities
#' @description Calculate the empirical allocation probability for each class
#' based on the sampled allocation probabilities.
#' @param mcmc_samples Output from ``callMDI``.
#' @param view Which view to calculate the allocation probabilities for.
#' @param burn The number of samples to discard.
#' @param method The point estimate to use. ``method = 'mean'`` or
#' ``method = 'median'``. ``'median'`` is the default.
#' @return An N x K matrix of class probabilities.
#' @export
calcAllocProb <- function(mcmc_samples, view, burn = 0, method = "median") {
  R <- mcmc_samples$R
  thin <- mcmc_samples$thin
  V <- mcmc_samples$V
  
  .alloc <- mcmc_samples$allocation_probabilities[[view]]
  
  if(burn > 0) {
    if(burn > R) {
      stop("Burn in exceeds number of iterations run.")
    }
    
    eff_burn <- floor(burn / thin)
    dropped_samples <- seq(1, eff_burn)
    .alloc <- .alloc[, , -dropped_samples]
  }
  
  if(view > V) {
    .err <- paste(
      "Requested view not in ``mcmc_samples``. Please check that the requested",
      "view is less than or equal to V, the number of views modelled."
      )
    stop(.err)
  }
  probs <- NULL
  
  if(method == "median") {
    probs <- apply(.alloc, c(1, 2), median)
  }
  if(method == "mean") {
    probs <- rowSums(.alloc, dims = 2) / dim(.alloc)[3]
  }
  if(length(probs) == 1) {
    if(is.null(probs)) {
      stop("``method`` must be one of 'mean' or 'median'")
    }
  }
  probs
}