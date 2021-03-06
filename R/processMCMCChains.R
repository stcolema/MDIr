#' @title Process MCMC chains
#' @description Applies a burn in to and finds a point estimate for each of the
#' chains outputted from ``runMCMCChains``.
#' @param mcmc_lst Output from ``runMCMCChains``
#' @param burn The number of MCMC samples to drop as part of a burn in.
#' @returns A named list similar to the output of
#' ``runMCMCChains`` with some additional entries:
#' 
#'  * ``allocation_probability``: $(N x K)$ matrix. The point estimate of
#'  the allocation probabilities for each data point to each class.
#'
#'  * ``prob``: $N$ vector. The point estimate of the probability of being
#'  allocated to the class with the highest probability.
#'
#'  * ``pred``: $N$ vector. The predicted class for each sample.
#'  
#'  @importFrom parallel parLapply
#' @export
processMCMCChains <- function(mcmc_lst, burn,
                              point_estimate_method = "median") {
  new_output <- lapply(
    mcmc_lst,
    processMCMCChain,
    burn,
    point_estimate_method
  )

  # Return the MCMC object with burn in applied and point estimates found
  new_output
}
