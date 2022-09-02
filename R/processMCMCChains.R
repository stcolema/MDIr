#' @title Process MCMC chains
#' @description Applies a burn in to and finds a point estimate for each of the
#' chains outputted from ``runMCMCChains``.
#' @param mcmc_lst Output from ``runMCMCChains``
#' @param burn The number of MCMC samples to drop as part of a burn in.
#' @param point_estimate_method Summary statistic used to define the point
#' estimate. Must be ``'mean'`` or ``'median'``. ``'median'`` is the default.
#' @param construct_psm Logical indicating if PSMs be constructed in the
#' unsupervised views. Defaults to FALSE. If TRUE the PSM is constructed and
#' this is used to infer the point estimate rather than the sampled partitions.
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
#' @export
processMCMCChains <- function(mcmc_lst, burn,
                              point_estimate_method = "median",
                              construct_psm = FALSE) {
  new_output <- lapply(
    mcmc_lst,
    processMCMCChain,
    burn,
    point_estimate_method,
    construct_psm
  )

  # Return the MCMC object with burn in applied and point estimates found
  new_output
}
