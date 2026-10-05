#' @title Get likelihood
#' @description Extracts the model fit traces from the output of ``callMDI`` (or
#' a processed chain): the complete-data log-likelihood (the likelihood of the
#' data given the sampled allocations, summed over views) and the observed-data
#' log-likelihood (each view's items marginalised over their component under the
#' view's own normalised weights; the cross-view coupling is not part of a
#' single view's likelihood).
#' @param mcmc_output The output from ``callMDI`` or ``processMCMCChain``.
#' @return A data.frame with the log-likelihood (``log_likelihood``), its
#' ``type`` and the MCMC ``iteration`` at which it was recorded.
#' @export
getLikelihood <- function(mcmc_output) {
  R <- mcmc_output$R
  thin <- mcmc_output$thin

  traces <- list(complete = mcmc_output$complete_likelihood)
  if (!is.null(mcmc_output$observed_likelihood)) {
    traces$observed <- mcmc_output$observed_likelihood
  }

  do.call(rbind, lapply(names(traces), function(type) {
    trace <- as.numeric(traces[[type]])
    n <- length(trace)
    # Saved samples are at iterations 0, thin, 2 thin, ..., floor(R / thin) thin
    # (R itself only if it is a multiple of thin); a burn in removes the earliest
    last_saved <- floor(R / thin) * thin
    data.frame(
      log_likelihood = trace,
      type = type,
      iteration = seq(last_saved - (n - 1) * thin, last_saved, by = thin)
    )
  }))
}
