#' @title Check number of samples
#' @description Internal function to check the number of MCMC samples.
#' @param R The number of iterations in the sampler.
#' @param thin The factor by which the samples generated are thinned, e.g. if
#' ``thin=50`` only every 50th sample is kept.
#' @param verbose Logical inficating if warnings should be printed.
#' @return No return value, called for side effects.
#' @examples
#' R <- 10000
#' thin <- 50
#' checkNumberOfSamples(R, thin)
checkNumberOfSamples <- function(R, thin, verbose = FALSE) {
  if (R < thin) {
    stop("Iterations to run less than thinning factor. No samples recorded.")
  }

  number_iterations_saved <- floor(R / thin)
  if ((number_iterations_saved < 20) & verbose) {
    warning("Number of saved iterations (before applying a burn in) is less than 20.")
  }
  NULL
}
