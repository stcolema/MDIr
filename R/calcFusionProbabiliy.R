#' @title Caluclate fusion probability
#' @description The mean of the number of MCMC samples for which each item has 
#' the same label in views v and w.
#' @param mcmc Output from ``runMCMCChain``
#' @param v First view considered.
#' @param w Second view considered.
#' @returns A vector of probabilities
#' @export
calcFusionProbabiliy <- function(mcmc, v, w) {
  V <- mcmc$V
  views <- seq(1, V)
  v_not_in_views <- ! (v %in% views)
  w_not_in_views <- ! (w %in% views)
  v_too_large <- v > V
  w_too_large <- w > V
  v_equals_w <- v == w
  
  bad_input <- (v_not_in_views 
    || w_not_in_views
    || v_too_large
    || w_too_large
    || v_equals_w
  )
  if(bad_input) {
    stop("`v` and `w` must be different natual numbers in the range [1, V].")
  }
  
  colMeans(mcmc$allocation_probabilities[[v]] == mcmc$allocation_probabilities[[w]])
}