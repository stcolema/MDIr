#' @title Calculate fusion probability
#' @description The mean of the number of MCMC samples for which each item has
#' the same label in views v and w.
#' @param mcmc Output from ``runMCMCChain``
#' @param v First view considered.
#' @param w Second view considered.
#' @param processed Has the chain been processed already (defaults to FALSE).
#' @returns A vector of probabilities
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
#' mcmc_out <- callMDI(data_modelled, R, thin, types, K = K)
#' calcFusionProbabiliy(mcmc_out, 1, 2)
#'
#' @export
calcFusionProbabiliy <- function(mcmc, v, w, processed = FALSE) {
  V <- mcmc$V
  views <- seq(1, V)
  v_not_in_views <- !(v %in% views)
  w_not_in_views <- !(w %in% views)
  v_too_large <- v > V
  w_too_large <- w > V
  v_equals_w <- v == w

  bad_input <- (v_not_in_views ||
    w_not_in_views ||
    v_too_large ||
    w_too_large ||
    v_equals_w
  )
  if (bad_input) {
    stop("`v` and `w` must be different natual numbers in the range [1, V].")
  }
  if (processed) {
    out <- colMeans(mcmc$allocations[[v]] == mcmc$allocations[[w]])
  } else {
    out <- colMeans(mcmc$allocations[, , v] == mcmc$allocations[, , w])
  }
  out
}
