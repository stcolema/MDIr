#' @title Read in Saved Samples
#' @description Reads in the saved output of `callMDIWritingToFile` and compiles
#' it into a more easily usable format.
#' @param run_details Output from `callMDIWritingToFile` - a list describing the
#' model run, saved sample location, etc.
#' @return A named list containing the sampled partitions, component weights,
#' phi and mass parameters and some details on the model call.
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
#' model_call <- callMDIWritingToFile(data_modelled, R, thin, types, K = K)
#' mcmc_out <- readInSavedSamples(model_call)
#' @export
readInSavedSamples <- function(run_details) {
  N <- run_details$N
  V <- run_details$V
  K <- run_details$K
  n_samples <- run_details$n_samples
  n_param <- run_details$n_param
  dir_path <- run_details$Save_dir
  .output <- readMCMCsamples(n_samples, n_param, dir_path)

  mcmc_output <- run_details
  mcmc_output$allocations <- array(dim = c(n_samples, N, V))
  mcmc_output$weights <- array(0, dim = c(n_samples, max(K), V))

  # Layout of a saved sample: labels (view by view), weights (view by view),
  # masses, phis, complete likelihood, observed likelihood.
  weight_start <- N * V
  for (v in seq_len(V)) {
    mcmc_output$allocations[, , v] <- .output[, N * (v - 1) + seq_len(N), drop = FALSE]
    mcmc_output$weights[, seq_len(K[v]), v] <- .output[, weight_start + seq_len(K[v]), drop = FALSE]
    weight_start <- weight_start + K[v]
  }

  n_phi <- choose(V, 2)
  mcmc_output$mass <- .output[, weight_start + seq_len(V), drop = FALSE]
  mcmc_output$phis <- .output[, weight_start + V + seq_len(n_phi), drop = FALSE]
  mcmc_output$complete_likelihood <- .output[, weight_start + V + n_phi + 1]
  mcmc_output$observed_likelihood <- .output[, weight_start + V + n_phi + 2]

  mcmc_output
}
