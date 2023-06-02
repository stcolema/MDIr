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
  
  n_samples <- run_details$n_samples
  n_param <- run_details$n_param
  dir_path <- run_details$Save_dir
  .output <- readMCMCsamples(n_samples, n_param, dir_path)
  
  mcmc_output <- run_details
  mcmc_output$allocations <- array(dim = c(n_samples, N, V))
  mcmc_output$weights <- array(dim = c(n_samples, max(K), V))
  
  for(v in seq(1, V)) {
    mcmc_output$allocations[, , v] <- .output[, seq(N * (v - 1) + 1, N * v)]
    if(v == 1) {
      mcmc_output$weights[, seq(1, K[v]), v] <- .output[, seq(N * V + 1, N * V + K[v])]
    } else {
      mcmc_output$weights[, seq(1, K[v]), v] <- .output[, seq(N * V + 1 + K[v - 1], N * V + K[v])]
    }
  }
  
  mcmc_output$mass <- .output[, seq(N * V + sum(K), N * V + sum(K) + V)]
  mcmc_output$phis <- .output[, seq(N * V + sum(K) + V, n_param)]
  
  mcmc_output
}
