#' @title Run MCMC Chains
#' @description Run multiple chains of Multiple Dataset Integration (MDI) using
#' the same inputs in each model run.
#' @param X Data to cluster. List of $L$ matrices with the $N$ items to cluster
#' held in rows.
#' @param n_chains Integer. Number of MCMC chains to run.
#' @param R The number of iterations in the sampler.
#' @param thin The factor by which the samples generated are thinned, e.g. if
#' ``thin=50`` only every 50th sample is kept.
#' @param types Character vector indicating density type to use. 'MVN'
#' (multivariate normal), 'TAGM' (t-adjust Gaussian mixture) or 'C' (categorical).
#' @param K Vector indicating the number of components to include (the upper
#' bound on the number of clusters in each dataset).
#' @param initial_labels Initial clustering. $N x L$ matrix.
#' @param fixed Which items are fixed in their initial label. $N x L$ matrix.
#' @param alpha The concentration parameter for the stick-breaking prior and the
#' weights in the model.
#' @param initial_labels_as_intended Logical indicating if the passed initial
#' labels are as intended or should ``generateInitialLabels`` be called.
#' @param proposal_windows List of the proposal windows for the Metropolis-Hastings
#' sampling of Gaussian process hyperparameters. Each entry corresponds to a
#' view. For views modelled using a Gaussian process, the first entry is the
#' proposal window for the ampltiude, the second is for the length-scale and the
#' third is for the noise. These are not used in other mixture types.
#' @return A named list containing the sampled partitions, component weights and
#' phi parameters, model fit measures and some details on the model call.
#' @examples
#' \donttest{
#' N <- 100
#' X <- matrix(c(rnorm(N, 0, 1), rnorm(N, 3, 1)), ncol = 2, byrow = TRUE)
#' Y <- matrix(c(rnorm(N, 0, 1), rnorm(N, 3, 1)), ncol = 2, byrow = TRUE)
#'
#' truth <- c(rep(1, N / 2), rep(2, N / 2))
#' data_modelled <- list(X, Y)
#'
#' V <- length(data_modelled)
#'
#' # MCMC parameters
#' R <- 5000
#' thin <- 50
#' burn <- 1000
#'
#' K_max <- 10
#' K <- rep(K_max, V)
#' types <- rep("G", V)
#'
#' n_chains <- 3
#' mcmc_out <- runMCMCChains(data_modelled, n_chains, R, thin, types, K = K)
#' }
#' @export
runMCMCChains <- function(X,
                          n_chains,
                          R,
                          thin,
                          types,
                          K = NULL,
                          initial_labels = NULL,
                          fixed = NULL,
                          alpha = NULL,
                          initial_labels_as_intended = FALSE,
                          proposal_windows = NULL) {
  mcmc_lst <- vector("list", n_chains)

  mcmc_lst <- lapply(mcmc_lst, function(x) {
    callMDI(X,
      R,
      thin,
      types,
      K = K,
      initial_labels = initial_labels,
      fixed = fixed,
      alpha = alpha,
      initial_labels_as_intended = initial_labels_as_intended,
      proposal_windows = proposal_windows
    )
  })

  # Record chain number
  for (ii in seq(n_chains)) {
    mcmc_lst[[ii]]$Chain <- ii
  }

  mcmc_lst
}
