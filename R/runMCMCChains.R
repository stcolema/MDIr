#' @title Run MCMC Chains
#' @description Run multiple chains of Multiple Dataset Integration (MDI) using
#' the same inputs in each model run.
#' @param X Data to cluster. List of $L$ matrices with the $N$ items to cluster 
#' held in rows.
#' @param n_chains Integer. Number of MCMC chains to run.
#' @param initial_labels Initial clustering. $N x L$ matrix.
#' @param fixed Which items are fixed in their initial label. $N x L$ matrix.
#' @param R The number of iterations in the sampler.
#' @param thin The factor by which the samples generated are thinned, e.g. if
#' ``thin=50`` only every 50th sample is kept.
#' @param type Character vector indicating density type to use. 'MVN'
#' (multivariate normal), 'TAGM' (t-adjust Gaussian mixture) or 'C' (categorical).
#' @param K_max Vector indicating the number of components to include (the upper
#' bound on the number of clusters in each dataset).
#' @param alpha The concentration parameter for the stick-breaking prior and the
#' weights in the model.
#' @param initial_labels_as_intended Logical indicating if the passed initial 
#' labels are as intended or should ``generateInitialLabels`` be called.
#' @return A named list containing the sampled partitions, component weights and
#' phi parameters, model fit measures and some details on the model call.
#' @export
runMCMCChains <- function(X,
                          n_chains,
                          R,
                          thin,
                          initial_labels,
                          fixed,
                          types,
                          K,
                          alpha = NULL,
                          initial_labels_as_intended = FALSE) {
  mcmc_lst <- vector("list", n_chains)
  
  mcmc_lst <- lapply(mcmc_lst, function(x) {
    callMDI(X,
            R,
            thin,
            initial_labels,
            fixed,
            types,
            K = K,
            alpha = alpha,
            initial_labels_as_intended = initial_labels_as_intended
    )
  })
  
  # Record chain number 
  for(ii in seq(n_chains)) {
    mcmc_lst[[ii]]$Chain <- ii
  }
  
  mcmc_lst
}
