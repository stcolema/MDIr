#' @title Run MCMC Chains Saving to File
#' @description Run multiple chains of Multiple Dataset Integration (MDI) using
#' the same inputs in each model run.
#' @param X Data to cluster. List of $L$ matrices with the $N$ items to cluster
#' held in rows.
#' @param n_chains Integer. Number of MCMC chains to run.
#' @param n_iter The number of iterations in the sampler.
#' @param thin The factor by which the samples generated are thinned, e.g. if
#' ``thin=50`` only every 50th sample is kept.
#' @param types Character vector indicating density type to use. 'MVN'
#' (multivariate normal), 'TAGM' (t-adjust Gaussian mixture) or 'C' (categorical).
#' @param dir_path A string for the directory to save samples to.
#' @param K_max Vector indicating the number of components to include (the upper
#' bound on the number of clusters in each dataset).
#' @param control A \code{\link{mdirControl}} object: the proposal windows of the
#' Gaussian process hyperparameters and whether and for how long they are adapted
#' during burn in.
#' @param R Deprecated; use \code{n_iter}.
#' @param K Deprecated; use \code{K_max}.
#' @param initial_labels Initial clustering. $N x L$ matrix.
#' @param fixed Which items are fixed in their initial label. $N x L$ matrix.
#' @param alpha The concentration parameter for the stick-breaking prior and the
#' weights in the model.
#' @param initial_labels_as_intended Logical indicating if the passed initial
#' labels are as intended or should ``generateInitialLabels`` be called.
#' @param proposal_windows Deprecated; use \code{control = mdirControl(proposal_windows = ...)}.
#' @param time_points See \code{\link{callMDI}}.
runMCMCChainsSavingToFile <- function(X,
                          n_chains,
                          n_iter,
                          thin,
                          types,
                          dir_path,
                          K_max = NULL,
                          initial_labels = NULL,
                          fixed = NULL,
                          alpha = NULL,
                          initial_labels_as_intended = FALSE,
                          proposal_windows = NULL,
                          time_points = NULL,
                          control = mdirControl(),
                          R = NULL,
                          K = NULL) {
  n_iter <- .mdirRenamedArg(if (missing(n_iter)) NULL else n_iter, R, "n_iter", "R", "runMCMCChainsSavingToFile", required = TRUE)
  K_max <- .mdirRenamedArg(if (missing(K_max)) NULL else K_max, K, "K_max", "K", "runMCMCChainsSavingToFile")
  K <- K_max
  control <- .mdirResolveControl(
    missing(control), control,
    if (!missing(proposal_windows)) list(proposal_windows = proposal_windows) else list(),
    "runMCMCChainsSavingToFile"
  )
  proposal_windows <- control$proposal_windows
  auto_tune <- control$auto_tune
  n_burn <- .mdirResolveAutoTune(auto_tune, control$n_burn, n_iter)
  model_call <- vector("list", n_chains)
  
  chain_dirs <- paste0(dir_path, "/Chain_", seq(1, n_chains))
  
  for(ii in seq(1, n_chains)) {
    directory_missing <- ! dir.exists(chain_dirs[ii])
    if(directory_missing) {
      dir.create(chain_dirs[ii])
    }
  model_call[[ii]] <- callMDIWritingToFile(X,
            n_iter,
            thin,
            types,
            K_max = K_max,
            initial_labels = initial_labels,
            fixed = fixed,
            alpha = alpha,
            initial_labels_as_intended = initial_labels_as_intended,
            dir_path = chain_dirs[ii],
            time_points = time_points,
            control = control
    )
  model_call[[ii]]$Chain <- ii
  
  }
  
  model_call
}
