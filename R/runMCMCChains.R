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
#' @param save_parameters,save_imputed,prior,density_prior See ``callMDI``.
#' @param save_pointwise,phi_update See ``callMDI``.
#' @param betas,swap_scheme,swap_every Parallel tempering within each chain; see
#' \code{\link{callMDI}}.
#' @param n_cores Number of cores on which to run the chains. The default, one,
#' runs them in turn exactly as before; \code{options(mdir.cores = )} changes the
#' default. With more than one, each chain has its own L'Ecuyer-CMRG random
#' number stream (Unix-alikes fork the workers, Windows starts a socket cluster
#' that loads \pkg{mdir}), so the result depends on the seed but not on the
#' number of cores; it is not the same as the serial result for the same seed.
#' The number is capped at the number of chains and of cores available, and at
#' two while a package is being checked (\code{_R_CHECK_LIMIT_CORES_} is set), as
#' CRAN requires. Progress messages are given when the chains have finished.
#' @param verbose Logical. Report the start and end of each chain (with run
#' time) as a message. \code{FALSE} by default here; \code{\link{fitMDI}}, which
#' also assesses convergence, defaults to \code{TRUE}.
#' @return An object of class \code{mdir_fit_list}: a list with one
#' \code{\link{callMDI}} result per chain, each with the chain number in
#' \code{Chain}. It prints as a short report; \code{summary()} gives per-chain
#' detail.
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
                          proposal_windows = NULL,
                          save_parameters = TRUE,
                          save_imputed = FALSE,
                          prior = mdiPrior(),
                          density_prior = densityPrior(),
                          verbose = FALSE,
                          save_pointwise = FALSE,
                          phi_update = c("slice", "gibbs"),
                          n_cores = NULL,
                          betas = 1,
                          swap_scheme = c("deo", "seo"),
                          swap_every = 1L) {
  phi_update <- match.arg(phi_update)
  swap_scheme <- match.arg(swap_scheme)
  if (!is.numeric(n_chains) || length(n_chains) != 1 || is.na(n_chains) || n_chains < 1) {
    stop("`n_chains` must be a single positive integer.", call. = FALSE)
  }
  n_cores <- .mdirResolveCores(n_cores, n_chains)

  # report prior warnings once rather than for every chain
  if (is.null(K)) K_used <- rep(floor(nrow(X[[1]]) / 2), length(X)) else K_used <- K
  for (m in .checkSparsity(X, types, K_used, prior)) message(m)

  fit_one <- function() {
    callMDI(X,
      R,
      thin,
      types,
      K = K,
      initial_labels = initial_labels,
      fixed = fixed,
      alpha = alpha,
      initial_labels_as_intended = initial_labels_as_intended,
      proposal_windows = proposal_windows,
      save_parameters = save_parameters,
      save_imputed = save_imputed,
      prior = prior,
      density_prior = density_prior,
      check_prior = FALSE,
      save_pointwise = save_pointwise,
      phi_update = phi_update,
      betas = betas,
      swap_scheme = swap_scheme,
      swap_every = swap_every
    )
  }

  if (n_cores > 1) {
    if (verbose) {
      message(sprintf("Running %d chains of %d iterations on %d cores...", n_chains, R, n_cores))
    }
    mcmc_lst <- .mdirRunChainsParallel(n_chains, n_cores, fit_one)
    for (ii in seq_len(n_chains)) {
      mcmc_lst[[ii]]$Chain <- ii
      if (verbose) {
        message(sprintf(
          "Chain %d/%d: finished in %s.", ii, n_chains, .mdirFormatTime(mcmc_lst[[ii]]$Time)
        ))
      }
    }
  } else {
    mcmc_lst <- vector("list", n_chains)
    for (ii in seq_len(n_chains)) {
      if (verbose) {
        message(sprintf("Chain %d/%d: running %d iterations...", ii, n_chains, R))
      }
      mcmc_lst[[ii]] <- fit_one()

      # Record chain number
      mcmc_lst[[ii]]$Chain <- ii

      if (verbose) {
        message(sprintf(
          "Chain %d/%d: finished in %s.", ii, n_chains, .mdirFormatTime(mcmc_lst[[ii]]$Time)
        ))
      }
    }
  }

  # A classed list (a plain list of chains underneath), see R/mdirFitMethods.R
  class(mcmc_lst) <- c("mdir_fit_list", "list")
  mcmc_lst
}
