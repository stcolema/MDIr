#' @title Fit the MDI model
#' @description The recommended entry point: runs several chains of Multiple
#' Dataset Integration (MDI), reports progress as they finish and checks
#' convergence before returning. It wraps \code{\link{runMCMCChains}} (same
#' arguments, same order) and adds
#' \code{\link{assessConvergence}} on the result.
#'
#' Several chains are the default workflow because a single chain cannot show
#' that it has reached the posterior: a mixture model chain can look stable and
#' still be stuck in a poor mode (Vehtari et al., 2021, recommend at least four
#' chains). The convergence check uses the rank-normalised split-\eqn{\hat{R}}
#' and effective sample sizes of Vehtari et al. (2021), computed on quantities
#' that do not depend on how the clusters are labelled.
#'
#' The sampler is compiled C++ and does not report progress within a chain, so
#' messages appear when each chain starts and finishes.
#' @inheritParams runMCMCChains
#' @param burn Number of iterations treated as burn in when assessing
#' convergence. Defaults to half of \code{R}. It only affects the diagnostics
#' attached to the result, not the chains themselves; choose the burn in for
#' point estimates in \code{\link{processMCMCChains}}.
#' @param save_pointwise,phi_update,n_cores,betas,swap_scheme,swap_every,split_merge See
#' \code{\link{runMCMCChains}}.
#' @param verbose Logical. Report progress and the convergence verdict as
#' messages. Defaults to \code{TRUE} unless \code{options(mdir.quiet = TRUE)}
#' is set.
#' @return An object of class \code{mdir_fit_list} (see
#' \code{\link{runMCMCChains}}) with the diagnostics from
#' \code{\link{assessConvergence}} attached as \code{attr(., "convergence")}.
#' Printing it gives a short report; \code{summary()} gives per-chain detail.
#' If the diagnostics cannot be computed (for example, too few saved
#' iterations) a warning is given and the attribute is absent.
#'
#' Pass the result to \code{\link{processMCMCChains}} or
#' \code{\link{predictFromMultipleChains}} for point estimates.
#' @references Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and
#' Burkner, P.-C. (2021). Rank-normalization, folding, and localization: an
#' improved Rhat for assessing convergence of MCMC. \emph{Bayesian Analysis},
#' 16(2), 667-718.
#' @seealso \code{\link{summary.mdir_fit_list}}, \code{\link{assessConvergence}},
#' \code{\link{predictFromMultipleChains}}
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' X <- lapply(1:2, function(v) {
#'   m <- matrix(rnorm(60 * 2, rep(c(0, 3), each = 30)), 60, 2)
#'   rownames(m) <- 1:60
#'   m
#' })
#'
#' # Far too few iterations for real use
#' fit <- fitMDI(X, n_chains = 3, R = 1000, thin = 5, types = c("MVN", "MVN"), K = c(4, 4))
#' fit
#' summary(fit)
#' }
fitMDI <- function(X,
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
                   burn = NULL,
                   verbose = !isTRUE(getOption("mdir.quiet")),
                   save_pointwise = FALSE,
                   phi_update = c("slice", "gibbs"),
                   n_cores = NULL,
                   betas = 1,
                   swap_scheme = c("deo", "seo"),
                   swap_every = 1L,
                   split_merge = 0L) {
  phi_update <- match.arg(phi_update)
  swap_scheme <- match.arg(swap_scheme)
  if (verbose && is.list(X) && length(X) > 0 && is.matrix(X[[1]])) {
    message(sprintf(
      "Fitting MDI to %d items in %s (%s): %s of %d iterations, thin = %d (%d samples saved per chain).",
      nrow(X[[1]]), .mdirPlural(length(X), "view"), paste(types, collapse = ", "),
      .mdirPlural(n_chains, "chain"), R, thin, floor(R / thin) + 1
    ))
  }

  chains <- runMCMCChains(
    X, n_chains, R, thin, types,
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
    verbose = verbose,
    save_pointwise = save_pointwise,
    phi_update = phi_update,
    n_cores = n_cores,
    betas = betas,
    swap_scheme = swap_scheme,
    swap_every = swap_every,
    split_merge = split_merge
  )

  convergence <- tryCatch(
    assessConvergence(chains, burn = burn),
    error = function(e) {
      warning(
        "Convergence diagnostics could not be computed: ", conditionMessage(e),
        call. = FALSE
      )
      NULL
    }
  )
  if (!is.null(convergence)) {
    attr(chains, "convergence") <- convergence
    if (verbose) {
      message(.mdirWrap(format(convergence)))
    }
  }

  chains
}
