# Internal: standardise a data argument to a list of matrices with row names
.asViewList <- function(X) {
  if (is.matrix(X) || is.data.frame(X)) {
    X <- list(as.matrix(X))
  }
  if (!is.list(X) || !all(vapply(X, is.matrix, logical(1)))) {
    stop("`X` must be a matrix or a list of matrices.")
  }
  X
}

# Internal: the C++ codes for a set of types
.typeCodes <- function(types) {
  list(density = translateTypes(types), outlier = setupOutlierComponents(types))
}

#' @title Simulate datasets from the prior predictive distribution
#' @description Draws the MDI parameters (concentrations, weights,
#' association parameters), the component parameters and the item allocations
#' from the model's priors and then simulates every view, for a prior
#' predictive check: does the model, before seeing the data, generate datasets
#' that look like the data you expect to see? (Gelman et al., 2013, Section 6.1;
#' Gelman et al., 2020, Section 2.) If the simulated datasets or cluster
#' structures are implausible, the prior needs reconsidering
#' (see \code{\link{mdiPrior}}).
#'
#' The simulation uses the same C++ prior code as the sampler, so it cannot drift
#' from the model that is fitted. The data-driven prior hyperparameters (means
#' and scales) are computed from \code{X}, so this checks whether that
#' empirical-Bayes construction produces sensible data at the chosen \code{K},
#' not a data-independent prior. \code{X} is used only for these hyperparameters
#' and its dimensions.
#' @param X A matrix or list of matrices (items in rows) as passed to
#' \code{\link{callMDI}}. \code{NA} entries are ignored when setting the
#' hyperparameters.
#' @param types Density type of each view, as for \code{\link{callMDI}}.
#' @param K Number of components in each view. Defaults to \code{min(10, N/2)}
#' in every view.
#' @param n_datasets Number of datasets to simulate.
#' @param prior MDI-level priors, see \code{\link{mdiPrior}}.
#' @param mimic_missingness If \code{TRUE} the entries that are \code{NA} in
#' \code{X} are set to \code{NA} in every simulated dataset so that summaries
#' compare like with like (this assumes missingness is independent of the
#' values, as the model does).
#' @return A list of class \code{mdir_predictive} with \code{replicates} (a list
#' with one entry per simulated dataset, each a list of matrices, one per view),
#' \code{parameters} (the simulated concentrations, phis, weights, labels and
#' outlier indicators of each dataset) and \code{source = "prior"}.
#' @references Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari,
#' A. and Rubin, D. B. (2013). \emph{Bayesian Data Analysis}, 3rd edition.
#' CRC Press.
#'
#' Gelman, A., Vehtari, A., Simpson, D., et al. (2020). Bayesian workflow.
#' arXiv:2011.01808.
#' @seealso \code{\link{simulatePosteriorPredictive}}, \code{\link{plotPredictiveCheck}}
#' @export
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(80 * 2, rep(c(0, 3), each = 40)), 80, 2)
#' rownames(X) <- seq_len(nrow(X))
#' prior_sims <- simulatePriorPredictive(X, "MVN", K = 3, n_datasets = 5)
#' plotPredictiveCheck(X, prior_sims, column = 1)
simulatePriorPredictive <- function(X,
                                    types,
                                    K = NULL,
                                    n_datasets = 20,
                                    prior = mdiPrior(),
                                    mimic_missingness = TRUE) {
  X <- .asViewList(X)
  checkDataCorrectInput(X, types)
  V <- length(X)
  N <- nrow(X[[1]])
  if (is.null(K)) {
    K <- rep(min(10, floor(N / 2)), V)
  }
  if (length(K) == 1) {
    K <- rep(K, V)
  }
  codes <- .typeCodes(types)
  raw <- simulatePriorPredictiveCpp(X, as.integer(K), codes$density, codes$outlier,
                                    as.integer(n_datasets), as.numeric(prior))

  replicates <- lapply(raw, function(r) {
    lapply(seq_len(V), function(v) {
      m <- r$data[[v]]
      if (mimic_missingness) {
        m[is.na(X[[v]])] <- NA_real_
      }
      dimnames(m) <- dimnames(X[[v]])
      m
    })
  })
  parameters <- lapply(raw, function(r) r[setdiff(names(r), "data")])
  structure(list(replicates = replicates, parameters = parameters, source = "prior"),
            class = "mdir_predictive")
}

#' @title Simulate replicates from the posterior predictive distribution
#' @description For a sample of saved MCMC iterations, simulates a replicate of
#' every item from that iteration's component parameters and allocations (or from
#' the outlier distribution for items sampled as outliers), for a posterior
#' predictive check: does the fitted model generate data that look like the data
#' we observed? (Gelman et al., 2013, Section 6.3; Gabry et al., 2019.) The
#' replicates condition on the sampled allocations, so they check the
#' within-cluster distributions (shape, scale, correlation) but not the
#' clustering itself. Cells that are missing in the observed data are simulated
#' like any other and, by default, set to \code{NA} so that summaries compare
#' like with like.
#'
#' Posterior predictive p-values are conservative because the data are used
#' twice (Bayarri and Berger, 2000): treat values near 0 or 1 as evidence of
#' misfit, and do not read values in the middle as proof of fit.
#' @param mcmc_output A single chain from \code{\link{callMDI}} or a list of
#' chains from \code{\link{runMCMCChains}}. Component parameters must have been
#' saved (\code{save_parameters = TRUE}, the default).
#' @param X The data that were modelled (matrix or list of matrices).
#' @param n_draws Number of iterations to replicate from, sampled without
#' replacement from those retained after the burn in (pooled across chains).
#' @param burn Number of iterations to discard as burn in.
#' @param mimic_missingness Set the entries that are \code{NA} in \code{X} to
#' \code{NA} in the replicates.
#' @return A list of class \code{mdir_predictive} with \code{replicates} (a list
#' with one entry per draw, each a list of matrices, one per view),
#' \code{draws} (the chain and saved-sample index of each replicate) and
#' \code{source = "posterior"}.
#' @references Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari,
#' A. and Rubin, D. B. (2013). \emph{Bayesian Data Analysis}, 3rd edition.
#' CRC Press.
#'
#' Gabry, J., Simpson, D., Vehtari, A., Betancourt, M. and Gelman, A. (2019).
#' Visualization in Bayesian workflow. \emph{Journal of the Royal Statistical
#' Society A}, 182(2), 389-402.
#'
#' Bayarri, M. J. and Berger, J. O. (2000). P values for composite null models.
#' \emph{Journal of the American Statistical Association}, 95(452), 1127-1142.
#' @seealso \code{\link{simulatePriorPredictive}}, \code{\link{plotPredictiveCheck}}
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' X <- matrix(rnorm(80 * 2, rep(c(0, 3), each = 40)), 80, 2)
#' rownames(X) <- seq_len(nrow(X))
#' fit <- callMDI(list(X), R = 500, thin = 5, types = "MVN", K = 4)
#' post_sims <- simulatePosteriorPredictive(fit, X, n_draws = 20, burn = 250)
#' plotPredictiveCheck(X, post_sims, column = 1)
#' }
simulatePosteriorPredictive <- function(mcmc_output,
                                        X,
                                        n_draws = 50,
                                        burn = 0,
                                        mimic_missingness = TRUE) {
  chains <- if (!is.null(mcmc_output$allocations)) list(mcmc_output) else mcmc_output
  X <- .asViewList(X)
  first <- chains[[1]]
  V <- first$V
  if (is.null(V)) V <- length(X)
  if (is.null(first$parameters) || any(vapply(first$parameters, function(p) is.null(p) || length(p) == 0, logical(1)))) {
    stop("Component parameters were not saved. Re-run with `save_parameters = TRUE`.")
  }
  if (length(X) != V || nrow(X[[1]]) != first$N) {
    stop("`X` does not match the data the chains were fitted to.")
  }

  first_kept <- floor(burn / first$thin) + 2
  pool <- do.call(rbind, lapply(seq_along(chains), function(i) {
    n_saved <- dim(chains[[i]]$allocations)[1]
    if (first_kept > n_saved) return(NULL)
    cbind(chain = i, sample = seq(first_kept, n_saved))
  }))
  if (is.null(pool) || nrow(pool) == 0) {
    stop("No saved iterations remain after applying `burn`.")
  }
  pick <- pool[sample.int(nrow(pool), min(n_draws, nrow(pool))), , drop = FALSE]

  types <- first$types
  codes <- .typeCodes(types)
  K <- first$K
  N <- first$N

  # Assemble the selected draws
  n_pick <- nrow(pick)
  allocations <- array(0, dim = c(n_pick, N, V))
  outliers <- array(0, dim = c(n_pick, N, V))
  parameters <- lapply(seq_len(V), function(v) {
    matrix(0, n_pick, ncol(first$parameters[[v]]))
  })
  for (i in seq_len(n_pick)) {
    ch <- chains[[pick[i, "chain"]]]
    s <- pick[i, "sample"]
    for (v in seq_len(V)) {
      allocations[i, , v] <- ch$allocations[s, , v]
      outliers[i, , v] <- ch$outliers[s, , v]
      parameters[[v]][i, ] <- ch$parameters[[v]][s, ]
    }
  }

  prior <- if (is.null(first$prior)) mdiPrior() else first$prior
  cubes <- simulatePosteriorPredictiveCpp(X, as.integer(K), codes$density, codes$outlier,
                                          parameters, allocations, outliers, as.numeric(prior))

  replicates <- lapply(seq_len(n_pick), function(i) {
    lapply(seq_len(V), function(v) {
      m <- matrix(cubes[[v]][i, , ], nrow = N)
      if (mimic_missingness) {
        m[is.na(X[[v]])] <- NA_real_
      }
      dimnames(m) <- dimnames(X[[v]])
      m
    })
  })
  structure(list(replicates = replicates, draws = pick, source = "posterior"),
            class = "mdir_predictive")
}

#' @title Predictive check statistics and p-values
#' @description Applies a scalar test statistic to the observed data and to
#' every replicate from \code{\link{simulatePriorPredictive}} or
#' \code{\link{simulatePosteriorPredictive}}, and returns the statistic values
#' with the tail-area probability \eqn{P(T(y^{rep}) \ge T(y))} (Gelman et al.,
#' 2013, Section 6.3). A value close to 0 or 1 flags a feature of the data
#' that the model does not reproduce. For a two-sided reading use
#' \code{min(p, 1 - p) * 2}.
#' @param X Observed data (matrix or list of matrices).
#' @param predictive Output of \code{\link{simulatePriorPredictive}} or
#' \code{\link{simulatePosteriorPredictive}}.
#' @param statistic A function of a numeric matrix (the observed items by
#' measurements of one view, possibly with \code{NA}s) returning one number.
#' @param view Which view to check.
#' @return A list with \code{observed}, \code{replicated} and \code{p_value}.
#' @export
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(60 * 2, rep(c(0, 3), each = 30)), 60, 2)
#' rownames(X) <- seq_len(60)
#' sims <- simulatePriorPredictive(X, "MVN", K = 3, n_datasets = 20)
#' predictiveCheck(X, sims, function(m) cor(m[, 1], m[, 2], use = "complete.obs"))
predictiveCheck <- function(X, predictive, statistic, view = 1) {
  X <- .asViewList(X)
  observed <- statistic(X[[view]])
  replicated <- vapply(predictive$replicates, function(r) statistic(r[[view]]), numeric(1))
  list(observed = observed, replicated = replicated, p_value = mean(replicated >= observed))
}

#' @title Plot a prior or posterior predictive check
#' @description Compares the observed data with replicates from
#' \code{\link{simulatePriorPredictive}} or
#' \code{\link{simulatePosteriorPredictive}} (Gabry et al., 2019). Two styles:
#' \itemize{
#'   \item \code{"density"}: the density of one measurement in the observed data
#'   (black) over the same density in each replicate (blue);
#'   \item \code{"statistic"}: the histogram of a scalar statistic across
#'   replicates with the observed value marked, and the tail-area probability
#'   \eqn{P(T(y^{rep}) \ge T(y))} in the subtitle.
#' }
#' \code{NA} entries are excluded from the observed data, and, if the
#' replicates were made with \code{mimic_missingness = TRUE}, from the
#' replicates too.
#' @param X Observed data (matrix or list of matrices).
#' @param predictive Output of \code{\link{simulatePriorPredictive}} or
#' \code{\link{simulatePosteriorPredictive}}.
#' @param style \code{"density"} or \code{"statistic"}.
#' @param view,column Which view and measurement to check.
#' @param statistic For \code{style = "statistic"}, a function of a numeric
#' vector (the column, without \code{NA}s) returning one number.
#' @param statistic_name Label for the statistic.
#' @param max_replicates Cap on the number of replicate densities drawn.
#' @return A \code{ggplot} object.
#' @references Gabry, J., Simpson, D., Vehtari, A., Betancourt, M. and Gelman,
#' A. (2019). Visualization in Bayesian workflow. \emph{Journal of the Royal
#' Statistical Society A}, 182(2), 389-402.
#' @export
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(60 * 2, rep(c(0, 3), each = 30)), 60, 2)
#' rownames(X) <- seq_len(60)
#' sims <- simulatePriorPredictive(X, "MVN", K = 3, n_datasets = 20)
#' plotPredictiveCheck(X, sims, style = "density", column = 1)
#' plotPredictiveCheck(X, sims, style = "statistic", column = 1,
#'   statistic = stats::sd, statistic_name = "SD")
plotPredictiveCheck <- function(X,
                                predictive,
                                style = c("density", "statistic"),
                                view = 1,
                                column = 1,
                                statistic = mean,
                                statistic_name = "statistic",
                                max_replicates = 50) {
  style <- match.arg(style)
  X <- .asViewList(X)
  obs <- X[[view]][, column]
  obs <- obs[!is.na(obs)]
  reps <- lapply(predictive$replicates, function(r) {
    v <- r[[view]][, column]
    v[!is.na(v)]
  })
  title <- paste0(if (identical(predictive$source, "prior")) "Prior" else "Posterior",
                  " predictive check: view ", view, ", column ", column)

  if (style == "density") {
    if (all(obs %in% c(0, 1))) {
      warning("The column looks binary; use style = 'statistic' with statistic = mean.")
    }
    reps <- reps[seq_len(min(length(reps), max_replicates))]
    rep_df <- do.call(rbind, lapply(seq_along(reps), function(i) data.frame(value = reps[[i]], draw = i)))
    return(
      ggplot2::ggplot() +
        ggplot2::geom_density(data = rep_df, ggplot2::aes(x = .data$value, group = .data$draw),
                              colour = "steelblue", linewidth = 0.3) +
        ggplot2::geom_density(data = data.frame(value = obs), ggplot2::aes(x = .data$value),
                              colour = "black", linewidth = 1) +
        ggplot2::labs(title = title,
                      subtitle = paste0(length(reps), " replicates (blue) and observed data (black)"),
                      x = paste0("Column ", column), y = "Density") +
        ggplot2::theme_minimal()
    )
  }

  obs_stat <- statistic(obs)
  rep_stats <- vapply(reps, statistic, numeric(1))
  p_value <- mean(rep_stats >= obs_stat)
  ggplot2::ggplot(data.frame(value = rep_stats), ggplot2::aes(x = .data$value)) +
    ggplot2::geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
    ggplot2::geom_vline(xintercept = obs_stat, colour = "black", linewidth = 1) +
    ggplot2::labs(title = title,
                  subtitle = paste0("Observed ", statistic_name, " = ", signif(obs_stat, 4),
                                    "; P(replicate >= observed) = ", signif(p_value, 3)),
                  x = statistic_name, y = "Replicates") +
    ggplot2::theme_minimal()
}
