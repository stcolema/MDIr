#' @title Density-level prior options
#' @description Options for the hierarchical (partially pooled) priors of the
#' component densities, for use with \code{\link{callMDI}},
#' \code{\link{runMCMCChains}} and \code{\link{simulatePriorPredictive}}.
#'
#' \strong{Pooling of the variance scale (\code{"MVN"}, \code{"TAGM"},
#' \code{"G"}).} The component covariances (variances for \code{"G"}) have an
#' inverse-Wishart (inverse-gamma) prior whose scale is diagonal, \eqn{\mathrm{diag}(s)}, with a per-measurement
#' hyperprior \eqn{s_p \sim \mathrm{Gamma}(a, a / c_p)}. Its mean \eqn{c_p} is set from the data: the average marginal variance divided by
#' \eqn{K^{2/P}}, and \eqn{a =} \code{scale_pool_shape} controls how tightly
#' \eqn{s_p} is held to it (the coefficient of variation is \eqn{1/\sqrt{a}}; the
#' prior on the scale of each component's spread is informed by the other
#' components, so small or empty components borrow strength). Richardson and Green
#' (1997) also give the corresponding scale hyperparameter a prior in their
#' univariate mixture. The update conditions on the occupied components only,
#' which is the exact conditional with the empty components integrated out.
#' \code{scale_pool_shape = 0} fixes the scale at \eqn{c_p}.
#'
#' \strong{Precision of the component-mean prior (\code{"MVN"}, \code{"TAGM"},
#' \code{"G"}).} The component means have the conjugate prior
#' \eqn{\mu_k \mid \Sigma_k \sim N(\xi, \Sigma_k / \kappa)}, so the prior
#' standard deviation of a mean is \eqn{\sigma_k / \sqrt{\kappa}}. The default
#' \eqn{\kappa = 0.01} (a small fraction of one observation of information,
#' as in the regularising prior of Fraley and Raftery, 2007) therefore says that
#' component means lie within about ten of their own standard deviations of the
#' grand mean. That is harmless for estimating the means of known components, but
#' the marginal likelihood of a partition depends on it: for well separated
#' groups the penalty on a component far from \eqn{\xi} can outweigh the
#' evidence for splitting it, so the posterior odds that a group splits can
#' change when unrelated groups elsewhere in the data are added or removed.
#' With \code{kappa_range = c(lo, hi)} \eqn{\kappa} is shared by the
#' components and given a log-uniform hyperprior on \eqn{[lo, hi]}, and is
#' sampled from its full conditional given the occupied components,
#' \eqn{\mathrm{Gamma}(P n_{occ} / 2, \sum_k (\mu_k - \xi)^\top \Sigma_k^{-1}
#' (\mu_k - \xi) / 2)} truncated to the range. It is then included among the
#' pooled hyperparameters (the last column of \code{pooled_hyperparameters}),
#' and \code{kappa} is ignored. \code{kappa_range = NULL} (the default) keeps
#' \eqn{\kappa} fixed at \code{kappa}.
#'
#' \strong{Gaussian process views (\code{"GP"}, \code{"TAGPM"}).} Component \eqn{k}
#' has mean function \eqn{\mu_k = \xi + f_k} with \eqn{\xi} the column means and
#' \eqn{f_k} a Gaussian process with squared exponential kernel
#' \eqn{a_k \exp(-\Delta^2 / (2 \lambda_k^2))} over the positions of the
#' measurements (one unit apart unless \code{time_points} is given to the fitting
#' function), plus measurement noise variance \eqn{\sigma_k^2}.
#' \emph{Length scale.} Fitting the length scale freely tends to drive it below
#' what the data can resolve; a length scale shorter than the spacing makes the
#' kernel diagonal, so the mean function is white noise, cannot be told apart
#' from the measurement noise, and fits noise (overfitting). Two safeguards are
#' used. \eqn{\lambda_k} cannot fall below \code{gp_min_length} times the smallest
#' spacing between measurements (by default one spacing), and it has an
#' inverse-gamma prior calibrated so that 1\% of its mass lies below that floor and
#' 1\% above the extent of the measurements, the distance from the first to the last
#' (\eqn{P - 1} for the default positions; the recipe of Betancourt, 2017, "Robust
#' Gaussian Process Modeling"). With unequally spaced measurements the smallest
#' spacing is a choice, not a result: it is the largest floor that never excludes a
#' length scale that the closest pair of measurements can resolve, but a very close
#' pair makes the floor low for the rest, in which case raise \code{gp_min_length}. \emph{Amplitude and noise.} \eqn{\log a_k} and
#' \eqn{\log \sigma_k^2} each have a normal population shared by the components,
#' \eqn{N(m, s^2)}, with \eqn{m \sim N(\log v, } \code{gp_center_sd}\eqn{^2)}
#' centred on the average data variance \eqn{v} (so the data need not be
#' standardised) and \eqn{s \sim} half-normal(\code{gp_pool_sd_scale})
#' (Gelman, 2006). With \code{pool_gp = FALSE} the population is fixed at
#' \eqn{N(\log v, 1)}.
#' @param scale_pool_shape Shape \eqn{a} of the hyperprior on the variance
#' scale; 0 disables pooling.
#' @param kappa Positive number. The precision of the component-mean prior when
#' it is not sampled.
#' @param kappa_range \code{NULL} (the default) or a numeric vector
#' \code{c(lo, hi)} with \code{0 < lo < hi}: the range of a log-uniform
#' hyperprior on \eqn{\kappa}, which is then sampled. For example
#' \code{c(1e-6, 1)}.
#' @param pool_gp Logical. Pool the GP amplitude and noise across components.
#' @param gp_min_length Smallest permitted GP length scale, in multiples of the
#' smallest spacing between measurements (the spacing itself for equally spaced
#' columns).
#' @param gp_center_sd Prior sd of the centre of the log amplitude and log noise
#' populations around the log average data variance.
#' @param gp_pool_sd_scale Scale of the half-normal prior on the sd of those
#' populations.
#' @return A numeric vector of class \code{mdir_density_prior}.
#' @references Fraley, C. and Raftery, A. E. (2007). Bayesian regularization for
#' normal mixture estimation and model-based clustering. \emph{Journal of
#' Classification}, 24, 155-181.
#'
#' Richardson, S. and Green, P. J. (1997). On Bayesian analysis of
#' mixtures with an unknown number of components. \emph{Journal of the Royal
#' Statistical Society B}, 59(4), 731-792.
#'
#' Gelman, A. (2006). Prior distributions for variance parameters in
#' hierarchical models. \emph{Bayesian Analysis}, 1(3), 515-534.
#'
#' Betancourt, M. (2017). Robust Gaussian process modeling. Stan case study,
#' \url{https://betanalpha.github.io/assets/case_studies/gp_part3/part3.html}.
#' @export
#' @examples
#' densityPrior()
#'
#' # No pooling of the variance scale, and a longer minimum GP length scale
#' densityPrior(scale_pool_shape = 0, gp_min_length = 2)
densityPrior <- function(scale_pool_shape = 2,
                         pool_gp = TRUE,
                         gp_min_length = 1,
                         gp_center_sd = 2,
                         gp_pool_sd_scale = 1,
                         kappa = 0.01,
                         kappa_range = NULL) {
  ok <- function(x) is.numeric(x) && length(x) == 1 && is.finite(x)
  if (!ok(scale_pool_shape) || scale_pool_shape < 0) {
    stop("`scale_pool_shape` must be a single non-negative number.")
  }
  if (!is.logical(pool_gp) || length(pool_gp) != 1 || is.na(pool_gp)) {
    stop("`pool_gp` must be TRUE or FALSE.")
  }
  for (nm in c("gp_min_length", "gp_center_sd", "gp_pool_sd_scale")) {
    if (!ok(get(nm)) || get(nm) <= 0) {
      stop("`", nm, "` must be a single positive number.")
    }
  }
  if (!ok(kappa) || kappa <= 0) {
    stop("`kappa` must be a single positive number.")
  }
  if (is.null(kappa_range)) {
    kappa_range <- c(0, 0)
  } else if (!is.numeric(kappa_range) || length(kappa_range) != 2 ||
             !all(is.finite(kappa_range)) || kappa_range[1] <= 0 ||
             kappa_range[2] <= kappa_range[1]) {
    stop("`kappa_range` must be NULL or c(lo, hi) with 0 < lo < hi.")
  }
  structure(
    c(scale_pool_shape = scale_pool_shape, pool_gp = as.numeric(pool_gp),
      gp_min_length = gp_min_length, gp_center_sd = gp_center_sd,
      gp_pool_sd_scale = gp_pool_sd_scale, kappa = kappa,
      kappa_min = kappa_range[1], kappa_max = kappa_range[2]),
    class = "mdir_density_prior"
  )
}

#' @rdname densityPrior
#' @param x An \code{mdir_density_prior} object.
#' @param ... Unused.
#' @export
print.mdir_density_prior <- function(x, ...) {
  cat("Density-level prior options\n")
  if (x[1] > 0) {
    cat(sprintf("  variance scale pooled across components, Gamma shape %g\n", x[1]))
  } else {
    cat("  variance scale fixed (no pooling)\n")
  }
  if (length(x) >= 8 && x[8] > x[7]) {
    cat(sprintf("  component-mean precision kappa sampled, log-uniform on [%g, %g]\n", x[7], x[8]))
  } else {
    cat(sprintf("  component-mean precision kappa fixed at %g\n", if (length(x) >= 6) x[6] else 0.01))
  }
  cat(sprintf("  GP amplitude / noise %s; min length scale %g grid units; centre sd %g; population sd scale %g\n",
              if (x[2] > 0) "pooled" else "not pooled", x[3], x[4], x[5]))
  invisible(x)
}
