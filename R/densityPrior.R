#' @title Density-level prior options
#' @description Options for the hierarchical (partially pooled) priors of the
#' component densities, for use with \code{\link{callMDI}},
#' \code{\link{runMCMCChains}} and \code{\link{simulatePriorPredictive}}.
#'
#' \strong{Pooling of the variance scale (\code{"MVN"}, \code{"TAGM"},
#' \code{"G"}).} The component covariances (variances for \code{"G"}) have an
#' inverse-Wishart (inverse-gamma) prior whose scale was fixed from the data. The
#' scale is now diagonal, \eqn{\mathrm{diag}(s)}, with a per-measurement
#' hyperprior \eqn{s_p \sim \mathrm{Gamma}(a, a / c_p)}. Its mean \eqn{c_p} is the
#' former data-driven value, the average marginal variance divided by
#' \eqn{K^{2/P}}, and \eqn{a =} \code{scale_pool_shape} controls how tightly
#' \eqn{s_p} is held to it (the coefficient of variation is \eqn{1/\sqrt{a}}; the
#' prior on the scale of each component's spread is informed by the other
#' components, so small or empty components borrow strength). Richardson and Green
#' (1997) also give the corresponding scale hyperparameter a prior in their
#' univariate mixture. The update conditions on the occupied components only,
#' which is the exact conditional with the empty components integrated out.
#' \code{scale_pool_shape = 0} restores the fixed scale.
#'
#' \strong{Gaussian process views (\code{"GP"}, \code{"TAGPM"}).} Component \eqn{k}
#' has mean function \eqn{\mu_k = \xi + f_k} with \eqn{\xi} the column means and
#' \eqn{f_k} a Gaussian process with squared exponential kernel
#' \eqn{a_k \exp(-\Delta^2 / (2 \lambda_k^2))} over measurements one unit
#' apart, plus measurement noise variance \eqn{\sigma_k^2}.
#' \emph{Length scale.} Fitting the length scale freely tends to drive it below
#' what the data can resolve; a length scale shorter than the spacing makes the
#' kernel diagonal, so the mean function is white noise, cannot be told apart
#' from the measurement noise, and fits noise (overfitting). Two safeguards are
#' used. \eqn{\lambda_k} cannot fall below \code{gp_min_length} (default one grid
#' spacing), and it has an inverse-gamma prior calibrated so that 1\% of its mass
#' lies below \code{gp_min_length} and 1\% above the extent of the grid,
#' \eqn{P - 1} (the recipe of Betancourt, 2017, "Robust Gaussian Process
#' Modeling"). \emph{Amplitude and noise.} \eqn{\log a_k} and
#' \eqn{\log \sigma_k^2} each have a normal population shared by the components,
#' \eqn{N(m, s^2)}, with \eqn{m \sim N(\log v, } \code{gp_center_sd}\eqn{^2)}
#' centred on the average data variance \eqn{v} (so the data need not be
#' standardised) and \eqn{s \sim} half-normal(\code{gp_pool_sd_scale})
#' (Gelman, 2006). With \code{pool_gp = FALSE} the population is fixed at
#' \eqn{N(\log v, 1)}.
#' @param scale_pool_shape Shape \eqn{a} of the hyperprior on the variance
#' scale; 0 disables pooling.
#' @param pool_gp Logical. Pool the GP amplitude and noise across components.
#' @param gp_min_length Smallest permitted GP length scale, in units of the
#' spacing between adjacent measurements.
#' @param gp_center_sd Prior sd of the centre of the log amplitude and log noise
#' populations around the log average data variance.
#' @param gp_pool_sd_scale Scale of the half-normal prior on the sd of those
#' populations.
#' @return A numeric vector of class \code{mdir_density_prior}.
#' @references Richardson, S. and Green, P. J. (1997). On Bayesian analysis of
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
                         gp_pool_sd_scale = 1) {
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
  structure(
    c(scale_pool_shape = scale_pool_shape, pool_gp = as.numeric(pool_gp),
      gp_min_length = gp_min_length, gp_center_sd = gp_center_sd,
      gp_pool_sd_scale = gp_pool_sd_scale),
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
  cat(sprintf("  GP amplitude / noise %s; min length scale %g grid units; centre sd %g; population sd scale %g\n",
              if (x[2] > 0) "pooled" else "not pooled", x[3], x[4], x[5]))
  invisible(x)
}
