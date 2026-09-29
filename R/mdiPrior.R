#' @title Prior hyperparameters of the MDI model
#' @description Creates the MDI-level prior hyperparameters used by
#' \code{\link{callMDI}}, \code{\link{runMCMCChains}} and
#' \code{\link{simulatePriorPredictive}}.
#'
#' \strong{The MDI-level priors.} For view \eqn{l} with \eqn{K_l} components the
#' unnormalised component weights are
#' \eqn{w_{lk} \sim \mathrm{Gamma}(a_l / K_l, b_w)} with concentration
#' \eqn{a_l \sim \mathrm{Gamma}(}\code{mass_shape}\eqn{,}\code{mass_rate}\eqn{)}, and
#' the dataset association parameters are
#' \eqn{\phi_{lm} \sim \mathrm{Gamma}(}\code{phi_shape}\eqn{,}\code{phi_rate}\eqn{)}.
#' Normalised, the weights of one view are a symmetric Dirichlet with total
#' concentration \eqn{a_l}. Items are allocated jointly across views with
#' probability proportional to
#' \eqn{\prod_l w_{l c_l} \prod_{l<m} (1 + \phi_{lm} 1[c_l = c_m])} (Kirk et
#' al., 2012).
#'
#' \strong{Choosing them.} (i) With \eqn{K} components and concentration
#' \eqn{a}, each component has Dirichlet parameter \eqn{e_0 = a / K}. In an
#' overfitted mixture (\eqn{K} larger than the number of clusters) the extra
#' components empty out if \eqn{e_0 < d / 2}, where \eqn{d} is the number of free
#' parameters of one component, and instead duplicate occupied components if
#' \eqn{e_0 > d / 2} (Rousseau and Mengersen, 2011). The defaults give a prior
#' mean of \eqn{a = 20}, so with a large \code{K} the extra components empty
#' out, but with \code{K} of about 5 or less the prior is not sparse. Adjust
#' \code{mass_shape} and \code{mass_rate} to your \code{K}. (ii) \eqn{\phi}
#' multiplies the probability that an item takes the same component in two
#' views by \eqn{1 + \phi}; the default prior has median around 8 and mean
#' 10, which expresses an expectation of strong agreement between views. If you
#' expect weaker sharing, lower the prior mean. (iii) Use
#' \code{\link{simulatePriorPredictive}} to look at datasets and cluster
#' agreement that these choices imply before fitting.
#'
#' \strong{Priors that are not set here.} Each density has conjugate priors
#' whose location and scale are set from the data (empirical Bayes), as is
#' standard for finite mixtures (Fraley and Raftery, 2007): component means
#' \eqn{N(\bar{x}, \Sigma_k / 0.01)}, covariances inverse-Wishart with \eqn{P + 2}
#' degrees of freedom and scale \eqn{\mathrm{tr}(S)/(P K^{2/P}) I} for
#' \code{"MVN"}/\code{"TAGM"} (independent normal-inverse-gamma with 3 degrees
#' of freedom for \code{"G"}), and Dirichlet priors with total concentration 1
#' equal to the observed category frequencies for \code{"C"}. Because the data
#' set these hyperparameters, the posterior uses the data twice (weakly, through
#' two location/scale summaries). \code{"GP"}/\code{"TAGPM"} views use
#' log-normal(0, 1) priors on the amplitude, length scale and noise variance,
#' which assume the data are on a scale of order one: \strong{standardise GP
#' views}. Outliers (\code{"TAGM"}, \code{"TAGPM"}) have weight
#' \eqn{\mathrm{Beta}(2, 10)} and a multivariate t with 4 degrees of freedom,
#' the data mean and half of the data covariance (Crook et al., 2018).
#' @param mass_shape,mass_rate Shape and rate of the Gamma prior on each view's
#' concentration \eqn{a_l}.
#' @param weight_rate Rate \eqn{b_w} of the Gamma prior on the component weights.
#' It sets the scale of the unnormalised weights only, and does not affect the
#' normalised weights of a view.
#' @param phi_shape,phi_rate Shape and rate of the Gamma prior on each
#' \eqn{\phi_{lm}}.
#' @return A numeric vector of class \code{mdir_prior}.
#' @references Kirk, P., Griffin, J. E., Savage, R. S., Ghahramani, Z. and
#' Wild, D. L. (2012). Bayesian correlated clustering to integrate multiple
#' datasets. \emph{Bioinformatics}, 28(24), 3290-3297.
#'
#' Rousseau, J. and Mengersen, K. (2011). Asymptotic behaviour of the posterior
#' distribution in overfitted mixture models. \emph{Journal of the Royal
#' Statistical Society B}, 73(5), 689-710.
#'
#' Fraley, C. and Raftery, A. E. (2007). Bayesian regularization for normal
#' mixture estimation and model-based clustering. \emph{Journal of
#' Classification}, 24, 155-181.
#'
#' Crook, O. M., Mulvey, C. M., Kirk, P. D. W., Lilley, K. S. and Gatto, L.
#' (2018). A Bayesian mixture modelling approach for spatial proteomics.
#' \emph{PLoS Computational Biology}, 14(11), e1006516.
#' @export
#' @examples
#' mdiPrior()
#'
#' # A sparser prior on the concentration and weaker view association
#' mdiPrior(mass_shape = 1, mass_rate = 0.5, phi_shape = 2, phi_rate = 0.5)
mdiPrior <- function(mass_shape = 2,
                     mass_rate = 0.1,
                     weight_rate = 2,
                     phi_shape = 2,
                     phi_rate = 0.2) {
  prior <- c(
    mass_shape = mass_shape, mass_rate = mass_rate, weight_rate = weight_rate,
    phi_shape = phi_shape, phi_rate = phi_rate
  )
  if (!is.numeric(prior) || length(prior) != 5 || any(!is.finite(prior)) || any(prior <= 0)) {
    stop("All prior hyperparameters must be single positive finite numbers.")
  }
  structure(prior, class = "mdir_prior")
}

#' @export
print.mdir_prior <- function(x, ...) {
  cat("MDI prior hyperparameters\n")
  cat(sprintf("  concentration a_l ~ Gamma(%g, %g)  [mean %g, sd %.3g]\n", x[1], x[2], x[1] / x[2], sqrt(x[1]) / x[2]))
  cat(sprintf("  weights w_lk ~ Gamma(a_l / K_l, %g)\n", x[3]))
  cat(sprintf("  phi_lm ~ Gamma(%g, %g)  [mean %g, sd %.3g]\n", x[4], x[5], x[4] / x[5], sqrt(x[4]) / x[5]))
  invisible(x)
}
