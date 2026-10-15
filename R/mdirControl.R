#' @title Sampler-tuning control for an mdir fit
#' @description Bundles the Metropolis-Hastings proposal windows and the
#' auto-tuning schedule, the arguments that control how hard the sampler works to
#' explore the posterior and not what is being fitted, into one object, the way
#' \code{stats::glm(control = glm.control(...))} does and as
#' \code{batchmix::batchmixControl()} does for batchmix. The windows are starting
#' values for the adaptation (\code{auto_tune = TRUE}, the default), so they need
#' not be right.
#'
#' Choices that change what is fitted, such as \code{prior}, \code{density_prior},
#' \code{K_max} and \code{alpha}, stay as ordinary arguments of the fitting
#' functions. So do the choices of the sampler's moves (\code{phi_update},
#' \code{split_merge}, \code{joint_allocation} and the tempering ladder).
#'
#' Only the Gaussian process hyperparameters have a user-set proposal window. Two
#' other Metropolis-Hastings proposals have fixed scales that are not adapted:
#' the random walk on the log of the mass parameter (standard deviation 0.3) and
#' the one on the standard deviation of the population of Gaussian process
#' amplitudes and noises (0.5).
#' @param proposal_windows \code{NULL}, or a list with an entry for each view. For a
#' view modelled with a Gaussian process (\code{"GP"} or \code{"TAGPM"}), the entry is
#' \code{NULL} or a vector of three standard deviations of the log-scale random-walk
#' proposals for the amplitude, the length scale and the noise (the default is
#' \code{c(0.125, 0.15, 0.10)}); the amplitude proposal is a mixture of three random
#' walks, with standard deviations 0.4, 1 and 2.5 times its window. The entries of
#' other views are not used.
#' @param auto_tune Logical. Adapt the windows of the Gaussian process
#' hyperparameters while the sampler burns in, towards the acceptance rate of about
#' 0.44 of a one-dimensional random-walk Metropolis update (Gelman, Roberts and
#' Gilks, 1996), by a Robbins-Monro recursion with decaying steps (Garthwaite, Fan
#' and Sisson, 2016), and then keep them fixed. Without a Gaussian process view
#' there is nothing to adapt, and no random numbers are used, so the chain is the
#' same whatever its value. The windows at the end are returned in
#' \code{final_proposal_windows}.
#' @param n_burn Number of sweeps over which the windows are adapted. The default
#' (\code{NULL}) is half of \code{n_iter}.
#' @return An object of class \code{mdir_control}, a list of the arguments.
#' @references Gelman, A., Roberts, G. O. and Gilks, W. R. (1996). Efficient
#' Metropolis jumping rules. In \emph{Bayesian Statistics 5}, 599-607. Oxford
#' University Press.
#'
#' Garthwaite, P. H., Fan, Y. and Sisson, S. A. (2016). Adaptive optimal scaling of
#' Metropolis-Hastings algorithms using the Robbins-Monro process.
#' \emph{Communications in Statistics - Theory and Methods}, 45(17), 5098-5111.
#' @seealso \code{\link{callMDI}}, \code{\link{fitMDI}}, \code{\link{densityPrior}}
#' @export
#' @examples
#' mdirControl()
#'
#' # Start the amplitude, length and noise windows of the first view lower, and adapt
#' # them over the first 500 sweeps only
#' mdirControl(proposal_windows = list(c(0.05, 0.05, 0.05)), n_burn = 500)
#'
#' # Fixed windows
#' mdirControl(auto_tune = FALSE)
mdirControl <- function(proposal_windows = NULL, auto_tune = TRUE, n_burn = NULL) {
  if (!is.null(proposal_windows)) {
    ok <- is.list(proposal_windows) && all(vapply(proposal_windows, function(w) {
      is.null(w) || (is.numeric(w) && length(w) >= 1 && all(is.finite(w)) && all(w >= 0))
    }, logical(1)))
    if (!ok) {
      stop("mdirControl(): `proposal_windows` must be NULL or a list with, for each view, NULL ",
        "or a vector of non-negative numbers.",
        call. = FALSE
      )
    }
  }
  if (!is.logical(auto_tune) || length(auto_tune) != 1L || is.na(auto_tune)) {
    stop("mdirControl(): `auto_tune` must be a single logical value (TRUE or FALSE).", call. = FALSE)
  }
  if (!is.null(n_burn) && (!is.numeric(n_burn) || length(n_burn) != 1L || is.na(n_burn) || n_burn < 0)) {
    stop("mdirControl(): `n_burn` must be NULL or a single non-negative number.", call. = FALSE)
  }
  structure(
    list(proposal_windows = proposal_windows, auto_tune = auto_tune, n_burn = n_burn),
    class = "mdir_control"
  )
}

#' @export
print.mdir_control <- function(x, ...) {
  cat("<mdir control>\n")
  cat(sprintf(
    "  auto_tune = %s, n_burn = %s\n",
    x$auto_tune, if (is.null(x$n_burn)) "NULL (half of n_iter)" else x$n_burn
  ))
  if (is.null(x$proposal_windows)) {
    cat("  proposal_windows: defaults (amplitude, length, noise = 0.125, 0.15, 0.10)\n")
  } else {
    cat("  proposal_windows (starting values if auto_tune = TRUE, fixed otherwise):\n")
    for (v in seq_along(x$proposal_windows)) {
      w <- x$proposal_windows[[v]]
      cat(sprintf("    view %d: %s\n", v, if (is.null(w)) "default" else paste(format(w), collapse = ", ")))
    }
  }
  invisible(x)
}

# Internal: the `control` of a fitting function, with the deprecated individual
# argument `proposal_windows` merged in as batchmix does: an explicit `control` wins
# over it (with a warning), and on its own it is moved into the default control (also
# with a warning).
# control_missing: missing(control) in the caller; control_value: its `control`;
# deprecated: a list of the deprecated arguments the caller was given, possibly empty;
# fn_name: the caller's name, for the message.
.mdirResolveControl <- function(control_missing, control_value, deprecated, fn_name) {
  if (control_missing) {
    control <- mdirControl()
  } else {
    if (!inherits(control_value, "mdir_control")) {
      stop("`control` must be created with mdirControl().", call. = FALSE)
    }
    control <- control_value
  }
  if (length(deprecated) > 0) {
    fields <- paste(sprintf("`%s`", names(deprecated)), collapse = ", ")
    if (control_missing) {
      warning(
        sprintf(
          "The %s argument of %s() is deprecated; use `control = mdirControl(%s = ...)` instead.",
          fields, fn_name, names(deprecated)[1]
        ),
        call. = FALSE
      )
      control[names(deprecated)] <- deprecated
    } else {
      warning(
        sprintf(
          "Both `control` and %s were given to %s(); using `control` and ignoring %s.",
          fields, fn_name, fields
        ),
        call. = FALSE
      )
    }
  }
  control
}
