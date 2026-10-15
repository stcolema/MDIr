# Internal: the positions of the columns of the Gaussian process views.
#
# `time_points` is NULL (every GP view's columns lie at 0, 1, ..., P - 1), or a list
# with, for each view, NULL or a numeric vector with one position per column (a plain
# vector is accepted when there is a single view). Positions need to be finite and
# distinct, not sorted. Returns a list with one numeric vector per view, `numeric(0)`
# where the default is used.
.mdirTimePoints <- function(time_points, X, types) {
  V <- length(types)
  none <- rep(list(numeric(0)), V)
  if (is.null(time_points)) {
    return(none)
  }
  if (is.numeric(time_points) && V == 1L) {
    time_points <- list(time_points)
  }
  if (!is.list(time_points) || length(time_points) != V) {
    stop("`time_points` must be NULL or a list with an entry for each view (NULL or one position per column).",
      call. = FALSE
    )
  }
  for (v in seq_len(V)) {
    t_v <- time_points[[v]]
    if (is.null(t_v)) next
    if (!types[v] %in% c("GP", "TAGPM")) {
      stop(sprintf(
        "`time_points` was given for view %d, which has type \"%s\"; it applies to \"GP\" and \"TAGPM\" views only.",
        v, types[v]
      ), call. = FALSE)
    }
    if (!is.numeric(t_v) || length(t_v) != ncol(X[[v]]) || any(!is.finite(t_v))) {
      stop(sprintf(
        "`time_points[[%d]]` must be a finite numeric vector with one entry for each of the %d columns of view %d.",
        v, ncol(X[[v]]), v
      ), call. = FALSE)
    }
    if (anyDuplicated(t_v) > 0) {
      stop(sprintf("`time_points[[%d]]` must not repeat a position.", v), call. = FALSE)
    }
    none[[v]] <- as.numeric(t_v)
  }
  none
}

# The time points of a fit, as stored by the fitting functions (older fits have none)
.mdirFitTimePoints <- function(fit, V) {
  tp <- fit$time_points
  if (is.null(tp)) rep(list(numeric(0)), V) else tp
}
