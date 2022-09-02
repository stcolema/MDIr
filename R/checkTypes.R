#' @title Check types
#' @description Internal function that checks the types passed to ``callMDI``
#' are the correct format.
#' @param types Density types passed to ``callMDI``. Should be a vector of
#' strings. Viable options are ``'MVN'``, ``'TAGM'`` and ``'C'``.
#' @return No return value, called for side effects.
#' @examples
#' types <- c("GP", "MVN", "C")
#' checkTypes(types)
checkTypes <- function(types) {
  allowed_types <- c("MVN", "TAGM", "C", "G", "GP", "TAGPM")
  not_a_vector <- !is.vector(types)

  if (not_a_vector) {
    stop("``types`` must be a vector of strings.")
  }

  V <- length(types)
  view_indices <- seq(1, V)

  # iterate over the types checking they are viable types
  for (v in view_indices) {
    type <- types[v]
    wrong_type <- !type %in% allowed_types
    if (wrong_type) {
      stop("Type not recognised. Please use 'MVN', 'TAGM', 'G', 'GP', 'TAGPM' or 'C'.")
    }
  }

  NULL
}
