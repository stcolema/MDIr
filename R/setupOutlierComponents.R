#' @title Setup outlier components
#' @description Internal function that translates the types passed to ``callMDI``
#' into the outlier component type vector.
#' @param types Density types passed to ``callMDI``.
#' @return Vector of  numbers to pass to the C++ MDI function indicating what
#' outlier component types should be used. 0, the default, means no outliers are
#' modelled, 1 uses a global MVT distribution to model outliers a la Crook et al.
setupOutlierComponents <- function(types) {
  V <- length(types)
  view_indices <- seq(1, V)
  outlier_types <- rep(0, V)

  # iterate over the types checking they are viable types
  for (v in view_indices) {
    type <- types[v]
    tagm_used <- (type == "TAGM" || type == "TAGPM")
    if (tagm_used) {
      outlier_types[v] <- 1
    }
  }


  outlier_types
}
