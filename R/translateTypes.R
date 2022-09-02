#' @title Translate types
#' @description Converts types from strings to integers passed to C++.
#' @param types Density types passed to ``callMDI``. Should be a vector of
#' strings. Viable options are ``'MVN'``, ``'TAGM'`` and ``'C'``.
#' @return Vector of integers
translateTypes <- function(types) {
  checkTypes(types)

  V <- length(types)
  view_indices <- seq(1, V)
  density_types <- rep(0, V)

  # The map from the inputted type to the form passed to the C++
  # Note that TAGM differs from MVN in the outlier component used, not the
  # density
  type_map <- data.frame(
    "Input_type" = c("MVN", "TAGM", "C", "G", "GP", "TAGPM"),
    "Used_type" = c(1, 1, 2, 0, 3, 3)
  )

  # iterate over the types checking they are viable types
  for (v in view_indices) {
    type <- types[v]
    map_row <- which(type_map$Input_type == type)
    translated_type <- type_map$Used_type[map_row]

    density_types[v] <- translated_type
  }

  density_types
}
