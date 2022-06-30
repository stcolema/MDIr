#' @title Prepare similarity matrices for ggplot
#' @description Given a list of symmatric similarity matrices, converts them
#' to the correct formatting for ggplot2, optionally giving all a common
#' ordering.
#' @param similarity_matrices List of similarity matrices
#' @param matrix_setting_order The matrix that defines the ordering of all
#' others, if a common ordering is used. Defaults to `1`.
#' @param use_common_ordering Logical indicating if a common ordering should be
#' used. If false all matrices are ordered individually, defaults to `TRUE`.
#' @returns A long data.frame containing columns `x` (the x-axis position of the
#' entry for geom_tile()), `y` (the y-axis position of the entry for
#' geom_tile()), `Entry` (value in similarity  matrix) and `Chain` (assumes
#' chains are ordered from one to the number of chains present).
#' @examples
#' plot_df <- prepSimilarityMatricesForGGplot(my_psms)
#' plot_df |>
#'   ggplot(aes(x = x, y = y, fill = Entry)) +
#'   geom_tile() +
#'   facet_wrap(~Chain) +
#'   scale_fill_gradient(low = "#FFFFFF", high = "#146EB4")
#' @export
prepSimilarityMatricesForGGplot <- function(similarity_matrices,
                                            matrix_setting_order = 1,
                                            use_common_ordering = TRUE,
                                            ignore_checks = FALSE) {
  not_list <- !is.list(similarity_matrices)
  if (not_list) {
    stop("`similarity_matrices` must be a list of matrices.")
  }
  n_matrices <- length(similarity_matrices)

  all_symmetric_matrices <- TRUE
  mismatched_dimensions <- FALSE
  if (!ignore_checks) {
    for (ii in seq(1, n_matrices)) {
      all_symmetric_matrices <- isSymmetric.matrix(similarity_matrices[[ii]])
      if (all_symmetric_matrices) {
        if (ii == 1) {
          N <- nrow(similarity_matrices[[ii]])
        }
        N_i <- nrow(similarity_matrices[[ii]])
        mismatched_dimensions <- (N_i != N)
      }
      if ((!all_symmetric_matrices) | mismatched_dimensions) {
        stop("Matrices must all be symmetric and of the same dimension.")
      }
    }
  }

  row_order <- col_order <- findOrder(similarity_matrices[[matrix_setting_order]])

  for (ii in seq(1, n_matrices)) {
    first_iteration <- ii == 1
    .df <- prepDataForggHeatmap(similarity_matrices[[ii]], row_order, col_order)
    .df$Chain <- ii
    if (first_iteration) {
      sim_df <- .df
    } else {
      sim_df <- rbind(sim_df, .df)
    }
  }
  sim_df$Chain <- factor(sim_df$Chain)
  sim_df
}
