#' @title Prepare similarity matrices for ggplot
#' @description Given a list of symmatric similarity matrices, converts them
#' to the correct formatting for ggplot2, optionally giving all a common
#' ordering.
#' @param similarity_matrices List of similarity matrices
#' @param matrix_setting_order The matrix that defines the ordering of all
#' others, if a common ordering is used. Defaults to `1`.
#' @param use_common_ordering Logical indicating if a common ordering should be
#' used. If false all matrices are ordered individually, defaults to `TRUE`.
#' @param ignore_checks Logical indicating if checks for matrix symmetry and
#' a commmon size should be ignored. Defaults to TRUE as these checks seem
#' to break when they should not.
#' @param x_axis The name for the column defining the x-axis. Defaults to 
#' ``"x"``.
#' @param y_axis The name for the column defining the y-axis. Defaults to 
#' ``"y"``.
#' @returns A long data.frame containing columns `x` (the x-axis position of the
#' entry for geom_tile()), `y` (the y-axis position of the entry for
#' geom_tile()), `Entry` (value in similarity  matrix) and `Chain` (assumes
#' chains are ordered from one to the number of chains present).
#' @examples
#' \donttest{
#' N <- 100
#' X <- matrix(c(rnorm(N, 0, 1), rnorm(N, 3, 1)), ncol = 2, byrow = TRUE)
#' Y <- matrix(c(rnorm(N, 0, 1), rnorm(N, 3, 1)), ncol = 2, byrow = TRUE)
#'
#' truth <- c(rep(1, N / 2), rep(2, N / 2))
#' data_modelled <- list(X, Y)
#'
#' V <- length(data_modelled)
#'
#' # MCMC parameters
#' R <- 5000
#' thin <- 50
#' burn <- 1000
#'
#' K_max <- 10
#' K <- rep(K_max, V)
#' types <- rep("G", V)
#'
#' n_chains <- 3
#' mcmc_out <- runMCMCChains(data_modelled, n_chains, R, thin, types, K = K)
#' mcmc_out <- processMCMCChains(mcmc_out, burn, construct_psm = TRUE)
#'
#' psms_v1 <- list()
#' for (ii in seq(1, n_chains)) {
#'   psms_v1[[ii]] <- mcmc_out[[ii]]$psms[[1]]
#' }
#'
#' plot_df <- prepSimilarityMatricesForGGplot(psms_v1)
#' plot_df |>
#'   ggplot2::ggplot(ggplot2::aes(x = x, y = y, fill = Entry)) +
#'   ggplot2::geom_tile() +
#'   ggplot2::facet_wrap(~Chain) +
#'   ggplot2::scale_fill_gradient(low = "#FFFFFF", high = "#146EB4")
#' }
#' @export
prepSimilarityMatricesForGGplot <- function(similarity_matrices,
                                            matrix_setting_order = 1,
                                            use_common_ordering = TRUE,
                                            ignore_checks = TRUE,
                                            y_axis = "y", 
                                            x_axis = "x") {
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
    if (!use_common_ordering) {
      row_order <- col_order <- findOrder(similarity_matrices[[ii]])
    }
    .df <- prepDataForggHeatmap(similarity_matrices[[ii]], row_order, col_order, x_axis, y_axis)
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
