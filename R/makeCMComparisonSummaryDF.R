#' @title Make CM Comparison Summary DF
#' @description Converts a list of consensus matrices into a long data.frame
#' ready for ggplot2. The matrices are given a common ordering defined by the
#' list entry in index ``matrix_setting_order``, this is the ``x`` and ``y``
#' entry in the output and is used in the ``ggplot2`` geom, ``geom_tile``, as
#' the ``x`` and ``y`` aesthetics.  The description of the ensembles generating
#' the matrices is given in ``model_description_df`` and is used to annotate the
#' matrices.
#' @param cms List of similarity matrices of common size.
#' @param model_description_df A data.frame with as many rows as ``cms`` has
#' entries with a common ordering. Each row is expected to contain the
#' width/number of chains used in constructing the corresponding consensus
#' matrix and the depth/iteration from each chain that is used. The column names
#' should be ``Depth`` and ``Width``.
#' @return A data frame of the mean absolute difference between consecutive 
#' ensembles from the ``model_description_df`` input.
#' @examples
#' \donttest{
#' model_description_df <- data.frame(
#'   "Depth" = rep(c(100, 500, 1000), 3),
#'   "Width" = c(rep(50, 3), rep(100, 3), rep(200, 3))
#' )
#' n_models <- nrow(model_description_df)
#' cms <- vector("list", n_models)
#' for (ii in seq(1, n_models)) {
#'   cms[[ii]] <- makeSimilarityMat(mcmc[[ii]]$allocations)
#' }
#' mean_abs_diff_df <- makeCMComparisonSummaryDF(cms, model_description_df)
#' }
#' @export
makeCMComparisonSummaryDF <- function(cms, model_description_df) {
  D <- model_description_df$Depth
  W <- model_description_df$Width

  D_considered <- unique(D)
  W_considered <- unique(W)

  number_chains <- length(W_considered)
  number_depths <- length(D_considered)

  mean_absolute_difference_df <- NULL

  # Make CMs and compare across depths
  for (ii in seq(1, number_chains)) {
    curr_w <- W_considered[ii]
    chains_used <- seq(1, curr_w)

    for (jj in seq(2, number_depths)) {
      curr_d <- D_considered[jj]

      curr_cm_index <- which(
        (model_description_df$Depth == curr_d) &
          (model_description_df$Width == curr_w)
      )

      comparison_d <- D_considered[jj - 1]
      comparison_cm_index <- which(
        (model_description_df$Depth == comparison_d) &
          (model_description_df$Width == curr_w)
      )

      # mean_absolute_difference <-
      mad <- mean(abs(cms[[curr_cm_index]] - cms[[comparison_cm_index]]))
      mad_entry <- data.frame(
        "Depth" = curr_d,
        "Width" = curr_w,
        "Mean_absolute_difference" = mad,
        "Depth_compared" = comparison_d,
        "Width_compared" = curr_w,
        "Quantity_varied" = "Depth"
      )

      if (is.null(mean_absolute_difference_df)) {
        mean_absolute_difference_df <- mad_entry
      } else {
        mean_absolute_difference_df <- rbind(mean_absolute_difference_df, mad_entry)
      }
    }
  }

  # Compare across widths
  for (jj in seq(1, number_depths)) {
    curr_d <- D_considered[jj]
    for (ii in seq(2, number_chains)) {
      curr_w <- W_considered[ii]
      chains_used <- seq(1, curr_w)
      curr_cm_index <- which((model_description_df$Depth == curr_d & model_description_df$Width == curr_w))
      comparison_w <- W_considered[ii - 1]
      comparison_cm_index <- which((model_description_df$Depth == curr_d) & (model_description_df$Width == comparison_w))


      mad <- mean(abs(cms[[curr_cm_index]] - cms[[comparison_cm_index]]))
      mad_entry <- data.frame(
        "Depth" = curr_d,
        "Width" = curr_w,
        "Mean_absolute_difference" = mad,
        "Depth_compared" = curr_d,
        "Width_compared" = comparison_w,
        "Quantity_varied" = "Width"
      )

      mean_absolute_difference_df <- rbind(mean_absolute_difference_df, mad_entry)
    }
  }
  # mean_absolute_difference_df$Depth <- factor(mean_absolute_difference_df$Depth)
  # mean_absolute_difference_df$Wdith <- factor(mean_absolute_difference_df$Width)

  mean_absolute_difference_df
}
