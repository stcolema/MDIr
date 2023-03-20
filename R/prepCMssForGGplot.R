#' @title Prepare consensus matrices for ggplot
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
#' @param matrix_setting_order The index for the consensus matrix that defines
#' the ordering of entries for the heatmap. Defaults to the final entry in
#' ``cms``.
#' @return Returns a data frame of the entries from consensus matrices (``cms``)
#' annotated by the details of the ensemble run (i.e., number of chains, chain 
#' depth, given in ``model_description_df``), in the appropriate format to be 
#' plotted as a set of heatmap facets in a ggplot2 object.
#' @examples
#' \donttest{
#' model_description_df <- data.frame(
#'   "Depth" = rep(c(100, 500, 1000), 3),
#'   "Width" = c(rep(50, 3), rep(100, 3), rep(200, 3))
#' )
#' n_models <- nrow(model_description_df)
#' cms <- vector("list", n_models)
#' for (ii in seq(1, n_models)) {
#'   cms[[ii]] <- createSimilarityMat(mcmc[[ii]]$allocations)
#' }
#' plot_df <- prepCMssForGGplot(cms, model_description_df)
#' plot_df |>
#'   ggplot2::ggplot(ggplot2::aes(x = x, y = y, fill = Entry)) +
#'   ggplot2::geom_tile() +
#'   ggplot2::facet_grid(Depth ~ Width, labeller = label_both()) +
#'   ggplot2::scale_fill_gradient(low = "#FFFFFF", high = "#146EB4")
#' }
#'
#' @export
prepCMssForGGplot <- function(cms, model_description_df,
                              matrix_setting_order = NULL) {
  if (is.null(matrix_setting_order)) {
    matrix_setting_order <- length(cms)
  }

  not_list <- !is.list(cms)
  if (not_list) {
    stop("`similarity_matrices` must be a list of matrices.")
  }

  n_matrices <- length(cms)
  n_models <- nrow(model_description_df)
  if (n_matrices != n_models) {
    .err <- paste(
      "Number of consensus matrices and number of models described",
      "in ``model_description_df`` are not matching."
    )
    stop(.err)
  }
  row_order <- col_order <- findOrder(cms[[matrix_setting_order]])

  depths <- model_description_df$Depth
  widths <- model_description_df$Width

  for (ii in seq(1, n_matrices)) {
    first_iteration <- ii == 1

    d <- depths[ii]
    w <- widths[ii]
    .df <- prepDataForggHeatmap(cms[[ii]], row_order, col_order)
    .df$Depth <- d
    .df$Width <- w

    if (first_iteration) {
      cm_df <- .df
    } else {
      cm_df <- rbind(cm_df, .df)
    }
  }
  cm_df$Depth <- factor(cm_df$Depth)
  cm_df$Width <- factor(cm_df$Width)
  cm_df
}
