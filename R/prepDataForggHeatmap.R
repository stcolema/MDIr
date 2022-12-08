#' @title Prepare data for ggplot heatmap
#' @description Converts a matrix to the correct format for heatmapping in 
#' ``ggplot2`` with an option to order rows and columns.
#' @param X Matrix.
#' @param row_order Order for rows to be represented in heatmap along the 
#' y-axis. Can be ``NULL`` (order defined on data using `hierarchical 
#' clustering), ``FALSE`` (no ordering applied), or a vector of indices.
#' @param col_order Order for columns to be represented in heatmap along the 
#' x-axis. Can be ``NULL`` (order defined on data using `hierarchical 
#' clustering), ``FALSE`` (no ordering applied), or a vector of indices.
#' @param x_axis The name for the column defining the x-axis. Defaults to 
#' ``"Feature"``.
#' @param y_axis The name for the column defining the y-axis. Defaults to 
#' ``"Item"``.
#' @returns A long data.frame containing columns `Feature` (the x-axis position
#' of the entry for geom_tile()), `Item` (the y-axis position of the entry for 
#' geom_tile()), and `Entry` (value in similarity  matrix).
#' @importFrom tidyr pivot_longer
#' @export
prepDataForggHeatmap <- function(X, 
                                 row_order = NULL, 
                                 col_order = NULL,
                                 x_axis = "Feature",
                                 y_axis = "Item") {
  
  N <- nrow(X)
  P <- ncol(X)
  
  X_not_matrix <- !is.matrix(X)
  if (X_not_matrix) {
    stop("X should be a matrix please.")
  }
  
  X_not_DF <- !is.data.frame(X)
  
  # Flag indicating if the function should cluster rows and columns
  cluster_rows <- is.null(row_order)
  cluster_cols <- is.null(col_order)
  
  # Flag indicating if no ordering is desired
  no_row_ordering <- isFALSE(row_order)
  no_col_ordering <- isFALSE(col_order)
  
  if (no_row_ordering) {
    row_order <- seq(1, N)
  }
  if (cluster_rows) {
    row_order <- findOrder(X)
  }
  
  if (no_col_ordering) {
    col_order <- seq(1, P)
  }
  if (cluster_cols) {
    col_order <- findOrder(t(X))
  }
  
  if (X_not_DF) {
    X <- data.frame(X)
  }
  
  # Order the data
  Y <- X[row_order, col_order]
  
  # The order of column names in the ordered data
  col_names <- colnames(Y)
  
  # Add the column names as a feature and set as a factor to ensure ordering is 
  # preserved
  Y[[y_axis]] <- factor(row.names(X), levels = row.names(X))
  
  # Pivot longer, from wide format to long
  Y_long <- Y |> 
    tidyr::pivot_longer(-tidyr::any_of(y_axis), values_to = "Entry", names_to = x_axis)
  
  # Make sure the features have the correct ordering
  Y_long[[x_axis]] <- factor(Y_long[[x_axis]], levels = col_names)
  
  Y_long
}

