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
#' @returns A long data.frame containing columns `x` (the x-axis position of the
#' entry for geom_tile()), `y` (the y-axis position of the entry for
#' geom_tile()), and `Entry` (value in similarity  matrix).
#' @importFrom tidyr pivot_longer any_of
#' @importFrom stringr str_extract
#' @export
prepDataForggHeatmap <- function(X, row_order = NULL, col_order = NULL) {
  N <- nrow(X)
  P <- ncol(X)

  X_not_matrix <- !is.matrix(X)
  if (X_not_matrix) {
    stop("X should be a matrix please.")
  }

  X_not_DF <- !is.data.frame(X)

  cluster_rows <- is.null(row_order)
  cluster_cols <- is.null(col_order)

  no_row_ordering <- isFALSE(row_order)
  no_col_ordering <- isFALSE(col_order)

  if (X_not_DF) {
    X <- data.frame(X)
  }
  Y <- X

  Y$y <- seq(1, N, by = 1)
  # Y$y <- seq(N, 1, by = -1)
  if (no_row_ordering) {
    row_order <- seq(1, N)
  }
  if (cluster_rows) {
    row_order <- findOrder(X)
  }

  Y$y <- match(Y$y, row_order)
  feature_name_order <- colnames(Y)[-ncol(Y)]
  Z <- tidyr::pivot_longer(Y, -tidyr::any_of("y"), names_to = "Feature", values_to = "Entry")
  # Z$x <- as.numeric(stringr::str_extract(Z$Feature, "[:digit:]+"))
  Z$x <- match(Z$Feature, feature_name_order)

  if (no_col_ordering) {
    col_order <- Z$x
  }
  if (cluster_cols) {
    col_order <- findOrder(t(X))
  }

  Z$x <- match(Z$x, col_order)
  Z$Feature <- NULL
  Z
}
