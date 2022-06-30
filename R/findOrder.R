#' @title Prepare data for ggplot heatmap
#' @description Find the ordering of the rows of X based on a hierarchical 
#' clustering of X for linkage ``linkage`` and distance ``dist``.
#' @param X Matrix.
#' @param dist Distance used when defining a distance matrix of X.
#' @param linkage Type of linkage used in the hierarchical clustering of X.
#' @returns The order of the rows of X based on a hierarchical clustering.
#' @importFrom stats hclust dist
#' @export
findOrder <- function(X, dist = "euclidean", linkage = "complete") {
  stats::hclust(stats::dist(X, method = dist), method =  linkage)$order
}
