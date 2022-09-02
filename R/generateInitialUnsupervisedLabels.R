#' @title Generate initial unsupervised labels
#' @description For a dataset with no observed labels, generate a vector of
#' initial labels based on the stick-breaking prior for component weights.
#' @param N The number of items modelled in the dataset.
#' @param alpha The concentration parameter.
#' @param K The number of components modelled.
#' @return An N vector of labels.
#' @export
generateInitialUnsupervisedLabels <- function(N, alpha, K) {

  # Breakdown of class proportions
  weights <- sampleStickBreakingPrior(alpha, K)
  components_modelled <- seq(1, K)

  labels <- sample(components_modelled, N,
    replace = TRUE,
    prob = weights
  )

  labels
}
