#' @title Sample view prior labels
#' @description Generate labels from the stick-breaking prior.
#' @param alpha The concentration parameter for the stick-breaking prior.
#' @param K The number of components to include (the upper bound on the number of unique labels generated).
#' @param N The number of labels to generate.
#' @return A vector of labels.
#' @examples
#' initial_labels <- sampleViewPriorLabels(1, 50, 100)
sampleViewPriorLabels <- function(alpha, K, N) {
  w <- sampleStickBreakingPrior(alpha, K)
  initial_labels <- sample(seq(1, K), N, replace = TRUE, prob = w)
}
