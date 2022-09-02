#' @title Generate initial semi-sueprvised labels
#' @description For a dataset with some observed labels, generate a vector of
#' initial labels based on the proportion of each class in the observed subset.
#' @param labels The vector of labels, represented by integers, for the initial
#' labelling in a dataset.
#' @param fixed The vector of 0s and 1s indicating which labels are observed /
#' to be held fixed.
#' @return An N vector of labels.
#' @export
generateInitialSemiSupervisedLabels <- function(labels, fixed) {
  N <- length(labels)
  N_fixed <- sum(fixed)
  N_unfixed <- N - N_fixed

  observed_indices <- which(fixed == 1)
  unobserved_indices <- which(fixed == 0)

  observed_labels <- labels[observed_indices]
  observed_classes <- unique(observed_labels)

  unobserved_labels <- rep(0, N_unfixed)

  # Breakdown of class proportions
  ratio <- table(observed_labels) / N_fixed

  unobserved_labels <- sample(observed_classes, N_unfixed,
    replace = TRUE,
    prob = ratio
  )

  labels[unobserved_indices] <- unobserved_labels

  labels
}
