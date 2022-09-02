#' @title Check labels
#' @description Internal function that checks the initial labels are the correct
#' type.
#' @param labels Vector of initial labels for a view in ``callMDI``.
#' @param K Number of components modelled.
#' @return No return value, called for side effects.
#' @examples
#' N <- 20
#' K <- 3
#' labels <- sample(seq(1, K), replace = TRUE, size = N)
#' checkLabels(labels, K)
checkLabels <- function(labels, K) {
  not_numeric <- any(!is.numeric(labels))

  if (not_numeric) {
    err_message <- paste(
      "Labels are not numeric. They must be contiguous",
      "integers, not factors or characters."
    )
    stop(err_message)
  }

  nas_present <- any(is.na(labels))

  if (nas_present) {
    err_message <- paste(
      "NAs present in labels. All labels must be integer values between 1 and",
      "K, the number of components modelled."
    )
    stop(err_message)
  }

  integer_labels <- as.integer(labels)
  distance_integers_to_passed <- sum(abs(labels - integer_labels))
  labels_not_integers <- (distance_integers_to_passed > .Machine$double.eps)

  if (labels_not_integers) {
    stop("Labels must be integers, not doubles.")
  }

  n_initial_components <- length(unique(labels))
  highest_labels <- max(labels)

  too_many_components <- (n_initial_components > K)

  if (too_many_components) {
    err_message <- paste(
      "There are more unique labels in the initial allocations than components",
      "modelled."
    )
    stop(err_message)
  }

  components_are_mislabelled <- (highest_labels > K)
  if (components_are_mislabelled) {
    err <- paste0(
      "The initial labels appear to be wrong. Please check that ",
      "the initial labels are a contiguous sequence. Currently there are less ",
      "unique labels than requested to be modelled (appropriately), but the ",
      "largest value used to represent a component is greater than the ",
      "the number of components modelled."
    )
    stop(err)
  }

  NULL
}
