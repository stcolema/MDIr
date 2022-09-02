#' @title Check fixed input
#' @description Checks if the ``fixed`` input for ``callMDI`` has the expected
#' properties.
#' @param fixed $N x V$ matrix indicating which labels are observed and hence
#' fixed. If the vth column has no 1's then this view is unsupervised.
#' @param N The number of items being clustered.
#' @param V The number of views being modelled.
#' @return No return value, called for side effects.
#' @examples
#' N <- 100
#' V <- 3
#' fixed <- matrix(0, N, V)
#' checkFixedInput(fixed, N, V)
#'
checkFixedInput <- function(fixed, N, V) {
  not_a_matrix <- !is.matrix(fixed)

  if (not_a_matrix) {
    stop("``fixed`` must be a binary matrix.")
  }

  non_binary <- !all(fixed %in% c(0, 1))
  if (non_binary) {
    stop("``fixed`` must be a binary matrix.")
  }

  V_in_fixed <- ncol(fixed)
  N_in_fixed <- nrow(fixed)

  wrong_number_of_columns <- V != V_in_fixed
  wrong_number_of_rows <- N != N_in_fixed
  wrong_dimensions <- wrong_number_of_columns & wrong_number_of_rows
  if (wrong_dimensions) {
    err_message <- paste(
      "``fixed`` must have a column for each dataset in ``X`` and the same",
      "number of rows as each dataset."
    )
    stop(err_message)
  }

  NULL
}
