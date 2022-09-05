#' @title Check data correct input
#' @description Internal function that checks the data passed to ``callMDI`` is
#' the correct format.
#' @param X Data passed to ``callMDI``. Should be a list of matrices each with
#' N items held in rows.
#' @return No return value, called for side effects.
#' @examples
#' N <- 100
#' X <- matrix(c(rnorm(N, 0, 1), rnorm(N, 3, 1)), ncol = 2, byrow = TRUE)
#' Y <- matrix(c(rnorm(N, 0, 1), rnorm(N, 3, 1)), ncol = 2, byrow = TRUE)
#' data_modelled <- list(X, Y)
#' checkDataCorrectInput(data_modelled)
checkDataCorrectInput <- function(X) {
  data_not_in_list <- !is.list(X)

  if (data_not_in_list) {
    stop("X is not a list. Data should be a list of matrices.")
  }

  V <- length(X)
  view_indices <- seq(1, V)
  N <- rep(0, V)
  different_numbers_of_items <- rep(T, V)

  for (v in view_indices) {
    vth_view_is_not_a_matrix <- !is.matrix(X[[v]])

    if (vth_view_is_not_a_matrix) {
      err_message <- paste0(
        "View ",
        v,
        " is not a matrix. Each dataset should be in matrix format."
      )
      stop(err_message)
    }

    N[v] <- nrow(X[[v]])
    different_numbers_of_items[v] <- (N[v] != N[1])
  }

  mismatching_number_of_rows <- any(different_numbers_of_items)

  if (mismatching_number_of_rows) {
    stop(
      paste(
        "Mismatch in number of rows across datasets. All datasets must have",
        "the same number of samples/rows."
      )
    )
  }

  row_names <- row.names(X[[1]])
  for (v in view_indices) {
    .x <- X[[v]]
    .row_names_not_matching <- !all(row.names(.x) == row_names)
    if (.row_names_not_matching) {
      err <- paste0(
        "All datasets must have the same order of row names. Dataset ",
        v,
        " has different row names to the first dataset, please check this."
      )
      stop(err)
    }
  }

  NULL
}
