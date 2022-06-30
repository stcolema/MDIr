#' @title Generate initial labels
#' @description For V views, generate initial labels allowing for both 
#' semi-supervised and unsupervised views.
#' @param labels $N x V$ matrix of initial labels. The actual values passed only
#' matter for semi-supervised views (i.e. the views for which some labels are 
#' observed).
#' @param fixed $N x V$ matrix indicating which labels are observed and hence 
#' fixed. If the vth column has no 1's then this view is unsupervised.
#' @param alpha The concentration parameter (vector).
#' @param K The number of components modelled in each  view.
#' @return NULL
checkFixedInput <- function(fixed, N, V) {

  not_a_matrix <- ! is.matrix(fixed)
  
  if (not_a_matrix) {
    stop("``fixed`` must be a binary matrix.")
  }
  
  non_binary <- ! all(fixed %in% c(0, 1))
  if (non_binary) {
    stop("``fixed`` must be a binary matrix.")
  }
  
  V_in_fixed <- ncol(fixed)
  N_in_fixed <- nrow(fixed)
  
  wrong_number_of_columns <- V != V_in_fixed
  wrong_number_of_rows <- N != N_in_fixed
  wrong_dimensions <- wrong_number_of_columns & wrong_number_of_rows
  if(wrong_dimensions) {
    err_message <- paste(
      "``fixed`` must have a column for each dataset in ``X`` and the same",
      "number of rows as each dataset."
    )
    stop(err_message)
  }
  
  NULL

}
