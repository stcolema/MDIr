#' @title Prepare Data for MDI or Mixture Model
#' @description Validates and formats data for input to callMDI or callMixtureModel.
#' Ensures consistent dimensions, ordering, and data types.
#' @param ... Matrix objects to be combined into a list for MDI, or a single 
#' matrix for mixture model
#' @param row_names Optional character vector of row names to apply consistently
#' @param ensure_numeric Logical. Convert all data to numeric type
#' @returns List of matrices (for MDI) or single matrix (for mixture model)
#' @export
#' @examples
#' X1 <- matrix(rnorm(200), ncol = 2)
#' X2 <- matrix(rnorm(300), ncol = 3)
#' data_list <- prepareData(X1, X2, row_names = paste0("Item", 1:100))
prepareData <- function(..., row_names = NULL, ensure_numeric = TRUE) {
  
  data_list <- list(...)
  n_datasets <- length(data_list)
  
  if (n_datasets == 0) {
    stop("No data provided")
  }
  
  # Validate all are matrices
  not_matrices <- !sapply(data_list, is.matrix)
  if (any(not_matrices)) {
    which_not <- which(not_matrices)
    stop(paste(
      "Inputs", paste(which_not, collapse = ", "), "are not matrices"
    ))
  }
  
  # Check dimensions
  N <- nrow(data_list[[1]])
  for (i in seq_along(data_list)) {
    if (nrow(data_list[[i]]) != N) {
      stop(paste(
        "Dataset", i, "has", nrow(data_list[[i]]),
        "rows, expected", N
      ))
    }
  }
  
  # Apply row names if provided
  if (!is.null(row_names)) {
    if (length(row_names) != N) {
      stop(paste(
        "Length of row_names (", length(row_names), 
        ") does not match number of rows (", N, ")"
      ))
    }
    data_list <- lapply(data_list, function(x) {
      rownames(x) <- row_names
      x
    })
  }
  
  # Ensure numeric if requested
  if (ensure_numeric) {
    data_list <- lapply(data_list, function(x) {
      if (!is.numeric(x)) {
        warning("Converting non-numeric data to numeric")
        x <- apply(x, 2, as.numeric)
      }
      x
    })
  }
  
  # Return single matrix if only one dataset, otherwise list
  if (n_datasets == 1) {
    data_list[[1]]
  } else {
    data_list
  }
}