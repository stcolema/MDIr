#' @title Check data correct input
#' @description Internal function that checks the data passed to ``callMDI`` is
#' the correct format.
#' @param X Data passed to ``callMDI``. Should be a list of matrices each with
#' N items held in rows. ``NA`` and ``NaN`` entries are treated as missing.
#' @param types Optional character vector of the density type of each view; if
#' given, type-specific checks (e.g. categorical coding) are applied.
#' @return No return value, called for side effects.
#' @examples
#' N <- 100
#' X <- matrix(c(rnorm(N, 0, 1), rnorm(N, 3, 1)), ncol = 2, byrow = TRUE)
#' Y <- matrix(c(rnorm(N, 0, 1), rnorm(N, 3, 1)), ncol = 2, byrow = TRUE)
#' data_modelled <- list(X, Y)
#' checkDataCorrectInput(data_modelled)
checkDataCorrectInput <- function(X, types = NULL) {
  data_not_in_list <- !is.list(X)

  if (data_not_in_list) {
    stop("X is not a list. Data should be a list of matrices.")
  }

  V <- length(X)
  if (V < 1) {
    stop("X must contain at least one dataset.")
  }
  if (!is.null(types) && length(types) != V) {
    stop("``types`` must have one entry for each dataset in ``X``.")
  }
  N <- rep(0, V)

  for (v in seq_len(V)) {
    if (!is.matrix(X[[v]])) {
      stop(paste0("View ", v, " is not a matrix. Each dataset should be in matrix format."))
    }
    if (!is.numeric(X[[v]])) {
      stop(paste0("View ", v, " is not numeric."))
    }
    N[v] <- nrow(X[[v]])

    # NA and NaN are missing data. Infinite values are ambiguous (an error in
    # the data, or a censored value) and are not accepted.
    if (any(is.infinite(X[[v]]))) {
      stop(paste0("View ", v, " contains infinite values. Use NA to mark missing entries."))
    }

    # A column with no observed values carries no information and its
    # data-driven prior hyperparameters cannot be defined.
    no_observed <- which(colSums(!is.na(X[[v]])) == 0)
    if (length(no_observed) > 0) {
      stop(paste0(
        "View ", v, ": column(s) ", paste(no_observed, collapse = ", "),
        " have no observed values. Remove them before modelling."
      ))
    }
    if (nrow(X[[v]]) < 2) {
      stop(paste0("View ", v, " must have at least two items."))
    }

    if (!is.null(types) && types[v] == "C") {
      x_obs <- X[[v]][!is.na(X[[v]])]
      if (any(x_obs < 0) || any(x_obs != round(x_obs))) {
        stop(paste0(
          "View ", v, " is categorical, so its entries must be non-negative ",
          "integers coding the categories 0, 1, 2, ..."
        ))
      }
    }
  }

  if (any(N != N[1])) {
    stop(
      paste(
        "Mismatch in number of rows across datasets. All datasets must have",
        "the same number of samples/rows."
      )
    )
  }

  row_names <- row.names(X[[1]])
  for (v in seq_len(V)) {
    if (!all(row.names(X[[v]]) == row_names)) {
      stop(paste0(
        "All datasets must have the same order of row names. Dataset ", v,
        " has different row names to the first dataset, please check this."
      ))
    }
  }

  # Items that are missing in every view are informed only by the prior
  all_missing <- rowSums(sapply(X, function(x) rowSums(!is.na(x)) > 0)) == 0
  if (any(all_missing)) {
    warning(sum(all_missing), " item(s) have no observed values in any view; ",
            "their allocations are determined by the prior alone.")
  }

  invisible(NULL)
}
