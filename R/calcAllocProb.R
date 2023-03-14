#' @title Calculate allocation probabilities
#' @description Calculate the empirical allocation probability for each class
#' based on the sampled allocation probabilities. Only makes sense in
#' semi-supervised views.
#' @param mcmc_samples Output from ``callMDI``.
#' @param view The view for which to calculate the allocation probabilities.
#' @param burn The number of samples to discard.
#' @param method The point estimate to use. ``method = 'mean'`` or
#' ``method = 'median'``. ``'median'`` is the default.
#' @return An N x K matrix of class probabilities.
#' @examples
#'
#' N <- 100
#' X <- matrix(c(rnorm(N, 0, 1), rnorm(N, 3, 1)), ncol = 2, byrow = TRUE)
#' Y <- matrix(c(rnorm(N, 0, 1), rnorm(N, 3, 1)), ncol = 2, byrow = TRUE)
#'
#' truth <- c(rep(1, N / 2), rep(2, N / 2))
#' data_modelled <- list(X, Y)
#'
#' V <- length(data_modelled)
#'
#' # These values are much too low for real applications
#' R <- 50
#' thin <- 5
#' burn <- 10
#'
#' alpha <- rep(1, V)
#' K_max <- 10
#' K <- rep(K_max, V)
#' labels <- matrix(1, nrow = N, ncol = V)
#' fixed <- matrix(0, nrow = N, ncol = V)
#'
#' # A random quarter of labels are known in view 1
#' fixed[, 1] <- sample(c(0, 1), N, replace = TRUE, prob = c(3, 1))
#' labels[, 1] <- generateInitialSemiSupervisedLabels(truth, fixed = fixed[, 1])
#' types <- rep("G", V)
#'
#' mcmc_out <- callMDI(data_modelled, R, thin, types, K, labels, fixed, alpha)
#' calcAllocProb(mcmc_out, 1)
#'
#' @export
calcAllocProb <- function(mcmc_samples, view, burn = 0, method = "mean") {
  R <- mcmc_samples$R
  thin <- mcmc_samples$thin
  V <- mcmc_samples$V
  
  .alloc <- mcmc_samples$allocation_probabilities[[view]]
  
  if (burn > 0) {
    if (burn > R) {
      stop("Burn in exceeds number of iterations run.")
    }
    
    eff_burn <- floor(burn / thin)
    dropped_samples <- seq(1, eff_burn)
    .alloc <- .alloc[, , -dropped_samples]
  }
  
  if (view > V) {
    .err <- paste(
      "Requested view not in ``mcmc_samples``. Please check that the requested",
      "view is less than or equal to V, the number of views modelled."
    )
    stop(.err)
  }
  probs <- NULL
  
  if (method == "median") {
    probs <- apply(.alloc, c(1, 2), median)
  }
  if (method == "mean") {
    probs <- rowSums(.alloc, dims = 2) / dim(.alloc)[3]
  }
  if (length(probs) == 1) {
    if (is.null(probs)) {
      stop("``method`` must be one of 'mean' or 'median'")
    }
  }
  probs
}
