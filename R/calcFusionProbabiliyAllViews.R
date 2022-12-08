#' @title Calculate fusion probability all views
#' @description The mean of the number of MCMC samples for which each item has 
#' the same label in each pair of views.
#' @param mcmc Output from ``runMCMCChain``
#' @param processed Has the chain been processed already (defaults to FALSE).
#' @returns A list of vector of probabilities for each item fusing across each 
#' pair of views.
#' @examples 
#' 
#' N <- 100
#' X <- matrix(c(rnorm(N, 0, 1), rnorm(N, 3, 1)), ncol = 2, byrow = TRUE)
#' Y <- matrix(c(rnorm(N, 0, 1), rnorm(N, 3, 1)), ncol = 2, byrow = TRUE)
#' data_modelled <- list(X, Y)
#' 
#' truth <- c(rep(1, N / 2), rep(2, N / 2))
#' 
#' V <- length(data_modelled)
#' 
#' # This R is much too low for real applications
#' R <- 100
#' thin <- 5
#' burn <- 10
#' 
#' K_max <- 10
#' K <- rep(K_max, V) 
#' types <- rep("G", V)
#' 
#' mcmc_out <- callMDI(data_modelled, R, thin, types, K = K)
#' calcFusionProbabiliyAllViews(mcmc_out)
#' 
#' @export
calcFusionProbabiliyAllViews <- function(mcmc, processed = FALSE) {
  V <- mcmc$V
  views <- seq(1, V)
  VC2 <- choose(V, 2)
  fusion_probabilities <- vector("list", VC2)
  
  entry <- 0
  names <- c()
  for(v in seq(1, V - 1)) {
    for(w in seq(v + 1, V)) {
      name <- paste0("fused_probs_", v, w)
      names <- c(names, name)
      entry <- entry + 1
      fusion_probabilities[[entry]] <- calcFusionProbabiliy(mcmc, v, w, processed)
    }
  }
  names(fusion_probabilities) <- names
  fusion_probabilities
}