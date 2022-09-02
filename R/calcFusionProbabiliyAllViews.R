#' @title Caluclate fusion probability all views
#' @description The mean of the number of MCMC samples for which each item has 
#' the same label in each pair of views.
#' @param mcmc Output from ``runMCMCChain``
#' @returns A list of vector of probabilities for each item fusing across each 
#' pair of views.
#' @export
calcFusionProbabiliyAllViews <- function(mcmc) {
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
      fusion_probabilities[[entry]] <- calcFusionProbabiliy(mcmc, v, w)
    }
  }
  names(fusion_probabilities) <- names
  fusion_probabilities
}