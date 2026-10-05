#' @title Stick breaking prior
#' @description Draw weights from the stick-breaking prior.
#' @param alpha The concentration parameter.
#' @param K The number of weights to generate.
#' @return A vector of weights.
#' @examples
#' weights <- stickBreakingPrior(1, 50)
#' @importFrom stats rbeta
stickBreakingPrior <- function(alpha, K) {
  sampleStickBreakingPrior(alpha, K)
}
