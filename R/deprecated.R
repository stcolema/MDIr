#' @title Deprecated functions
#' @description These functions keep their earlier names and forward to the
#' current ones. They will be removed in a future release.
#' \itemize{
#'   \item \code{calcFusionProbabiliy()} is now \code{\link{calcFusionProbability}}.
#'   \item \code{calcFusionProbabiliyAllViews()} is now
#'   \code{\link{calcFusionProbabilityAllViews}}.
#' }
#' @param ... Arguments passed to the current function.
#' @return As the current function.
#' @name mdir-deprecated
#' @keywords internal
NULL

#' @rdname mdir-deprecated
#' @export
calcFusionProbabiliy <- function(...) {
  .Deprecated("calcFusionProbability", package = "mdir")
  calcFusionProbability(...)
}

#' @rdname mdir-deprecated
#' @export
calcFusionProbabiliyAllViews <- function(...) {
  .Deprecated("calcFusionProbabilityAllViews", package = "mdir")
  calcFusionProbabilityAllViews(...)
}
