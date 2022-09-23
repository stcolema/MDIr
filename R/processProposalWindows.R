#' @title Process proposal windows
#' @description Internal function to check that the proposal windows passed to
#' MDI are of the expected type. If the default is used then the GP
#' hyper-parameters have proposal windows 0.125 (ampltidue), 0.15 (length-scale),
#' and 0.1 (noise).
#' @param proposal_windows List of the proposal windows for the Metropolis-Hastings
#' sampling of Gaussian process hyperparameters. Each entry corresponds to a
#' view. For views modelled using a Gaussian process, the first entry is the
#' proposal window for the ampltiude, the second is for the length-scale and the
#' third is for the noise. These are not used in other mixture types.
#' @param types Character vector indicating density types to use. 'G' (Gaussian
#' with diagonal covariance matrix) 'MVN' (multivariate normal), 'TAGM'
#' (t-adjust Gaussian mixture), 'GP' (MVN with Gaussian process prior on the
#' mean), 'TAGPM' (TAGM with GP prior on the mean), 'C' (categorical).
#' @return List of vectors of proposal windwos (for GP or TAGPM mixtures) or
#' NULLs (all other types).
#' @examples
#' V <- 3
#' proposal_windows <- vector("list", V)
#' types <- c("GP", "MVN", "C")
#' proposal_windows[[1]] <- c(0.15, 0.20, 0.10)
#' processProposalWindows(proposal_windows, types)
#'
processProposalWindows <- function(proposal_windows, types) {
  V <- length(types)
  gp_used <- types %in% c("GP", "TAGPM")

  null_passed <- is.null(proposal_windows)
  if (null_passed) {
    proposal_windows <- vector("list", V)
  }
  wrong_length <- length(proposal_windows) != V
  if (wrong_length) {
    err <- paste(
      "Please check that proposal windows is either a list with an",
      "entry for each view or else ``NULL``"
    )
    stop(err)
  }
  for (v in seq(1, V)) {
    if (gp_used[v]) {
      default_used <- is.null(proposal_windows[[v]])
      if (default_used) {
        proposal_windows[[v]] <- c(0.125, 0.15, 0.10)
      } else {
        wrong_number_of_windows_passed <- length(proposal_windows[[v]]) != 3
        if (wrong_number_of_windows_passed) {
          err_message <- paste(
            "Wrong number of proposal windows passed for GP",
            "hyper parameters. Please pass either ``NULL`` or a list with an",
            "entry for each view containing either ``NULL`` or a 3-vector for",
            "the ampltidue, length-scale and noise hyperparameters of the GP",
            "respectively."
          )
          stop(err_message)
        }
      }
    } else {
      window_passed_for_non_gp <- !(
        is.null(proposal_windows[[v]]) ||
          proposal_windows[[v]] == 0
      )

      if (window_passed_for_non_gp) {
        warning_message <- paste(
          "Proposal window passed for a non-GP mixture",
          "type. Please check arguments are correct."
        )
        warning(warning_message)
      }
      proposal_windows[[v]] <- 0
    }
  }
  proposal_windows
}
