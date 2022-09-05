#' @title Get likelihood
#' @description Extracts the model fit score from the mixture model output.
#' @param mcmc_output The output from the mixture model.
#' @return A data.frame containing the model log-likelihood and the associated
#' iteration.
#' @export
getLikelihood <- function(mcmc_output) {
  R <- mcmc_output$R
  thin <- mcmc_output$thin
  burn <- mcmc_output$burn
  first_recorded_iter <- burn

  iters <- seq(first_recorded_iter, R, by = thin)

  V <- mcmc_output$V
  view_indices <- seq(1, V)

  for (v in view_indices) {
    .l_df <- data.frame(
      "log_likelihood" = mcmc_output$complete_likelihood[, v],
      "view" = v,
      "iteration" = iters
    )

    if (v == 1) {
      lkl_df <- .l_df
    } else {
      lkl_df <- rbind(lkl_df, .l_df)
    }
  }

  lkl_df
}
