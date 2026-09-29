#' @title Plot likelihoods
#' @description Plots the model fit for multiple chains.
#' @param mcmc_outputs The output from ``runMCMCChains`` (or a list of
#' processed chains).
#' @param colour_by_chain Logical indcating if plots should be coloured by chain
#' or all the same colour. Defaults to ``TRUE``.
#' @return A ggplot2 object. Line plot of likelihood across iteration.
#' @importFrom ggplot2 ggplot aes geom_line facet_wrap label_both
#' @export
plotLikelihoods <- function(mcmc_outputs,
                            colour_by_chain = TRUE) {
  lkl_lst <- lapply(mcmc_outputs, getLikelihood)

  n_chains <- length(lkl_lst)
  for (ii in seq(1, n_chains)) {
    lkl_lst[[ii]]$chain <- if (is.null(mcmc_outputs[[ii]]$Chain)) ii else mcmc_outputs[[ii]]$Chain
  }

  lkl_df <- do.call(rbind, lkl_lst)
  lkl_df$chain <- factor(lkl_df$chain)

  if (colour_by_chain) {
    p <- ggplot2::ggplot(
      data = lkl_df,
      mapping = ggplot2::aes(
        x = .data$iteration,
        y = .data$log_likelihood,
        colour = .data$chain
      )
    ) +
      ggplot2::geom_line()
  } else {
    p <- ggplot2::ggplot(
      data = lkl_df,
      mapping = ggplot2::aes(
        x = .data$iteration,
        y = .data$log_likelihood,
        group = .data$chain
      )
    ) +
      ggplot2::geom_line()
  }
  p <- p +
    ggplot2::facet_wrap(~type, scales = "free_y", labeller = ggplot2::label_both)

  p
}
