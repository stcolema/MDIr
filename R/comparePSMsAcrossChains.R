#' @title Compare Posterior Similarity Matrices Across Chains
#' @description Takes the PSMs from the output of ``processMCMCChains`` and
#' converts them to a data.frame prepared for ``ggplot2``.
#' @param mcmc Output of ``processMCMCChains`` with the argument
#' ``construct_psm = TRUE``.
#' @param matrix_setting_order The matrix that defines the ordering of all
#' others, if a common ordering is used. Defaults to `1`.
#' @param use_common_ordering Logical indicating if a common ordering should be
#' used. If false all matrices are ordered individually, defaults to `TRUE`.
#' @param ignore_checks Logical indicating if checks for matrix symmetry and
#' a commmon size should be ignored. Defaults to TRUE as these checks seem
#' to break when they should not.
#' @returns A long data.frame containing columns `x` (the x-axis position of the
#' entry for geom_tile()), `y` (the y-axis position of the entry for
#' geom_tile()), `Entry` (value in similarity  matrix), `Chain` (assumes
#' chains are ordered from one to the number of chains present) and `View`.
#' @examples
#' \donttest{
#' N <- 500
#' K <- 4
#' P <- 10
#'
#' X <- generateSimulationDataset(K, N, P)
#' Y <- generateSimulationDataset(K, N, P)
#' Z <- generateSimulationDataset(K, N, P)
#'
#' row.names(Z$data) <- row.names(Y$data) <- row.names(X$data)
#'
#' data_modelled <- list(X$data, Y$data, Z$data)
#'
#' n_chains <- 4
#' R <- 10000
#' thin <- 25
#' types <- c("G", "G", "G")
#' K <- c(10, 10, 10)
#' mcmc <- runMCMCChains(data_modelled, n_chains, R, thin, types, K = K)
#'
#' burn <- 2500
#' mcmc <- processMCMCChains(mcmc, burn, construct_psm = TRUE)
#' psm_df <- comparePSMsAcrossChains(mcmc)
#'
#' psm_df |>
#'   ggplot2::ggplot(ggplot2::aes(x = x, y = y, fill = Entry)) +
#'   ggplot2::geom_tile() +
#'   ggplot2::facet_grid(View ~ Chain, labeller = ggplot2::label_both) +
#'   ggplot2::scale_fill_gradient(low = "#FFFFFF", high = "#146EB4") +
#'   ggplot2::labs(x = "Item", y = "Item", fill = "Coclustering\nproportion") +
#'   ggplot2::theme(
#'     axis.text = ggplot2::element_blank(),
#'     axis.ticks = ggplot2::element_blank(),
#'     panel.grid = ggplot2::element_blank(),
#'     axis.title.y = ggplot2::element_text(size = 10.5),
#'     axis.title.x = ggplot2::element_text(size = 10.5),
#'     plot.title = ggplot2::element_text(size = 18, face = "bold"),
#'     plot.subtitle = ggplot2::element_text(size = 14),
#'     strip.text.x = ggplot2::element_text(size = 10.5),
#'     legend.text = ggplot2::element_text(size = 10.5)
#'   )
#' }
#' @export
comparePSMsAcrossChains <- function(mcmc,
                                    matrix_setting_order = 1,
                                    use_common_ordering = TRUE,
                                    ignore_checks = TRUE) {
  n_chains <- length(mcmc)
  V <- mcmc[[1]]$V

  chain_inds <- seq(1, n_chains)
  view_inds <- seq(1, V)

  psm_df <- NULL
  first_iteration <- is.null(psm_df)
  psm_lst <- list()
  for (v in view_inds) {
    psm_lst[[v]] <- list()
    for (ii in chain_inds) {
      psm_lst[[v]][[ii]] <- mcmc[[ii]]$psm[[v]]
    }
  }
  for (v in view_inds) {
    .psm_df <- prepSimilarityMatricesForGGplot(psm_lst[[v]],
      matrix_setting_order = matrix_setting_order,
      use_common_ordering = use_common_ordering,
      ignore_checks = ignore_checks
    )
    .psm_df$View <- v
    if (first_iteration) {
      psm_df <- .psm_df
      first_iteration <- FALSE
    } else {
      psm_df <- rbind(psm_df, .psm_df)
    }
  }
  psm_df$Chain <- factor(psm_df$Chain, levels = chain_inds)
  psm_df$View <- factor(psm_df$View, levels = view_inds)
  psm_df
}
