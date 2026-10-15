# Internal: the monitored scalar quantities of one chain, one column each, over
# every stored sample, with the iteration at which each was saved.
.mdirTraceMatrix <- function(x) {
  n <- nrow(x$mass)
  cols <- list(complete_likelihood = as.numeric(x$complete_likelihood))
  if (!is.null(x$joint_likelihood)) {
    cols$joint_likelihood <- as.numeric(x$joint_likelihood)
  }
  for (v in seq_len(x$V)) {
    cols[[sprintf("mass[%d]", v)]] <- x$mass[, v]
    cols[[sprintf("occupied[%d]", v)]] <- .mdirOccupied(x, v, seq_len(n))
  }
  if (x$V > 1 && !is.null(x$phis)) {
    pairs <- utils::combn(x$V, 2)
    for (i in seq_len(ncol(pairs))) {
      cols[[sprintf("phi[%d,%d]", pairs[1, i], pairs[2, i])]] <- x$phis[, i]
    }
  }
  # Sample s of an unprocessed chain was saved at iteration (s - 1) thin; a
  # processed chain has dropped its first floor(burn / thin) + 1 samples.
  dropped <- if (.mdirIsProcessed(x)) floor(x$burn / x$thin) + 1 else 0
  list(values = do.call(cbind, cols), iteration = (dropped + seq_len(n) - 1) * x$thin)
}

# Internal: long data frame (chain, iteration, quantity, value) of the traces
.mdirTraceTable <- function(chains, pars = NULL) {
  out <- lapply(seq_along(chains), function(i) {
    tm <- .mdirTraceMatrix(chains[[i]])
    data.frame(
      chain = factor(i, levels = seq_along(chains)),
      sample = rep(seq_len(nrow(tm$values)), ncol(tm$values)),
      iteration = rep(tm$iteration, ncol(tm$values)),
      quantity = factor(rep(colnames(tm$values), each = nrow(tm$values)), levels = colnames(tm$values)),
      value = as.numeric(tm$values)
    )
  })
  out <- do.call(rbind, out)
  if (!is.null(pars)) {
    unknown <- setdiff(pars, levels(out$quantity))
    if (length(unknown) > 0) {
      stop("Unknown quantit", if (length(unknown) > 1) "ies" else "y", ": ", paste(unknown, collapse = ", "),
           ". Available: ", paste(levels(out$quantity), collapse = ", "), ".", call. = FALSE)
    }
    out <- out[out$quantity %in% pars, , drop = FALSE]
    out$quantity <- factor(out$quantity, levels = pars)
  }
  out
}

# Internal: the retained (post burn in) allocations of one view, all chains stacked
.mdirPooledAllocations <- function(chains, view, burn) {
  do.call(rbind, lapply(chains, function(ch) {
    idx <- .mdirRetained(ch, burn)$index
    matrix(ch$allocations[idx, , view], nrow = length(idx))
  }))
}

#' @title Plot an mdir fit
#' @description Graphical summaries of one or several chains from
#' \code{\link{callMDI}}, \code{\link{runMCMCChains}} or \code{\link{fitMDI}},
#' as ggplot2 objects. The default, \code{type = "trace"}, is the first thing to
#' look at: the sampled value of every monitored quantity against the iteration,
#' one line per chain, with the burn in marked (the initial state is left out). Chains that overlap and show no
#' trend are consistent with convergence; chains in different places are not.
#' None of the plots proves convergence (see \code{\link{assessConvergence}}).
#'
#' \describe{
#'   \item{\code{"trace"}}{Trace plots of the complete- and joint-data
#'   log-likelihoods, the concentrations (\code{mass}), the numbers of occupied
#'   components, and the view-association parameters \code{phi}. Quantities that
#'   do not depend on how clusters are labelled are used throughout.}
#'   \item{\code{"density"}}{Density of each quantity after the burn in, one
#'   curve per chain. Curves that differ show chains that disagree.}
#'   \item{\code{"rhat"}}{Rank-normalised split-Rhat and effective sample sizes
#'   (Vehtari et al., 2021) of every monitored quantity against their
#'   thresholds; see \code{\link{assessConvergence}}.}
#'   \item{\code{"psm"}}{Posterior similarity matrix of one view: the proportion
#'   of retained draws (all chains pooled) in which two items share a cluster,
#'   with the items ordered by hierarchical clustering. Blocks on the diagonal
#'   are clusters; off-diagonal grey is uncertainty about the partition.}
#'   \item{\code{"fusion"}}{For each pair of views, the fusion probability of each
#'   item (see \code{\link{calcFusionProbability}}) sorted from lowest to highest.
#'   It compares component indices across views, so it measures fusion of
#'   clusters only if the indices are aligned between views.}
#' }
#' @param x Output of \code{\link{runMCMCChains}} or \code{\link{fitMDI}}
#' (\code{mdir_fit_list}), or of \code{\link{callMDI}} (\code{mdir_fit}).
#' @param type One of \code{"trace"}, \code{"density"}, \code{"rhat"}, \code{"psm"} and
#' \code{"fusion"}.
#' @param pars Optional character vector naming the quantities for
#' \code{"trace"} and \code{"density"}, for example \code{"phi[1,2]"}. The
#' default is all of them.
#' @param burn Number of iterations treated as burn in. The default is half of
#' \code{R}, as in \code{\link{assessConvergence}}. Ignored for chains that have
#' been through \code{\link{processMCMCChain}}.
#' @param view The view shown by \code{type = "psm"}.
#' @param ... Unused.
#' @return A ggplot2 object.
#' @references Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and
#' Burkner, P.-C. (2021). Rank-normalization, folding, and localization: an
#' improved Rhat for assessing convergence of MCMC. \emph{Bayesian Analysis},
#' 16(2), 667-718.
#' @seealso \code{\link{summary.mdir_fit_list}}, \code{\link{assessConvergence}}
#' @export
#' @examples
#' set.seed(1)
#' X <- lapply(1:2, function(v) {
#'   m <- matrix(rnorm(40 * 2, rep(c(0, 3), each = 20)), 40, 2)
#'   rownames(m) <- 1:40
#'   m
#' })
#' fit <- runMCMCChains(X, n_chains = 2, R = 200, thin = 5, types = c("G", "G"), K = c(4, 4))
#' plot(fit)
#' plot(fit, type = "psm", view = 1)
plot.mdir_fit_list <- function(x,
                               type = c("trace", "density", "rhat", "psm", "fusion"),
                               pars = NULL,
                               burn = NULL,
                               view = 1,
                               ...) {
  type <- match.arg(type)
  chains <- unclass(x)[seq_along(x)]
  first <- chains[[1]]
  if (type %in% c("psm") && (!is.numeric(view) || length(view) != 1 || view < 1 || view > first$V)) {
    stop("`view` must be a single view between 1 and ", first$V, ".", call. = FALSE)
  }
  if (type == "fusion" && first$V < 2) {
    stop("Fusion probabilities need at least two views.", call. = FALSE)
  }
  switch(type,
    trace = .mdirPlotTrace(chains, pars, burn),
    density = .mdirPlotDensity(chains, pars, burn),
    rhat = plot(attr(x, "convergence") %||% assessConvergence(chains, burn = burn)),
    psm = .mdirPlotPSM(chains, view, burn),
    fusion = .mdirPlotFusion(chains, burn)
  )
}

#' @rdname plot.mdir_fit_list
#' @export
plot.mdir_fit <- function(x, ...) {
  plot.mdir_fit_list(structure(list(x), class = c("mdir_fit_list", "list")), ...)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

.mdirPlotTrace <- function(chains, pars, burn) {
  d <- .mdirTraceTable(chains, pars)
  # The initial state is a random draw far from the posterior and would set the scale
  d <- d[d$iteration > 0, , drop = FALSE]
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$iteration, y = .data$value, colour = .data$chain)) +
    ggplot2::geom_line(linewidth = 0.3, alpha = 0.8) +
    ggplot2::facet_wrap(~quantity, scales = "free_y") +
    ggplot2::labs(x = "Iteration", y = NULL, colour = "Chain", title = "Trace plots") +
    ggplot2::theme_minimal()
  if (!.mdirIsProcessed(chains[[1]])) {
    burn_used <- .mdirRetained(chains[[1]], burn)$burn
    p <- p +
      ggplot2::geom_vline(xintercept = burn_used, linetype = "dashed", colour = "grey40") +
      ggplot2::labs(subtitle = paste0("Dashed line: end of burn in (", burn_used, " iterations)"))
  }
  p
}

.mdirPlotDensity <- function(chains, pars, burn) {
  d <- .mdirTraceTable(chains, pars)
  retained <- lapply(chains, function(ch) .mdirRetained(ch, burn)$index)
  d <- d[mapply(function(ch, smp) smp %in% retained[[ch]], as.integer(d$chain), d$sample), , drop = FALSE]
  # A quantity that never changes (for example one occupied component) has no density
  varies <- tapply(d$value, d$quantity, function(v) length(unique(v)) > 1)
  d <- d[varies[as.character(d$quantity)], , drop = FALSE]
  ggplot2::ggplot(d, ggplot2::aes(x = .data$value, colour = .data$chain)) +
    ggplot2::geom_density() +
    ggplot2::facet_wrap(~quantity, scales = "free") +
    ggplot2::labs(x = NULL, y = "Density", colour = "Chain", title = "Posterior densities by chain") +
    ggplot2::theme_minimal()
}

.mdirPlotPSM <- function(chains, view, burn) {
  a <- .mdirPooledAllocations(chains, view, burn)
  psm <- createSimilarityMat(a)
  ord <- stats::hclust(stats::as.dist(1 - psm), method = "average")$order
  ids <- chains[[1]]$sample_ids %||% seq_len(ncol(a))
  N <- ncol(a)
  d <- data.frame(
    item_x = rep(seq_len(N), times = N),
    item_y = rep(seq_len(N), each = N),
    probability = as.numeric(psm[ord, ord])
  )
  ggplot2::ggplot(d, ggplot2::aes(x = .data$item_x, y = .data$item_y, fill = .data$probability)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_gradient(low = "white", high = "navy", limits = c(0, 1)) +
    ggplot2::coord_equal(expand = FALSE) +
    ggplot2::scale_y_reverse() +
    ggplot2::labs(
      x = "Item (ordered by clustering)", y = NULL, fill = "P(same cluster)",
      title = paste0("Posterior similarity matrix, view ", view),
      subtitle = paste0(nrow(a), " retained draws from ", length(chains), " chain", if (length(chains) > 1) "s" else "")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text = ggplot2::element_blank(), panel.grid = ggplot2::element_blank())
}

.mdirPlotFusion <- function(chains, burn) {
  V <- chains[[1]]$V
  pairs <- utils::combn(V, 2)
  d <- do.call(rbind, lapply(seq_len(ncol(pairs)), function(i) {
    a <- .mdirPooledAllocations(chains, pairs[1, i], burn)
    b <- .mdirPooledAllocations(chains, pairs[2, i], burn)
    p <- sort(colMeans(a == b))
    data.frame(
      pair = sprintf("views %d and %d", pairs[1, i], pairs[2, i]),
      rank = seq_along(p), probability = p
    )
  }))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$rank, y = .data$probability)) +
    ggplot2::geom_point(size = 0.8) +
    ggplot2::facet_wrap(~pair) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Item (sorted by fusion probability)", y = "P(same component in both views)",
      title = "Fusion probabilities"
    ) +
    ggplot2::theme_minimal()
}

#' @title Plot convergence diagnostics
#' @description Rank-normalised split-Rhat (left) and bulk and tail effective
#' sample sizes (right) of every monitored quantity, with the thresholds used by
#' \code{\link{assessConvergence}} as dashed lines. Points beyond a threshold
#' are coloured.
#' @param x Output of \code{\link{assessConvergence}}.
#' @param ... Unused.
#' @return A ggplot2 object.
#' @export
#' @examples
#' set.seed(1)
#' X <- list(matrix(rnorm(80, rep(c(0, 3), each = 20)), 40, 2))
#' rownames(X[[1]]) <- 1:40
#' fit <- runMCMCChains(X, n_chains = 2, R = 400, thin = 2, types = "G", K = 4)
#' plot(assessConvergence(fit))
plot.mdir_convergence <- function(x, ...) {
  threshold <- attr(x, "threshold")
  min_ess <- attr(x, "min_ess")
  q <- factor(x$quantity, levels = rev(x$quantity))
  d <- rbind(
    data.frame(quantity = q, metric = "Rhat (max of bulk and folded)", value = x$rhat,
               series = "Rhat", flagged = !is.na(x$rhat) & x$rhat >= threshold),
    data.frame(quantity = q, metric = "Effective sample size", value = x$ess_bulk,
               series = "bulk", flagged = !is.na(x$ess_bulk) & x$ess_bulk < min_ess),
    data.frame(quantity = q, metric = "Effective sample size", value = x$ess_tail,
               series = "tail", flagged = !is.na(x$ess_tail) & x$ess_tail < min_ess)
  )
  d$metric <- factor(d$metric, levels = c("Rhat (max of bulk and folded)", "Effective sample size"))
  limits <- data.frame(
    metric = factor(c("Rhat (max of bulk and folded)", "Effective sample size"), levels = levels(d$metric)),
    x = c(threshold, min_ess)
  )
  ggplot2::ggplot(d, ggplot2::aes(x = .data$value, y = .data$quantity, shape = .data$series, colour = .data$flagged)) +
    ggplot2::geom_vline(data = limits, ggplot2::aes(xintercept = .data$x), linetype = "dashed", colour = "grey40") +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~metric, scales = "free_x") +
    ggplot2::scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "firebrick"), guide = "none") +
    ggplot2::labs(
      x = NULL, y = NULL, shape = NULL,
      title = "Convergence diagnostics",
      subtitle = paste0(attr(x, "n_chains"), " chain(s); dashed: Rhat = ", threshold, " and ESS = ", min_ess)
    ) +
    ggplot2::theme_minimal()
}

#' @title Plot a weighted ensemble from smcMDI
#' @description The annealing path of \code{\link{smcMDI}}: the inverse
#' temperature, the effective sample size of the weights and the conditional
#' effective sample size used to choose the next temperature, at each step, with
#' the steps at which the particles were resampled marked. A path on which the
#' effective sample size collapses at a sharp step is the usual sign of an
#' unreliable ensemble (see \code{\link{smcDiagnostics}}).
#' @param x Output of \code{\link{smcMDI}}.
#' @param ... Unused.
#' @return A ggplot2 object.
#' @export
plot.mdir_smc <- function(x, ...) {
  tr <- x$trace
  d <- rbind(
    data.frame(step = tr$step, panel = "Inverse temperature", value = tr$beta, resampled = tr$resampled),
    data.frame(step = tr$step, panel = "Effective sample size (fraction of particles)",
               value = tr$ess / x$n_particles, resampled = tr$resampled),
    data.frame(step = tr$step, panel = "Conditional ESS fraction", value = tr$cess, resampled = tr$resampled)
  )
  d$panel <- factor(d$panel, levels = unique(d$panel))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$step, y = .data$value)) +
    ggplot2::geom_line(linewidth = 0.4) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$resampled), size = 1) +
    ggplot2::scale_colour_manual(values = c("FALSE" = "grey60", "TRUE" = "firebrick"), name = "Resampled") +
    ggplot2::facet_wrap(~panel, scales = "free_y", ncol = 1) +
    ggplot2::labs(x = "Step", y = NULL, title = "Annealing path",
                  subtitle = paste0(x$n_particles, " particles, ", x$schedule_type, " schedule")) +
    ggplot2::theme_minimal()
}
