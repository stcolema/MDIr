# S3 print(), summary() and format() methods for the objects the fitting
# functions return: "mdir_fit" (one chain, from callMDI()), "mdir_fit_list"
# (several chains, from runMCMCChains() and fitMDI()) and, in
# convergenceDiagnostics.R, "mdir_convergence".
#
# The objects are still plain lists underneath (x$allocations, chains[[i]]$phis
# and so on keep working); the class only changes what the console shows. Before
# these methods, printing a fit dumped every sampled array (thousands of lines).
# The layout follows mclust (short print(), rule-delimited summary() with a
# clustering table) and Stan (a header stating chains/iterations/draws, then a
# table of posterior summaries).

# Rule used to delimit summary() output, as in mclust
.mdirRule <- function(width = 52) strrep("-", width)

.mdirPlural <- function(n, word) paste0(n, " ", word, if (n == 1) "" else "s")

# Wrap a paragraph to the console width
.mdirWrap <- function(text, indent = 0) {
  paste(strwrap(text, width = getOption("width") - 2, exdent = indent), collapse = "\n")
}

# Time as a short string, e.g. "12.3 secs"
.mdirFormatTime <- function(x) {
  if (is.null(x)) {
    return("not recorded")
  }
  format(signif(x, 3))
}

# One row per view: what was modelled. A fit that has had a burn in applied
# has fewer saved samples than R / thin + 1 but the same view information.
.mdirViewTable <- function(x) {
  data.frame(
    view = seq_len(x$V),
    type = x$types,
    P = vapply(x$P, function(p) as.integer(p), integer(1)),
    K = as.integer(x$K),
    supervision = ifelse(x$Semisupervised, "semi-supervised", "unsupervised"),
    missing = sprintf("%.0f%%", 100 * x$missing_proportion),
    stringsAsFactors = FALSE
  )
}

# Has processMCMCChain() been applied? (It adds point estimates and drops the
# burn in from the stored arrays.)
.mdirIsProcessed <- function(x) !is.null(x$pred)

# Indices into the stored samples that are used for summaries. For an
# unprocessed fit saved sample 1 is the initial state, and sample s is
# iteration (s - 1) * thin, as in assessConvergence().
.mdirRetained <- function(x, burn = NULL) {
  n_stored <- nrow(x$mass)
  if (.mdirIsProcessed(x)) {
    return(list(index = seq_len(n_stored), burn = x$burn, burn_default = FALSE))
  }
  burn_default <- is.null(burn)
  if (burn_default) {
    burn <- floor(x$R / 2)
  }
  first <- floor(burn / x$thin) + 2
  if (first > n_stored) {
    stop(
      "The burn in (", burn, ") leaves no saved samples; ", n_stored - 1,
      " iterations were saved (R = ", x$R, ", thin = ", x$thin, ").",
      call. = FALSE
    )
  }
  list(index = seq(first, n_stored), burn = burn, burn_default = burn_default)
}

# Number of occupied components in each retained sample of a view
.mdirOccupied <- function(x, v, index) {
  z <- matrix(x$allocations[index, , v], nrow = length(index))
  apply(z, 1, function(row) length(unique(row)))
}

# Posterior summaries of the phi parameters (view-association), Stan style
.mdirPhiTable <- function(x, index) {
  if (x$V < 2 || is.null(x$phis)) {
    return(NULL)
  }
  pairs <- utils::combn(x$V, 2)
  phis <- x$phis[index, , drop = FALSE]
  out <- t(apply(phis, 2, function(p) {
    c(
      mean = mean(p), sd = stats::sd(p),
      stats::quantile(p, c(0.025, 0.5, 0.975), names = FALSE)
    )
  }))
  colnames(out) <- c("mean", "sd", "2.5%", "50%", "97.5%")
  rownames(out) <- sprintf("phi[%d,%d]", pairs[1, ], pairs[2, ])
  out
}

#' @title Print an mdir fit
#' @description Prints a short description of a fitted chain: the views
#' modelled, the sampler settings and whether a burn in has been applied. It
#' replaces the default printing of a list, which dumps every sampled array.
#' Use \code{summary()} for posterior summaries. \code{x} is still a plain list
#' underneath (\code{x$allocations}, \code{x$phis} and so on keep working).
#' @param x Output of \code{\link{callMDI}} or \code{\link{processMCMCChain}}
#' (one chain).
#' @param ... Unused.
#' @return \code{x}, invisibly.
#' @export
#' @examples
#' set.seed(1)
#' X <- list(matrix(rnorm(80, rep(c(0, 3), each = 20)), 40, 2))
#' fit <- callMDI(X, R = 100, thin = 5, types = "MVN", K = 4)
#' fit
print.mdir_fit <- function(x, ...) {
  cat(sprintf(
    "<mdir fit> %s, %d items\n",
    .mdirPlural(x$V, "view"), x$N
  ))
  views <- .mdirViewTable(x)
  for (v in seq_len(x$V)) {
    cat(sprintf(
      "  view %d: %s, %d measurements, K = %d, %s%s\n",
      v, views$type[v], views$P[v], views$K[v], views$supervision[v],
      if (x$missing_proportion[v] > 0) paste0(", ", views$missing[v], " missing") else ""
    ))
  }
  if (!is.null(x$tempering)) {
    cat(sprintf(
      "  parallel tempering: %d replicas, beta from %s to 1 (see ptDiagnostics())\n",
      length(x$tempering$betas), format(signif(min(x$tempering$betas), 3))
    ))
  }
  if (.mdirIsProcessed(x)) {
    cat(sprintf(
      "  R = %d, thin = %d; burn = %d applied, %d samples retained\n",
      x$R, x$thin, x$burn, nrow(x$mass)
    ))
  } else {
    cat(sprintf(
      "  R = %d, thin = %d; %d samples saved (burn in not applied)\n",
      x$R, x$thin, nrow(x$mass)
    ))
  }
  cat("  Run time: ", .mdirFormatTime(x$Time), "\n", sep = "")
  if (.mdirIsProcessed(x)) {
    cat("\nPoint estimates are in x$pred; summary() gives details.\n")
  } else {
    cat("\nUse summary() for posterior summaries, processMCMCChain() for point estimates.\n")
  }
  invisible(x)
}

#' @title Summarise an mdir fit
#' @description Posterior summaries for one chain, laid out as in
#' \pkg{mclust}'s \code{summary()}: the views modelled, the number of occupied
#' components in each, posterior summaries of the dataset association
#' parameters \eqn{\phi} (mean, standard deviation and quantiles, as in Stan),
#' mean log-likelihoods and, if \code{\link{processMCMCChain}} has been
#' applied, a clustering table for each view.
#' @param object Output of \code{\link{callMDI}} or
#' \code{\link{processMCMCChain}} (one chain).
#' @param burn Number of iterations to discard before summarising. Defaults to
#' half of \code{R}, as in \code{\link{assessConvergence}}. Ignored (the applied
#' burn in is used) if the chain has already been through
#' \code{\link{processMCMCChain}}.
#' @param ... Unused.
#' @return An object of class \code{summary.mdir_fit}.
#' @export
#' @examples
#' set.seed(1)
#' X <- lapply(1:2, function(v) matrix(rnorm(80, rep(c(0, 3), each = 20)), 40, 2))
#' fit <- callMDI(X, R = 200, thin = 5, types = c("MVN", "MVN"), K = c(4, 4))
#' summary(fit)
summary.mdir_fit <- function(object, burn = NULL, ...) {
  kept <- .mdirRetained(object, burn)
  index <- kept$index
  V <- object$V

  views <- .mdirViewTable(object)
  occupied <- lapply(seq_len(V), function(v) .mdirOccupied(object, v, index))
  views$occupied <- vapply(occupied, function(o) {
    r <- range(o)
    if (r[1] == r[2]) {
      sprintf("%d", r[1])
    } else {
      sprintf("%.0f (%d-%d)", stats::median(o), r[1], r[2])
    }
  }, character(1))

  clustering <- NULL
  if (.mdirIsProcessed(object)) {
    clustering <- lapply(object$pred, table)
  }

  structure(
    list(
      N = object$N, V = V, R = object$R, thin = object$thin,
      burn = kept$burn, burn_default = kept$burn_default,
      n_samples = length(index),
      processed = .mdirIsProcessed(object),
      views = views,
      phi = .mdirPhiTable(object, index),
      complete_likelihood = mean(object$complete_likelihood[index]),
      observed_likelihood = if (is.null(object$observed_likelihood)) {
        NA_real_
      } else {
        mean(object$observed_likelihood[index])
      },
      joint_likelihood = if (is.null(object$joint_likelihood)) {
        NA_real_
      } else {
        mean(object$joint_likelihood[index])
      },
      mass_acceptance_rate = object$mass_acceptance_rate,
      time = object$Time,
      clustering = clustering
    ),
    class = "summary.mdir_fit"
  )
}

#' @rdname summary.mdir_fit
#' @param x Output of \code{summary.mdir_fit}.
#' @param digits Significant digits to print.
#' @export
print.summary.mdir_fit <- function(x, digits = 3, ...) {
  cat(.mdirRule(), "\n", sep = "")
  cat("MDI model fitted by MCMC\n")
  cat(.mdirRule(), "\n\n")

  cat(sprintf(
    "%d items, %s; R = %d, thin = %d.\n",
    x$N, .mdirPlural(x$V, "view"), x$R, x$thin
  ))
  cat(sprintf(
    "Burn in = %d%s; %d samples used.\n\n",
    x$burn, if (x$burn_default) " (default, R / 2)" else "", x$n_samples
  ))

  v <- x$views
  tab <- data.frame(
    View = v$view, Type = v$type, P = v$P, K = v$K,
    Supervision = v$supervision, Missing = v$missing,
    "Occupied components" = v$occupied, check.names = FALSE
  )
  print(tab, row.names = FALSE)

  if (!is.null(x$phi)) {
    cat("\nDataset association (phi), posterior summary:\n")
    print(signif(x$phi, digits))
  }

  cat(sprintf(
    "\nMean log-likelihood: complete-data = %s, observed-data = %s",
    format(signif(x$complete_likelihood, 5)),
    if (is.na(x$observed_likelihood)) "not recorded" else format(signif(x$observed_likelihood, 5))
  ))
  if (!is.null(x$joint_likelihood) && !is.na(x$joint_likelihood)) {
    cat(sprintf(", joint = %s", format(signif(x$joint_likelihood, 5))))
  }
  cat("\n")
  cat("Run time: ", .mdirFormatTime(x$time), "\n", sep = "")

  if (x$processed) {
    for (i in seq_along(x$clustering)) {
      cat(sprintf("\nClustering table, view %d:\n", i))
      print(x$clustering[[i]])
    }
  } else {
    cat("\nNo point estimates yet: use processMCMCChain() to obtain clusterings.\n")
  }
  invisible(x)
}

#' @title Print a list of mdir fits
#' @description Prints a short description of a set of chains: the views
#' modelled, the sampler settings and, if computed, a one-paragraph verdict on
#' convergence. It replaces the default printing of a list, which dumps every
#' sampled array of every chain. The object is still a plain list of chains
#' underneath (\code{x[[i]]$phis} and so on keep working).
#' @param x Output of \code{\link{runMCMCChains}}, \code{\link{fitMDI}} or
#' \code{\link{processMCMCChains}}.
#' @param ... Unused.
#' @return \code{x}, invisibly.
#' @export
print.mdir_fit_list <- function(x, ...) {
  n_chains <- length(x)
  first <- x[[1]]
  cat(sprintf(
    "<mdir fit list> %s, %s, %d items\n",
    .mdirPlural(n_chains, "chain"), .mdirPlural(first$V, "view"), first$N
  ))
  views <- .mdirViewTable(first)
  for (v in seq_len(first$V)) {
    cat(sprintf(
      "  view %d: %s, %d measurements, K = %d, %s\n",
      v, views$type[v], views$P[v], views$K[v], views$supervision[v]
    ))
  }
  if (.mdirIsProcessed(first)) {
    cat(sprintf(
      "  R = %d, thin = %d per chain; burn = %d applied, %d samples retained per chain\n",
      first$R, first$thin, first$burn, nrow(first$mass)
    ))
  } else {
    cat(sprintf(
      "  R = %d, thin = %d per chain; %d samples saved per chain (burn in not applied)\n",
      first$R, first$thin, nrow(first$mass)
    ))
  }
  total_time <- sum(vapply(x, function(ch) as.numeric(ch$Time, units = "secs"), numeric(1)))
  cat("  Total run time: ", .mdirFormatTime(round(total_time, 1)), " secs\n", sep = "")

  convergence <- attr(x, "convergence")
  if (!is.null(convergence)) {
    cat("\n", .mdirWrap(format(convergence)), "\n", sep = "")
    cat("\nprint(attr(x, \"convergence\")) gives the full table.\n")
  } else {
    cat("\nConvergence has not been assessed: see assessConvergence().\n")
  }
  cat("Use summary() for per-chain detail, predictFromMultipleChains() to pool chains.\n")
  invisible(x)
}

#' @title Subset a list of mdir fits
#' @description Keeps the \code{mdir_fit_list} class so a subset of chains still
#' prints as a summary. The convergence diagnostics attached by
#' \code{\link{fitMDI}} refer to the full set of chains, so they are dropped;
#' recompute them with \code{\link{assessConvergence}}.
#' @param x Output of \code{\link{runMCMCChains}} or \code{\link{fitMDI}}.
#' @param i Chains to keep.
#' @return An \code{mdir_fit_list}.
#' @export
`[.mdir_fit_list` <- function(x, i) {
  structure(unclass(x)[i], class = class(x))
}

#' @title Summarise a list of mdir fits
#' @description A per-chain table (run time, mean log-likelihood, mean
#' \eqn{\phi} for each view pair) and, if attached, the convergence
#' diagnostics. Differences between chains in the mean log-likelihood are worth
#' looking at: a chain sitting well below the others may be stuck in a poorer
#' mode.
#' @param object Output of \code{\link{runMCMCChains}}, \code{\link{fitMDI}} or
#' \code{\link{processMCMCChains}}.
#' @param burn Number of iterations to discard before summarising; see
#' \code{\link{summary.mdir_fit}}.
#' @param ... Unused.
#' @return An object of class \code{summary.mdir_fit_list}.
#' @export
summary.mdir_fit_list <- function(object, burn = NULL, ...) {
  per_chain <- lapply(object, summary, burn = burn)
  structure(
    list(
      n_chains = length(object),
      per_chain = per_chain,
      convergence = attr(object, "convergence")
    ),
    class = "summary.mdir_fit_list"
  )
}

#' @rdname summary.mdir_fit_list
#' @param x Output of \code{summary.mdir_fit_list}.
#' @param digits Significant digits to print.
#' @export
print.summary.mdir_fit_list <- function(x, digits = 3, ...) {
  first <- x$per_chain[[1]]
  cat(.mdirRule(), "\n", sep = "")
  cat("MDI model fitted by MCMC: multiple chains\n")
  cat(.mdirRule(), "\n\n")
  cat(sprintf(
    "%s; %d items, %s; R = %d, thin = %d.\n",
    .mdirPlural(x$n_chains, "chain"), first$N, .mdirPlural(first$V, "view"),
    first$R, first$thin
  ))
  cat(sprintf(
    "Burn in = %d%s; %d samples used per chain.\n\n",
    first$burn, if (first$burn_default) " (default, R / 2)" else "", first$n_samples
  ))

  tab <- data.frame(
    Chain = seq_len(x$n_chains),
    "Time (s)" = vapply(x$per_chain, function(s) round(as.numeric(s$time, units = "secs"), 1), numeric(1)),
    "Mean complete log-lik" = vapply(x$per_chain, function(s) signif(s$complete_likelihood, 6), numeric(1)),
    check.names = FALSE
  )
  if (!is.null(first$phi)) {
    phi_means <- vapply(x$per_chain, function(s) s$phi[, "mean"], numeric(nrow(first$phi)))
    phi_means <- matrix(phi_means, nrow = nrow(first$phi))
    for (i in seq_len(nrow(first$phi))) {
      tab[[paste0("Mean ", rownames(first$phi)[i])]] <- signif(phi_means[i, ], digits)
    }
  }
  print(tab, row.names = FALSE)

  if (!is.null(x$convergence)) {
    cat("\n")
    print(x$convergence, digits = digits)
  }
  invisible(x)
}
