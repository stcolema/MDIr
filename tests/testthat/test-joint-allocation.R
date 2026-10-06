# The joint draw of the labels of one item across views, against enumeration of
# its exact conditional, and the move inside the sampler.

# Exact joint conditional of the labels of the views in `block` (1-based), the other
# views held at `current` (1-based), by enumeration
ja_exact <- function(log_g, w, K, phi, block, current) {
  grid <- as.matrix(expand.grid(lapply(K[block], seq_len)))
  lp <- numeric(nrow(grid))
  for (r in seq_len(nrow(grid))) {
    cc <- current
    cc[block] <- grid[r, ]
    v <- 0
    for (b in block) v <- v + log(w[cc[b], b]) + log_g[cc[b], b]
    for (l in seq_len(length(K) - 1)) for (m in (l + 1):length(K)) {
      if (l %in% block || m %in% block) v <- v + log1p(phi[l, m] * (cc[l] == cc[m]))
    }
    lp[r] <- v
  }
  p <- exp(lp - max(lp))
  list(grid = grid, prob = p / sum(p))
}

ja_state <- function(K, phi_scale = 8, seed = 1, excluded = 0) {
  set.seed(seed)
  L <- length(K)
  Kmax <- max(K)
  w <- matrix(0, Kmax, L)
  log_g <- matrix(-Inf, Kmax, L)
  for (l in seq_len(L)) {
    w[seq_len(K[l]), l] <- rgamma(K[l], 1.5, 1)
    log_g[seq_len(K[l]), l] <- rnorm(K[l], 0, 2)
  }
  if (excluded > 0) log_g[sample(K[1], excluded), 1] <- -Inf     # components the item cannot take
  phi <- matrix(0, L, L)
  phi[upper.tri(phi)] <- rgamma(choose(L, 2), 2, 2 / phi_scale)
  phi <- phi + t(phi)
  list(w = w, log_g = log_g, K = K, phi = phi)
}

# Frequencies of n draws against exact probabilities: no cell more than 4.5 standard errors
# away, and the chi-square statistic of the cells with expected count >= 5 not extreme
ja_check <- function(st, block, current, n = 60000, label = "") {
  ex <- ja_exact(st$log_g, st$w, st$K, st$phi, block, current)
  d <- mdir:::jointAllocationDrawCpp(st$log_g, st$w, as.integer(st$K), st$phi,
                                     as.integer(block - 1), as.integer(current - 1), n)
  code <- function(m) as.numeric((m) %*% cumprod(c(1, st$K[block][-length(block)])))
  counts <- tabulate(code(d) + 1, nrow(ex$grid))
  e <- n * ex$prob
  z <- (counts - e) / sqrt(pmax(e * (1 - ex$prob), 1e-12))
  expect_lt(max(abs(z[e > 1])), 4.5, label = paste("max |z|", label))
  big <- e >= 5
  chi <- sum((counts[big] - e[big])^2 / e[big])
  expect_gt(pchisq(chi, sum(big) - 1, lower.tail = FALSE), 1e-4, label = paste("chi-square", label))
  # no draw outside the support
  expect_true(all(counts[ex$prob == 0] == 0), label = paste("support", label))
  invisible(list(counts = counts, exact = ex))
}

test_that("the joint draw matches enumeration of its exact conditional", {
  # ragged numbers of components, strong association, an excluded component
  st <- ja_state(c(3, 4, 2), phi_scale = 15, seed = 11, excluded = 1)
  ja_check(st, 1:3, c(1, 1, 1), label = "all views, K = (3, 4, 2)")
  # a block of two views, the third held at its current label
  ja_check(st, c(1, 3), c(1, 2, 2), label = "block (1, 3)")
  ja_check(st, c(2, 3), c(3, 1, 1), label = "block (2, 3)")
  # four views, large phi: the sequential conditioning must keep the coupling
  st4 <- ja_state(c(3, 3, 4, 3), phi_scale = 40, seed = 12)
  ja_check(st4, 1:4, c(1, 1, 1, 1), label = "four views, large phi")
  ja_check(st4, c(1, 2, 4), c(1, 1, 3, 1), label = "three of four views")
})

test_that("the joint draw reduces to the prior of the labels when the likelihood is flat", {
  # the case used for the prior draws beyond the enumeration limit
  st <- ja_state(c(4, 3, 4), phi_scale = 10, seed = 13)
  st$log_g[is.finite(st$log_g)] <- 0
  ja_check(st, 1:3, c(1, 1, 1), label = "flat likelihood")
})

test_that("the joint draw has the marginals of mdiClassProbabilities", {
  st <- ja_state(c(3, 4, 3), phi_scale = 20, seed = 14)
  probs <- mdir:::mdiClassProbabilitiesCpp(st$log_g, st$w, as.integer(st$K), st$phi)
  d <- mdir:::jointAllocationDrawCpp(st$log_g, st$w, as.integer(st$K), st$phi, 0:2, c(0L, 0L, 0L), 80000)
  for (l in 1:3) {
    freq <- tabulate(d[, l] + 1, st$K[l]) / nrow(d)
    se <- sqrt(probs[seq_len(st$K[l]), l] * (1 - probs[seq_len(st$K[l]), l]) / nrow(d))
    expect_true(all(abs(freq - probs[seq_len(st$K[l]), l]) < 4.5 * pmax(se, 1e-4)), info = paste("view", l))
  }
})

test_that("the joint draw copes with very large and very small scores", {
  st <- ja_state(c(3, 3), phi_scale = 50, seed = 15)
  st$log_g <- st$log_g * 300                       # likelihoods spanning hundreds of log units
  ja_check(st, 1:2, c(1, 1), n = 20000, label = "extreme scores")
})

test_that("prior labels are drawn beyond the enumeration limit", {
  skip_on_cran()
  X <- lapply(1:4, function(v) { m <- matrix(rnorm(20), 10, 2); rownames(m) <- 1:10; m })
  # 50^4 combinations: enumerating them was refused
  sims <- simulatePriorPredictive(X, rep("G", 4), K = rep(50, 4), n_datasets = 3)
  expect_equal(dim(sims$parameters[[1]]$labels), c(10, 4))
  expect_true(all(sims$parameters[[1]]$labels >= 0 & sims$parameters[[1]]$labels < 50))
})

raw_joint <- function(X, types, K, labels, fixed, R = 200, joint = 2L, split_merge = 0L, betas = numeric(0)) {
  V <- length(X)
  mdir:::runMDI(
    R, 1L, X, as.integer(K), mdir:::translateTypes(types), mdir:::setupOutlierComponents(types),
    labels, fixed, rep(list(0), V), FALSE, FALSE, as.numeric(mdiPrior()), as.numeric(densityPrior()),
    rep(0L, V), FALSE, TRUE, betas, 0L, 1L, as.integer(split_merge), as.integer(joint)
  )
}

ja_data <- function(N = 45, V = 3, seed = 7) {
  set.seed(seed)
  m <- rep(c(-4, 0, 4), each = N / 3)
  X <- lapply(seq_len(V), function(v) {
    x <- matrix(rnorm(N * 2, m), N, 2)
    rownames(x) <- seq_len(N)
    x
  })
  list(X = X, truth = rep(1:3, each = N / 3))
}

test_that("observed labels are never changed by the joint update", {
  d <- ja_data()
  N <- 45
  fixed <- matrix(0, N, 3)
  fixed[c(1:4, 31:34), 1] <- 1
  fixed[c(1:3, 16:18), 2] <- 1
  labels <- matrix(1L, N, 3)
  labels[1:4, 1] <- 0L; labels[31:34, 1] <- 2L        # classes {0, 2}: a gap
  labels[1:3, 2] <- 4L; labels[16:18, 2] <- 1L
  for (tp in list(c("G", "MVN", "G"), c("TAGM", "G", "MVN"), c("MVN", "TAGM", "TAGM"))) {
    for (jb in c(2L, 3L)) {
      out <- raw_joint(d$X, tp, c(5, 5, 5), labels, fixed, joint = jb)
      for (v in 1:2) {
        idx <- which(fixed[, v] == 1)
        expect_true(all(apply(out$allocations[, idx, v], 1, function(r) all(r == labels[idx, v]))),
                    info = paste(paste(tp, collapse = "/"), "block", jb, "view", v))
      }
      expect_true(all(out$outliers[, which(fixed[, 1] == 1), 1] == 0))
    }
  }
})

test_that("the joint update keeps the recorded likelihood and the labels consistent", {
  d <- ja_data()
  N <- 45
  out <- raw_joint(d$X, c("G", "MVN", "G"), c(4, 4, 4), matrix(0L, N, 3), matrix(0, N, 3), R = 60, joint = 3L)
  # The complete-data log-likelihood is recomputed after the joint redraw: it equals the
  # sum over items of the log density of each item at its (current) component, to rounding.
  # Compare with the likelihood of the joint labels to be finite and move with the labels.
  expect_true(all(is.finite(out$complete_likelihood)))
  expect_gt(length(unique(round(out$complete_likelihood, 6))), 30)
  expect_true(all(out$allocations >= 0 & out$allocations < 4))
  # counts of each view add up
  for (s in c(1, 30, 61)) for (v in 1:3) expect_equal(sum(out$N_k[1:4, v, s]), N)
})

test_that("a block smaller than the number of views is accepted", {
  d <- ja_data(V = 4)
  N <- 45
  out <- raw_joint(d$X, rep("G", 4), rep(4, 4), matrix(0L, N, 4), matrix(0, N, 4), R = 60, joint = 2L)
  expect_true(all(out$allocations >= 0 & out$allocations < 4))
  expect_error(raw_joint(d$X[1], "G", 4, matrix(0L, N, 1), matrix(0, N, 1), R = 5, joint = 2L), "two views")
  expect_error(callMDI(d$X, R = 10, thin = 5, types = rep("G", 4), K = rep(4, 4), joint_allocation = 1L),
               "at least 2")
  expect_error(callMDI(d$X, R = 10, thin = 5, types = rep("G", 4), K = rep(4, 4), joint_allocation = 13L),
               "at most 12")
})

test_that("the sampler with joint allocation recovers the prior under a constant likelihood", {
  skip_on_cran()
  set.seed(5)
  K <- c(3L, 3L); N <- 10; L <- 2
  X <- lapply(1:2, function(l) { m <- matrix(0, N, 1); rownames(m) <- 1:N; m })
  stats_of <- function(cc, w) {
    rho <- sweep(w, 2, colSums(w), "/")
    top <- apply(w, 2, which.max)
    c(own1 = mean(rho[cbind(cc[, 1], 1)]), own2 = mean(rho[cbind(cc[, 2], 2)]),
      top1 = mean(cc[, 1] == top[1]), top2 = mean(cc[, 2] == top[2]),
      agree = mean(cc[, 1] == cc[, 2]))
  }
  forward <- function() {
    mass <- rgamma(L, 2, 0.1); phi <- rgamma(1, 2, 0.2)
    w <- sapply(seq_len(L), function(l) rgamma(K[l], mass[l] / K[l], 2))
    grid <- as.matrix(expand.grid(lapply(K, seq_len)))
    lp <- rowSums(sapply(seq_len(L), function(l) log(w[grid[, l], l]))) + log1p(phi * (grid[, 1] == grid[, 2]))
    p <- exp(lp - max(lp)); p <- p / sum(p)
    stats_of(grid[sample.int(nrow(grid), N, TRUE, p), , drop = FALSE], w)
  }
  ref <- t(replicate(20000, forward()))
  out <- callMDI(X, R = 200000, thin = 10, types = c("C", "C"), K = K, joint_allocation = 2L,
                 check_prior = FALSE, save_parameters = FALSE)
  keep <- -(1:50)
  al <- out$allocations[keep, , , drop = FALSE] + 1L
  wt <- out$weights[keep, , , drop = FALSE]
  mc <- t(vapply(seq_len(dim(al)[1]), function(s) stats_of(al[s, , ], wt[s, , ]), numeric(5)))
  batch_se <- function(x, nb = 40) sd(vapply(split(x, cut(seq_along(x), nb, labels = FALSE)), mean, numeric(1))) / sqrt(nb)
  for (j in colnames(ref)) {
    se <- sqrt(batch_se(mc[, j])^2 + var(ref[, j]) / nrow(ref))
    expect_lt(abs(mean(mc[, j]) - mean(ref[, j])), 5 * se, label = paste("prior mean of", j))
  }
})
