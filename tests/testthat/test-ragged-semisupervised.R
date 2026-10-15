# Views with different numbers of components, some semi-supervised and some
# unsupervised. The label swap exchanges the labels, weights and parameters of
# two components of one view; if it also moved the weights of the other views,
# real weights were exchanged with the unused (zero) slots of views with fewer
# components, which left non-zero values in the unused slots, put zeros in real
# slots and, in some runs, stopped the sampler.

ragged_run <- function(K, types, supervised, R = 4000, N = 36, seed = 1) {
  set.seed(seed)
  L <- length(K)
  truth <- sample(3, N, TRUE)
  X <- lapply(seq_len(L), function(l) {
    m <- if (types[l] == "C") matrix(rbinom(N * 4, 1, c(.2, .5, .8)[truth]), N, 4) else matrix(rnorm(N * 2, 1.5 * truth), N, 2)
    rownames(m) <- seq_len(N); m
  })
  fixed <- matrix(0, N, L); init <- matrix(1, N, L)
  for (l in supervised) {
    obs <- sort(sample(N, 12)); K_fix <- min(K[l], 2)
    fixed[obs, l] <- 1; init[obs, l] <- ((truth[obs] - 1) %% K_fix) + 1
    init[obs[seq_len(K_fix)], l] <- seq_len(K_fix)       # every observed class is present
  }
  out <- callMDI(X, R = R, thin = 10, types = types, K = K, initial_labels = init, fixed = fixed)
  list(out = out, fixed = fixed, init = init, K = K)
}

check_ragged_run <- function(run, label) {
  out <- run$out; K <- run$K; L <- length(K); KM <- max(K)
  al <- out$allocations + 1L; w <- out$weights
  for (l in seq_len(L)) {
    expect_true(all(al[, , l] >= 1 & al[, , l] <= K[l]), label = paste(label, "labels in range, view", l))
    expect_true(all(is.finite(w[, seq_len(K[l]), l]) & w[, seq_len(K[l]), l] > 0), label = paste(label, "real weights positive, view", l))
    if (K[l] < KM) expect_true(all(w[, (K[l] + 1):KM, l] == 0), label = paste(label, "unused weight slots are zero, view", l))
  }
  expect_true(all(sapply(seq_len(dim(al)[1]), function(s) all(al[s, , ][run$fixed == 1] == run$init[run$fixed == 1]))),
              label = paste(label, "observed labels are kept"))
  expect_true(all(is.finite(out$phis) & out$phis > 0) && all(is.finite(out$mass) & out$mass > 0), label = paste(label, "phi and mass positive"))
}

test_that("five views with different numbers of components and mixed supervision stay consistent", {
  for (seed in 1:3) {
    run <- ragged_run(K = c(2L, 4L, 3L, 5L, 3L), types = c("MVN", "G", "C", "G", "MVN"), supervised = c(1L, 3L), seed = seed)
    check_ragged_run(run, paste("seed", seed))
  }
})

test_that("a view with a single component survives label swaps in the other views", {
  # A swap in another view must leave the single weight of this view in its slot
  for (seed in 1:3) {
    run <- ragged_run(K = c(1L, 3L, 2L, 4L, 2L), types = c("G", "G", "C", "MVN", "G"), supervised = c(2L, 4L), seed = seed)
    check_ragged_run(run, paste("seed", seed))
  }
})

test_that("unsupervised views with different numbers of components stay consistent", {
  run <- ragged_run(K = c(5L, 2L, 4L, 3L), types = c("G", "C", "G", "MVN"), supervised = integer(0), seed = 7)
  check_ragged_run(run, "unsupervised")
})
