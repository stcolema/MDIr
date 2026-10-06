# Simulated multi-view data sets for the joint-allocation experiments.

ari <- function(a, b) {
  t <- table(a, b)
  n <- sum(t)
  s <- function(x) sum(choose(x, 2))
  expected <- s(rowSums(t)) * s(colSums(t)) / choose(n, 2)
  maximum <- 0.5 * (s(rowSums(t)) + s(colSums(t)))
  if (maximum == expected) return(NA_real_)
  (s(t) - expected) / (maximum - expected)
}

# Three views: views 1 and 2 share three clusters except for 18 re-drawn items, view 3 is
# unrelated (the data of the multimodal vignette)
make_small <- function() {
  set.seed(2026)
  N <- 120
  z1 <- rep(1:3, each = N / 3)
  z2 <- z1
  broken <- sort(sample(N, 18))
  z2[broken] <- sample(1:3, 18, replace = TRUE)
  z3 <- sample(1:2, N, replace = TRUE)
  mu1 <- rbind(c(-3, 0), c(0, 3), c(3, 0))
  mu2 <- rbind(c(-2, -2, 0), c(2, -2, 2), c(0, 3, -2))
  p3 <- rbind(c(0.8, 0.2, 0.2, 0.7), c(0.2, 0.8, 0.7, 0.2))
  X <- list(
    mu1[z1, ] + matrix(rnorm(2 * N), N, 2),
    mu2[z2, ] + matrix(rnorm(3 * N), N, 3),
    matrix(rbinom(4 * N, 1, p3[z3, ]), N, 4)
  )
  for (v in 1:3) rownames(X[[v]]) <- sprintf("item%03d", 1:N)
  list(X = X, types = c("MVN", "G", "C"), K = c(6, 6, 6), truth = list(z1, z2, z3), broken = broken,
       prior = mdiPrior(mass_shape = 2, mass_rate = 0.5))
}

# Many clusters: `n_clusters` true clusters shared by the first two views (a fraction of
# the items re-drawn in view 2), a third view sharing half of them. Features are diagonal
# Gaussian around cluster means drawn at random.
make_large <- function(n_clusters = 30, per = 10, P = 4, redraw = 0.1, K_model = 45, seed = 2027) {
  set.seed(seed)
  N <- n_clusters * per
  z1 <- rep(seq_len(n_clusters), each = per)
  z2 <- z1
  broken <- sort(sample(N, round(redraw * N)))
  z2[broken] <- sample(n_clusters, length(broken), replace = TRUE)
  z3 <- ifelse(seq_len(N) %in% sample(N, N / 2), z1, sample(n_clusters, N, replace = TRUE))
  centres <- function(P) matrix(rnorm(n_clusters * P, 0, 2.5), n_clusters, P)
  m1 <- centres(P); m2 <- centres(P); m3 <- centres(P)
  X <- list(m1[z1, ] + matrix(rnorm(N * P), N, P), m2[z2, ] + matrix(rnorm(N * P), N, P),
            m3[z3, ] + matrix(rnorm(N * P), N, P))
  for (v in 1:3) rownames(X[[v]]) <- sprintf("item%04d", 1:N)
  list(X = X, types = rep("G", 3), K = rep(K_model, 3), truth = list(z1, z2, z3), broken = broken,
       prior = mdiPrior(mass_shape = 2, mass_rate = 0.5))
}

key_quantities <- function(nm, V) {
  keys <- c("joint_likelihood", paste0("occupied_components[", seq_len(V), "]"),
            "phi[1,2]", "agreement[1,2]")
  intersect(nm, keys)
}
