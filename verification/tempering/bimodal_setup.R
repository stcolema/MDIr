# Four well-separated clusters, fitted with K = 3 components: one pair must merge
# and which pair is a mode of the partition posterior that a Gibbs chain rarely leaves.
make_data <- function(seed = 1, n_per = 20, sep = 4, sd = 1) {
  set.seed(seed)
  centres <- rbind(c(-sep, -sep), c(-sep, sep), c(sep, -sep), c(sep, sep))
  truth <- rep(1:4, each = n_per)
  X <- centres[truth, ] + matrix(rnorm(length(truth) * 2, 0, sd), ncol = 2)
  rownames(X) <- seq_len(nrow(X))
  list(X = X, truth = truth)
}
# canonical set partition of the four true groups induced by one draw's labels
group_partition <- function(labels, truth) {
  lab <- vapply(1:4, function(g) as.integer(names(which.max(table(labels[truth == g])))), integer(1))
  canon <- match(lab, unique(lab))
  paste(vapply(split(1:4, canon), paste, character(1), collapse = ""), collapse = "|")
}
