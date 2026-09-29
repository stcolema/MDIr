# Internal: dimension d of the component-specific parameter of a density, as in
# Rousseau and Mengersen (2011). NA where the result is not applicable (GP).
.componentDimension <- function(x, type) {
  P <- ncol(x)
  switch(type,
    "G" = 2 * P,
    "MVN" = P + P * (P + 1) / 2,
    "TAGM" = P + P * (P + 1) / 2,
    "C" = sum(apply(x, 2, function(col) max(col, na.rm = TRUE)))   # sum of (categories - 1)
    ,
    NA_real_
  )
}

# Internal: messages for views where the prior on the weights sits on the side of
# d / 2 that duplicates rather than empties superfluous components
.checkSparsity <- function(X, types, K, prior) {
  msgs <- character(0)
  if (isTRUE(getOption("mdir.quiet"))) {
    return(msgs)
  }
  a_median <- stats::qgamma(0.5, shape = prior[["mass_shape"]], rate = prior[["mass_rate"]])
  for (v in seq_along(X)) {
    d <- .componentDimension(X[[v]], types[v])
    if (is.na(d) || d <= 0) next
    e0 <- a_median / K[v]
    if (e0 > d / 2) {
      msgs <- c(msgs, sprintf(
        paste0("View %d: the prior median of mass / K is %.2f, above d / 2 = %.2f. ",
               "Superfluous components are then expected to duplicate occupied ones ",
               "rather than empty out (Rousseau and Mengersen, 2011; asymptotic guidance). ",
               "See ?mdiPrior."), v, e0, d / 2))
    }
  }
  msgs
}
