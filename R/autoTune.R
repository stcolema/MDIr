# Internal: validate the arguments that switch on and bound the adaptation of the
# proposal windows. Returns the number of sweeps over which they are adapted.
.mdirResolveAutoTune <- function(auto_tune, n_burn, n_iter) {
  if (!is.logical(auto_tune) || length(auto_tune) != 1 || is.na(auto_tune)) {
    stop("`auto_tune` must be TRUE or FALSE.", call. = FALSE)
  }
  if (is.null(n_burn)) {
    n_burn <- floor(n_iter / 2)
  }
  if (!is.numeric(n_burn) || length(n_burn) != 1 || is.na(n_burn) || n_burn < 0) {
    stop("`n_burn` must be NULL or a single non-negative number.", call. = FALSE)
  }
  as.integer(min(floor(n_burn), n_iter))
}
