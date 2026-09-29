# mdir 0.10.0

This release corrects several errors in the sampler, removes the parallel
execution that made the package unbuildable on CRAN, and adds tools for a
Bayesian workflow. Results from earlier versions should not be reused: the
posterior of the weights, `phi` and mass parameters was wrong (see
"Corrections").

## Missing data

* `NA`/`NaN` entries in any view are treated as missing at random and imputed
  within the sampler (partially collapsed Gibbs: the allocation uses the
  marginal likelihood of the observed entries; the missing entries are then drawn
  from their exact conditional given the new allocation and before the
  parameters are updated). Imputation is exact for multivariate normal,
  diagonal Gaussian, GP, categorical and, for outliers, the multivariate t
  (conditional t with `df + p_obs` degrees of freedom).
* Previously imputed values were drawn under the old allocation, which
  biased clusters together (in one simulation with 40% of entries missing the
  separation of two clusters, truly 2.5, was recovered as 1.6).
* Data-driven prior hyperparameters use the observed entries only
  (pairwise-complete covariance, repaired to be positive definite).
* Input checks: infinite values, columns with no observed value and malformed
  categorical coding are rejected with informative messages.
* `save_imputed = TRUE` records the imputed values at each saved iteration.

## Bayesian workflow

* `simulatePriorPredictive()`, `simulatePosteriorPredictive()`,
  `predictiveCheck()` and `plotPredictiveCheck()`. The simulators live in C++
  next to the priors and samplers so they cannot drift from the fitted model.
* `rankNormalizedRhat()` and `assessConvergence()`: rank-normalised, folded,
  split R-hat and bulk/tail ESS (Vehtari et al., 2021; agrees with
  `posterior::rhat()`/`ess_bulk()`/`ess_tail()` to under 0.2%), applied to
  label-switching-invariant quantities.
* `mdiPrior()` exposes the MDI-level priors, with guidance on choosing them.
* The sampler now records the observed-data log-likelihood, outlier weights and
  (by default) the component parameters.

## Corrections

* Component weights had posterior shape `mass / K + N_k + 1`; the `+ 1` is
  removed (it stopped empty components from being pruned).
* The mass update used rate 1 for the weight prior instead of the rate used to
  draw the weights and lacked the Jacobian of its log-scale proposal.
* `logChoose()` dropped a term for every `k >= 1`, mis-weighting the mixture
  from which the `phi` posterior is sampled.
* The `phi` rate used the wrong pairs of views when there were three or more.
* Label swaps now exchange the component parameters as well as labels and
  weights, and their acceptance ratio includes the change in the normalising
  constant.
* Diagonal Gaussian: variances were used as standard deviations when drawing the
  component means and imputing values; the empirical scale used a test that
  selected incomplete rows.
* GP: the squared-exponential kernel had a sign error (introduced with the
  missing-data work); the hyperparameter Metropolis target was the posterior of
  `mu` rather than the GP prior density of `mu`; empty components drew `mu`
  before their hyperparameters.
* Outliers: a Beta hyperparameter was overwritten with a random number during
  sampling; the outlier likelihood was computed before missing patterns were
  known; labels and outlier status are now drawn jointly.
* The observed-data likelihood summed log terms over components (now a
  log-sum-exp); label sampling could read one past the last component;
  `invWishartLogLikelihood()` had the wrong sign on the scale determinant.
* `predictFromMultipleChains()`/`processMCMCChain()` no longer mis-handle a burn
  in shorter than `thin` or longer than the run.
* `callMDIWritingToFile()`/`readInSavedSamples()` used inconsistent layouts and
  undefined variables.

## Performance and build

* The normalising constant and the rates of the weight and `phi` conditionals
  are computed exactly by a set-partition recursion in `O(3^L + 2^L K L)`
  instead of enumerating `K^L` assignments (2.7 s to 0.02 s for 10 iterations
  with `L = 6`, `K = 10`). In benchmarks with complete data, serial runs were 2 to 3.5 times
  faster than the previous four-threaded runs.
* All `std::execution::par` loops, RcppParallel, TBB and OpenMP flags are
  removed. They shared the random number generator across threads, so results
  were neither reproducible with `set.seed()` nor CRAN-compatible. C++17 is no
  longer forced. `R CMD check` reports no package-related notes or warnings.
* The multivariate normal caches the inverse and log-determinant of each
  covariance rather than recomputing them for every item and component.

## Tests

* A `testthat` suite checks the exact normalising constant and rates against
  brute force, that the sampler recovers its prior under a constant likelihood,
  closed-form conjugate posteriors, the GP hyperparameter posterior against a
  grid, imputation exactness, and the diagnostics against `posterior`.
