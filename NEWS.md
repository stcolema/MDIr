# mdir 0.10.2

## CRAN readiness

* The package was archived on CRAN (2023-05-31) for installation failures on
  macOS and Fedora/clang and a GNU make `SystemRequirements` note. Both stemmed
  from the parallel build (C++17 `<execution>`, RcppParallel/TBB, `$(shell ...)`
  in `Makevars`), which is removed. See `cran-comments.md`.
* Vignettes: the Quarto tutorial, which needs Bioconductor data, moved to
  `articles/` (not built by CRAN). A new knitr vignette, "A Bayesian workflow
  with mdir", runs on simulated data.
* Examples: several were broken (undefined objects, a missing function, a wrong
  column mapping) or very slow (500 chains); all 42 now run, in under five
  seconds each.
* `DESCRIPTION`: institutional proxy DOI replaced, typo fixed. `aes_string()`
  (deprecated) and the native pipe (which needs R >= 4.1) are no longer used.
* The package help page is now `?mdir-package`.

# mdir 0.10.1

## Partial pooling and guarded Gaussian process priors

* **Variance scale pooled across components** (`"MVN"`, `"TAGM"`, `"G"`). The
  inverse-Wishart / inverse-gamma scale, previously a fixed data-driven value, is
  now diagonal with a Gamma hyperprior centred on that value, updated from the
  occupied components (the exact conditional with empty components integrated
  out). `densityPrior(scale_pool_shape = 0)` restores the fixed scale. The
  pooled scales are saved as `pooled_hyperparameters` and monitored by
  `assessConvergence()`.
* **Gaussian process views** (`"GP"`, `"TAGPM"`) have new priors:
  * The length scale is now an actual length scale, `lambda` in
    `a * exp(-d^2 / (2 * lambda^2))` (previously `lambda^2` was stored under the
    name `length`, so saved values change meaning). It has a hard floor at
    `gp_min_length` (default one measurement spacing) and an inverse-gamma prior
    calibrated to put 1% of its mass below the floor and 1% above the extent of
    the grid. On data with no structure, an unguarded prior put 94% of the
    length-scale draws below one grid unit (median 0.09), where the kernel is
    diagonal and the mean function is white noise; with the guard none fall below
    the floor.
  * Amplitude and noise variance have a log-normal population shared by the
    components, centred on the average data variance, with a half-normal prior on
    its sd. GP data no longer need to be standardised, and the mean function is
    centred on the column means rather than zero.
  * An interweaved non-centred amplitude move and a mixture of random-walk scales
    improve mixing of the amplitude (effective sample size 7 to about 400 in a
    test where the data barely constrain it).
  * At least three measurements per item are required.
* `densityPrior()` collects these options; `mdiPrior()` documentation now states
  the Rousseau and Mengersen (2011) condition (`mass / K` below `d / 2` empties
  superfluous components) with its assumptions (asymptotic, fixed weight prior,
  regular kernel), and `callMDI()` reports a message when the prior median of
  `mass / K` is above `d / 2` (`options(mdir.quiet = TRUE)` silences it).

## Tests

* Marginal posteriors of the pooled scales (one occupied component, with and
  without an empty second component) and of the GP hyperparameters are compared
  with exact numerical posteriors; the GP population update is compared with its
  analytic conditional.

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
