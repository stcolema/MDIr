# mdir (development version)

## New: weighted ensembles, and a map of what is guaranteed

* **`smcMDI()`** runs many particles from exact draws of the prior to the posterior along the
  likelihood-tempered path and weights them (annealed importance sampling / sequential Monte
  Carlo, with a conditional-ESS adaptive or a fixed schedule, systematic or multinomial
  resampling). **`smcReplicates()`** fixes the schedule from a pilot run and runs independent
  ensembles; **`combineSMC()`** pools them by their evidence estimates.
* Guarantees, each with its status in `verification/tempering/theory.md`: with a fixed schedule
  the weighted sum is exactly unbiased for the unnormalised posterior and the evidence estimate is
  unbiased, for any number of particles and without any mixing assumption (Lean 4 for annealed
  importance sampling; exact rational enumeration for resampling schemes; Del Moral, 2004); the
  ratio estimators are consistent as the number of independent runs or particles grows. An
  adaptive schedule is consistent but exactly biased (shown by enumeration). Starting above the
  prior (`beta_start`) is conditional on an exact start, with a proved bound on the bias in
  terms of the start's total-variation error (Lean 4). No bound on the variance is available.
* In the four-cluster example the variance is the problem: with a sharp transition on the
  path from prior to posterior, 20,000 particles had an effective sample size of 1.6, and with
  resampling the effective sample size looked healthy while the mode masses were wrong. Pooling
  independent plain chains converges to the basin weights (total variation distance 0.31 from the
  exact masses with 160 chains).
* Helpers: `smcDiagnostics()`, `smcWeightDiagnostic()` (Pareto k-hat of Vehtari et al., 2024, for
  unresampled runs), `smcWeights()`, `weightedPSM()`, `weightedConsensus()`, `smcPosterior()`,
  `resampleSMC()`, `smcAsChain()` (feeds `processMCMCChain()`), `smcSE()` (jackknife over
  independent runs), `compareRuns()` (agreement of independent runs, for tempered chains too).
  Diagnostics are labelled as having no guarantee.
* Supported: complete data, `G`/`MVN`/`C` views, no outlier component, unsupervised views.

## New: parallel tempering

* **`betas`** (in `callMDI()`, `runMCMCChains()` and `fitMDI()`) runs one replica of the
  sampler at each inverse temperature of a ladder ending at 1 and exchanges neighbouring
  replicas. The replica at `beta` targets `L^beta P`, with `L` the likelihood of the data
  given the labels and component parameters and `P` the rest of the model (priors, the
  coupling of the views through the weights and phis, the strategic latent variable); the
  recorded draws are those of the replica at `beta = 1`. Exchanges are accepted with
  probability `min(1, exp((beta_i - beta_j)(l_j - l_i)))`, `l` the data log-likelihood, and use the
  deterministic even-odd schedule of Syed et al. (2022) by default (`swap_scheme = "deo"`;
  `"seo"` is its reversible, stochastic-parity counterpart). `betas = 1` (the default)
  is the previous sampler, draw for draw (checked on seeded chains of every density,
  with missing data and outliers).
* `ptLadder()` builds a ladder, `ptDiagnostics()` reports exchange rates, the
  estimated communication barrier and round trips, `tuneLadder()` and `adaptLadder()`
  equalise neighbouring rejection rates (Syed et al., 2022, Algorithm 2).
* Supported: complete data, `"G"`, `"MVN"` and `"C"` views, no outlier component; other
  models stop with an error (the tempered conditionals of the missing-data, t-outlier and
  Gaussian process updates are not implemented).
* What was checked (`verification/tempering/`, `tests/testthat/test-tempering.R`): the
  tempered conjugate updates against the closed forms, symbolically (SymPy) and by
  simulation; the exchange acceptance, detailed balance, invariance of every exchange round
  and the cold marginal in Lean 4 with Mathlib (finite state spaces, no `sorry`); the same
  in exact rational arithmetic for a small chain, including that the even-odd scheme is
  invariant but not reversible; and the sampler's label posterior against exact references
  for tiny categorical models (one and two views) and Gaussian models, with negative
  controls that must fail.
* What parallel tempering does not guarantee: finite-time convergence. The
  ergodicity of the product chain and the bounds on mixing are not proved here. In a
  four-cluster example the model's tempering path has a sharp transition, hot replicas
  made 0 to 3 round trips in 20,000 rounds although the exchange rates predicted about
  900, and the cold chain nevertheless recovered the exact mode weights when 12 runs were
  pooled (single runs were still noisy); see the README.

## New: likelihood of the whole model, and prediction

* **`joint_likelihood`** is recorded at every saved draw: the log-likelihood of the
  data under the MDI model with every item's component assignments in all views
  summed out, `log Z(w g) - log Z(w)` per item, where `g` holds the per-view
  component likelihoods (outliers marginalised). It is exact, reuses the
  recursion that gives `Z`, and is checked against enumeration of the joint
  assignments, including missing data and ragged `K`. The existing
  `observed_likelihood` sums over each view separately with that view's own
  normalised weights, which ignores the coupling between views; it equals the
  joint likelihood only when every `phi` is zero (one view, or `phi = 0`), and is
  left unchanged. In a semi-supervised view the observed labels are treated as
  data. The joint likelihood is also monitored by `assessConvergence()` and shown
  by `summary()`. Recording it added no measurable run time (a 3-view, 500-item
  chain took 0.44 s for 200 sweeps recorded every sweep, before and after).
* **`pointwiseLogLik()`** returns the log-likelihood of each item at each draw
  (`save_pointwise = TRUE`), as the draws x items matrix that `loo::waic()` and
  `loo::loo()` take, with the chain of each row.
* **`predictMDI()`** evaluates new items at every saved draw without refitting:
  the log posterior predictive density of each, the probability that it shares a
  component with each fitted item, and, for semi-supervised views, class
  probabilities. Items may have missing entries or whole views missing. The
  calculation is exact given a draw and agrees with enumeration. Class
  probabilities are not returned for unsupervised views, whose labels are not
  identified across draws.

## New: parallel chains

* `runMCMCChains()` and `fitMDI()` gain **`n_cores`** (default 1, or
  `options(mdir.cores)`). Each chain gets its own L'Ecuyer-CMRG stream drawn from
  the user's generator, so results follow `set.seed()` and do not depend on the
  number of cores, nor on whether workers are forked (Unix-alikes) or sockets
  (Windows); the user's generator kind and state are restored. Parallel results
  differ from the serial results for the same seed. As `data.table` does, at most
  two cores are used when `_R_CHECK_LIMIT_CORES_` is set (as `R CMD check
  --as-cran` does); `parallel` (a base package) is the only addition to
  `Imports`. Parallel execution still never starts unless asked.

## Sampling

* **The phis are updated with the strategic latent variable integrated out**
  (`phi_update = "slice"`, the default). `Z` is linear in each `phi`, `Z = A + B
  phi`, so the conditional of `phi` given the weights and labels is explicit and
  one-dimensional, and a slice-sampling update (Neal, 2003) draws from it with
  nothing to tune: the width only affects the number of evaluations, each O(1)
  once `A` and `B` are known. It is a joint draw of `(phi, v)` given the
  weights and labels, so the target is unchanged. The previous update given
  `v` remains available as `phi_update = "gibbs"`, and reproduces the previous
  version's chains exactly for the same seed (checked with three views, missing
  data and outliers).
  * Checked: the draws match the exact conditional (quantiles, L = 3 and 5) with
    lag-1 autocorrelation below 0.05 at very different scales; the sampler
    recovers the prior of every `phi` and mass under a constant likelihood with
    four and five views, as it did with two and three; the posterior of the phis
    agrees with the Gibbs update with three views (16 independent chains each).
    Changing the power of `Z` in the update on purpose makes the full-sampler
    tests fail (checked by hand).
  * Efficiency, measured on simulated data: with two strongly associated views
    (600 items) the effective sample size of `phi` rose about 15-fold (124 to
    1800 per 4000 draws); with three to five moderately associated views (300
    items) the median `phi` effective sample size rose between 1.2 and 2 times
    (three replicates, noisy), with the same run time. Mixing of `mass` and of the
    cluster structure, which limit most fits, is unchanged. With four views and a
    small data set both updates leave a few chains in a slower, higher-`phi`
    mode, so several chains remain necessary.
  * Cost: each `phi` needs two evaluations of the partition recursion (`A` and
    `B`), so for up to eight views the run time is unchanged; with ten views a
    sweep took 11 ms instead of 4.5 ms (100 items, K = 4). Use
    `phi_update = "gibbs"` if the number of views is large and the densities are cheap.
* The recorded `evidence` is unchanged: `Z` at the start of the sweep.

## Output size

* **Allocation probabilities are recorded for semi-supervised views only.** The
  N x K x draws array of `allocation_probabilities` was recorded for every view
  but used only by semi-supervised ones (92 MB against 0.9 MB for the allocations
  in a two-view run with K = 100 and 200 draws; about 1.2 GB per view per chain for
  1000 items, K = 500 and 300 draws). It is now `NULL` for other views.
  `calcAllocProb()` says so if asked for one.

## Console output and workflow

* **Fits print as reports.** `callMDI()` returns an `mdir_fit` and
  `runMCMCChains()` an `mdir_fit_list`. Both are still plain lists (`$` and `[[`
  work as before), but `print()` now gives a short description (views, sampler
  settings, run time, whether a burn in has been applied) instead of dumping
  every sampled array; printing one chain used to produce thousands of lines.
* **`summary()` methods**, laid out after `mclust` and Stan: occupied
  components per view, posterior summaries of the view association `phi` (mean,
  sd, quantiles), mean log-likelihoods and, after `processMCMCChain()`, a
  clustering table per view. For several chains, a per-chain table and the
  convergence diagnostics. The burn in defaults to half of `R` and is stated in
  the output; `summary(fit, burn = )` overrides it.
* **`fitMDI()`** is the recommended entry point. It is `runMCMCChains()` plus
  `assessConvergence()` (same arguments, same order, plus `burn` and
  `verbose`): it reports each chain as it finishes, then a verdict on
  convergence, and attaches the diagnostics as `attr(., "convergence")`. It
  warns if the diagnostics cannot be computed. The sampler is compiled C++ and
  does not report progress within a chain.
* **`runMCMCChains()`** gains `verbose` (default `FALSE`, so existing behaviour
  is unchanged) and rejects an invalid `n_chains`.
* **Convergence output**: `print()` on an `assessConvergence()` result is now a
  fixed-width table (the `converged` column used to wrap onto a separate block),
  shows each chain's mean log-likelihood, and ends with a verdict that separates
  quantities with a high Rhat (chains disagree) from those with only a low ESS
  (chains agree, too few effective samples). `format()` returns the verdict.
  The returned data frame gains attributes `min_ess` and `chain_loglik`; its
  columns are unchanged.
* `processMCMCChains()` keeps the class and the convergence attribute; `[` on a
  list of chains keeps the class and drops the convergence attribute (it
  described the full set).
* Progress and the convergence verdict respect `options(mdir.quiet = TRUE)`,
  which already silenced the prior-sparsity message.

## Performance

* **Faster sampling with many views.** The sums over set partitions that give
  the MDI normalising constant and the rates of the weight conditionals now
  reuse the part that depends on the `phi`s alone, and take the rates of all
  the weights of a view from one pass instead of one pass per weight. On one
  machine a 100-iteration chain took about 0.5 s for 8 views and 2.4 s for 10
  views before, and 0.15 s and 0.55 s after. Fewer than five views are
  unaffected, because the densities dominate there.
* **Faster multivariate normal likelihood** for complete items (no per-item,
  per-component allocation): about 1.7 times faster for a two-view MVN chain
  with 500 items and ten components.
* Apart from the label-swap correction below, the targets and the random number
  streams are unchanged. Compared with the previous version before that
  correction, allocations and all other discrete output were identical and
  continuous output agreed to about 1e-13 relative (bit-identical for one and
  two views); the last digits can differ with more views because the sums are
  accumulated in a different order.

## Corrections

* **Label swaps now exchange the weights of one view only.** The move that
  exchanges two components of a view (labels, weights and component parameters)
  also exchanged the weights of every other view, without their labels. The
  same permutation of all views leaves the normalising constant `Z` unchanged, so
  the `-v (Z' - Z)` term of the acceptance ratio was always zero, and the other
  views' weights no longer matched their labels; the ratio was not the ratio of
  the model's target. The acceptance ratio is now `sum_m log(1 + phi_m,l) (A' -
  A) - v (Z' - Z)` with `Z'` from exchanging the two weights of view `l` alone,
  and equals the exact log ratio of the target (checked against enumeration and
  symbolically). Effects of the old move: with equal numbers of components the
  sampler's stationary distribution was slightly wrong (with a constant
  likelihood, 3 components in each of two views and 10 items, the share of items in
  the largest-weight component of a view was about 0.011 too low, against a
  forward simulation of the prior); with different numbers of components the
  recorded weights had non-zero values in unused slots, and a view's weights
  could be moved entirely out of its real slots, making `Z` zero and stopping the
  sampler (a `randg()` error or "Non-finite allocation probabilities"). In 240
  random models of three to five views with unequal numbers of components and
  30,000 sweeps each, the old move stopped 6 of 112 runs without a
  single-component view and 95 of 128 with one, and left non-zero weights in
  unused slots in nearly all runs; none of these occurred after the correction,
  with semi-supervised and unsupervised views mixed. Chains with label swaps (every tenth sweep)
  change; results from earlier versions of models with several views are not
  reproduced exactly.

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
