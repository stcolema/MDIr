# Guarantees for tempering and weighting: statements, assumptions, proofs

Every claim is tagged with the kind of statement it is:

* **[G1] exact, finite sample**: holds for any number of particles/iterations, no mixing assumed.
* **[G2] asymptotic**: holds in a stated limit under stated assumptions.
* **[C] conditional**: true if an assumption holds that the sampler's output cannot verify.
* **[D] diagnostic, no guarantee**: reported to help a human judge, not an estimator.

"Proved" means a Lean 4 proof (Mathlib v4.33.0, finite state spaces, standard axioms only) in
`lean/Tempering/`, or a proof below. "Checked" means exact rational arithmetic or simulation
against an exact reference; a check is evidence, not a proof. "Cited" means taken from a
source that has not been re-derived.

## 0. Setting

State `x = (labels c, component parameters theta, weights w, phis, masses, pooled scales, v)`.
For `beta in [0, 1]`, `pi_beta(x) ∝ L(x)^beta P(x)`, `L(x) = prod_{l,n} f_l(x_nl | theta_{l,c_nl})`
(complete data), `P` the remaining factors (priors, label coupling, the strategic latent variable);
`pi_0` is the prior (exactly samplable here, by ancestral sampling) and `pi_1` the posterior.
`Z_beta = ∫ L^beta P`. Only complete data, `G`/`MVN`/`C` views, no outlier component and (for
the weighted ensembles) no observed labels are covered.

## 1. Tempered conditionals [G1 given the state; checked]

The Gibbs conditionals of `theta` and of the labels under `pi_beta` are the untempered ones
with `n` replaced by `beta n` in the counts and `beta` times the scatter (NIG/NIW), or
`beta n_k` added to Dirichlet counts, or `beta` times the component log-likelihood in the
label draw; the weights, phis, masses, pooled scales and `v` have the same conditionals
as at `beta = 1` because `L` does not involve them.
*Evidence*: Lean `tempered_square`, `tempered_power` (scalar identities); SymPy (scalar NIG and
`P = 2` NIW, symbolic); quadrature for the marginal likelihoods; simulation of the C++ updates
against the closed forms; the label posterior against exact enumeration (`run_L1.R`,
`run_L2.R`, `run_L1_gauss.R`, `run_L1_fixed.R`). A draw from a full conditional leaves the
target invariant (Lean `gibbsK_invariant`, finite blocks).

## 2. Replica exchange

**Theorem 1 [G1, proved for finite state spaces]**. Let `Pi(s) = prod_t pi_{beta_t}(s_t)` and let the
exchange of temperatures `i, j` be accepted with probability `min(1, Pi(s∘(i j))/Pi(s))`. Then
(a) the ratio is `exp((beta_i - beta_j)(l_j - l_i))`, `l = log L`: the priors cancel
(`prodTarget_exch_ratio`); (b) the exchange kernel is `Pi`-reversible and `Pi`-invariant
(`swapK_detailedBalance`, `swapK_invariant`); (c) composition and mixtures of invariant
kernels are invariant, so any list of exchanges, the even/odd rounds and one full PT iteration
(independent local kernels followed by exchanges) are invariant (`invariant_comp`,
`invariant_mix`, `roundK_invariant`, `pt_iteration_invariant`); (d) if each `pi_{beta_t}` is
normalised the `t`-th marginal of `Pi` is `pi_{beta_t}` (`marginal_of_product`). The
non-reversibility of the even-odd scheme is not a violation: it is invariant and not reversible
(exact rational check, `sympy_checks.py`).
For general state spaces the same algebra holds with densities (Geyer, 1991; Syed et al., 2022);
This is not formalised.

**Proposition 2 [G2, argued, assumption A1]**. Assume (A1): every block update of the sweep has a
strictly positive conditional density with respect to a dominating measure on the support of
`pi_beta` (true by inspection for the supported updates: categorical draws from finite positive
logits, Gamma, inverse-Wishart, Dirichlet, normal; the mass update is a Metropolis step with a
positive log-normal proposal; the phi update is a slice sampler). Then the sweep kernel is
`pi_beta`-irreducible and aperiodic, with invariant law `pi_beta`; by the general theory of
Markov chains for posterior exploration (Tierney, 1994, Section 3 on irreducibility,
recurrence and the Gibbs/Metropolis hybrids; the theorem numbers were not re-checked) it is positive Harris
recurrent, so ergodic averages converge to `E_{pi_beta} f` from every starting point, and the same
holds for the product kernel with exchanges (positive density on the product support). This
proposition is argued here from the form of the conditionals, not proved formally. Hence the cold chain of PT
targets the posterior asymptotically **for any ladder**.
*What this does not say*: the plain chain (`betas = 1`) satisfies the same proposition.
Asymptotic correctness does not separate PT from a single chain; only the finite-time
behaviour does, and **there is no finite-time guarantee here**. Quantitative bounds
exist under conditions (Woodard, Schmidler and Huber, 2009; Syed et al., 2022, assume local
equilibrium and exact sampling of the reference); none has been verified for this model,
and the example shows they fail in practice (round trips 0 to 3 against about 900 predicted).
`ptDiagnostics()` is therefore [D]. Choosing the ladder (`ptLadder`, `tuneLadder`,
`adaptLadder`) affects efficiency only: validity does not depend on it (Theorem 1), provided
a tuning pilot is independent of the final run [G1].

## 3. Annealed importance sampling and sequential Monte Carlo from the prior

Particles start at exact draws from `pi_0`, each is moved by a `pi_{beta_k}`-invariant
kernel at step `k`, and carries the weight `prod_k exp((beta_k - beta_{k-1}) l(x_{k-1}))`.

**Theorem 3 [G1]**. For a **fixed** schedule, with or without resampling that depends only on the
weights, `E[Zhat] = Z_1/Z_0` and `E[Zhat · fhat(f)] = (1/Z_0) ∫ f L P` for every integrable `f`,
for every number of particles, whether or not the kernels mix.
*Proof/evidence*: without resampling, the forward-measure recursion is exactly
`annealMeasure_eq_last` and `ais_unbiased` (Lean, finite spaces; the identification of the
recursion with the expectation of "weight times indicator of state" is the definition of the
forward measure). With resampling: Del Moral (2004, Theorem 7.4.2, as stated by Naesseth et al.,
2019, Theorem 2.3.1). Checked in exact rational arithmetic (`sympy_smc_exact.py`):
`E[Zhat] = Z` and `E[Zhat fhat(s)] = Z pi(s)` hold to the last digit for never/always/ESS-triggered
resampling, multinomial and systematic, `N = 2, 3`.
Proof sketch of the resampling case: given the weighted system before resampling and the decision
to resample (a function of the weights), resampling replaces the weighted measure by an equally
weighted one with the same conditional expectation; the factor already absorbed in `Zhat`
carries the lost normalisation, and the next move is invariant, so the induction of the
no-resampling case goes through.

**Corollary 4 [G2]**. The ratio estimators (weighted frequency of a mode, weighted similarity
matrix, any `fhat(f)/Zhat`) are biased at finite size (exact value shown in
`sympy_smc_exact.py`, e.g. 2e-2 at `N = 2`) but consistent. Two limits are available:
(i) fixed `N`, `J -> infinity` independent runs combined by their evidence estimates:
`sum_j Zhat_j fhat_j / sum_j Zhat_j -> E_pi f` almost surely by the strong law of large numbers
for the i.i.d. pairs `(Zhat_j fhat_j, Zhat_j)`, whose means are `(Z E_pi f, Z)/Z_0`, with `Z>0`
(**no mixing assumption, no assumption on `N`**; `combineSMC`, `smcReplicates`); with finite
variance of `Zhat_j`, a central limit theorem (so a jackknife standard error over runs is
asymptotically valid: `smcSE`; cited, Shao and Tu, 1995, from memory);
(ii) `N -> infinity` in one run (Del Moral, Doucet and Jasra, 2006; Chopin, 2004; cited).

**Proposition 5 [G2, not G1]**. With an adaptive schedule (`schedule = "adaptive"`) the
temperatures depend on the particles; consistency as `N -> infinity` is known
(Beskos, Jasra, Kantas and Thiery, 2016; cited) but `E[Zhat] = Z` is false: exact biases `-2.2e-3`
and `-1.2e-3` in `sympy_smc_exact.py`. For unbiased estimates fix the schedule
(independently of the runs), as `smcReplicates` does.

**Proposition 6 [C, proved bound]**. If the particles are started at `beta_start > 0` from a law
`mu_0` instead of `pi_{beta_start}`, then for `|f| <= B` and incremental ratios bounded by
`R_k`: `|E[W f(x_K)] - E_{pi_start}[W f(x_K)]| <= B (prod_k R_k) sum_x |mu_0(x) - pi_start(x)|`
(Lean `bias_bound`; finite spaces). Under exact start it is Theorem 3. `sum_x|mu_0 - pi_start|`,
the total-variation error of the start, cannot be estimated from the output when the start is
made by running the sampler, so `beta_start > 0` carries **no guarantee** unless the start
is known to be exact. The experiment moved toward the exact values as the number of
start sweeps grew, which is an illustration.

**What Theorem 3 does not give**: a bound on the variance, which can be astronomically
large. In the four-cluster example the path from the prior to the posterior has a sharp
transition and the weights degenerate there: effective sample size 1.6 of 20,000 particles
without resampling; with resampling, wrong and unstable mode masses with a large ESS. The
finite-variance condition for the CLT is not checkable in general; `smcWeightDiagnostic()`
(Pareto `k̂`, Vehtari et al., 2024) estimates the tail of independent weights and is [D].
If the likelihood is bounded above (for example `l <= 0` for categorical data), the weights
are bounded and Hoeffding bounds hold for independent runs, but they are vacuous in
practice because `Z` is far below the maximum weight.

## 4. Pooling independent chains

**Proposition 7 [G1 for the identity, equal-weight pooling is not consistent]**. Let the posterior
be `pi = sum_m pi(A_m) pi_m` with `pi_m = pi(. | A_m)` mutually singular, and let independent chains be
confined to the basin in which they start, `b_m` the probability of landing in `A_m`. The pooled
equal-weight estimate converges to `sum_m b_m pi_m`, which equals `pi` iff `b_m = pi(A_m)` for all
`m`. (Immediate: the `pi_m` are linearly independent.) This is the basin-weight limit; more
chains do not change it (`pooled_plain.R`: TV 0.31 against the exact masses with 160 chains).
The weights of Section 3 are what correct it, with the guarantees stated there.

## 5. Weighted summaries

`weightedPSM`, `smcPosterior`: the ratio estimators of Corollary 4 [G2]. `weightedConsensus`:
the resampled weighted draws give a Monte Carlo estimate of the expected loss [G2], and
`salso` searches partitions to minimise it; **the search has no guarantee** of reaching the
minimiser (as in `processMCMCChain`) [D]. `smcDiagnostics`, `compareRuns`: [D]; agreement
of independent runs is necessary and not sufficient.

## 5b. Split-merge (sequentially allocated reallocation of a component pair)

**The move.** Pick two distinct components `a, b` uniformly. Let `F` be the non-fixed items currently
in `a` or `b`; the fixed (observed-label) items in `a, b` supply base statistics. Draw a uniform random
visiting order of `F`, allocate the items one at a time to `a` or `b` with probability proportional to
`w_c u_{c,n} m(X_c + x_n) / m(X_c)` (`w` weight, `u` the MDI coupling factor of item `n`, `m` the
collapsed marginal likelihood at `beta`, `X_c` the items already in `c` including fixed ones), run the
same walk along the current labels to get the reverse density, accept with
`min(1, prod_t Z_t(proposed path) / prod_t Z_t(current path))` (`Z_t` the step normaliser). On
acceptance the component parameters of `a, b` are redrawn from their conditionals, so the sweep
remains a valid Gibbs-type scheme. This is the sequentially allocated *style* of Dahl (2005) without
anchor items; it is **not** his algorithm verbatim, and it is a restricted reallocation within a pair,
so splits (one of `a, b` empty), merges and re-shuffles all arise from the one rule.

**Proposition SM1 [G1, finite state space; proved in Lean].**
(i) For any strictly positive proposal `Q` on a finite space, the Metropolis-Hastings kernel
`mhK` has `pi` invariant (`mhK_detailedBalance`, `mhK_invariant`). (ii) If the target factorises along
the visiting order as `prod_t s_t` and the proposal is `prod_t s_t / Z_t`, the MH ratio equals
`prod_t Z_t(y) / prod_t Z_t(x)` (`seq_ratio`), the quantity coded in `mixtureModel.cpp`.
**Application.** Condition on `(a, b, order)`, the cell of all labels outside `F`, and the items in
`F`: the move never changes which items are in `a ∪ b`, so it stays inside a finite cell on which
`pi` and `Q` are strictly positive; each cell kernel is invariant (i), the choice of `(a, b, order)` is
independent of the state, so the mixture is invariant (`invariant_mix`). The target is the collapsed
conditional of the labels of the view given the weights, the coupling factors and the fixed items;
`u` is held fixed during the move.
**Outlier component (TAGM).** Each item's state is a pair (component, flag). A non-outlier
score is as above times `(1 - eps)`; an outlier score is `w_c u_{c,n} eps l_out(n)`, where `l_out(n)` is the
fixed global-t density of the item (its parameters never change in this sampler; only `eps` is updated
and it is held fixed during the move) and the block statistics are unchanged by the item. Flagged items
are excluded from the collapsed statistics, as in the parameter update. The step normaliser sums four
terms; (i) and (ii) apply unchanged, so Proposition SM1 covers it. Missing values: the sampler imputes a missing cell of a flagged item
from the outlier law (conditional t) and of any other item from its component. The state of the move
includes the imputations (held fixed during the move): a non-outlier's complete-data density is part of its
component's collapsed marginal, an outlier's is the t density of the *complete* (observed + imputed) vector, so
the move uses that, not the observed-only `l_out` of the allocation step. The joint that results,
`p(c, o, theta, x_mis | x_obs)`, is the one the imputation draws target, so the move is a block update given the
imputation and the usual allocation step (observed-only, marginalising the imputations) followed by
re-imputation remains valid. [C: the argument; checked below.] Tempering still refuses outlier
components (the outlier density is not raised to `beta` in the allocation step).

**Conditions not proved here [C]:** that the collapsed marginals equal the integral of the tempered
likelihood against the prior (checked against quadrature and closed forms to 1e-8 to 1e-14, not
proved in Lean); that the C++ implements the stated ratio (checked below).

**Checks [evidence, not proof].**
* `sympy_splitmerge_exact.py`: exact rational enumeration of the kernel with fixed anchor items,
  averaged over all visiting orders (8, 16, 32 states): rows sum to one, detailed balance and
  invariance hold exactly; a mutant that drops the proposal terms violates detailed balance.
* `run_splitmerge.R` (move alone, weights fixed): 16 chains x 120,000 moves against exact
  enumeration for C, G, MVN, `beta = 0.5`, and observed labels: total variation at the chain-noise
  level in every case; control (`beta = 0.5` chain against the `beta = 1` target) TV 0.64.
* `run_L1_sm.R`, `run_L1_semisup.R`: full sampler with `split_merge = 2` (plain, and PT cold chain)
  against the exact label posterior, unsupervised and with two observed labels: consistent with
  exact; controls rejected. One cell (plain, `beta = 1`, unsupervised) had mean z^2 1.45 against 1.12
  expected for 60 states (about 1.8 standard deviations), the others 0.68 to 1.08.
* `run_splitmerge_outlier.R` (outlier component, move alone, `eps` fixed): exact enumeration over
  (label, flag) of the free items (N = 5, K = 3, 72 to 87 states with non-negligible mass) for `eps` 0.1 and 0.4 and with two
  observed labels: TV at chain-noise level, mean z^2 0.96 to 1.19; control (chain `eps` 0.1 against exact `eps` 0.4) TV 0.77.
  A mutant whose reverse walk ignores the current flags gives TV 0.12 and mean z^2 up to 676.
* `run_splitmerge_outlier_missing.R`: as above with three items having missing cells, imputed values held fixed
  (hook returns them and every subset's marginal): TV 0.017 (noise 0.016), mean z^2 1.20; the mutant that uses the observed-only outlier
  density in the move gives TV 0.71. The observed-only and complete outlier log-densities differ by up to 3.1 for these items.
* `run_splitmerge_tagm_missing.R`: full sampler, TAGM with four missing cells, alone and coupled to an MVN view, 48 chains each:
  max |z| 1.92 and 2.45, mean z^2 0.89 and 1.02 (agreement of two samplers).
* `run_splitmerge_tagm.R`: full sampler, TAGM view alone and TAGM + MVN (MDI), plain chain vs split-merge chain,
  co-clustering and per-item outlier probabilities (45 quantities): max |z| 2.64 and 2.29, mean z^2 1.97 and 1.15; the 1.97
  became 0.57 with 48 independent chains, so it is read as noise. Agreement of two samplers, not an exact reference.
* `run_L2_sm.R`: two views (MDI coupling), N = 4, K = 2, split-merge on, plain / `beta = 0.5` / PT cold
  chain against the exact posterior (MC prior reference, TV contribution of its error <= 1e-4): TV at chain-noise
  level; control rejected (TV 0.22). K = 2 makes the pair choice trivial.
* Missing data (`run_splitmerge_missing.R`): no closed-form reference, so the check is agreement
  of co-clustering probabilities between a plain chain and a split-merge chain (28 pairs, 24 chains
  each): max |z| 2.34, mean z^2 1.17. This is agreement of two samplers sharing the data-
  augmentation step, with limited power against a subtle bias. The move conditions on the current
  imputation; the argument is that the imputation is part of the state and the move is a block
  update given it [C].
* `tests/testthat/test-splitmerge.R` includes a mutant (always accept) that the exact test rejects
  in two of its three exact cases (the `beta = 0.5` case did not flag it).

**What it does to the four-cluster example [empirical, no guarantee of mixing].** 160 independent
plain chains with `split_merge = 2` pooled with equal weights: TV 0.015 to the group-unit
reference masses (0.31 without), every chain's dominant pattern the modal one. At separation 4
(where the reference is accurate) 80 chains x 20,000 sweeps gave TV 0.001. At separation 2.5 a gap of
about 0.015 on the second pattern persists with 5 times longer chains and an across-chain standard error
of 0.001, so it is not Monte Carlo noise; the reference treats each true group as an indivisible unit
and ignores item-level misallocation, which matters at that separation, and the discrepancy vanishes
at separation 4. This is attributed to the reference but not proved: computing the exact
mass without the group-unit approximation is not feasible here.
Mixing of the move is not guaranteed by Prop. SM1: invariance only. That it removed the basin
weights in this example is an observation about this example.

## 5c. Joint allocation of an item across views

**The move.** For each item, redraw its labels in a block of views from their exact joint conditional given the
weights, phis, component parameters and the item's labels in the other views, then (with outlier components) its
outlier flags given the components, then redraw any missing values given the new allocation. It follows the ordinary
per-view allocation step in the sweep (`joint_allocation`).

**Proposition JA1 [G1; the conditional is derived below, the invariance is Lean `gibbsK_invariant`].** Given
`(w, phi, theta, v)` and the labels of the other items, the items are independent (Z does not involve the labels),
and the conditional of the labels `c_l, l in B` of one item in a block `B` of views is proportional to
`prod_{l in B} a_l(c_l) prod_{l<m; l or m in B} (1 + phi_lm 1[c_l = c_m])`, with `a_l(k) = w_lk g_l(k)`, `g_l(k)`
the likelihood of the item in component `k` raised to `beta`, marginalised over its outlier flag where the view has an outlier
component, and a single point mass at the observed component for an observed label. Drawing a block from its full
conditional leaves the target invariant. Outlier components: the (component, flag) pair has weight
`w_k [(1 - eps) f_k + eps f_out]` summed over the flag, so drawing the component from the marginal and then the flag
given the component is an exact draw of the pair. Missing values: the allocation uses the likelihood of the observed
entries (the imputed entries marginalised) and the imputation is then redrawn, the same partially collapsed
ordering as the ordinary allocation step.

**The sequential draw.** The first view of the block is drawn from its marginal `a_1(k) dZ(a)/da_1(k) / Z(a)`, where `Z`
is the polynomial of `mdiPartition.h` in the block's scores (it is multilinear in each view's column, so the sum over the
other views' labels is `dZ/da_1(k)`, which the recursion of `mdiWeightRates()` computes). Given `c_1 = k`, the terms
that involve view 1 are `(1 + phi_1m 1[c_m = k])`, which multiply `a_m(k)` by `1 + phi_1m` and leave the other
terms of the form of the same polynomial over the remaining views with the sub-matrix of `phi`. Repeating this draws
the block's labels exactly, at a cost `O(b 3^b + b K)` per item for a block of `b` views, so that the number of
components enters linearly (a 125-component view costs no more than a 5-component one). Views outside the block
enter through `(1 + phi)` factors at the component their label indicates. The same draw gives exact prior labels beyond the
5e6 combinations that enumeration is limited to (`samplePriorLabels`, used by `smcMDI()` and
`simulatePriorPredictive()`).

**Checks [evidence, not proof].**
* `tests/testthat/test-joint-allocation.R`: the draw against enumeration of the exact conditional for ragged numbers
  of components, excluded components, blocks of two or three of four views, large `phi` (40), a flat likelihood (the prior draw), and
  scores spanning hundreds of log units (max |z| < 4.5 over all cells, chi-square p > 1e-4); the marginals of the draw against
  `mdiClassProbabilities()`; a mutant that drops the `(1 + phi)` conditioning fails 14 of the 21 enumeration and marginal expectations and three of the five of the prior-recovery test.
* `run_L2_joint.R` (two views, N = 4, K = 2, exact posterior with the Monte Carlo prior of the cell counts): the full
  sampler with the joint allocation, single chain `beta = 1`, `beta = 0.5` and the cold chain of parallel tempering, total variation
  0.0037, 0.0031 and 0.0034 against chain noise 0.0037, 0.0033 and 0.0032; the control (chain at `beta = 1`
  against the exact `beta = 0.5`) is rejected (0.22). Largest |z| over about 80 states 3.0, 3.2 and 4.5 (the last with mean z^2 1.30 against 1.15 expected).
* `test-joint-allocation.R`: the sampler recovers the prior under a constant likelihood; observed labels (with a gap in
  the classes) are never changed, for `"G"`, `"MVN"` and `"TAGM"` views.
* Not checked exactly: outlier components with missing values inside the joint step (only the observed-label invariance
  and the prior-recovery tests apply); the move with a block smaller than the number of views inside a full sampler.

**What it does for mixing [empirical, `verification/joint_allocation`].** On the larger scenario (three views, 30 true
clusters, `K = 45`, N = 300, four chains per configuration, three seeds, R = 4000) the mean Rhat of the key
log-likelihoods is 2.83 for the plain sampler, 2.42 with joint allocation, 2.75 with split-merge and 2.28 with both; the
number of occupied components in view 1 is 19 to 21 against 30 true clusters in all four, and the across-chain sd of the
joint log-likelihood is 114, 64, 87 and 71. Every configuration is far from converged, so the differences are not
established (three seeds, no interval). Cost per sweep rises by about 1.8 times with joint allocation. The mechanism
that the experiments point to is that aligned splits are held by the rich-get-richer prior across views, which no
single-item move removes. A collapsed joint split-merge across aligned components (Proposition SM1 extends to it) is the
move that would address this and is not implemented.

## 6. What would give a certificate, and is not built

A finite-sample certificate of mixing would need a bound on the total-variation distance of a
chain to stationarity, for example from L-lag couplings (Biswas, Jacob and Vanetti, 2019; cited,
not re-read), which needs coupled versions of every update of the sweep. Not implemented; in the
multimodal regime the meeting times would be the quantity that is exponentially long.

## 7. Reference status

Read in the source text: Syed et al. (2022) eq. 6, eq. 8, Cor. 1, Thm 3, Algs 1-3 (arXiv 1905.02939);
Zhou, Johansen and Aston (2016) eq. 3.16 (CESS) (arXiv 1303.3123); Naesseth, Lindsten and Schon
(2019) Thm 2.3.1 and Remark 2.2.1 (arXiv 1903.04797); Vehtari et al. (2024) the `k̂` thresholds and
the finite-variance equivalence (arXiv 1507.02646).
Existence and details confirmed by search only: Neal (2001, Stat. Comput. 11, 125-139); Del Moral,
Doucet and Jasra (2006, JRSSB 68, 411-436) and (2012, Bernoulli 18, 252-278); Chopin (2002,
Biometrika 89, 539-551) and (2004, Ann. Statist. 32, 2385-2411); Beskos et al. (2016, Ann. Appl.
Probab. 26, 1111-1146); Tierney (1994, Ann. Statist. 22, 1701-1728); Woodard, Schmidler and Huber
(2009, Ann. Appl. Probab. 19, 617-640); Jasra, Holmes and Stephens (2005); Yao, Vehtari and
Gelman (2022); Coleman, Kirk and Wallace (2022). Split-merge: Jain and Neal (2004, JCGS 13, 158-182), Dahl (2005, sequentially allocated merge-split; existence and
description via search results only, the paper itself not read), Bouchard-Cote, Doucet and Roth (2017, JMLR 18(28),
15-397; read in abstract only), Nguyen, Trippe and Broderick (2022, AISTATS, PMLR 151, 3483-3514; abstract only),
Monteiller et al. (2019, NeurIPS; abstract only). From memory, not checked: Del Moral (2004,
*Feynman-Kac Formulae*, Springer, Thm 7.4.2 as a numbered statement), Shao and Tu (1995, *The
Jackknife and Bootstrap*), Geyer (1991), Biswas, Jacob and Vanetti (2019), Meyn and Tweedie.
