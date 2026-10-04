# Guarantees for tempering and weighting: statements, assumptions, proofs

Every claim is tagged with the kind of statement it is:

* **[G1] exact, finite sample**: holds for any number of particles/iterations, no mixing assumed.
* **[G2] asymptotic**: holds in a stated limit under stated assumptions.
* **[C] conditional**: true if an assumption holds that the sampler's output cannot verify.
* **[D] diagnostic, no guarantee**: reported to help a human judge, not an estimator.

"Proved" means a Lean 4 proof (Mathlib v4.33.0, finite state spaces, standard axioms only) in
`lean/Tempering/`, or a proof below. "Checked" means exact rational arithmetic or simulation
against an exact reference; a check is evidence, not a proof. "Cited" means taken from a
source I have not re-derived.

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
I did not formalise that.

**Proposition 2 [G2, argued, assumption A1]**. Assume (A1): every block update of the sweep has a
strictly positive conditional density with respect to a dominating measure on the support of
`pi_beta` (true by inspection for the supported updates: categorical draws from finite positive
logits, Gamma, inverse-Wishart, Dirichlet, normal; the mass update is a Metropolis step with a
positive log-normal proposal; the phi update is a slice sampler). Then the sweep kernel is
`pi_beta`-irreducible and aperiodic, with invariant law `pi_beta`; by the general theory of
Markov chains for posterior exploration (Tierney, 1994, Section 3 on irreducibility,
recurrence and the Gibbs/Metropolis hybrids; I did not re-check the theorem numbers) it is positive Harris
recurrent, so ergodic averages converge to `E_{pi_beta} f` from every starting point, and the same
holds for the product kernel with exchanges (positive density on the product support). This
proposition is argued here from the form of the conditionals, not proved formally. Hence the cold chain of PT
targets the posterior asymptotically **for any ladder**.
*What this does not say*: the plain chain (`betas = 1`) satisfies the same proposition.
Asymptotic correctness does not separate PT from a single chain; only the finite-time
behaviour does, and **there is no finite-time guarantee here**. Quantitative bounds
exist under conditions (Woodard, Schmidler and Huber, 2009; Syed et al., 2022, assume local
equilibrium and exact sampling of the reference); I have not verified any of them for this model,
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
Gelman (2022); Coleman, Kirk and Wallace (2022). From memory, not checked: Del Moral (2004,
*Feynman-Kac Formulae*, Springer, Thm 7.4.2 as a numbered statement), Shao and Tu (1995, *The
Jackknife and Bootstrap*), Geyer (1991), Biswas, Jacob and Vanetti (2019), Meyn and Tweedie.
