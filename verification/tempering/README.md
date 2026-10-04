# Parallel tempering: derivations and checks

What is claimed, what was checked, how, and what was not. Nothing here is a proof
of convergence in finite time; none can be given for a general multimodal posterior.

## The algorithm

State of one replica: labels, component parameters, weights, phis, masses, pooled
scales, strategic latent variable `v`. Replica `t` has inverse temperature `beta_t`
(increasing, `beta_T = 1`) and targets

    pi_beta(state)  proportional to  L(state)^beta * P(state)

* `L(state) = prod_{l,n} f_l(x_nl | theta_{l, c_nl})`, the likelihood of the data given the
  labels and the component parameters;
* `P(state)` is everything else: the priors on parameters, weights, phis and masses, and the
  coupling `prod_n prod_l w_{l c_nl} prod_{l<m}(1 + phi_lm 1[c_nl = c_nm]) / Z^N`
  (with the strategic latent variable). `P` is not tempered.

Updates that involve `L`: the component parameters (conjugate; `beta` multiplies the counts
and the sufficient statistics) and the label draw (`beta` multiplies the component
log-likelihood). Weights, phis, masses, pooled scales and `v` have the same conditional
at every `beta`. The label-swap move exchanges component parameters with labels, which
leaves `L` unchanged, so it is valid at every `beta` (its acceptance does not involve `L`).

Exchange between neighbouring temperatures `i` and `j`: accept with probability
`min(1, exp((beta_i - beta_j) (l_j - l_i)))` where `l = log L(state)` of the state now at `i`
and `j` respectively; the priors cancel. Rounds alternate between even and odd neighbour
pairs (`"deo"`, Syed et al., 2022) or choose the parity at random (`"seo"`).

The `beta = 1` replica's draws are returned.

## What was checked, and by what

| Claim | Check | File |
|---|---|---|
| Exchange ratio is `exp((beta_i - beta_j)(l_j - l_i))`; priors cancel | SymPy; Lean (`prodTarget_exch_ratio`) | `sympy_checks.py`, `lean/Tempering/ProductSwap.lean` |
| Same as the form `exp min{0, (beta_{i+1}-beta_i)(V_{i+1} - V_i)}` of Syed et al. (2022, eq. 6), `V = -l` | Lean (`acc_exch_eq`) | `ProductSwap.lean` |
| Metropolis kernel with an involutive proposal is reversible, hence invariant | Lean (`swapK_detailedBalance`, `swapK_invariant`), finite state space | `SwapKernel.lean` |
| Composition and mixtures of invariant kernels are invariant (so any list of exchanges, even/odd rounds, SEO) | Lean (`invariant_comp`, `invariant_mix`, `roundK_invariant`) | `SwapKernel.lean`, `PTKernel.lean` |
| Independent local kernels, each invariant for its tempered target, give an invariant product kernel; one full PT iteration is invariant | Lean (`prodKernel_invariant`, `pt_iteration_invariant`) | `ProductSwap.lean`, `PTKernel.lean` |
| The cold coordinate's marginal of the product target is the (normalised) cold target | Lean (`marginal_of_product`) | `ProductSwap.lean` |
| Drawing a block from its full conditional leaves the target invariant (finite blocks) | Lean (`gibbsK_invariant`) | `Conjugate.lean` |
| Tempered update = untempered update with `beta n` observations and `beta` times the centred sum of squares (scalar) | Lean (`tempered_square`); SymPy (scalar NIG and `P = 2` NIW, symbolic); quadrature for the marginal likelihood | `Conjugate.lean`, `sympy_checks.py` |
| Tempered categorical: exponent `alpha - 1 + beta n` | Lean (`tempered_power`); quadrature | `Conjugate.lean`, `sympy_checks.py` |
| On a small finite product space, every swap, even round, odd round, DEO and SEO iteration and a local-sweep + exchange iteration leave the product target invariant; DEO is invariant but **not** reversible; SEO is reversible | exact rational arithmetic | `sympy_checks.py` |
| C++ conjugate updates at `beta < 1` draw from the tempered posterior (G, MVN, C) | simulation against the closed forms | `tests/testthat/test-tempering.R` |
| The quantity the exchange uses (`l`) equals an independent R recomputation from the saved allocation and parameters | testthat (G, MVN, C) | `test-tempering.R` |
| At `beta = 1` and one replica the new code is the old sampler, draw for draw | seeded comparison of the two builds (G, MVN, C, TAGM, missing data, 1 to 3 views) | run when the feature was added, see the report |
| Sampler (single tempered chain, and PT cold chain) reproduces the exact label posterior, one and two views | exact enumeration (categorical; Gaussian with fixed scale) and Monte Carlo prior for two views; negative controls with the wrong `beta` are rejected | `run_L1.R`, `run_L2.R`, `run_L1_gauss.R`, `exact_reference.R` |
| Benefit and exactness in a multimodal case | exact masses of the merge patterns against plain chains and PT | `mode_mass.R`, `pooled_plain.R` |
| Ladder update equalises rejection rates | known Gaussian path | `test-tempering.R` |

Tests were also run against two deliberately broken versions (wrong sign in the exchange,
counts not tempered in the MVN update); each was caught.

## What was not checked

* Lean covers finite state spaces and real arithmetic. It does not cover general state
  spaces, the continuous full conditionals (checked symbolically and by simulation
  instead), irreducibility, ergodicity or any rate of convergence. It also does not
  verify that the C++ implements the algorithm; that link is the simulation and the
  independent recomputation of `l`.
* Models other than complete-data `G`, `MVN`, `C` views without an outlier component. The
  missing-data, t-outlier and Gaussian-process updates are not tempered, and a ladder is
  refused for them.
* The mixing time of PT for any model. `mode_mass.R` is one dataset with one model.
* The matrix (NIW) identity is symbolic for `P = 2`, not for general `P` and not in Lean.

## Running

```sh
python3 sympy_checks.py                       # needs sympy, mpmath
# simulations (install the package first; arguments: library, chains, sweeps, cores)
Rscript run_L1.R <lib> 20 150000 3
Rscript run_L2.R <lib> 20 150000 3
Rscript run_L1_gauss.R <lib> 20 120000 2
Rscript make_ladder.R <lib> 2.5 mm.rds.ladder.rds
Rscript mode_mass.R <lib> 2.5 12 20000 1 mm.rds
Rscript pooled_plain.R <lib> mm.rds 160 6000 1
# Lean: Mathlib v4.33.0 with its cache fetched (lake exe cache get), Lean 4.33.0
cd lean && ./check.sh /path/to/mathlib4 SwapKernel ProductSwap Conjugate PTKernel Check
```

`lean/lakefile.toml` pins Mathlib `v4.33.0`; it was not used to build (the files were
compiled with `check.sh` against a local checkout). The Lean toolchain used here came from
conda-forge (`lean4` 4.33.0), whose bundled `leantar` (0.1.19) could not read the current
Mathlib cache; `leantar` 0.1.20 built from `digama0/leangz` was substituted.

## References (checked against the source text where noted)

* Syed, Bouchard-Cote, Deligiannidis, Doucet (2022), Non-reversible parallel tempering:
  a scalable highly parallel MCMC scheme, JRSSB 84(2), 321-350 (arXiv 1905.02939; the swap
  acceptance, the even/odd alternation, Corollary 1, Theorem 3 and Algorithms 2-3 were
  read in the arXiv text). Their round-trip results assume local equilibrium at each
  temperature and exact sampling of the reference; neither holds for this sampler.
* Woodard, Schmidler, Huber (2009), Conditions for rapid mixing of parallel and simulated
  tempering on multimodal distributions, Ann. Appl. Probab. 19(2), 617-640 (exists;
  cited for the limitation that tempering can mix slowly under phase transitions; the
  full text was not read for this work).
* Jasra, Holmes, Stephens (2005), Markov chain Monte Carlo methods and the label switching
  problem in Bayesian mixture modeling, Statist. Sci. 20(1), 50-67 (exists; relevance to
  tempering for mixtures not re-read).
* Coleman, Kirk, Wallace (2022), Consensus clustering for Bayesian mixture models, BMC
  Bioinformatics 23, 290. Its abstract frames the method as applying consensus clustering to
  heuristic clustering algorithms obtained from Bayesian mixture models by early stopping.
* Yao, Vehtari, Gelman (2022), Stacking for non-mixing Bayesian computations: the curse
  and blessing of multimodal posteriors, JMLR 23(79), 1-45 (exists; weights parallel runs
  by stacking for predictive performance, a different aim from recovering the posterior).
