// runSMC.h
// =============================================================================
// include guard
#ifndef RUNSMC_H
#define RUNSMC_H

// =============================================================================
// included dependencies
# include "logLikelihoods.h"
# include "genericFunctions.h"
# include "mdi.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp ;
using namespace arma ;

//' @title Sequential Monte Carlo sampler for MDI (C++)
//' @description Annealed importance sampling / sequential Monte Carlo over the
//' likelihood-tempered targets \eqn{\pi_\beta \propto L^\beta P}, from the prior
//' (\eqn{\beta = 0}, sampled exactly) to the posterior (\eqn{\beta = 1}). Each
//' particle is moved by the tempered Gibbs sweep, which leaves \eqn{\pi_\beta}
//' invariant.
//' @param n_particles Number of particles.
//' @param Y,K,mixture_types,outlier_types,fixed,prior,density_prior,phi_slice As in `runMDI`.
//' `fixed` must be all zero (unsupervised views only).
//' @param betas Fixed schedule: strictly increasing values in (0, 1] ending at 1.
//' Ignored if `adaptive`.
//' @param adaptive Choose each next inverse temperature so that the conditional
//' effective sample size of the incremental weights is `cess_target` times
//' the number of particles (Zhou, Johansen and Aston, 2016, eq. 3.16).
//' @param cess_target Target conditional ESS fraction, in (0, 1).
//' @param resample_threshold Resample when the ESS of the weights falls below
//' this fraction of the number of particles (0 never, 1 always).
//' @param resample_scheme 0 systematic, 1 multinomial.
//' @param sweeps_per_step Gibbs sweeps of every particle after each temperature change.
//' @param max_steps Upper bound on the number of temperatures; if it is reached the
//' next temperature is set to 1.
//' @param final_sweeps Extra sweeps at beta = 1 after the last temperature, recorded
//' every `final_thin` sweeps (the state at the end of annealing is always recorded).
//' @param final_thin Thinning of the extra sweeps.
//' @param beta_start Inverse temperature at which the population starts, in [0, 1).
//' At 0 the particles are exact draws from the prior. Above 0 they are drawn from the
//' prior and then moved by `start_sweeps` sweeps at `beta_start`, so they follow
//' \eqn{\pi_{\beta_{start}}} only to the extent that those sweeps mix.
//' @param start_sweeps Sweeps at `beta_start` before annealing (ignored if `beta_start = 0`).
//' @param split_merge Split-merge attempts per view and sweep (see `runMDI`).
//' @return A list with the particles' recorded draws and weights, the evidence
//' estimate and the trace of the run.
// [[Rcpp::export]]
Rcpp::List runMDISMC(
  arma::uword n_particles,
  arma::field<arma::mat> Y,
  arma::uvec K,
  arma::uvec mixture_types,
  arma::uvec outlier_types,
  arma::umat fixed,
  arma::vec prior,
  arma::vec density_prior,
  bool phi_slice,
  arma::vec betas,
  bool adaptive,
  double cess_target,
  double resample_threshold,
  arma::uword resample_scheme,
  arma::uword sweeps_per_step,
  arma::uword max_steps,
  arma::uword final_sweeps,
  arma::uword final_thin,
  double beta_start,
  arma::uword start_sweeps,
  arma::uword split_merge
);

#endif /* RUNSMC_H */
