// runMDI.h
// =============================================================================
// include guard
#ifndef RUNMDI_H
#define RUNMDI_H

// =============================================================================
// included dependencies
# include "logLikelihoods.h"
# include "genericFunctions.h"
# include "mdi.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp ;
using namespace arma ;

// =============================================================================
// runMDI function header

//' @title Call Multiple Dataset Integration
//' @description C++ function to perform MCMC sampling for MDI.
//' @param R The number of iterations to run for.
//' @param thin thinning factor for samples recorded.
//' @param Y The list of data matrices to perform integrative clustering upon 
//' with items to cluster in rows. Non-finite entries are treated as missing.
//' @param K Vector of the number of components to model in each view. This is the upper 
//' limit on the number of clusters that can be found.
//' @param mixture_types Character vector of densities used in each view
//' @param outlier_types Character vector of outlier components used in each 
//' view ('MVT' or 'None').
//' @param labels Matrix item labels to initialise from. Rows correspond to the
//' items being clustered, columns to views.
//' @param fixed Binary matrix of the items that are fixed in their initial
//' label.
//' @param proposal_windows List/field of vectors
//' @param save_parameters Record the component parameters at each saved 
//' iteration (needed for posterior predictive checks).
//' @param save_imputed Record the imputed values of missing entries at each 
//' saved iteration.
//' @param prior Optional vector of MDI-level prior hyperparameters (mass shape, 
//' mass rate, weight rate, phi shape, phi rate). Empty for the defaults.
//' @param density_prior Options of the density-level priors (variance scale 
//' pooling and Gaussian process priors), see `densityPrior()` in R.
//' @param save_allocation_probabilities For each view, 1 to record the allocation
//' probabilities of every item at every saved iteration (an N x K x draws array, 
//' needed only for semi-supervised views), 0 to leave them out. A single value is
//' recycled over the views.
//' @param save_pointwise Record the log-likelihood of every item at every saved 
//' iteration (see `pointwiseLogLikelihood`); the total is always recorded.
//' @param phi_slice Update the phis with the strategic latent variable 
//' integrated out (slice sampling; TRUE) or by Gibbs sampling given it (FALSE).
//' @param betas Inverse temperatures of the replicas for parallel tempering,
//' strictly increasing, the last equal to one (the posterior). Empty or a single
//' one is an ordinary run. A single value other than one samples that tempered
//' target and is for testing only.
//' @param swap_scheme How replicas are paired for exchange: 0, deterministic
//' even-odd (non-reversible); 1, stochastic even-odd (reversible).
//' @param swap_every Attempt replica exchanges after every `swap_every` sweeps.
//' @param split_merge Attempts per view and sweep of the sequentially-allocated 
//' re-partition of two randomly chosen components (0 turns it off).
//' @param joint_allocation Size of the blocks of views whose labels are redrawn jointly for each
//' item, from their exact joint conditional (0 turns it off; see `mdi::joint_block`).
//' @return Named list of the different quantities drawn by the sampler.
// [[Rcpp::export]]
Rcpp::List runMDI(
  arma::uword R,
  arma::uword thin,
  arma::field<arma::mat> Y,
  arma::uvec K,
  arma::uvec mixture_types,
  arma::uvec outlier_types,
  arma::umat labels,
  arma::umat fixed,
  arma::field< arma::vec > proposal_windows,
  bool save_parameters,
  bool save_imputed,
  arma::vec prior,
  arma::vec density_prior,
  arma::uvec save_allocation_probabilities,
  bool save_pointwise,
  bool phi_slice,
  arma::vec betas,
  arma::uword swap_scheme,
  arma::uword swap_every,
  arma::uword split_merge,
  arma::uword joint_allocation
);

#endif /* RUNMDI_H */
