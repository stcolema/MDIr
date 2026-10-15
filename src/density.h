// density.h
// =============================================================================
// include guard
#ifndef DENSITY_H
#define DENSITY_H

// =============================================================================
// included dependencies
# include <RcppArmadillo.h>

# include "logLikelihoods.h"
# include "genericFunctions.h"

// Validates and fills the density-level prior options (empty gives the defaults)
arma::vec resolveDensityPrior(const arma::vec& prior);

using namespace arma ;

// =============================================================================
// virtual density class
//
// Missing data. Entries of X that are not finite (NA / NaN / Inf in R) are
// treated as missing at random and handled by data augmentation: they are
// imputed from their full conditional distribution given the item's current
// component, and the imputed X is then used for the parameter updates. The
// allocation step uses the likelihood of the observed entries only (missing
// entries are marginalised), and the imputation follows immediately after the
// allocation step (a partially collapsed Gibbs sampler).

class density {

public:

  uword
    // The number of components modelled
    K,

    // The number of components occupied (i.e. clusters/groups)
    K_occ,

    // The dimensions of the dataset, samples and columns respectively
    N,
    P,

    // The number of parameters in the model
    n_param = 0;

  double complete_likelihood = 0.0, observed_likelihood = 0.0, BIC = 0.0;

  uvec
    // The cluster/class labels
    labels,

    // The number of items in each class
    N_k,

    // Sequence of integers that we iterate over
    K_inds,

    // Acceptance count of MH sampled parameters
    acceptance_count = zeros< uvec >(0);

  vec
    // Used in recording GP hyperparameters
    hypers = zeros< vec >(0);

  // The data, with any missing entries replaced by their current imputation
  mat X;

  // Missing value storage (common to all densities)
  arma::field<arma::uvec> missing_indices;
  arma::field<arma::uvec> observed_indices;
  arma::umat has_missing;

  // Options for the density-level (hierarchical) priors, ordered as
  //  (scale_pool_shape, gp_pool, gp_min_length, gp_center_sd, gp_pool_sd_scale);
  // see resolveDensityPrior()
  arma::vec density_prior;
  
  density(
    arma::uword _K,
    arma::uvec _labels,
    arma::mat _X,
    arma::vec _density_prior = arma::vec());

  virtual ~density() { };

  // === Sampling ==============================================================
  virtual void sampleFromPriors() = 0;
  virtual void sampleKthComponentParameters(
    uword k,
    const umat& members,
    const uvec& non_outliers
  ) = 0;
  
  // Updates the component parameters. Occupied components (those with at least 
  // one non-outlier member) are updated first, then the pooled hyperparameters 
  // given the occupied components only, then the empty components are drawn 
  // from their prior given the new hyperparameters. Conditioning the 
  // hyperparameters on the occupied components alone is the exact conditional
  // with the empty components integrated out; using the empty ones too would 
  // make the hyperparameters a random walk driven by their own prior draws.
  virtual void sampleParameters(const arma::umat& members, const arma::uvec& non_outliers);
  
  // Hyperparameters shared across components (partial pooling)
  virtual void updatePooledHyperparameters(const arma::uvec& occupied) { }
  virtual arma::vec pooledHyperparameters() const { return arma::vec(); }

  // === Likelihood ============================================================
  // Log-likelihood of the *observed* entries of item n in each / one component
  virtual arma::vec itemLogLikelihood(arma::uword n) = 0;
  virtual double logLikelihood(arma::uword n, arma::uword k) = 0;

  virtual void receiveHyperParametersProposalWindows(vec proposal_windows) {};

  // === Component relabelling =================================================
  // Exchange every component-specific parameter of components k and k'
  virtual void swapComponents(uword k, uword kprime) = 0;

  // === Predictive checks =====================================================
  // The component-specific parameters flattened into a vector, the inverse,
  // and a draw of one full observation from component k
  virtual arma::vec parameters() const = 0;
  virtual void setParameters(const arma::vec& theta) = 0;
  virtual arma::vec simulate(arma::uword k) const = 0;

  // The data-driven hyperparameters of the priors, for reporting
  virtual Rcpp::List hyperparameterList() const { return Rcpp::List::create(); }

  // === New data ==============================================================
  // Replace the data (and its missing-value patterns) by new items with the same
  // columns, so the likelihood of the new items can be evaluated at saved
  // parameters. The hyperparameters set from the original data are kept.
  virtual void replaceData(const arma::mat& X_new);

  // === Missing data ==========================================================
  void identifyMissingValues();

  // Starting values for the missing entries (column mean plus noise by default)
  virtual void initializeMissingValues();

  // Sample the missing entries of item n given its label
  virtual void sampleMissingForObservation(arma::uword n) = 0;

};

#endif /* DENSITY_H */
