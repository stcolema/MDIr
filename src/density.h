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

  density(
    arma::uword _K,
    arma::uvec _labels,
    arma::mat _X);

  virtual ~density() { };

  // === Sampling ==============================================================
  virtual void sampleFromPriors() = 0;
  virtual void sampleKthComponentParameters(
    uword k,
    const umat& members,
    const uvec& non_outliers
  ) = 0;
  virtual void sampleParameters(const arma::umat& members, const arma::uvec& non_outliers);

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

  // === Missing data ==========================================================
  void identifyMissingValues();

  // Starting values for the missing entries (column mean plus noise by default)
  virtual void initializeMissingValues();

  // Sample the missing entries of item n given its label
  virtual void sampleMissingForObservation(arma::uword n) = 0;

};

#endif /* DENSITY_H */
