// density.h
// =============================================================================
// include guard
#ifndef DENSITY_H
#define DENSITY_H

// =============================================================================
// included dependencies
// #define ARMA_WARN_LEVEL 0 // Turn off warnings that occur due to point errors.

# include <RcppArmadillo.h>
# include <RcppParallel.h>
# include <execution>

# include "logLikelihoods.h"
# include "genericFunctions.h"

using namespace arma ;

// =============================================================================
// virtual density class

class density {
  
private:
  
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
    N_inds,
    P_inds,
    K_inds,
    
    // Acceptance count of MH sampled parameters
    acceptance_count = zeros< uvec >(0);
  
  vec ll, likelihood,
  
    // Used in recording GP hyperparameters
    proposal_windows = zeros< vec >(0),
    hypers = zeros< vec >(0);
  
  mat X, X_t;
  
  // Parametrised class
  density(
    arma::uword _K,
    arma::uvec _labels,
    arma::mat _X);
  
  
  // Destructor
  virtual ~density() { };
  
  // The virtual functions that will be defined in every subclasses
  virtual void sampleFromPriors() = 0;
  virtual void sampleKthComponentParameters(
    uword k, 
    umat members, 
    uvec non_outliers
  ) = 0;
  
  virtual void sampleParameters(arma::umat members, arma::uvec non_outliers);
  virtual arma::vec itemLogLikelihood(arma::vec x) = 0;
  virtual double logLikelihood(arma::vec x, arma::uword k) = 0;

  virtual void receiveHyperParametersProposalWindows(vec proosal_windows) {};
  
};

#endif /* DENSITY_H */