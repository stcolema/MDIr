// mixture.h
// =============================================================================
// include guard
#ifndef MVNMIXTURE_H
#define MVNMIXTURE_H

// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "mixture.h"
# include "genericFunctions.h"

using namespace arma ;

// =============================================================================
// mvnMixture class

//' @name mvnMixture
//' @title Multivariate Normal mixture type
//' @description The sampler for the Multivariate Normal mixture model for batch effects.
//' @field new Constructor \itemize{
  //' \item Parameter: K - the number of components to model
//' \item Parameter: labels - the initial clustering of the data
  //' \item Parameter: concentration - the vector for the prior concentration of
//' the Dirichlet distribution of the component weights
  //' \item Parameter: X - the data to model
//' }
//' @field printType Print the sampler type called.
//' @field updateWeights Update the weights of each component based on current
//' clustering.
//' @field updateAllocation Sample a new clustering.
//' @field sampleFromPrior Sample from the priors for the multivariate normal
//' density.
//' @field calcBIC Calculate the BIC of the model.
//' @field logLikelihood Calculate the likelihood of a given data point in each
//' component. \itemize{
//' \item Parameter: point - a data point.
//' }
class mvnMixture: virtual public mixture {

public:

  // Each component has a weight, a mean vector and a symmetric covariance matrix.
  arma::uword n_param = 0;

  double kappa, nu;
  arma::vec xi, cov_log_det;
  arma::mat scale, mu, cov_comb_log_det;
  arma::cube cov, cov_inv;

  using mixture::mixture;

  mvnMixture(
    arma::uword _K,
    arma::uvec _labels,
    // arma::vec _concentration,
    arma::mat _X
  ) ;

  // Destructor
  virtual ~mvnMixture() { };

  // Sampling from priors
  virtual void sampleCovPrior();
  virtual void sampleMuPrior();
  virtual void sampleFromPriors();

  // Update the common matrix manipulations to avoid recalculating N times
  virtual void matrixCombinations();

  // The log likelihood of a item belonging to each cluster
  virtual arma::vec itemLogLikelihood(arma::vec item);
  
  // The log likelihood of a item belonging to a specific cluster
  virtual double logLikelihood(arma::vec item, arma::uword k);

  virtual void calcBIC();

  virtual void sampleParameters();
  
  double posteriorPredictive(arma::vec x, arma::uvec indices);
};

#endif /* MVNMIXTURE_H */