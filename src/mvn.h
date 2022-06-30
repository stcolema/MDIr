// mvn.h
// =============================================================================
// include guard
#ifndef MVN_H
#define MVN_H

// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "density.h"
# include "genericFunctions.h"

using namespace arma ;

// =============================================================================
// virtual mvn class

//' @name mvn
//' @title Multivariate Normal density
//' @description Class for the MVN density.
//' @field new Constructor \itemize{
//' \item Parameter: K - the number of components to model
//' \item Parameter: labels - the initial clustering of the data
//' \item Parameter: X - the data to model
//' }
//' @field sampleFromPrior Sample from the priors for the multivariate normal
//' density.
//' @field calcBIC Calculate the BIC of the model.
//' @field logLikelihood Calculate the likelihood of a given data point in each
//' component. \itemize{
//' \item Parameter: point - a data point.
//' }
class mvn : virtual public density
{
  
public:
  
  // Parameters and hyperparameters
  double kappa, nu;
  
  arma::vec xi, cov_log_det;
  arma::mat scale, mu, cov_comb_log_det;
  arma::cube cov, cov_inv;
  
  using density::density;
  
  mvn(arma::uword _K, arma::uvec _labels, arma::mat _X);
  
  // Destructor
  virtual ~mvn() { };
  
  // Calculate the empirical hyperparameters 
  arma::vec empiricalMean();
  arma::mat empiricalScaleMatrix();
  void empiricalBayesHyperparameters();
  
  // Sampling from priors
  void sampleCovPrior();
  void sampleMuPrior();
  void sampleFromPriors();
  
  void sampleKthComponentParameters(uword k, umat members, uvec non_outliers);
  void sampleParameters(arma::umat members, arma::uvec non_outliers);
  double posteriorPredictive(arma::vec x, arma::uvec indices);
  
  // Update the common matrix manipulations to avoid recalculating N times
  void matrixCombinations();
  
  // The log likelihood of a item belonging to each cluster
  arma::vec itemLogLikelihood(arma::vec item);
  
  // The log likelihood of a item belonging to a specific cluster
  double logLikelihood(arma::vec item, arma::uword k);
  
  
};


#endif /* MVN_H */