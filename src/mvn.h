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
private:
  // Pre-allocated working matrices
  mutable arma::mat temp_cov_obs, temp_cov_miss, temp_cov_cross, temp_L;
  mutable arma::vec temp_mu_miss, temp_mu_obs, temp_residual;
  mutable arma::mat temp_solve_matrix;
  mutable arma::vec temp_solve_vector;
  
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
  
  // Missing value methods
  void initializeMissingValues() override;
  void sampleMissingForObservation(arma::uword n) override;
  
  // Modified likelihood functions
  arma::vec itemLogLikelihood(arma::uword n) override;
  double logLikelihood(arma::uword n, arma::uword k) override;
  
};


#endif /* MVN_H */