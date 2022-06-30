// gaussian.h
// =============================================================================
// include guard
#ifndef GAUSSIAN_H
#define GAUSSIAN_H

// =============================================================================
// included dependencies
# include "density.h"

using namespace arma ;

// =============================================================================
// virtual gaussian class

//' @name gaussian
//' @title Gaussian density
//' @description Class for the MVN density with covariance matrices restricted 
//' to be diagonal.
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
class gaussian : virtual public density
{
  
public:
  
  // Parameters and hyperparameters
  double kappa = 0.01, nu = 3.0;
  
  arma::vec xi, scale;
  arma::mat mu, std_devs, precisions, log_std_devs;
  
  using density::density;
  
  gaussian(arma::uword _K, arma::uvec _labels, arma::mat _X);
  
  // Destructor
  virtual ~gaussian() { };
  
  // Calculate the empirical hyperparameters 
  arma::vec empiricalMean();
  arma::mat empiricalScaleVector();
  void empiricalBayesHyperparameters();
  
  // Sampling from priors
  void sampleStdDevPrior();
  void sampleMuPrior();
  void sampleFromPriors();
  
  void sampleKthComponentParameters(uword k, umat members, uvec non_outliers);
  // void sampleParameters(arma::umat members, arma::uvec non_outliers);
  double posteriorPredictive(arma::vec x, arma::uvec indices);
  
  // The log likelihood of a item belonging to each cluster
  arma::vec itemLogLikelihood(arma::vec item);
  
  // The log likelihood of a item belonging to a specific cluster
  double logLikelihood(arma::vec item, arma::uword k);
  
};


#endif /* GAUSSIAN_H */