// categorical.h
// =============================================================================
// include guard
#ifndef CATEGORICAL_H
#define CATEGORICAL_H

// =============================================================================
// included dependencies
# include "density.h"

using namespace arma ;

// =============================================================================
// categorical class

//' @name categorical
//' @title Categorical density
//' @description The Categorical density for the mixture model.
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
class categorical : virtual public density
{
  
public:
  
  // The number of categories in each measurement
  uvec n_cat;
  
  // This will hold the data converted to a matrix of integers
  umat Y;
  
  // The prior probability of being in a given category in each measurement
  field<vec> cat_prior_probability;
  
  // The probability of each class within category; it will be a N_cat x K x P
  // array
  arma::field<arma::mat> category_probabilities;
  
  using density::density;
  
  categorical(arma::uword _K, arma::uvec _labels, arma::mat _X);
  
  // Destructor
  virtual ~categorical() { };
  
  // Sampling from priors
  void sampleFromPriors();
  void sampleKthComponentParameters(uword k, umat members, uvec non_outliers);
  // void sampleParameters(arma::umat members, arma::uvec non_outliers);
  void initialiseParameters();
  
  // The log likelihood of a item belonging to each cluster
  arma::vec itemLogLikelihood(arma::vec item);
  
  // The log likelihood of a item belonging to a specific cluster
  double logLikelihood(arma::vec item, arma::uword k);
  // double posteriorPredictive(arma::vec x, arma::uvec indices);
  
  
};


#endif /* CATEGORICAL_H */