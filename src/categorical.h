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
  
  // The number of categories in each measurement (categories are coded 0, 1, ...)
  uvec n_cat;
  
  // The data converted to a matrix of integers (missing entries hold their 
  // current imputation)
  umat Y;
  
  // The prior concentration of each category in each measurement: the observed 
  // category frequencies, so the Dirichlet prior on each component's category 
  // probabilities has total concentration 1 (floored for absent categories)
  field<vec> cat_prior_probability;
  
  // The probability of each category within each component; entry p is a 
  // n_cat(p) x K matrix
  arma::field<arma::mat> category_probabilities;
  
  categorical(arma::uword _K, arma::uvec _labels, arma::mat _X);
  
  virtual ~categorical() { };
  
  void sampleFromPriors() override;
  void sampleKthComponentParameters(uword k, const umat& members, const uvec& non_outliers) override;
  void initialiseParameters();
  
  void initializeMissingValues() override;
  void sampleMissingForObservation(arma::uword n) override;
  
  arma::vec itemLogLikelihood(arma::uword n) override;
  double logLikelihood(arma::uword n, arma::uword k) override;
  
  void replaceData(const arma::mat& X_new) override;
  
  void swapComponents(uword k, uword kprime) override;
  
  // Layout: category_probabilities(0), category_probabilities(1), ..., each 
  // n_cat(p) x K column-major
  arma::vec parameters() const override;
  void setParameters(const arma::vec& theta) override;
  arma::vec simulate(arma::uword k) const override;
  
  Rcpp::List hyperparameterList() const override;
  
private:
  void normaliseColumns(uword p);
};

#endif /* CATEGORICAL_H */
