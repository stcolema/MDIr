// tagmMixture.h
// =============================================================================
// include guard
#ifndef TAGMMIXTURE_H
#define TAGMMIXTURE_H

// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "mvnMixture.h"
# include "tAdjustedMixture.h"
# include "genericFunctions.h"

using namespace arma ;

// =============================================================================
// tagmMixture class

//' @name tagmMixture
//' @title T-ADjusted Gaussian Mixture (TAGM) type
//' @description The sampler for the TAGM mixture model.
//' @field new Constructor \itemize{
//' \item Parameter: K - the number of components to model
//' \item Parameter: labels - the initial clustering of the data
//' \item Parameter: concentration - the vector for the prior concentration of
//' the Dirichlet distribution of the component weights
//' \item Parameter: X - the data to model
//' }
//' @field printType Print the sampler type called.
//' @field calcBIC Calculate the BIC of the model.
//' @field calcTdistnLikelihood Calculate the likelihood of a given data point
//' the gloabl t-distirbution. \itemize{
//' \item Parameter: point - a data point.
//' }
class tagmMixture : public tAdjustedMixture, public mvnMixture {
  
private:
  
public:
  // for use in the outlier distribution
  // arma::uvec outlier;
  // arma::vec global_mean;
  // arma::mat global_cov;
  // double df = 4.0, u = 2.0, v = 10.0, b = 0.0, outlier_weight = 0.0;
  
  using mvnMixture::mvnMixture;
  
  tagmMixture(arma::uword _K,
              arma::uvec _labels,
              // arma::vec _concentration,
              arma::mat _X
  ) ;
  
  // Destructor
  virtual ~tagmMixture() { };
  
  virtual void calcBIC();
};

#endif /* TAGMMIXTURE_H */