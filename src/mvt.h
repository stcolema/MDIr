// mvt.h
// =============================================================================
// include guard
#ifndef MVT_H
#define MVT_H

# include "outlierComponent.h"
# include "genericFunctions.h"
# include "logLikelihoods.h"

using namespace arma ;

// =============================================================================
// mvt class of outlier component
//
// Multivariate t with fixed location (the column means of the observed data), 
// scale matrix (half the pairwise-complete data covariance) and df = 4 degrees 
// of freedom.
class mvt : virtual public outlierComponent {
  
public:
  
  double df = 4.0;
      
  // The dataset location and scale matrix
  vec global_mean;
  mat global_cov;
  
  mvt(
    arma::uvec _fixed, 
    arma::mat _X,
    const arma::field<arma::uvec>* miss_idx,
    const arma::field<arma::uvec>* obs_idx
  );
  
  virtual ~mvt() { };
  
  double calculateItemLogLikelihood(arma::uword n) override; 
  arma::vec sampleMissingValues(arma::uword n) const override;
  double completeLogDensity(const arma::vec& x) const override;
  arma::vec simulate() const override;
};

#endif /* MVT_H */
