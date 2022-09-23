// mvt.h
// =============================================================================
// include guard
#ifndef MVT_H
#define MVT_H

// =============================================================================
// included dependencies
# include "outlierComponent.h"
# include "genericFunctions.h"
# include "logLikelihoods.h"

using namespace arma ;

// =============================================================================
// mvt class of outlier component

class mvt : virtual public outlierComponent {
  
private:
  
public:
  
  double global_log_det = 0.0,
    t_likelihood_const = 0.0,
    
    // T degrees of freedom
    df = 4.0;
      
  // The dataset mean
  vec global_mean;
  
  // The dataset covariance
  mat global_cov, global_cov_inv;
  
  using outlierComponent::outlierComponent;
  
  // Parametrised class
  mvt(arma::uvec _fixed, arma::mat _X);
  
  // Destructor
  virtual ~mvt() { };
  
  // Calculate the likelihood of each item being an outlier
  double calculateItemLogLikelihood(arma::vec x) override;
  
  arma::mat findInvertibleGlobalCov(double threshold = DBL_EPSILON);
  
};

#endif /* MVT_H */