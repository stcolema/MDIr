// noOutliers.h
// =============================================================================
// include guard
#ifndef NOOUTLIERS_H
#define NOOUTLIERS_H

// =============================================================================
// included dependencies
# include "outlierComponent.h"

using namespace arma ;

// =============================================================================
// noOutliers class of outlier component

class noOutliers : virtual public outlierComponent {
  
private:
  
public:
  
  using outlierComponent::outlierComponent;
  
  // Parametrised class
  noOutliers(arma::uvec _fixed, arma::mat _X);
  
  // Destructor
  virtual ~noOutliers() { };
  
  // The likelihood of a given item
  void calculateAllLogLikelihoods();
  double calculateItemLogLikelihood(arma::vec x);
  
  // Update the outlier weights
  void updateWeights(uvec non_outliers, uvec outliers);
  
  // Sample if a given item is an outlier or not
  arma::uword sampleOutlier(double non_outlier_likelihood_n,
                            double outlier_likelihood_n);
  
};

#endif /* NOOUTLIERS_H */