// outlierComponent.h
// =============================================================================
// include guard
#ifndef OUTLIERCOMPONENT_H
#define OUTLIERCOMPONENT_H

# include <RcppArmadillo.h>

using namespace arma ;

// =============================================================================
// outlierComponent class
//
// A global (component-independent) distribution that absorbs items poorly
// described by any cluster, as in the t-augmented Gaussian mixture model of
// Crook et al. (2018). Each unlabelled item is an outlier with probability
// eps = outlier_weight, eps ~ Beta(outlier_prior_a, outlier_prior_b). Items 
// with an observed label are never outliers.
class outlierComponent {
  
public:
  
  uword 
    // The number of items in the dataset
    N, 
    
    // Number of measurements
    P;
  
  double
    // Outlier component weights
    non_outlier_weight = 1.0, outlier_weight = 0.0,
  
    // Beta prior on the outlier weight, Beta(a, b), mean a / (a + b)
    outlier_prior_a = 2.0, outlier_prior_b = 10.0;
    
  uvec outliers, non_outliers;
    
  // Log-likelihood of the observed entries of each item under the outlier 
  // distribution (the outlier distribution has fixed parameters)
  vec outlier_likelihood;

  // Data (used for the observed entries only)
  mat X;
  
  // Missing patterns, owned by the density
  const arma::field<arma::uvec>* missing_indices_ref = nullptr;   
  const arma::field<arma::uvec>* observed_indices_ref = nullptr;  
  
  outlierComponent(
    arma::uvec _fixed, 
    arma::mat _X,
    const arma::field<arma::uvec>* miss_idx,
    const arma::field<arma::uvec>* obs_idx
  );
  
  virtual ~outlierComponent() { };

  // True if the component can absorb items
  virtual bool active() const { return true; }
  
  // Replace the data by new items (the outlier distribution keeps its parameters)
  void replaceData(const arma::mat& X_new);
  
  void calculateAllLogLikelihoods();
  virtual double calculateItemLogLikelihood(arma::uword n) = 0;
  
  // Draw the missing entries of item n given its observed entries, assuming the
  // item is an outlier
  virtual arma::vec sampleMissingValues(arma::uword n) const = 0;
  
  // Draw a full observation from the outlier distribution
  virtual arma::vec simulate() const = 0;
  
  // Update the outlier weight from the current indicators
  virtual void updateWeights(const uvec& non_outliers, const uvec& outliers);
  
  // Draw the outlier weight from its Beta prior
  virtual void sampleFromPrior();
  
};

#endif /* OUTLIERCOMPONENT_H */
