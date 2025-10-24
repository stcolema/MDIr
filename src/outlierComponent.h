// outlierComponent.h
// =============================================================================
// include guard
#ifndef OUTLIERCOMPONENT_H
#define OUTLIERCOMPONENT_H

// =============================================================================
// included dependencies
// #define ARMA_WARN_LEVEL 0 // Turn off warnings that occur due to point errors.
# include <RcppArmadillo.h>

using namespace arma ;

// =============================================================================
// virtual outlierComponent class

class outlierComponent {
  
private:
  
public:
  
  uword 
    // The number of items in the dataset
    N, 
    
    // Number of measurements
    P;
  
  double
    // Outlier component weights
    non_outlier_weight = 1.0, outlier_weight = 0.0,
  
    // Hyperparameters for outlier weights
    u = 2.0, v = 10.0, tau_1 = 0.0, tau_2 = 0.0;
    
  uvec outliers, non_outliers;
    
  // Assume a global outlier likelihood with constant parameters
  vec outlier_likelihood;

  // The data and its transpose
  mat X, X_t;
  
  // Parametrised class
  outlierComponent(arma::uvec _fixed, arma::mat _X);
  
  // Destructor
  virtual ~outlierComponent() { };

  // Calculate the likelihood of each item being an outlier
  virtual void calculateAllLogLikelihoods();
  
  
  // Reference to missing patterns (shared with density)
  const arma::field<arma::uvec>* missing_indices_ref = nullptr;   
  const arma::field<arma::uvec>* observed_indices_ref = nullptr;  
  const arma::umat* has_missing_ref = nullptr;
  
  // Virtual missing value methods
  virtual void initializeMissingValues() = 0;
  virtual void sampleMissingForObservation(arma::uword n) = 0;
  virtual double calculateItemLogLikelihood(arma::uword n) = 0; // Modified signature
  
  // Set references to missing patterns (called by mixtureModel)
  void setMissingPatterns(const arma::field<arma::uvec>& miss_idx, 
                          const arma::field<arma::uvec>& obs_idx,
                          const arma::umat& has_miss) {
    missing_indices_ref = &miss_idx;
    observed_indices_ref = &obs_idx; 
    has_missing_ref = &has_miss;
  }
  
  // Update the outlier weights
  void updateWeights(uvec non_outliers, uvec outliers);
  
  // Sample if a given item is an outlier or not
  virtual arma::uword sampleOutlier(
    double non_outlier_likelihood_n,
    double outlier_likelihood_n
  );
  
};

#endif /* OUTLIERCOMPONENT_H */