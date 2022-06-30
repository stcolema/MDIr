// tAdjustedMixture.h
// =============================================================================
// include guard
#ifndef TADJMIXTURE_H
#define TADJMIXTURE_H

// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "mixture.h"
# include "genericFunctions.h"

using namespace arma ;

// =============================================================================
// tAdjustedMixture class

//' @name tAdjustedMixture
//' @title Base class for adding a t-distribution to sweep up outliers in the
//' model.
//' @description The class that the specific TAGM types inherit from.
//' @field new Constructor \itemize{
//' \item Parameter: K - the number of components to model
//' \item Parameter: labels - the initial clustering of the data
//' \item Parameter: concentration - the vector for the prior concentration of
//' the Dirichlet distribution of the component weights
//' \item Parameter: X - the data to model
//' }
//' @field updateOutlierWeights Updates the weight of the outlier component.
//' @field updateWeights Update the weights of each component based on current
//' clustering, excluding the outliers.
//' @field sampleOutlier Sample is the nth item is an outlier. \itemize{
//' \item Parameter n: the index of the individual in the data matrix and
//' allocation vector.
//' }
//' @field updateAllocation Sample a new clustering and simultaneously allocate
//' items to the outlier distribution.
//' @field calcTdistnLikelihood Virtual placeholder for the function that calculates
//' the likelihood of a given point in a t-distribution. \itemize{
//' \item Parameter: point - a data point.
//' }
class tAdjustedMixture : virtual public mixture {
  
private:
  
public:
  // for use in the outlier distribution
  double global_log_det = 0.0, t_likelihood_const = 0.0;
  
  // Vector indicating if the item is an outlier (value of 1) or not (value of 0)
  arma::uvec outlier;
  
  arma::vec 
    // The dataset mean
    global_mean, 
    
    // The outleier component likelihood
    t_ll;
  
  // The dataset covariance
  arma::mat global_cov_inv;
  
  double 
    // T degrees of freedom
    df = 4.0, 
    
    // Hyperparameters  
    u = 2.0, v = 10.0, b = 0.0, 
    
    // Outlier component weight
    outlier_weight = 0.0;
  
  using mixture::mixture;
  
  tAdjustedMixture(arma::uword _K,
                   arma::uvec _labels,
                   arma::mat _X
  ) ;
  
  // Destructor
  virtual ~tAdjustedMixture() { };
  
  // double calcTdistnLikelihood(arma::vec point) {
  double calcTdistnLikelihood(arma::uword n);
  
  void updateOutlierWeights();
  
  // void updateWeights(){
  //
  //   double a = 0.0;
  //
  //   for (arma::uword k = 0; k < K; k++) {
  //
  //     // Find how many labels have the value
  //     members.col(k) = (labels == k) % outliers;
  //     N_k(k) = arma::sum(members.col(k));
  //
  //     // Update weights by sampling from a Gamma distribution
  //     a  = concentration(k) + N_k(k);
  //     w(k) = arma::randg( arma::distr_param(a, 1.0) );
  //   }
  //
  //   // Convert the cluster weights (previously gamma distributed) to Beta
  //   // distributed by normalising
  //   w = w / arma::sum(w);
  // };
  
  arma::uword sampleOutlier(arma::uword n);
  
  
  
  virtual void updateAllocation(arma::vec weights, arma::mat upweigths);
};

#endif /* TADJMIXTURE_H */