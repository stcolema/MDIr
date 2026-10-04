// mixtureModel.h
// =============================================================================
// include guard
#ifndef MIXTUREMODEL_H
#define MIXTUREMODEL_H

// =============================================================================
// included dependencies
# include "density.h"
# include "outlierComponent.h"
# include "densityFactory.h"
# include "outlierComponentFactory.h"

using namespace arma ;

// =============================================================================
// mixtureModel class
//
// One view of the MDI model: a finite mixture of K components of a chosen
// density, optionally with a global outlier component. The MDI class supplies
// the component weights and the cross-view upweights that enter the allocation.
class mixtureModel {
  
public:
  
  uword 
    // The type of density modelled
    mixture_type,
    
    // The type of outlier component modelled, 0 (none) or 1 (multivariate t)
    outlier_type = 0,
    
    // The number of components modelled
    K, 
    
    // The number of components occupied (i.e. clusters/groups)
    K_occ, 
    
    // The dimensions of the dataset, samples and columns respectively
    N, P, 
    
    // The number of observed labels
    N_fixed = 0,
    
    n_param = 0;
  
  double complete_likelihood = 0.0, observed_likelihood = 0.0, BIC = 0.0;
  
  uvec 
    // The cluster/class labels
    labels, 
    
    // The number of items in each class
    N_k, 
    
    // Indicators of items with an observed label
    fixed,
    
    // Outlier indicators
    outliers,
    non_outliers;
  
  vec 
    // The contribution of each item to the complete-data and observed-data 
    // log-likelihoods
    complete_likelihood_vec, 
    observed_likelihood_vec;
  
  umat members;
  
  // Allocation probabilities of each item to each component
  mat alloc;
  
  std::unique_ptr<density> density_ptr;
  std::unique_ptr<outlierComponent> outlierComponent_ptr;
  
  mixtureModel(
    arma::uword _mixture_type,
    arma::uword _outlier_type,
    arma::uword _K,
    arma::uvec _labels,
    arma::uvec _fixed,
    arma::mat _X,
    arma::vec _density_prior);
  
  virtual ~mixtureModel() { };
  
  // One sweep: outlier weight, then (label, outlier status) for each item, then
  // the missing values given the new allocation.
  void updateAllocation(const arma::vec& log_weights, const arma::mat& log_upweights);
  void updateItemAllocation(uword n, const arma::vec& log_weights, const arma::vec& log_upweights);
  
  void updateOutlierWeights();
  
  void initialiseDensity(arma::uword type, const arma::mat& X, const arma::vec& density_prior);
  void initialiseOutlierComponent(arma::uword type, const arma::mat& X);
  void initialiseMixture(const arma::vec& log_weights, const arma::mat& log_upweights);
  
  void sampleFromPriors();
  void sampleParameters();
  void calcBIC();
  
  // Sample every missing value from its full conditional given the current 
  // allocation
  void sampleAllMissingValues();
  
  // Exchange components k and k'
  void swapComponents(uword k, uword kprime);
  
  // Set the inverse temperature of the tempered target (see mdi::beta). The
  // allocation uses beta times the component log-likelihood and the density
  // draws its parameters from the tempered conditional.
  void setBeta(double beta_new);
  double beta = 1.0;
  
  // The data with missing entries replaced by the current imputation
  const arma::mat& imputedData() const { return density_ptr->X; }
  
};

#endif /* MIXTUREMODEL_H */
