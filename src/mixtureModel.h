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

class mixtureModel {
  
private:
  
public:
  
  uword 
    
    // The type of density modelled, must be 1 or 2
    mixture_type,
    
    // The type of outlier component modelled, must be 0 or 1
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
    
    // Vector of ones used ocasionally
    vec_of_ones,
    fixed,
    fixed_ind,
    unfixed_ind,
    
    // Outlier vectors (not really used here, declared so can be accessed in MDI class)
    outliers,
    non_outliers,
    
    // Used for looping over data indices
    N_inds;
  
  vec 
    // Concentration hyperparameter for ocmponent weights
    concentration, 
    
    // Component weights
    w,
    
    // The log-likelihood of an item in each component
    ll, 
    
    // The contribution of each item to the complete log likelihood
    complete_likelihood_vec, 
    
    // The contribution of each item to the observed log likelihood
    observed_likelihood_vec, 
    
    likelihood, 
    
    // Log-likelihood of being non-outlier or outlier
    outlier_likelihood;
  umat members;
  mat X, X_t, alloc;
  
  // Create a unique_ptr 
  // std::unique_ptr<density> density_ptr = std::make_unique<density>();
  // std::unique_ptr<outlierComponent> outlierComponent_ptr = std::make_unique<outlierComponent>();
  std::unique_ptr<density> density_ptr;
  std::unique_ptr<outlierComponent> outlierComponent_ptr;
  
  // Parametrised class
  mixtureModel(
    arma::uword _mixture_type,
    arma::uword _outlier_type,
    arma::uword _K,
    arma::uvec _labels,
    arma::uvec _fixed,
    arma::mat _X);
  
  
  // Destructor
  virtual ~mixtureModel() { };
  
  void updateAllocation(arma::vec log_weights, arma::mat log_upweigths);
  void updateItemAllocation(uword n, vec log_weights, vec log_upweigths);
  
  arma::uword sampleOutlier(
    double non_outlier_likelihood_n,
    double outlier_likelihood_n
  );
  
  void updateOutlierWeights();
  
  // Initialise the density, outlier component and mixture model
  void initialiseDensity(arma::uword type);
  void initialiseOutlierComponent(arma::uword type);
  void initialiseMixture(arma::vec log_weights, arma::mat log_upweigths);
  
  // The functions collected from the density
  void sampleFromPriors();
  void sampleParameters();
  void calcBIC();
  arma::vec itemLogLikelihood(arma::vec x);
  double logLikelihood(arma::vec x, arma::uword k);
  
};

#endif /* MIXTUREMODEL_H */