// mixture.h
// =============================================================================
  // include guard
#ifndef MIXTURE_H
#define MIXTURE_H

// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
using namespace arma ;
  
// =============================================================================
// virtual mixture class

class mixture {
  
private:
  
public:
  
  uword K, N, P, K_occ;
  double complete_likelihood = 0.0, observed_likelihood = 0.0, BIC = 0.0;
  uvec labels, 
    N_k, 
    batch_vec, 
    N_b, 
    outliers, 
    non_outliers, 
    vec_of_ones,
    fixed,
    fixed_ind,
    unfixed_ind;
    
  vec concentration, w, ll, likelihood;
  umat members;
  mat X, X_t, alloc;

  // Parametrised class
  mixture(
    arma::uword _K,
    arma::uvec _labels,
    // arma::vec _concentration,
    arma::mat _X);


  // Destructor
  virtual ~mixture() { };
  
  virtual void updateAllocation(arma::vec weights, arma::mat upweigths);
    
  // The virtual functions that will be defined in every subclasses
  virtual void sampleFromPriors() = 0;
  virtual void sampleParameters() = 0;
  virtual void calcBIC() = 0;
  virtual arma::vec itemLogLikelihood(arma::vec x) = 0;
  virtual double logLikelihood(arma::vec x, arma::uword k) = 0;
  
  // // Not every class needs to save matrix combinations, so this is not purely
  // // virtual
  // virtual void matrixCombinations() {};
};
  
#endif /* MIXTURE_H */