// mixtureFactory.h
// =============================================================================
// include guard
#ifndef MIXTUREFACTORY_H
#define MIXTUREFACTORY_H

// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "mixture.h"
# include "mvnMixture.h"
# include "tagmMixture.h"

using namespace std ;

// =============================================================================
// mixtureFactory class


// Factory for creating instances of samplers
//' @name mixtureFactory
//' @title Factory for different types of mixtures.
//' @description The factory allows the type of mixture implemented to change
//' based upon the user input.
//' @field new Constructor \itemize{
//' \item Parameter: samplerType - the density type to be modelled
//' \item Parameter: K - the number of components to model
//' \item Parameter: labels - the initial clustering of the data
//' \item Parameter: concentration - the vector for the prior concentration of
//' the Dirichlet distribution of the component weights
//' \item Parameter: X - the data to model
//' }
class mixtureFactory
{
public:
  
  // empty contructor
  mixtureFactory(){ };
  
  // destructor
  virtual ~mixtureFactory(){ };
  
  // copy constructor
  mixtureFactory(const mixtureFactory &L);
  
  enum mixtureType {
    // G = 0,
    MVN = 1,
    // C = 2,
    TMVN = 3 //,
    // TG = 4
  };
  
  static unique_ptr<mixture> createMixture(mixtureType type,
                                                arma::uword K,
                                                arma::uvec labels,
                                                arma::mat X) ;
  
};


#endif /* MIXTUREFACTORY_H */