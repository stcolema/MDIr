
# include "mixtureFactory.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace std ;
  
unique_ptr<mixture> mixtureFactory::createMixture(mixtureType type,
                                              arma::uword K,
                                              arma::uvec labels,
                                              arma::mat X) {
  switch (type) {
  // case G: return std::make_unique<gaussianSampler>(K, labels, concentration, X);
  case MVN: return std::make_unique<mvnMixture>(K, labels, X);
    // case C: return std::make_unique<categoricalSampler>(K, labels, concentration, X);
  case TMVN: return std::make_unique<tagmMixture>(K, labels, X);
    // case TG: return std::make_unique<tagmGaussian>(K, labels, concentration, X);
  default: throw std::invalid_argument( "invalid sampler type." );
  }
};
