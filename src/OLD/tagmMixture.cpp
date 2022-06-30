
# include "tagmMixture.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;


tagmMixture::tagmMixture(arma::uword _K,
  arma::uvec _labels,
  // arma::vec _concentration,
  arma::mat _X
) : mvnMixture(_K, _labels, _X),
tAdjustedMixture(_K, _labels, _X),
mixture(_K, _labels, _X){
};
  
void tagmMixture::calcBIC(){
    
  // The outlier component introduces an additional n_param parameters (arguably 
  // + 1 for the d.o.f)
  BIC = 2 * complete_likelihood - n_param * (K_occ + 1) * std::log(N);
  
  // BIC = n_param * std::log(N) - 2 * model_likelihood;
  
};