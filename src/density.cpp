// density.cpp
// =============================================================================
// included dependencies
# include "density.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;

// =============================================================================
// virtual density class

// Parametrised class
density::density(
  arma::uword _K,
  arma::uvec _labels,
  arma::mat _X)
{
  
  K = _K;
  labels = _labels;
  
  X = _X;
  X_t = X.t();
  
  // Dimensions
  N = X.n_rows;
  P = X.n_cols;
  
  // Class populations
  N_k = zeros<uvec>(K);
  
  // N_inds.set_size(N);
  // P_inds.set_size(P);
  N_inds = linspace< uvec >(0, N - 1, N);
  P_inds = linspace< uvec >(0, P - 1, P);
  K_inds = linspace< uvec >(0, K - 1, K);
  
  // Log likelihood (individual and model)
  ll = zeros<vec>(K);
  likelihood = zeros<vec>(N);

};

void density::sampleParameters(arma::umat members, arma::uvec non_outliers) {
  // for(uword k = 0; k < K; k++) {
  std::for_each(
    std::execution::par,
    K_inds.begin(),
    K_inds.end(),
    [&](uword k) {
      sampleKthComponentParameters(k, members, non_outliers);
    }
  );
};
