// density.cpp
// =============================================================================
// included dependencies
# include "density.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;

// =============================================================================
// virtual density class

density::density(
  arma::uword _K,
  arma::uvec _labels,
  arma::mat _X)
{
  K = _K;
  labels = _labels;
  X = _X;

  N = X.n_rows;
  P = X.n_cols;

  N_k = zeros<uvec>(K);
  K_inds = linspace< uvec >(0, K - 1, K);
};

void density::sampleParameters(const arma::umat& members, const arma::uvec& non_outliers) {
  for(uword k = 0; k < K; k++) {
    sampleKthComponentParameters(k, members, non_outliers);
  }
};

void density::identifyMissingValues() {
  has_missing.set_size(N, P);
  has_missing.zeros();
  missing_indices.set_size(N);
  observed_indices.set_size(N);

  for(uword n = 0; n < N; n++) {
    for(uword p = 0; p < P; p++) {
      has_missing(n, p) = std::isfinite(X(n, p)) ? 0 : 1;
    }
    missing_indices(n) = arma::find(has_missing.row(n) == 1);
    observed_indices(n) = arma::find(has_missing.row(n) == 0);
  }
}

void density::initializeMissingValues() {
  for(uword p = 0; p < P; p++) {
    const arma::uvec finite_indices = arma::find_finite(X.col(p));
    double col_mean = 0.0, col_sd = 1.0;
    if(finite_indices.n_elem > 1) {
      const arma::vec col_p = X.col(p);
      const arma::vec observed = col_p.elem(finite_indices);
      col_mean = arma::mean(observed);
      col_sd = arma::stddev(observed);
      if(!(col_sd > 0.0)) {
        col_sd = 1.0;
      }
    } else if(finite_indices.n_elem == 1) {
      col_mean = X(finite_indices(0), p);
    }
    for(uword n = 0; n < N; n++) {
      if(has_missing(n, p) == 1) {
        X(n, p) = col_mean + 0.5 * col_sd * arma::randn();
      }
    }
  }
}
