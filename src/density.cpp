// density.cpp
// =============================================================================
// included dependencies
# include "density.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;

// =============================================================================
// virtual density class

arma::vec resolveDensityPrior(const arma::vec& prior) {
  // scale_pool_shape (0 switches pooling of variance scales off), gp_pool (0/1), 
  // gp_min_length (grid units), gp_center_sd, gp_pool_sd_scale
  arma::vec out = {2.0, 1.0, 1.0, 2.0, 1.0};
  if(prior.n_elem == 0) {
    return out;
  }
  if(prior.n_elem != 5 || !prior.is_finite() || prior(0) < 0.0 
     || (prior(1) != 0.0 && prior(1) != 1.0) || prior(2) <= 0.0 
     || prior(3) <= 0.0 || prior(4) <= 0.0) {
    Rcpp::stop("Invalid density prior options; see densityPrior().");
  }
  return prior;
}

density::density(
  arma::uword _K,
  arma::uvec _labels,
  arma::mat _X,
  arma::vec _density_prior)
{
  density_prior = resolveDensityPrior(_density_prior);
  K = _K;
  labels = _labels;
  X = _X;

  N = X.n_rows;
  P = X.n_cols;

  N_k = zeros<uvec>(K);
  K_inds = linspace< uvec >(0, K - 1, K);
};

void density::sampleParameters(const arma::umat& members, const arma::uvec& non_outliers) {
  arma::uvec is_occupied(K, arma::fill::zeros);
  for(uword k = 0; k < K; k++) {
    is_occupied(k) = (accu((members.col(k) == 1) && (non_outliers == 1)) > 0) ? 1 : 0;
  }
  const arma::uvec occupied = find(is_occupied == 1);
  
  for(uword i = 0; i < occupied.n_elem; i++) {
    sampleKthComponentParameters(occupied(i), members, non_outliers);
  }
  updatePooledHyperparameters(occupied);
  for(uword k = 0; k < K; k++) {
    if(is_occupied(k) == 0) {
      sampleKthComponentParameters(k, members, non_outliers);
    }
  }
};

void density::replaceData(const arma::mat& X_new) {
  if(X_new.n_cols != P) {
    Rcpp::stop("The new data must have %d columns.", (int) P);
  }
  X = X_new;
  N = X.n_rows;
  identifyMissingValues();
}

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
