// outlierComponent.cpp
// =============================================================================
# include <RcppArmadillo.h>
# include "outlierComponent.h"
# include "genericFunctions.h"

using namespace arma ;

outlierComponent::outlierComponent(
    arma::uvec _fixed, 
    arma::mat _X,
    const arma::field<arma::uvec>* miss_idx,
    const arma::field<arma::uvec>* obs_idx
) {
  X = _X;
  missing_indices_ref = miss_idx;
  observed_indices_ref = obs_idx;
  
  // Every item starts as a non-outlier so that the first parameter update sees
  // all of the data; outlier status is resampled at every allocation step.
  outliers = zeros<uvec>(_X.n_rows);
  non_outliers = ones<uvec>(_X.n_rows);
  
  N = X.n_rows;
  P = X.n_cols;
  
  outlier_likelihood = zeros< vec >(N);
  updateWeights(non_outliers, outliers);
};

void outlierComponent::calculateAllLogLikelihoods() {
  for(uword n = 0; n < N; n++) {
    outlier_likelihood(n) = calculateItemLogLikelihood(n);
  }
}

void outlierComponent::updateWeights(const uvec& non_outliers, const uvec& outliers) {
  const double tau_1 = (double) sum(non_outliers);
  const double tau_2 = (double) sum(outliers);
  non_outlier_weight = rBeta(tau_1 + outlier_prior_b, tau_2 + outlier_prior_a);
  outlier_weight = 1.0 - non_outlier_weight;
};
