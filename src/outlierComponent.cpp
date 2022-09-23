// outlierComponent.cpp
// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "outlierComponent.h"
# include "genericFunctions.h"

using namespace arma ;

// =============================================================================
// virtual outlierComponent class

// Parametrised class
outlierComponent::outlierComponent(arma::uvec _fixed, arma::mat _X) {
  
  X = _X;
  X_t = X.t();
  
  outliers = 1 - _fixed;
  non_outliers = _fixed;
  
  N = X.n_rows;
  P = X.n_cols;
  
  // Initialise the outlier log-likelihood vector
  outlier_likelihood = zeros< vec >(N);
  outlier_likelihood.fill(-DBL_MAX);
  
  updateWeights(non_outliers, outliers);
};


void outlierComponent::calculateAllLogLikelihoods() {
  
  for(uword n = 0; n < N; n++) {
    outlier_likelihood(n) = calculateItemLogLikelihood(X_t.col(n));
  }
}
  
void outlierComponent::updateWeights(uvec non_outliers, uvec outliers) {

  // Parameters are number of items allocated as outlier/non-outlier
  tau_1 = (double) sum(non_outliers);
  tau_2 = (double) sum(outliers);
  
  // Sample values for the weights
  // outlier_weight = rBeta(tau_2 + u, N + v - tau_2);
  non_outlier_weight = rBeta(tau_1 + v, tau_2 + u);
  outlier_weight = 1.0 - non_outlier_weight;
};

arma::uword outlierComponent::sampleOutlier(double non_outlier_likelihood_n,
                                            double outlier_likelihood_n) {
  
  uword pred_outlier = 0;
  // arma::uword k = labels(n);
  
  arma::vec outlier_prob(2);
  outlier_prob.zeros();

  // The likelihood of the point in the current cluster and in the outlier 
  // subcomponent
  outlier_prob(0) = non_outlier_likelihood_n + log(non_outlier_weight);
  outlier_prob(1) = outlier_likelihood_n + log(outlier_weight);

  // Convert from log-likelihoods to likelihoods while allowing for overflow
  outlier_prob = exp(outlier_prob - max(outlier_prob));

  // Normalise - convert form likelihoods to probabilities
  outlier_prob = outlier_prob / sum(outlier_prob);
  
  // Prediction and update
  u = arma::randu<double>( );
  pred_outlier = sum(u > cumsum(outlier_prob));
  
  return pred_outlier;
  
};
