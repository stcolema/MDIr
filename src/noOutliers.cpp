// noOutliers.cpp
// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "noOutliers.h"

using namespace arma ;

// =============================================================================
// noOutliers class of outlier component

// Parametrised class
noOutliers::noOutliers(arma::uvec _fixed, arma::mat _X) : outlierComponent(_fixed, _X) {
  
  // No outliers
  outliers = zeros<uvec>(N);
  non_outliers = ones<uvec>(N);
  
  outlier_likelihood.set_size(N);
  outlier_likelihood.zeros();
  
};


void noOutliers::initializeMissingValues() {
  // No-op: This component represents "no outliers", so no missing values to initialize
  // All missing value initialization is handled by the main density component
  return;
}

void noOutliers::sampleMissingForObservation(arma::uword n) {
  // No-op: This component represents "no outliers", so no missing values to sample
  // All missing value sampling is handled by the main density component
  return;
}

double noOutliers::calculateItemLogLikelihood(arma::uword n) {
  // Return extremely low likelihood for any observation being an outlier
  // This maintains the "no outliers" behavior
  return -DBL_MAX;
}

// double noOutliers::calculateItemLogLikelihood(arma::vec x) {
//   return 0.0;
// };

arma::uword noOutliers::sampleOutlier(double non_outlier_likelihood_n,
                          double outlier_likelihood_n) {
  return 0;
};

void noOutliers::updateWeights(uvec non_outliers, uvec outliers) {
};

void noOutliers::calculateAllLogLikelihoods() {
  outlier_likelihood.zeros(); 
};