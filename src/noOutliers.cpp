// noOutliers.cpp
// =============================================================================
# include <RcppArmadillo.h>
# include "noOutliers.h"

using namespace arma ;

noOutliers::noOutliers(
    arma::uvec _fixed, 
    arma::mat _X,
    const arma::field<arma::uvec>* miss_idx,
    const arma::field<arma::uvec>* obs_idx
) : outlierComponent(_fixed, _X, miss_idx, obs_idx) {
  outliers = zeros<uvec>(N);
  non_outliers = ones<uvec>(N);
  non_outlier_weight = 1.0;
  outlier_weight = 0.0;
  outlier_likelihood.fill(-arma::datum::inf);
};

double noOutliers::calculateItemLogLikelihood(arma::uword n) {
  return -arma::datum::inf;
}

arma::vec noOutliers::sampleMissingValues(arma::uword n) const {
  Rcpp::stop("noOutliers cannot impute values.");
  return arma::vec();
}

arma::vec noOutliers::simulate() const {
  Rcpp::stop("noOutliers cannot simulate values.");
  return arma::vec();
}

void noOutliers::updateWeights(const uvec& non_outliers, const uvec& outliers) { }

void noOutliers::sampleFromPrior() { }
