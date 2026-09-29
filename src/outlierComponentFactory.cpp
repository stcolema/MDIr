// outlierComponentFactory.cpp
// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "outlierComponentFactory.h"

using namespace arma ;

// =============================================================================
// virtual outlierComponentFactory class

// empty contructor
outlierComponentFactory::outlierComponentFactory() { };
outlierComponentFactory::outlierComponentFactory(const outlierComponentFactory &L) { };

std::unique_ptr<outlierComponent> outlierComponentFactory::createOutlierComponent(
    outlierType type, arma::uvec fixed, arma::mat X,
    const arma::field<arma::uvec>* miss_idx,
    const arma::field<arma::uvec>* obs_idx
) {
  // Rcpp::Rcout << "\nMake component.";
  
  switch (type) {
  case E: {
    // Rcpp::Rcout << "\nMake empty component.";
    return std::unique_ptr<noOutliers>(new noOutliers(fixed, X, miss_idx, obs_idx));
  }
  case MVT: {
    // Rcpp::Rcout << "\nMake MVT component.";
    return std::unique_ptr<mvt>(new mvt(fixed, X, miss_idx, obs_idx));
  }
  default : {
    // Rcpp::Rcout << "\nThrow an error.";
    Rcpp::stop("invalid outlier type.");
  }
  }
};

